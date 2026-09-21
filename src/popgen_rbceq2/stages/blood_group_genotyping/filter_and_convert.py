"""Turns a DRAGEN gVCF into the VCF rbceq2 reads, and extracts the sites the QC stage judges."""

import typing

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import constants, off_design, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import posthoc_genotype

# Which sites of a BED fall inside no span of another BED, as a BED.
#
# `awk -v spans=spans.bed <this> spans.bed sites.bed` prints the sites no span contains. A span
# is `[$2, $3)`, 0-based half-open. The spans are always the covered spans of the DRAGEN records
# overlapping a defining site: bcftools' %END is POS + rlen - 1, which is INFO/END for a
# reference block and POS + len(REF) - 1 for a variant, so `%POS0\t%END` is that span, and this
# containment test is the one `rbceq2_call_qc_job.GvcfRecord.covers` applies. The two
# definitions of "covered" have to agree — a site this calls uncovered is a site the QC would
# flag NOCOV — and tests/test_posthoc_merge.py holds them together.
#
# It also served the capture-design subtraction once, over the same containment test. That half
# is now `bedtools intersect -v` in scripts/gen_off_design_sites.py, run once per design and
# committed under resources/, so the vendor BED's extra columns and `track` lines are no longer
# this program's problem. What is left here reads only bcftools output, three columns.
#
# The spans file is named by `-v spans=`, and must be, rather than being detected with the
# usual `NR == FNR`. That idiom reads "still in the first file" only while the first file has
# produced records, so an *empty* spans file makes awk take the sites file for the span list
# and print nothing at all. Every caller below would then subtract the sites from themselves:
# an empty covered.bed, which is what a gVCF with no record at any defining site produces,
# would fill no holes and fail later blaming contig naming. `ARGIND` would say the same thing
# more directly but is a gawk extension, and this runs under mawk.
#
# Held as a plain string rather than inlined: the command below is an f-string, and every
# brace in an awk program would have to be doubled.
_SITES_OUTSIDE_SPANS_AWK = """
FILENAME == spans { n = ++c[$1]; lo[$1, n] = $2 + 0; hi[$1, n] = $3 + 0; next }
{ for (i = 1; i <= c[$1]; i++) if (lo[$1, i] <= $2 + 0 && $2 + 0 < hi[$1, i]) next; print }
"""

# Stamp INFO/POSTHOC on every record of the post-hoc supplement. `.` is the empty INFO column
# and has to be replaced rather than appended to, or the record grows a `.;POSTHOC=...` INFO
# that no parser accepts.
_TAG_POSTHOC_AWK = """
/^#/ { print; next }
{ $8 = ($8 == "." ? tag : $8 ";" tag); print }
"""

# The INFO field naming the caller that supplied a record where the primary gVCF was silent.
#
# Declared on the intermediate for *every* run, genome included, and not only where records
# carry it. `bcftools query -f '%INFO/POSTHOC'` fails outright on a tag the header does not
# define ("no such tag defined in the VCF header"); it renders `.` only for a declared tag a
# record happens to lack. Without the unconditional declaration the extract would abort on
# every genome run, and on any exome whose gVCF turned out to have no holes to fill.
_POSTHOC_HEADER_LINE = (
    '##INFO=<ID=POSTHOC,Number=1,Type=String,Description='
    '"Caller that supplied this record at a site the primary gVCF had no record for">'
)


# The INFO flag marking a post-hoc record whose span reaches a defining site it may not fill:
# one the primary caller already has a record for, or a hole inside the capture design. Set
# inside the merge and stripped again before the supplement is concatenated, so it reaches
# neither the merged file nor the extract.
_COVERED_HEADER_LINE = (
    '##INFO=<ID=COVERED,Number=0,Type=Flag,Description="Overlaps a defining site the post-hoc caller may not fill">'
)


# The `bcftools query` format the QC extract is written in. One field per entry of
# `rbceq2_call_qc_job.EXTRACT_COLUMNS`, in the same order; test_posthoc_merge holds the two to
# the same field count, since a mismatch would be caught only when a job ran and the parser
# rejected the file.
_EXTRACT_FORMAT = r'%CHROM\t%POS\t%REF\t%ALT\t%INFO/END\t[%GT\t%DP\t%GQ\t%MIN_DP]\t%INFO/POSTHOC\n'


def _primary_records_guard(sites_bed: str, genome: str) -> str:
    """Shell that fails the job unless a DRAGEN record overlaps some blood-group defining site.

    Reads `dragen.vcf.gz`, the DRAGEN-only intermediate, and nothing else. An empty answer means
    the input is wrong rather than the sample poor: a gVCF called against another build, or one
    naming its contigs `1` where the BEDs say `chr1`, matches no region and `bcftools norm -R`
    reports that as zero records with exit 0.

    It runs before the merge on purpose. The post-hoc caller takes its contigs from the
    reference fasta, not from the gVCF, so its records are well-formed whatever the gVCF looks
    like. Asked of the merged file, this question is answered by the supplement alone: a
    wrong-build exome gVCF yields an empty DRAGEN set, every off-design site becomes a hole,
    the trespass check passes with no DRAGEN blocks to trip it, and the post-hoc records fill
    the holes and populate the extract. That sample is then typed from the second caller with
    everything else NOCOV, and looks like a poor exome rather than a broken input.

    Args:
        sites_bed: Localised `bg_defining_sites.<genome>.bed`.
        genome: The configured genome build, for the message.

    Returns:
        Shell lines, indented for the job command.
    """
    return f"""
            bcftools query -T {sites_bed} --targets-overlap 2 -f '%POS\\n' dragen.vcf.gz > dragen_at_sites.txt
            if [ ! -s dragen_at_sites.txt ]; then
                echo "ERROR: no DRAGEN gVCF record overlaps any blood-group defining site." >&2
                echo "Check the gVCF contig naming, and that references.genome_build" >&2
                echo "({genome}) matches the build the gVCF was called against." >&2
                exit 1
            fi"""


def _sample_check_commands(posthoc_gvcf: str, gvcf: str) -> str:
    """Shell that fails the job if the post-hoc gVCF and the DRAGEN gVCF name different samples.

    The post-hoc caller takes its sample name from the CRAM's read group and DRAGEN named the
    gVCF from the same run, so a mismatch means the CRAM and the gVCF this sequencing group
    resolves to do not describe one individual. Merging them would splice another person's
    genotypes into this one's calls at exactly the sites nothing else covers, and the result
    would look like an ordinary recovery.

    Relabelling the supplement instead would satisfy `concat`, which requires identical sample
    sets, and would bury that. Some mackenzie DRAGEN 3.7.8 test CRAMs do trip this benignly,
    carrying a retired sequencing-group ID for the same individual from an upstream test-set
    reheadering bug that is fixed for newer additions. That is a reason to fix those inputs,
    not to weaken the check for every cohort: this is the only place the pipeline compares the
    two files it was handed, and a real swap and a stale header are indistinguishable from
    here.

    A separate fragment from the merge so the stage can run it first in the job, before the
    norm pass that builds the DRAGEN intermediate, and whether or not this sample has a hole
    to fill: the check is about the two inputs, not about the data in them, so it reads both
    from the inputs, depends on nothing computed, and fails in seconds.

    Args:
        posthoc_gvcf: Localised post-hoc gVCF from PosthocGenotypeOffTargetSites.
        gvcf: The localised DRAGEN gVCF, as handed to the job.

    Returns:
        The shell fragment, for interpolation at the top of the stage's command.
    """
    return f"""
        posthoc_sample=$(bcftools query -l {posthoc_gvcf})
        dragen_sample=$(bcftools query -l {gvcf})
        if [ "$posthoc_sample" != "$dragen_sample" ]; then
            echo "ERROR: the post-hoc calls and the gVCF name different samples." >&2
            echo "  CRAM/post-hoc: $posthoc_sample" >&2
            echo "  primary gVCF:  $dragen_sample" >&2
            echo "The CRAM and gVCF registered for this sequencing group are not from one" >&2
            echo "DRAGEN run of one individual. Either the CRAM is registered against the" >&2
            echo "wrong sequencing group, or its read group was never updated to the" >&2
            echo "current ID. Check somalier, then fix the input; do not merge." >&2
            exit 1
        fi
    """


def _convert_commands(out_vcf: str, cpu: int) -> str:
    """Shell turning `merged.vcf.gz` into the VCF rbceq2 reads, and indexing it.

    Two things happen and no more. `<NON_REF>` records are dropped, because the symbolic allele
    breaks rbceq2, which does not accept native gVCFs; and the ALT alleles left unused by the
    earlier split are trimmed. **Nothing touches FORMAT/GT.**

    Separate from the stage, and tested against real bcftools, so that last part is assertable:
    a one-token GT on single-copy chrX must reach rbceq2 as one token. Keeping the helper free
    of any ploidy rewrite is what `tests/test_haploid_passthrough.py` pins.

    The no-rewrite guarantee is this helper's alone. `merged.vcf.gz` may hold a second caller's
    records; reconciling their ploidy is `_merge_posthoc_commands`' job.

    Args:
        out_vcf: Path to write the bgzipped converted VCF to. Its `.tbi` goes beside it.
        cpu: Thread count for the BGZF deflation and for reading it back to index.

    Returns:
        Shell, indented to sit inside the stage's command block.
    """
    return f"""
            bcftools view \\
                    -e 'ALT="<NON_REF>"' \\
                    --trim-alt-alleles --threads {cpu} -Oz -o {out_vcf} merged.vcf.gz
            bcftools index -t --threads {cpu} {out_vcf}"""


def _merge_posthoc_commands(
    posthoc_gvcf: str,
    sites_bed: str,
    off_design_bed: str,
    design_key: str,
    cpu: int,
    genome: str,
) -> str:
    """Shell to fill the primary gVCF's blind spots from the post-hoc caller's gVCF.

    Reads `dragen.vcf.gz` and writes `merged.vcf.gz`, which is what the extract and the
    conversion then run over.

    The fill is empirical and per sample, restricted to the capture design: a post-hoc record
    survives only at a defining site that lies outside the design's intervals *and* that the
    DRAGEN gVCF has no record covering. DRAGEN wins wherever both speak, and a hole inside the
    design is left as it is, to reach the QC as NOCOV, because a site the design targeted and
    DRAGEN still said nothing about is a question about that sample's DRAGEN run, not one a
    second caller should answer.

    The design BED has to be the one the gVCF was called against, and this is checked: DRAGEN
    emits reference blocks over exactly the target BED (no padding), so a DRAGEN reference
    block reaching a defining site outside the configured design means the wrong file is
    configured — on the validation cohorts a target-regions file in place of the probe-footprint
    one leaves 309 defining sites "outside" the design with DRAGEN records, and would have gated
    off three quarters of the recoveries while the run looked fine. That case fails the job. A
    DRAGEN *variant* is not proof of anything: it is anchored inside the target, but its REF can
    run past the edge, so a deletion at a capture edge can legitimately cover an off-design
    site. Such a site is covered and so not filled, and the QC reports DEL from DRAGEN's record.

    A record is selected for reaching a fillable hole and is selected whole, so one anchored in
    a hole can extend over a neighbouring defining site it may not fill: one DRAGEN did call,
    or a hole inside the design. Blood-group defining sites are dense — most have another
    within 20bp — so this is the ordinary case near a capture edge, not a corner one, and what
    happens next depends on what the record is.

    A post-hoc **variant** that reaches a site it may not fill is dropped. At a called base,
    keeping it would hand rbceq2 two callers' alleles at one base with nothing to choose
    between them: a DRAGEN SNP at a site, and a post-hoc deletion removing it. The QC could not
    report the conflict either, because `resolve_coverage` prefers the primary record and so
    reads the site as an ordinary PASS. At an in-design hole, keeping it would let the second
    caller decide the genotype at a site the design targeted, which is exactly what the design
    bound exists to prevent. Either way the hole the dropped record would have filled goes back
    to reaching the QC as NOCOV, which is the honest answer.

    A post-hoc **reference block** that reaches a site it may not fill is kept whole. It
    asserts nothing rbceq2 ever sees, since the conversion drops every <NON_REF>-only record
    before rbceq2 reads the file, and dropping the block instead would throw away the hole it
    was kept for. That does leave the block covering the site in the extract. At a called base
    two records cover it and `resolve_coverage` prefers the one with no INFO/POSTHOC; at an
    in-design hole the block is the only record there, so the QC is handed the fillable sites
    (`FlagBloodGroupCallQc` reads the same off-design BED) and disregards a post-hoc record at
    any other site, which is what keeps such a hole NOCOV. Splitting blocks on the boundary
    would be the alternative and is not worth it.

    Assumes `_sample_check_commands` has already run: the two files name one sample.

    Args:
        posthoc_gvcf: Localised post-hoc gVCF from PosthocGenotypeOffTargetSites.
        sites_bed: The committed defining-sites BED.
        off_design_bed: Localised off-design defining sites for the configured design, the
            committed `off_design.resource_path`; the only sites a post-hoc record may fill.
        design_key: The `[references]` key the design came from, for the error message.
        cpu: Threads to give the BGZF steps.
        genome: The configured genome build, for the single-copy chrX bounds the fill is
            kept out of.

    Returns:
        The shell fragment, for interpolation into the stage's command.

    Raises:
        KeyError: `genome` has no recorded chrX PAR boundaries. See `constants.non_par_x`.
    """
    non_par_lo, non_par_hi = constants.non_par_x(genome)
    return f"""
        # --targets-overlap 1 here, not 2: this needs every record whose *span* reaches a
        # defining site, which is what %END reports and what the QC counts as covering. Mode 2
        # asks whether the *variant* overlaps, and drops a deletion anchored on the site
        # itself, which would make a covered site look like a hole and let a post-hoc record
        # displace a DRAGEN call.
        bcftools index -t --threads {cpu} dragen.vcf.gz

        bcftools query -T {sites_bed} --targets-overlap 1 \\
            -f '%CHROM\\t%POS0\\t%END\\n' dragen.vcf.gz > covered.bed

        # Which of the off-design sites the DRAGEN gVCF has no record at. Only those are
        # filled. The off-design set itself is not computed here: it depends on the configured
        # design and the committed sites and on nothing about this sample, so it is subtracted
        # once per design and committed under resources/.
        awk -v spans=covered.bed '{_SITES_OUTSIDE_SPANS_AWK}' \\
            covered.bed {off_design_bed} > uncovered.all.bed

        # Single-copy chrX is not filled. HaplotypeCaller runs at its default ploidy of 2 and
        # writes a two-token GT there, where DRAGEN writes one token for a male sample. Both in
        # one file tell rbceq2 the sample has one chromosome copy and two, and it reports the
        # blood group Undetermined; a het post-hoc call claims fewer copies than it has, passes
        # silently and flips the phenotype. See "Haploid genotypes on the sex chromosomes" in
        # the README.
        #
        # Drop the site, do not rewrite its genotype. A rewrite has to decide the sample's
        # ploidy from the DRAGEN side, and a wrong decision yields a confident wrong call where
        # this yields a NOCOV flag.
        #
        # Costs the off-design XK sites: 10 on Twist, 2 on Agilent CREv2. GATA1 and ATP11C have
        # no off-design site, so nothing else in single-copy chrX is affected.
        awk -v lo={non_par_lo} -v hi={non_par_hi} 'BEGIN{{FS=OFS="\\t"}} \\
            !($1 == "chrX" && $3 >= lo && $3 <= hi)' uncovered.all.bed > uncovered.bed
        n_haploid=$(($(wc -l < uncovered.all.bed) - $(wc -l < uncovered.bed)))
        if [ "$n_haploid" -gt 0 ]; then
            echo "post-hoc: $n_haploid off-design site(s) in single-copy chrX left unfilled, to" >&2
            echo "keep one ploidy in the merged VCF; they reach the QC as NOCOV" >&2
        fi

        # A DRAGEN reference block reaching a site outside the design means the configured BED
        # is not the one the gVCF was called against. The likely case is a target-regions file
        # where the gVCF used the probe footprint, which would quietly gate off most of the
        # recoveries. Blocks, not every record: a block asserts hom-ref over bases DRAGEN
        # evaluated, and those stop at the target edge, whereas a variant is anchored inside
        # the target and its REF can run past the edge. A deletion at a capture edge covering
        # an off-design site is a carrier, not a wrong file; it is covered above, so not
        # filled, and reaches the QC as DEL from the DRAGEN record. Nearly every site is
        # hom-ref, so a wrong design still trips this on almost all of its sites.
        bcftools query -T {off_design_bed} --targets-overlap 1 -i 'INFO/END!="."' \\
            -f '%CHROM\\t%POS0\\t%END\\n' dragen.vcf.gz > off_design_blocks.bed
        awk -v spans=off_design_blocks.bed '{_SITES_OUTSIDE_SPANS_AWK}' \\
            off_design_blocks.bed {off_design_bed} > off_design_outside_blocks.bed
        sort {off_design_bed} > off_design_sites.sorted.bed
        sort off_design_outside_blocks.bed > off_design_outside_blocks.sorted.bed
        comm -23 off_design_sites.sorted.bed off_design_outside_blocks.sorted.bed > off_design_in_blocks.bed
        if [ -s off_design_in_blocks.bed ]; then
            n_in_blocks=$(wc -l < off_design_in_blocks.bed | tr -d ' ')
            echo "ERROR: DRAGEN reference blocks cover $n_in_blocks defining site(s) outside the capture design." >&2
            echo "{stage_support.DESIGN_CONFIG_PATH} = {design_key} is not the BED this gVCF" >&2
            echo "was called against." >&2
            echo "For an Agilent design that usually means Regions configured where the gVCF" >&2
            echo "used Covered. First sites:" >&2
            head -n 10 off_design_in_blocks.bed >&2
            exit 1
        fi
        sort uncovered.bed > uncovered.sorted.bed

        awk -v spans=covered.bed '{_SITES_OUTSIDE_SPANS_AWK}' covered.bed {sites_bed} > holes.bed
        n_holes=$(wc -l < holes.bed | tr -d ' ')
        n_fill=$(wc -l < uncovered.bed | tr -d ' ')
        echo "post-hoc: $n_holes defining site(s) with no DRAGEN record; $n_fill outside the" >&2
        echo "capture design and filled, $((n_holes - n_fill)) inside it and left for the QC to flag NOCOV" >&2

        # An empty -T file is a hard error in bcftools ("Failed to read the targets"), so the
        # no-holes case has to branch rather than fall through the same pipeline.
        if [ -s uncovered.bed ]; then
            # The defining sites a post-hoc record may not fill: every site, less the fillable
            # holes. That is every site DRAGEN did cover, plus every hole inside the design.
            # Subtracting holes.bed instead would leave an in-design hole out of the mark, and
            # a record kept for the off-design hole beside it would fill both. The set is
            # deliberately narrower than covered.bed: marking against the DRAGEN records' whole
            # spans would also drop a post-hoc variant that merely overlaps the tail of a long
            # reference block, reaching no defining site and so contradicting nothing rbceq2
            # reads.
            sort {sites_bed} > sites.sorted.bed
            comm -23 sites.sorted.bed uncovered.sorted.bed \\
                | sort -k1,1 -k2,2n | bgzip -c --threads {cpu} > unfillable_sites.bed.gz
            tabix -p bed unfillable_sites.bed.gz
            echo '{_COVERED_HEADER_LINE}' > covered_hdr.txt

            # rbceq2 and the extract read GT, DP, GQ, MIN_DP and END; every other tag the
            # post-hoc caller emits is dropped here. That keeps the supplement to the fields
            # the pipeline actually reads, and keeps bcftools concat from having to reconcile
            # two callers' definitions of tags nothing downstream looks at.
            #
            # posthoc_hdr.txt is written by the caller, before dragen.vcf.gz is built, because
            # the primary intermediate needs the same declaration. Both sides carrying it is
            # what lets `concat` merge the two headers without a conflict.
            #
            # Zero-depth records are dropped. Given -L, the post-hoc caller emits a reference
            # block over the whole interval, including stretches with no reads at all, as
            # DP=0,GQ=0. Keeping those would put a record over every hole and so retire NOCOV
            # for exomes entirely: a site with no reads would read LOWQ(DP=0), which says
            # "poor data" where the truth is "no data". A record carrying no DP field at all
            # is kept, since absence of the field is not proof of absence of reads.
            #
            # `-T uncovered.bed --targets-overlap 2` selects a record whose span reaches a
            # hole, and selects it whole. Mode 2 asks whether the *variant* overlaps, but a
            # gVCF record is still multiallelic here — the real ALT plus <NON_REF> — and the
            # symbolic allele makes bcftools match on the whole record span, so a deletion
            # anchored on the hole is selected too. That is wanted; what it drags in is
            # handled by the COVERED mark below, not by the selection mode.
            #
            # INFO/COVERED marks every selected record whose span also reaches a defining
            # site it may not fill. `annotate -m` matches on the record's span, not on POS, so
            # a deletion anchored on a hole and reaching such a site one base away is marked,
            # which is the case this whole check exists for.
            #
            # FILTER is set here because HaplotypeCaller leaves it `.` on every record and
            # rbceq2 discards an allele whose defining variant is not literally PASS. Left as
            # `.`, every post-hoc alternate allele would be thrown away before genotyping and
            # each recovered site would be typed as reference by absence, while its QC flag
            # said the call rested on the recall. `filter -s` writes PASS on every record the
            # expression does not match, and the one filter DRAGEN applies to these cohorts
            # that means the same thing on both callers is its depth rule, LowDepth at DP<=1.
            # Its QUAL rule is not copied: DRAGEN's QUAL is ML-recalibrated and the two
            # callers' scales are not comparable, so a post-hoc variant is used at any QUAL
            # and its DP and GQ reach the QC flags, as a PASS DRAGEN variant's already do.
            bcftools view -T uncovered.bed --targets-overlap 2 -e 'FORMAT/DP=0' -Ou {posthoc_gvcf} \\
                | bcftools filter -s LowDepth -e 'FORMAT/DP<=1' -Ou - \\
                | bcftools annotate -x '^INFO/END,^FORMAT/GT,FORMAT/DP,FORMAT/GQ,FORMAT/MIN_DP' -Ou - \\
                | bcftools annotate -a unfillable_sites.bed.gz -h covered_hdr.txt -c CHROM,FROM,TO -m COVERED \\
                    -Ob -o posthoc_marked.bcf -

            # A marked variant is dropped: it asserts something about a defining site it may
            # not fill, in the file rbceq2 reads. A marked reference block is kept, see the
            # docstring. The two are told apart before `norm -m -any` splits anything, while a
            # gVCF variant is still `<real ALT>,<NON_REF>` and a block is `<NON_REF>` alone, so
            # the test is the alleles and nothing about INFO tags. Keying it on a missing
            # INFO/END would rest on HaplotypeCaller's habit of not writing END on a variant,
            # and a variant that carried one would survive; and after the split, a variant's
            # `<NON_REF>` twin inherits its REF span and any END, so it would pass as a block
            # and fill the hole with an apparent hom-ref call. Here there is no twin yet, and
            # the count is one per variant however many alleles `norm` later splits it into.
            n_trespass=$(bcftools view -H -i 'INFO/COVERED=1 && (N_ALT>1 || ALT!="<NON_REF>")' posthoc_marked.bcf \\
                | wc -l | tr -d ' ')
            bcftools view -e 'INFO/COVERED=1 && (N_ALT>1 || ALT!="<NON_REF>")' -Ou posthoc_marked.bcf \\
                | bcftools norm -m -any --threads {cpu} -Ou - \\
                | bcftools annotate -x INFO/COVERED -h posthoc_hdr.txt -Ov - \\
                | awk -v OFS='\\t' -v tag='POSTHOC={constants.POSTHOC_CALLER}' '{_TAG_POSTHOC_AWK}' \\
                | bgzip -c --threads {cpu} > posthoc_tagged.vcf.gz

            # What was actually merged, which is not the hole count reported above: most holes
            # yield no record at all (DP=0), and a kept variant contributes two records here,
            # itself and its <NON_REF> twin. Reporting only the holes hid the all-zero-depth
            # case entirely — a header-only supplement that `concat` merges silently.
            n_kept=$(bcftools view -H posthoc_tagged.vcf.gz | wc -l | tr -d ' ')
            echo "post-hoc: $n_kept record(s) kept over $n_fill hole(s); $n_trespass variant(s) dropped for" >&2
            echo "reaching a defining site they may not fill, the rest had no reads there (DP=0)" >&2

            bcftools index -t --threads {cpu} posthoc_tagged.vcf.gz
            bcftools concat -a --threads {cpu} -Oz -o merged.vcf.gz dragen.vcf.gz posthoc_tagged.vcf.gz
        else
            echo "post-hoc: no defining site outside the capture design lacks a DRAGEN record; nothing to fill" >&2
            mv dragen.vcf.gz merged.vcf.gz
        fi
    """


class FilterAndConvertGvcfsForRbceq2(cpg_flow.stage.SequencingGroupStage):
    """Convert a sequencing group's gVCF into the VCF rbceq2 reads, with bcftools.

    Also extracts the per-site DP/GQ that FlagBloodGroupCallQc turns into a QC flag.

    Emits two outputs from one pass over the gVCF: `vcf`, the rbceq2 input, and
    `defining_sites`, the DP/GQ at every allele-defining coordinate. Both derive from a
    blood-group-regions intermediate that retains <NON_REF>, and so retains the DRAGEN
    reference blocks the extract needs.

    Restriction to resources/bg_regions.<genome>.bed is unconditional, so the gVCF .tbi
    must exist. The BED must be a strict superset of every coordinate rbceq2 queries for
    the configured reference build, or blood-group calls are silently wrong.

    We split multiallelics, drop the <NON_REF> symbolic allele (which breaks rbceq2),
    then trim now-unused ALT alleles. A tabix index is written alongside the VCF
    because rbceq2 fetches blood-group regions by coordinate.

    Genotypes are NOT filtered on FORMAT/DP or FORMAT/GQ. rbceq2 reads any blood-group
    site absent from its input as a confident homozygous reference call, so dropping a
    borderline genotype does not produce a no-call — it manufactures a wild-type call at
    a site that defines a blood-group antigen. DRAGEN has already hard-filtered these
    gVCFs: a record that failed keeps its filter name (DRAGENSnpHardQUAL,
    DRAGENIndelHardQUAL or LowDepth on the recalibrated gVCFs), and rbceq2 excludes an
    allele whose defining variant is not PASS. DP and GQ are reported as a per-system QC
    flag instead of silently removing data. Post-hoc records are given the same FILTER
    semantics in the merge, since HaplotypeCaller leaves the column `.`, which rbceq2
    would read as a failure.

    For an exome sequencing group this stage also merges in the post-hoc calls from
    PosthocGenotypeOffTargetSites, which fill the defining sites the capture-target BED
    stopped DRAGEN emitting at. Only sites outside the capture design are filled, and which
    sites those are is the committed subtraction for the configured design (`off_design`); see
    `_merge_posthoc_commands` for the rest of the fill rule. A genome sequencing group has no
    post-hoc input and never reads the design key, so nothing is merged for it and the merge is
    a plain rename.

    Its command is not otherwise unchanged, though. Every run, genome included, now declares
    INFO/POSTHOC on the intermediate and extracts a trailing POSTHOC column, because
    `bcftools query` aborts on a tag the header does not declare rather than rendering `.`.
    That is what the release version bump records.

    The `norm -m -any` split must stay ahead of the <NON_REF> exclusion. In a gVCF a
    variant record carries <NON_REF> as a trailing ALT (A -> G,<NON_REF>) and
    ALT="<NON_REF>" matches if any ALT matches, so filtering before the split would
    delete every variant in the file. After the split only the symbolic-only record
    matches and the real variant survives.

    Genotypes are never rewritten. DRAGEN calls single-copy chrX at its real ploidy in a male
    sample (GT="1"/"0" outside PAR), rbceq2 scores a one-token GT as one copy, and a rewrite to
    "1|1" would make a hemizygous XK null render as `XK*N.16/XK*N.16`, indistinguishable in the
    genotype TSV from a female homozygote.

    **The file rbceq2 reads must carry one ploidy per region.** rbceq2 derives a single
    chromosome-copy count per blood group and refuses any record claiming more: the system
    reports `Undetermined/Undetermined` and empty phenotypes. A record claiming *fewer* passes
    silently and resolves the contradiction the wrong way, flipping the phenotype.

    That invariant belongs to `merged.vcf.gz`, not to this conversion. A genome's merged VCF is
    DRAGEN's records alone and satisfies it by construction. An exome's also holds
    HaplotypeCaller's, which run at ploidy 2, so `_merge_posthoc_commands` keeps the fill out of
    single-copy chrX. See the haploid-encoding section of the README.
    """

    def expected_outputs(
        self, sequencing_group: cpg_flow.targets.SequencingGroup
    ) -> stage_support.ExpectedOutputs | None:
        if not sequencing_group.gvcf:
            return None
        prefix = stage_support.get_sg_output_prefix(sequencing_group, stage_name=self.name, category='tmp')
        return {
            'vcf': prefix / f'{sequencing_group.id}.converted.vcf.gz',
            'index': prefix / f'{sequencing_group.id}.converted.vcf.gz.tbi',
            'defining_sites': prefix / f'{sequencing_group.id}.defining_sites.tsv',
        }

    def queue_jobs(
        self,
        sequencing_group: cpg_flow.targets.SequencingGroup,
        inputs: cpg_flow.stage.StageInput,
    ) -> cpg_flow.stage.StageOutput | None:
        outputs = self.expected_outputs(sequencing_group)
        if outputs is None:
            return None
        cfg = stage_support.config_section(self)
        cpu = cpg_utils.config.config_retrieve(['workflow', cfg, 'cpu'], 4)
        genome = cpg_utils.config.genome_build()

        b = cpg_utils.hail_batch.get_batch()
        j = b.new_bash_job(
            f'FilterAndConvertGvcfsForRbceq2/{sequencing_group.id}',
            self.get_job_attrs(sequencing_group) | {'tool': 'bcftools'},
        )
        j = stage_support.configure_job(
            j,
            self,
            cpu=cpu,
            memory='highmem',
            storage='40Gi',
            image=cpg_utils.config.image_path('bcftools', '1.24-1'),
        )

        # -R index-jumps to the blood-group regions, so the gVCF .tbi is required. Both
        # BEDs come from one gen_bg_resources.py pass over the db.tsv the pinned rbceq2
        # image uses, so the converted regions and the QC sites cannot drift apart.
        gvcf = b.read_input_group(
            **{'g.vcf.gz': str(sequencing_group.gvcf), 'g.vcf.gz.tbi': f'{sequencing_group.gvcf}.tbi'},
        )['g.vcf.gz']
        regions_bed = b.read_input(stage_support.blood_group_resource(f'bg_regions.{genome}.bed'))
        sites_bed = b.read_input(stage_support.blood_group_resource(f'bg_defining_sites.{genome}.bed'))
        j.declare_resource_group(out={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
        # declare_resource_group returns the job, so the group comes back through
        # Job.__getattr__, which is typed as a plain Resource. It is the group just declared.
        out = typing.cast('hailtop.batch.resource.ResourceGroup', j.out)

        # The same predicate the post-hoc stage gates its own outputs on, so this cannot ask
        # cpg_flow for an input that stage produced nothing for.
        sample_check = ''
        merge_posthoc = ''
        if posthoc_genotype.applies_to(sequencing_group):
            # Both keys come from the producer, rather than the index being spelled here as
            # the gVCF path plus '.tbi'. The merge reads the index, for the -T targeted read,
            # so a producer that renamed it would fail inside a running job on a missing file.
            # Reading the key it declared makes that a graph-build error instead.
            posthoc_paths = inputs.as_dict(sequencing_group, posthoc_genotype.PosthocGenotypeOffTargetSites)
            posthoc_gvcf = b.read_input_group(
                **{
                    'g.vcf.gz': str(posthoc_paths['gvcf']),
                    'g.vcf.gz.tbi': str(posthoc_paths['index']),
                },
            )['g.vcf.gz']
            # The design BED itself is never localised. Only the committed 167-line
            # subtraction of it is, which saves moving 5.5Mb of vendor intervals to every
            # sequencing group to re-derive one design-constant answer.
            design_key = stage_support.exome_design_bed()
            off_design_bed = b.read_input(off_design.resource_path(design_key))
            sample_check = _sample_check_commands(str(posthoc_gvcf), str(gvcf))
            merge_posthoc = _merge_posthoc_commands(
                str(posthoc_gvcf),
                str(sites_bed),
                str(off_design_bed),
                design_key,
                cpu,
                genome,
            )
        else:
            merge_posthoc = '        mv dragen.vcf.gz merged.vcf.gz'

        # The extract reads the intermediate, not the converted VCF: dropping <NON_REF>
        # deletes every reference block, and a reference block is exactly what covers a
        # defining site the caller saw no variant at.
        #
        # --threads only ever parallelises BGZF (de)compression, so it belongs on the steps
        # that do some: norm decompresses the bgzipped gVCF (the shared pool is attached to
        # input readers as well as the output, synced_bcf_reader.c bcf_sr_add_hreader), the
        # final view deflates the -Oz converted VCF, and index reads that back.
        #
        # --targets-overlap 2 is what makes the extract see a reference block that starts
        # before a defining site and spans it. Streamed targets default to `pos`, which
        # requires POS inside the region and so would silently miss every spanning block,
        # the case the QC flag exists to catch. Streaming is fine here: the intermediate is
        # only the blood-group regions, so there is nothing to gain from indexing it.
        #
        # GT is extracted for the records that reach a defining site from an earlier POS: a
        # deletion spanning the site removes the base its antigen is defined on, and whether
        # it does so on one haplotype or both is the difference between rbceq2's reference
        # call being half-supported and being unsupported.
        #
        # INFO/POSTHOC is extracted unconditionally, on genome runs as well as exome ones,
        # where it renders `.` for every record. One extract format everywhere means the
        # parser needs no per-sequencing-type branch to know how many columns to expect.
        #
        # The wrong-build guard reads the DRAGEN intermediate, before the merge, and not the
        # extract: see _primary_records_guard for why the extract cannot answer it on an exome.
        j.command(
            f"""
            set -euxo pipefail
{sample_check}
            echo '{_POSTHOC_HEADER_LINE}' > posthoc_hdr.txt
            bcftools norm -m -any --threads {cpu} -R {regions_bed} -Ou {gvcf} \\
                | bcftools annotate -h posthoc_hdr.txt --threads {cpu} -Oz -o dragen.vcf.gz -
{_primary_records_guard(str(sites_bed), genome)}
{merge_posthoc}
            bcftools query -T {sites_bed} --targets-overlap 2 \\
                -f '{_EXTRACT_FORMAT}' \\
                merged.vcf.gz > {j.sites}
{_convert_commands(str(out['vcf.gz']), cpu)}
            """,
        )
        # write_output base drops the suffix; the resource group re-adds .vcf.gz / .vcf.gz.tbi.
        b.write_output(out, str(outputs['vcf']).removesuffix('.vcf.gz'))
        b.write_output(j.sites, str(outputs['defining_sites']))
        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])
