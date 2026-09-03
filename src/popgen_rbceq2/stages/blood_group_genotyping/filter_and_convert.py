"""Turns a DRAGEN gVCF into the VCF rbceq2 reads, and extracts the sites the QC stage judges."""

import typing

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import posthoc_genotype

# Which sites of a BED fall inside no span of another BED, as a BED.
#
# `awk -v spans=spans.bed <this> spans.bed sites.bed` prints the sites no span contains. A span
# is `[$2, $3)`, 0-based half-open, so the program serves two inputs without change: the covered
# spans of every DRAGEN record overlapping a defining site (bcftools' %END is POS + rlen - 1,
# which is INFO/END for a reference block and POS + len(REF) - 1 for a variant, so `%POS0\t%END`
# is that span and this containment test is the one `rbceq2_call_qc_job.GvcfRecord.covers`
# applies), and the capture design's target intervals. The two definitions of "covered" have
# to agree — a site this calls uncovered is a site the QC would flag NOCOV — and
# tests/test_posthoc_merge.py holds them together.
#
# The spans file is named by `-v spans=`, and must be, rather than being detected with the
# usual `NR == FNR`. That idiom reads "still in the first file" only while the first file has
# produced records, so an *empty* spans file makes awk take the sites file for the span list
# and print nothing at all. Every caller below would then subtract the sites from themselves:
# an empty covered.bed, which is what a gVCF with no record at any defining site produces,
# would fill no holes and fail later blaming contig naming. `ARGIND` would say the same thing
# more directly but is a gawk extension, and this runs under mawk.
#
# Extra BED columns and `track`/`browser` header lines are harmless: only $1-$3 are read, and a
# header line becomes an empty span on a contig no site is on.
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


# The INFO flag marking a post-hoc record whose span reaches a defining site the primary caller
# already has a record for. Set from covered.bed inside the merge and stripped again before the
# supplement is concatenated, so it reaches neither the merged file nor the extract.
_COVERED_HEADER_LINE = (
    '##INFO=<ID=COVERED,Number=0,Type=Flag,Description='
    '"Overlaps a defining site the primary caller already has a record for">'
)


# The `bcftools query` format the QC extract is written in. One field per entry of
# `rbceq2_call_qc_job.EXTRACT_COLUMNS`, in the same order; test_posthoc_merge holds the two to
# the same field count, since a mismatch would be caught only when a job ran and the parser
# rejected the file.
_EXTRACT_FORMAT = r'%CHROM\t%POS\t%REF\t%ALT\t%INFO/END\t[%GT\t%DP\t%GQ\t%MIN_DP]\t%INFO/POSTHOC\n'


# The stage config key naming the capture design an exome cohort was called against, as a key
# into the `[references]` section, e.g.
# `exome_probesets_hg38/agilent_sureselect_clinical_research_exome_v2_covered_by_probes_bed`.
# Required for an exome run; a genome run never reads it.
EXOME_DESIGN_KEY = 'exome_design_bed'


def exome_design_bed(stage: cpg_flow.stage.Stage) -> tuple[str, str]:
    """The reference key and path of the capture design BED an exome run fills holes outside of.

    Read at graph-build time, so a run missing it fails before a job starts rather than on the
    first exome sequencing group's merge.

    Args:
        stage: The conversion stage, whose config section holds the key.

    Returns:
        The `[references]` key as configured, and the path it resolves to.

    Raises:
        cpg_utils.config.ConfigError: The key is not set, or names no reference.
    """
    section = stage_support.config_section(stage)
    try:
        key = cpg_utils.config.config_retrieve(['workflow', section, EXOME_DESIGN_KEY])
    except cpg_utils.config.ConfigError as e:
        raise cpg_utils.config.ConfigError(
            f'An exome run needs workflow.{section}.{EXOME_DESIGN_KEY}: the [references] key of the '
            'capture design BED the gVCFs were called against, e.g. '
            "'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'. Post-hoc calls fill "
            'defining sites only outside that design.'
        ) from e
    return key, cpg_utils.config.reference_path(key)


def _merge_posthoc_commands(
    posthoc_gvcf: str,
    sites_bed: str,
    design_bed: str,
    design_key: str,
    cpu: int,
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
    emits records over exactly the target BED (no padding), so a DRAGEN record at a defining
    site outside the configured design means the wrong file is configured — on the validation
    cohorts a target-regions file in place of the probe-footprint one leaves 309 defining sites
    "outside" the design with DRAGEN records, and would have gated off three quarters of the
    recoveries while the run looked fine. That case fails the job.

    A record is selected for reaching a hole and is selected whole, so one anchored in a hole
    can extend over a neighbouring defining site DRAGEN did call. Blood-group defining sites
    are dense — most have another within 20bp — so this is the ordinary case near a capture
    edge, not a corner one, and what happens next depends on what the record is.

    A post-hoc **variant** that reaches a called base is dropped. Keeping it would hand rbceq2
    two callers' alleles at one base with nothing to choose between them: a DRAGEN SNP at a
    site, and a post-hoc deletion removing it. The QC could not report the conflict either,
    because `resolve_coverage` prefers the primary record and so reads the site as an ordinary
    PASS. The hole the dropped record would have filled goes back to reaching the QC as NOCOV,
    which is the honest answer — DRAGEN wins wherever both speak, and here both spoke.

    A post-hoc **reference block** that reaches a called base is kept whole. It asserts nothing
    rbceq2 ever sees, since the conversion drops every <NON_REF>-only record before rbceq2
    reads the file, and dropping the block instead would throw away the hole it was kept for.
    That does leave two records covering the called site in the extract, which is why
    `resolve_coverage` prefers the record with no INFO/POSTHOC. Splitting blocks on the
    boundary would be the alternative and is not worth it.

    Fails rather than merging if the CRAM and the gVCF name different samples. They are two
    outputs of one DRAGEN run, so a disagreement means this sequencing group's inputs do not
    describe one individual, and merging would splice another person's genotypes into the
    calls at exactly the sites nothing else covers.

    Args:
        posthoc_gvcf: Localised post-hoc gVCF from PosthocGenotypeOffTargetSites.
        sites_bed: The committed defining-sites BED.
        design_bed: Localised capture design BED the gVCF was called against.
        design_key: The `[references]` key `design_bed` came from, for the error message.
        cpu: Threads to give the BGZF steps.

    Returns:
        The shell fragment, for interpolation into the stage's command.
    """
    return f"""
        # --targets-overlap 1 here, not 2: this needs every record whose *span* reaches a
        # defining site, which is what %END reports and what the QC counts as covering. Mode 2
        # asks whether the *variant* overlaps, and drops a deletion anchored on the site
        # itself, which would make a covered site look like a hole and let a post-hoc record
        # displace a DRAGEN call.
        bcftools index -t --threads {cpu} dragen.vcf.gz
        bcftools query -T {sites_bed} --targets-overlap 1 \\
            -f '%CHROM\\t%POS0\\t%END\\n' dragen.vcf.gz > covered.bed

        # The defining sites the capture design did not target, then those of them the DRAGEN
        # gVCF has no record at. Only that second set is filled.
        awk -v spans={design_bed} '{_SITES_OUTSIDE_SPANS_AWK}' \\
            {design_bed} {sites_bed} > off_design_sites.bed
        awk -v spans=covered.bed '{_SITES_OUTSIDE_SPANS_AWK}' \\
            covered.bed off_design_sites.bed > uncovered.bed

        # A DRAGEN record at a site outside the design means the configured BED is not the one
        # the gVCF was called against. The likely case is a target-regions file where the gVCF
        # used the probe footprint, which would quietly gate off most of the recoveries.
        sort off_design_sites.bed > off_design_sites.sorted.bed
        sort uncovered.bed > uncovered.sorted.bed
        comm -23 off_design_sites.sorted.bed uncovered.sorted.bed > off_design_called.bed
        if [ -s off_design_called.bed ]; then
            n_called=$(wc -l < off_design_called.bed | tr -d ' ')
            echo "ERROR: DRAGEN has records at $n_called defining site(s) outside the capture design." >&2
            echo "{EXOME_DESIGN_KEY} = {design_key} is not the BED this gVCF was called against." >&2
            echo "For an Agilent design that usually means Regions configured where the gVCF" >&2
            echo "used Covered. First sites:" >&2
            head -n 10 off_design_called.bed >&2
            exit 1
        fi

        awk -v spans=covered.bed '{_SITES_OUTSIDE_SPANS_AWK}' covered.bed {sites_bed} > holes.bed
        n_holes=$(wc -l < holes.bed | tr -d ' ')
        n_fill=$(wc -l < uncovered.bed | tr -d ' ')
        echo "post-hoc: $n_holes defining site(s) with no DRAGEN record; $n_fill outside the" >&2
        echo "capture design and filled, $((n_holes - n_fill)) inside it and left for the QC to flag NOCOV" >&2

        # An empty -T file is a hard error in bcftools ("Failed to read the targets"), so the
        # no-holes case has to branch rather than fall through the same pipeline.
        if [ -s uncovered.bed ]; then
            # The defining sites DRAGEN did cover: every site, less the holes. That is what a
            # post-hoc record may not reach, and it is deliberately narrower than covered.bed.
            # Marking against the DRAGEN records' whole spans would also drop a post-hoc
            # variant that merely overlaps the tail of a long reference block, reaching no
            # defining site DRAGEN called and so contradicting nothing rbceq2 reads.
            sort {sites_bed} > sites.sorted.bed
            sort holes.bed > holes.sorted.bed
            comm -23 sites.sorted.bed holes.sorted.bed \
                | sort -k1,1 -k2,2n | bgzip -c --threads {cpu} > called_sites.bed.gz
            tabix -p bed called_sites.bed.gz
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
            # site DRAGEN called. `annotate -m` matches on the record's span, not on POS, so
            # a deletion anchored on a hole and reaching a called site one base away is
            # marked, which is the case this whole check exists for.
            bcftools view -T uncovered.bed --targets-overlap 2 -e 'FORMAT/DP=0' -Ou {posthoc_gvcf} \\
                | bcftools annotate -x '^INFO/END,^FORMAT/GT,FORMAT/DP,FORMAT/GQ,FORMAT/MIN_DP' -Ou - \\
                | bcftools norm -m -any --threads {cpu} -Ou - \\
                | bcftools annotate -a called_sites.bed.gz -h covered_hdr.txt -c CHROM,FROM,TO -m COVERED \\
                    -Ob -o posthoc_marked.bcf -

            # A marked record with no INFO/END is dropped: it is a variant, or the <NON_REF>
            # twin `norm` split off one, and either way it asserts something about a defining
            # site DRAGEN already called, in the file rbceq2 reads. INFO/END is what tells a
            # real reference block from that twin, and a block is kept — see the docstring.
            n_trespass=$(bcftools view -H -i 'INFO/COVERED=1 && INFO/END="."' posthoc_marked.bcf | wc -l | tr -d ' ')
            bcftools view -e 'INFO/COVERED=1 && INFO/END="."' -Ou posthoc_marked.bcf \\
                | bcftools annotate -x INFO/COVERED -h posthoc_hdr.txt -Ov - \\
                | awk -v OFS='\\t' -v tag='POSTHOC={constants.POSTHOC_CALLER}' '{_TAG_POSTHOC_AWK}' \\
                | bgzip -c --threads {cpu} > posthoc_tagged.vcf.gz

            # What was actually merged, which is not the hole count reported above: most holes
            # yield no record at all (DP=0), and a kept variant contributes two records here,
            # itself and its <NON_REF> twin. Reporting only the holes hid the all-zero-depth
            # case entirely — a header-only supplement that `concat` merges silently.
            n_kept=$(bcftools view -H posthoc_tagged.vcf.gz | wc -l | tr -d ' ')
            echo "post-hoc: $n_kept record(s) kept over $n_fill hole(s); $n_trespass dropped for" >&2
            echo "reaching a defining site DRAGEN called, the rest had no reads there (DP=0)" >&2

            # The two callers must agree on whose sample this is, and disagreeing is fatal.
            #
            # The post-hoc caller takes its sample name from the CRAM's read group and DRAGEN
            # named the gVCF from the same run, so a mismatch means the CRAM and the gVCF this
            # sequencing group resolves to do not describe one individual. Merging them would
            # splice another person's genotypes into this one's calls at exactly the sites
            # nothing else covers, and the result would look like an ordinary recovery.
            #
            # Relabelling the supplement instead would satisfy `concat`, which requires
            # identical sample sets, and would bury that. Some mackenzie DRAGEN 3.7.8 test
            # CRAMs do trip this benignly, carrying a retired sequencing-group ID for the same
            # individual from an upstream test-set reheadering bug that is fixed for newer
            # additions. That is a reason to fix those inputs, not to weaken the check for
            # every cohort: this is the only place the pipeline compares the two files it was
            # handed, and a real swap and a stale header are indistinguishable from here.
            posthoc_sample=$(bcftools query -l posthoc_tagged.vcf.gz)
            dragen_sample=$(bcftools query -l dragen.vcf.gz)
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
    gVCFs (records are FILTER=PASS); DP and GQ are reported as a per-system QC flag
    instead of silently removing data.

    For an exome sequencing group this stage also merges in the post-hoc calls from
    PosthocGenotypeOffTargetSites, which fill the defining sites the capture-target BED
    stopped DRAGEN emitting at. The design BED is named by `exome_design_bed` in this stage's
    config section and only sites outside it are filled; see `_merge_posthoc_commands` for
    the fill rule. A genome sequencing group has no post-hoc input and never reads the key, so
    nothing is merged for it and the merge is a plain rename.

    Its command is not otherwise unchanged, though. Every run, genome included, now declares
    INFO/POSTHOC on the intermediate and extracts a trailing POSTHOC column, because
    `bcftools query` aborts on a tag the header does not declare rather than rendering `.`.
    That is what the release version bump records.

    The `norm -m -any` split must stay ahead of the <NON_REF> exclusion. In a gVCF a
    variant record carries <NON_REF> as a trailing ALT (A -> G,<NON_REF>) and
    ALT="<NON_REF>" matches if any ALT matches, so filtering before the split would
    delete every variant in the file. After the split only the symbolic-only record
    matches and the real variant survives.

    `bcftools +fixploidy` normalises DRAGEN's true haploid calls (GT="1"/"0") in
    non-PAR chrX/Y for male samples into pseudo-diploid (GT="1|1" or "0|0"). rbceq2
    assumes diploid GTs everywhere and crashes on haploid calls at the XK/GATA1/ATP11C
    blood-group loci (get_ref asserts len(GT) == 3); fixploidy only touches haploid
    genotypes, leaving already-diploid calls and their phasing untouched.
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
            design_key, design_path = exome_design_bed(self)
            design_bed = b.read_input(design_path)
            merge_posthoc = _merge_posthoc_commands(str(posthoc_gvcf), str(sites_bed), str(design_bed), design_key, cpu)
        else:
            merge_posthoc = '        mv dragen.vcf.gz merged.vcf.gz'

        # The extract reads the intermediate, not the converted VCF: dropping <NON_REF>
        # deletes every reference block, and a reference block is exactly what covers a
        # defining site the caller saw no variant at.
        #
        # --threads only ever parallelises BGZF (de)compression, so it belongs on the steps
        # that do some: norm decompresses the bgzipped gVCF (the shared pool is attached to
        # input readers as well as the output, synced_bcf_reader.c bcf_sr_add_hreader),
        # +fixploidy deflates the -Oz output, and index reads that back. The middle view has
        # an uncompressed BCF stream on both sides, so a thread count there does nothing.
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
        j.command(
            f"""
            set -euxo pipefail
            echo '{_POSTHOC_HEADER_LINE}' > posthoc_hdr.txt
            bcftools norm -m -any --threads {cpu} -R {regions_bed} -Ou {gvcf} \\
                | bcftools annotate -h posthoc_hdr.txt --threads {cpu} -Oz -o dragen.vcf.gz -
{merge_posthoc}
            bcftools query -T {sites_bed} --targets-overlap 2 \\
                -f '{_EXTRACT_FORMAT}' \\
                merged.vcf.gz > {j.sites}
            if [ ! -s {j.sites} ]; then
                echo "ERROR: no gVCF record overlaps any blood-group defining site." >&2
                echo "Check the gVCF contig naming, and that references.genome_build" >&2
                echo "({genome}) matches the build the gVCF was called against." >&2
                exit 1
            fi
            bcftools view \\
                    -e 'ALT="<NON_REF>"' \\
                    --trim-alt-alleles -Ou merged.vcf.gz \\
                | bcftools +fixploidy --threads {cpu} -Oz -o {out['vcf.gz']} -
            bcftools index -t --threads {cpu} {out['vcf.gz']}
            """,
        )
        # write_output base drops the suffix; the resource group re-adds .vcf.gz / .vcf.gz.tbi.
        b.write_output(out, str(outputs['vcf']).removesuffix('.vcf.gz'))
        b.write_output(j.sites, str(outputs['defining_sites']))
        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])
