"""Turns a DRAGEN gVCF into the VCF rbceq2 reads, and extracts the sites the QC stage judges."""

import typing

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import posthoc_genotype

# Which allele-defining sites the primary gVCF says nothing about, as a BED.
#
# Reads the covered spans of every record overlapping a defining site, then prints the sites
# no span contains. A span is `[POS0, END)`: bcftools' %END is POS + rlen - 1, which is
# INFO/END for a reference block and POS + len(REF) - 1 for a variant, so this containment
# test is the same one `rbceq2_call_qc_job.GvcfRecord.covers` applies. The two have to agree —
# a site this calls uncovered is a site the QC would flag NOCOV.
#
# Held as a plain string rather than inlined: the command below is an f-string, and every
# brace in an awk program would have to be doubled.
_UNCOVERED_SITES_AWK = """
NR == FNR { n = ++c[$1]; lo[$1, n] = $2 + 0; hi[$1, n] = $3 + 0; next }
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


# The `bcftools query` format the QC extract is written in. One field per entry of
# `rbceq2_call_qc_job.EXTRACT_COLUMNS`, in the same order; test_posthoc_merge holds the two to
# the same field count, since a mismatch would be caught only when a job ran and the parser
# rejected the file.
_EXTRACT_FORMAT = r'%CHROM\t%POS\t%REF\t%ALT\t%INFO/END\t[%GT\t%DP\t%GQ\t%MIN_DP]\t%INFO/POSTHOC\n'


def _merge_posthoc_commands(
    posthoc_gvcf: str,
    sites_bed: str,
    cpu: int,
) -> str:
    """Shell to fill the primary gVCF's blind spots from the post-hoc caller's gVCF.

    Reads `dragen.vcf.gz` and writes `merged.vcf.gz`, which is what the extract and the
    conversion then run over.

    The fill is empirical and per sample: a post-hoc record survives only where the DRAGEN
    gVCF has no record covering a defining site, so DRAGEN wins wherever both speak. No
    capture BED is consulted anywhere, which is deliberate — a BED that misdescribes the real
    footprint of a sample's capture could otherwise overwrite a DRAGEN call or miss a hole.

    A kept post-hoc reference block is kept whole, so one that straddles a capture edge can
    also cover a defining site DRAGEN called. That leaves two records covering that site in
    the merged file, which is why `resolve_coverage` prefers the record with no INFO/POSTHOC.
    Splitting blocks on the boundary would be the alternative and is not worth it.

    Fails rather than merging if the CRAM and the gVCF name different samples: they are two
    outputs of one DRAGEN run and disagreeing means they are not a matched pair. The failure
    mode being guarded against is a CRAM registered against the wrong sequencing group, which
    would write another individual's genotypes into this one's call at precisely the sites
    with no other evidence.

    Args:
        posthoc_gvcf: Localised post-hoc gVCF from PosthocGenotypeOffTargetSites.
        sites_bed: The committed defining-sites BED.
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
        awk '{_UNCOVERED_SITES_AWK}' covered.bed {sites_bed} > uncovered.bed

        # An empty -T file is a hard error in bcftools ("Failed to read the targets"), so the
        # no-holes case has to branch rather than fall through the same pipeline.
        if [ -s uncovered.bed ]; then
            echo "post-hoc: filling $(wc -l < uncovered.bed) defining site(s) with no DRAGEN record" >&2
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
            bcftools view -T uncovered.bed --targets-overlap 2 -e 'FORMAT/DP=0' -Ou {posthoc_gvcf} \\
                | bcftools annotate -x '^INFO/END,^FORMAT/GT,FORMAT/DP,FORMAT/GQ,FORMAT/MIN_DP' -Ou - \\
                | bcftools norm -m -any --threads {cpu} -Ou - \\
                | bcftools annotate -h posthoc_hdr.txt -Ov - \\
                | awk -v OFS='\\t' -v tag='POSTHOC={constants.POSTHOC_CALLER}' '{_TAG_POSTHOC_AWK}' \\
                | bgzip -c --threads {cpu} > posthoc_tagged.vcf.gz

            # The two callers must agree on whose sample this is. The post-hoc caller names
            # its sample from the CRAM's read group and DRAGEN named the gVCF from the same
            # run, so on matched inputs these are equal — verified across the mackenzie
            # DRAGEN 3.7.8 outputs, where the embedded name differs from the sequencing-group
            # ID for some samples but the CRAM and the gVCF always agree with each other.
            #
            # Renaming the supplement to match instead of checking would be the easy path and
            # the wrong one: a CRAM holding another individual's reads is exactly what this
            # catches, and the consequence of missing it is that individual's genotypes being
            # written into this one's blood-group call, at the only sites we have no other
            # evidence for. concat also requires identical sample sets, so this is what makes
            # the concat below legal rather than a second thing to keep in step.
            posthoc_sample=$(bcftools query -l posthoc_tagged.vcf.gz)
            dragen_sample=$(bcftools query -l dragen.vcf.gz)
            if [ "$posthoc_sample" != "$dragen_sample" ]; then
                echo "ERROR: the post-hoc calls and the gVCF name different samples." >&2
                echo "  CRAM/post-hoc: $posthoc_sample" >&2
                echo "  primary gVCF:  $dragen_sample" >&2
                echo "Both should come from one DRAGEN run on one individual, so this means" >&2
                echo "the CRAM and the gVCF for this sequencing group are not a matched pair." >&2
                echo "Check what Metamist has registered as its cram before re-running." >&2
                exit 1
            fi
            bcftools index -t --threads {cpu} posthoc_tagged.vcf.gz
            bcftools concat -a --threads {cpu} -Oz -o merged.vcf.gz dragen.vcf.gz posthoc_tagged.vcf.gz
        else
            echo "post-hoc: every defining site has a DRAGEN record; nothing to fill" >&2
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
    stopped DRAGEN emitting at. See `_merge_posthoc_commands` for the fill rule. A genome
    sequencing group has no post-hoc input and its command is unchanged by that stage
    existing.

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
            posthoc_path = inputs.as_path(
                sequencing_group,
                posthoc_genotype.PosthocGenotypeOffTargetSites,
                key='gvcf',
            )
            posthoc_gvcf = b.read_input_group(
                **{'g.vcf.gz': str(posthoc_path), 'g.vcf.gz.tbi': f'{posthoc_path}.tbi'},
            )['g.vcf.gz']
            merge_posthoc = _merge_posthoc_commands(str(posthoc_gvcf), str(sites_bed), cpu)
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
