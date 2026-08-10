"""Turns a DRAGEN gVCF into the VCF rbceq2 reads, and extracts the sites the QC stage judges."""

import typing

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import stage_support


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
        inputs: cpg_flow.stage.StageInput,  # noqa: ARG002
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
        j.command(
            f"""
            set -euxo pipefail
            bcftools norm -m -any --threads {cpu} -R {regions_bed} -Oz -o bg_regions.vcf.gz {gvcf}
            bcftools query -T {sites_bed} --targets-overlap 2 \\
                -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/END\\t[%GT\\t%DP\\t%GQ\\t%MIN_DP]\\n' \\
                bg_regions.vcf.gz > {j.sites}
            if [ ! -s {j.sites} ]; then
                echo "ERROR: no gVCF record overlaps any blood-group defining site." >&2
                echo "Check the gVCF contig naming, and that references.genome_build" >&2
                echo "({genome}) matches the build the gVCF was called against." >&2
                exit 1
            fi
            bcftools view \\
                    -e 'ALT="<NON_REF>"' \\
                    --trim-alt-alleles -Ou bg_regions.vcf.gz \\
                | bcftools +fixploidy --threads {cpu} -Oz -o {out['vcf.gz']} -
            bcftools index -t --threads {cpu} {out['vcf.gz']}
            """,
        )
        # write_output base drops the suffix; the resource group re-adds .vcf.gz / .vcf.gz.tbi.
        b.write_output(out, str(outputs['vcf']).removesuffix('.vcf.gz'))
        b.write_output(j.sites, str(outputs['defining_sites']))
        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])
