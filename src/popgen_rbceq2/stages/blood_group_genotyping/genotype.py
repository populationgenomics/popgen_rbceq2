"""Runs rbceq2 over a sequencing group's converted VCF."""

from typing import Literal

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import filter_and_convert


class GenotypeBloodGroupsWithRbceq2(cpg_flow.stage.SequencingGroupStage):
    """Run rbceq2 blood-group genotyping on a sequencing group.

    One bash job consumes the filtered VCF from FilterAndConvertGvcfsForRbceq2 and emits three
    TSVs (geno, pheno_numeric, pheno_alphanumeric).

    The Analysis this registers records the inferred blood type in its meta (see
    analysis_meta.blood_group_calls); it does NOT write the blood type onto the SequencingGroup
    record (sg.meta). If we later want it queryable directly on the SG (like sg.meta['qc']), add
    a SequencingGroup update mutation.
    """

    def expected_outputs(
        self, sequencing_group: cpg_flow.targets.SequencingGroup
    ) -> stage_support.ExpectedOutputs | None:
        if not sequencing_group.gvcf:
            return None
        prefix = stage_support.get_sg_output_prefix(sequencing_group, self.name)
        return {key: prefix / f'{sequencing_group.id}.{key}.tsv' for key in constants.RBCEQ2_TSV_KEYS}

    def queue_jobs(
        self,
        sequencing_group: cpg_flow.targets.SequencingGroup,
        inputs: cpg_flow.stage.StageInput,
    ) -> cpg_flow.stage.StageOutput | None:
        outputs = self.expected_outputs(sequencing_group)
        if outputs is None:
            return None
        prefix = stage_support.get_sg_output_prefix(sequencing_group, self.name)
        cfg = stage_support.config_section(self)
        output_pdfs = cpg_utils.config.config_retrieve(['workflow', cfg, 'output_pdfs'], False)

        b: cpg_utils.hail_batch.Batch = cpg_utils.hail_batch.get_batch()
        vcf_path: str = str(
            inputs.as_path(sequencing_group, filter_and_convert.FilterAndConvertGvcfsForRbceq2, key='vcf'),
        )
        vcf: hailtop.batch.resource.ResourceGroup = b.read_input_group(
            **{'vcf.gz': vcf_path, 'vcf.gz.tbi': f'{vcf_path}.tbi'},
        )

        j = b.new_bash_job(
            f'GenotypeBloodGroupsWithRbceq2/{sequencing_group.id}',
            self.get_job_attrs(sequencing_group) | {'tool': 'rbceq2'},
        )
        j = stage_support.configure_job(
            j,
            self,
            cpu=2,
            memory='standard',
            storage='10Gi',
            image=cpg_utils.config.image_path('rbceq2', constants.RBCEQ2_IMAGE_TAG),
        )
        # Resource-group keys are the dot form (<sg>.geno.tsv on write); the {root}_<key>.tsv
        # templates capture rbceq2's underscore-named local output.
        j.declare_resource_group(out={f'{k}.tsv': f'{{root}}_{k}.tsv' for k in constants.RBCEQ2_TSV_KEYS})
        pdf_flag: Literal['--PDFs', ''] = '--PDFs' if output_pdfs else ''
        cmd = f"""
            set -euxo pipefail
            rbceq2 --vcf {vcf['vcf.gz']} --out {j.out} --reference_genome {cpg_utils.config.genome_build()} {pdf_flag}
        """
        # Hail names resource-group outputs <dest>.<key>, so this writes <sg>.geno.tsv etc.
        b.write_output(j.out, dest=str(prefix / sequencing_group.id))
        if output_pdfs:
            # rbceq2 writes a "<out>_PDFs/" directory; capture it as a tarball.
            cmd += f'\n            tar -czf {j.pdfs} -C "$(dirname {j.out})" "$(basename {j.out})_PDFs"'
            b.write_output(j.pdfs, str(prefix / f'{sequencing_group.id}_PDFs.tar.gz'))
        j.command(cmd)

        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])
