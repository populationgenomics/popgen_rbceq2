"""Flags the blood-group systems whose call rests on a low-quality defining site."""

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.job

from popgen_rbceq2 import stage_support
from popgen_rbceq2.stages.blood_group_genotyping import filter_and_convert, genotype


class FlagBloodGroupCallQc(cpg_flow.stage.SequencingGroupStage):
    """Flag each blood-group system whose call rests on a low-quality defining site.

    Joins the per-site DP/GQ extracted by FilterAndConvertGvcfsForRbceq2 to the committed
    site -> system map and writes a `<sg>.qc.tsv` in the same wide layout as rbceq2's
    genotype and phenotype TSVs, so it joins to the calls 1:1 by column. The QC TSV's
    columns are taken from `<sg>.geno.tsv`, so it covers the systems rbceq2 actually
    emitted rather than every system in the db.

    DP and GQ thresholds are read from this stage's config section, since this is the
    stage that reports them rather than filtering on them. They are recorded in the
    Analysis meta from the same values the job is given, so the Analysis cannot record a
    threshold the run did not use.

    A flag names both the defining site and what the caller reported at it, which are not
    always the same thing: the DP and GQ can come from a reference block covering the
    coordinate, or from a deletion that removed the base the antigen is defined on (`DEL`).
    Systems whose only defining alleles are large structural variants have no assessable
    site and are reported `NA` rather than `PASS`.

    Reads only the small extract, not the gVCF.

    Registers a per-SG Analysis of its own (analysis_type='blood_group_qc', output = the QC
    TSV), separate from the one GenotypeBloodGroupsWithRbceq2 registers, with the geno TSV
    path in its meta to join the two records.
    """

    def expected_outputs(
        self, sequencing_group: cpg_flow.targets.SequencingGroup
    ) -> stage_support.ExpectedOutputs | None:
        if not sequencing_group.gvcf:
            return None
        return {'qc': stage_support.get_sg_output_prefix(sequencing_group, self.name) / f'{sequencing_group.id}.qc.tsv'}

    def queue_jobs(
        self,
        sequencing_group: cpg_flow.targets.SequencingGroup,
        inputs: cpg_flow.stage.StageInput,
    ) -> cpg_flow.stage.StageOutput | None:
        outputs = self.expected_outputs(sequencing_group)
        if outputs is None:
            return None
        cfg = stage_support.config_section(self)
        genome = cpg_utils.config.genome_build()
        min_depth: int = cpg_utils.config.config_retrieve(['workflow', cfg, 'min_depth'], 10)
        min_gq: int = cpg_utils.config.config_retrieve(['workflow', cfg, 'min_gq'], 20)

        b: cpg_utils.hail_batch.Batch = cpg_utils.hail_batch.get_batch()
        j: hailtop.batch.job.BashJob = b.new_bash_job(
            f'FlagBloodGroupCallQc/{sequencing_group.id}',
            self.get_job_attrs(sequencing_group) | {'tool': 'rbceq2'},
        )
        j = stage_support.configure_job(j, self, cpu=1, memory='standard', storage='10Gi')

        geno_tsv: str = inputs.as_str(sequencing_group, genotype.GenotypeBloodGroupsWithRbceq2, key='geno')
        args: dict[str, stage_support.JobArg] = {
            'geno-tsv': geno_tsv,
            'site-systems': stage_support.blood_group_resource(f'bg_site_systems.{genome}.tsv'),
            'defining-sites-extract': inputs.as_str(
                sequencing_group,
                filter_and_convert.FilterAndConvertGvcfsForRbceq2,
                key='defining_sites',
            ),
            'output': str(outputs['qc']),
            'min-depth': str(min_depth),
            'min-gq': str(min_gq),
        }
        j.command(stage_support.build_python_command('rbceq2_call_qc_job.py', args))

        # Static Analysis meta: the calls this QC describes live in a different stage's
        # Analysis, so record the path rather than leaving the two records unlinked; and the
        # thresholds are the values the job was actually given, not a re-read of config.
        # update_analysis_meta only receives the output path and could derive neither.
        return self.make_outputs(
            sequencing_group,
            data=outputs,
            jobs=[j],
            meta={'blood_group_genotypes_path': geno_tsv, 'min_depth': min_depth, 'min_gq': min_gq},
        )
