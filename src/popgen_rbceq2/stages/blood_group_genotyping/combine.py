"""Concatenates a cohort's per-sequencing-group TSVs into combined cohort TSVs."""

import json

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils
import cpg_utils.hail_batch
import hailtop.batch.job

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import genotype
from popgen_rbceq2.stages.blood_group_qc import call_qc


class CombineRbceq2OutputsPerCohort(cpg_flow.stage.CohortStage):
    """Concatenate a cohort's per-SG rbceq2 TSVs and QC TSVs into combined cohort TSVs.

    The per-SG TSV paths are resolved from the workflow dependency graph, restricted to this
    cohort's sequencing groups, and written to a JSON manifest the job reads. A manifest
    rather than repeated flags because the command line grows with the cohort: at four paths
    per sequencing group, argv would exceed ARG_MAX around ~3,000 of them.

    The QC TSV is a fourth output rather than a member of constants.RBCEQ2_TSV_KEYS: it comes
    from FlagBloodGroupCallQc, not from rbceq2, and the keys drive rbceq2's own resource group.
    """

    def expected_outputs(self, cohort: cpg_flow.targets.Cohort) -> stage_support.ExpectedOutputs:
        prefix = stage_support.get_output_prefix(cohort, self.name)
        return {key: prefix / f'combined.{cohort.id}.{key}.tsv' for key in (*constants.RBCEQ2_TSV_KEYS, 'qc')}

    def queue_jobs(
        self,
        cohort: cpg_flow.targets.Cohort,
        inputs: cpg_flow.stage.StageInput,
    ) -> cpg_flow.stage.StageOutput | None:
        outputs = self.expected_outputs(cohort)

        # Per-SG TSVs from the upstream SequencingGroupStages, restricted to this cohort.
        # Both stages skip a sequencing group with no gVCF — the same set excluded here — so
        # the two input lists cover the same samples. Every row carries its own sample id in
        # column 1, so the combined TSVs join per sample without relying on row order.
        per_sg: dict[str, dict[str, cpg_utils.Path]] = inputs.as_dict_by_target(
            genotype.GenotypeBloodGroupsWithRbceq2,
        )
        per_sg_qc: dict[str, dict[str, cpg_utils.Path]] = inputs.as_dict_by_target(call_qc.FlagBloodGroupCallQc)
        cohort_sg_ids: set[str] = {sg.id for sg in cohort.get_sequencing_groups(only_active=True) if sg.gvcf}
        grouped: dict[str, list[str]] = {
            key: [str(tsvs[key]) for sg_id, tsvs in per_sg.items() if sg_id in cohort_sg_ids]
            for key in constants.RBCEQ2_TSV_KEYS
        }
        grouped['qc'] = [str(tsvs['qc']) for sg_id, tsvs in per_sg_qc.items() if sg_id in cohort_sg_ids]

        # With the paths in a manifest, build_python_command's empty-list guard no longer
        # covers the zero-eligible-SGs case, so it is checked here, where the cause can be
        # named. The job re-checks each list on its own side.
        if not cohort_sg_ids:
            raise ValueError(f'Cohort {cohort.id} has no sequencing groups with a gVCF; there is nothing to combine')

        manifest = stage_support.get_output_prefix(cohort, self.name, category='tmp') / f'{cohort.id}.manifest.json'
        manifest.write_text(json.dumps(grouped, indent=2))

        b: cpg_utils.hail_batch.Batch = cpg_utils.hail_batch.get_batch()
        j: hailtop.batch.job.BashJob = b.new_bash_job(self.name, self.get_job_attrs(cohort) | {'tool': 'rbceq2'})
        j = stage_support.configure_job(j, self, cpu=2, memory='standard', storage='10Gi')

        args: dict[str, stage_support.JobArg] = {
            'output-geno': str(outputs['geno']),
            'output-pheno-numeric': str(outputs['pheno_numeric']),
            'output-pheno-alphanumeric': str(outputs['pheno_alphanumeric']),
            'output-qc': str(outputs['qc']),
            'manifest': str(manifest),
        }
        j.command(stage_support.build_python_command('rbceq2_gather_job.py', args))

        return self.make_outputs(cohort, data=outputs, jobs=[j])
