"""Which blood-group defining sites an exome cohort's capture design never targeted."""

import re
import typing

import cpg_flow.inputs
import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import posthoc_genotype

# The stage config key naming the capture design an exome cohort was called against, as a key
# into the `[references]` section, e.g.
# `exome_probesets_hg38/agilent_sureselect_clinical_research_exome_v2_covered_by_probes_bed`.
# Required for an exome run; a genome run never reads it.
EXOME_DESIGN_KEY = 'exome_design_bed'


def _design_segment(design_key: str) -> str:
    """The output-path segment identifying which design was subtracted.

    The result depends on the configured design and on nothing else about the run, and the
    release segment above it cannot see a config change. Without this, repointing
    EXOME_DESIGN_KEY at a different design would find the previous run's BED already written
    and reuse it, and every sample would then be gated on the wrong design without a word in
    any log.

    Args:
        design_key: The `[references]` key the design came from.

    Returns:
        The key with anything outside `[A-Za-z0-9._-]` replaced, so it is one path segment.
    """
    return re.sub(r'[^A-Za-z0-9._-]', '_', design_key)


class SelectOffDesignDefiningSites(cpg_flow.stage.MultiCohortStage):
    """Subtract the capture design's intervals from the defining sites, once for the run.

    The result is the set of defining coordinates the capture never targeted, which is what
    bounds the post-hoc fill: `FilterAndConvertGvcfsForRbceq2` may fill a hole only at a site in
    this set. See `_merge_posthoc_commands` for the rest of the rule, and `docs/PRODUCT.md` for
    why the design bounds the fill at all.

    Run once per workflow run, not per sequencing group, because the answer depends only on the
    configured design and the committed defining sites, both fixed for the whole run.

    A MultiCohortStage because the output is one file per run, and the run's MultiCohort is the
    one target a per-sequencing-group consumer can always name to read it back (see
    `off_design_bed`). Why not a CohortStage is resolved question 5 in
    `docs/rbceq2_posthoc_exome/SPEC.md`; why bedtools rather than awk, and why the per-sample
    subtraction in the conversion job did not follow it, is the bedtools entry in that spec's
    section 10.
    """

    def expected_outputs(self, multicohort: cpg_flow.targets.MultiCohort) -> stage_support.ExpectedOutputs | None:
        # Nothing to bound if no sequencing group in the run can be filled. Skipped rather than
        # failed, matching how the post-hoc stage treats a genome, and it is the same predicate
        # so the two cannot disagree about whether a run needs a design.
        if not any(posthoc_genotype.applies_to(sg) for sg in multicohort.get_sequencing_groups(only_active=True)):
            return None
        design_key, _ = exome_design_bed()
        genome = cpg_utils.config.genome_build()
        prefix = stage_support.get_multicohort_output_prefix(multicohort, stage_name=self.name, category='tmp')
        return {'bed': prefix / _design_segment(design_key) / f'off_design_defining_sites.{genome}.bed'}

    def queue_jobs(
        self,
        multicohort: cpg_flow.targets.MultiCohort,
        inputs: cpg_flow.stage.StageInput,  # noqa: ARG002
    ) -> cpg_flow.stage.StageOutput | None:
        outputs = self.expected_outputs(multicohort)
        if outputs is None:
            return None
        cfg = stage_support.config_section(self)
        cpu = cpg_utils.config.config_retrieve(['workflow', cfg, 'cpu'], 1)
        genome = cpg_utils.config.genome_build()
        design_key, design_path = exome_design_bed()

        b = cpg_utils.hail_batch.get_batch()
        j = b.new_bash_job(self.name, self.get_job_attrs(multicohort) | {'tool': 'bedtools'})
        j = stage_support.configure_job(
            j,
            self,
            cpu=cpu,
            memory='lowmem',
            storage='10Gi',
            image=cpg_utils.config.image_path('bedtools', constants.BEDTOOLS_IMAGE_TAG),
        )

        sites_bed = b.read_input(stage_support.blood_group_resource(f'bg_defining_sites.{genome}.bed'))
        design_bed = b.read_input(design_path)

        # An empty result is legitimate — a design that targets every defining site leaves
        # nothing to fill, and every hole then reaches the QC as NOCOV. An empty *design* is
        # not: it would report every defining site as off-design, which the per-sample check on
        # DRAGEN's records would then fail on, one sample at a time, naming the design rather
        # than the fact that it was empty. Caught here instead, once, before any of that.
        j.command(
            f"""
            set -euxo pipefail
            if [ ! -s {design_bed} ]; then
                echo "ERROR: the capture design BED is empty." >&2
                echo "{DESIGN_CONFIG_PATH} = {design_key}" >&2
                echo "resolved to {design_path}" >&2
                exit 1
            fi
            bedtools intersect -v -a {sites_bed} -b {design_bed} > {j.bed}
            n_off=$(wc -l < {j.bed} | tr -d ' ')
            n_sites=$(wc -l < {sites_bed} | tr -d ' ')
            echo "off-design: $n_off of $n_sites defining site(s) are outside {design_key}," >&2
            echo "and so are the only sites a post-hoc record may fill" >&2
            """,
        )
        b.write_output(j.bed, str(outputs['bed']))
        return self.make_outputs(multicohort, data=outputs, jobs=[j])


# The config section EXOME_DESIGN_KEY is read from, and its fully-qualified path for the error
# messages that have to name it. Both come from the class above rather than a repeated string
# literal, so renaming the stage moves its section, its messages and its docs together — and
# they sit below the class for that reason, not above it with the key they complete.
DESIGN_CONFIG_SECTION = stage_support.camel_to_snake(SelectOffDesignDefiningSites.__name__)
DESIGN_CONFIG_PATH = f'workflow.{DESIGN_CONFIG_SECTION}.{EXOME_DESIGN_KEY}'


def exome_design_bed() -> tuple[str, str]:
    """The reference key and path of the capture design BED an exome run fills holes outside of.

    Read at graph-build time, so a run missing it fails before a job starts rather than on the
    first exome sequencing group's merge.

    Returns:
        The `[references]` key as configured, and the path it resolves to.

    Raises:
        cpg_utils.config.ConfigError: The key is not set, or names no reference.
    """
    try:
        key = cpg_utils.config.config_retrieve(['workflow', DESIGN_CONFIG_SECTION, EXOME_DESIGN_KEY])
    except cpg_utils.config.ConfigError as e:
        raise cpg_utils.config.ConfigError(
            f'An exome run needs {DESIGN_CONFIG_PATH}: the [references] key of the capture design '
            'BED the gVCFs were called against, e.g. '
            "'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'. Post-hoc calls "
            'fill defining sites only outside that design.'
        ) from e
    return key, cpg_utils.config.reference_path(key)


def off_design_bed(
    inputs: cpg_flow.stage.StageInput,
) -> hailtop.batch.resource.Resource:
    """Localise this run's off-design defining sites, for a consumer's job.

    Args:
        inputs: The consuming stage's StageInput. It must list SelectOffDesignDefiningSites in
            its `requires`.

    Returns:
        The localised BED.
    """
    # A MultiCohortStage's output is keyed by the multicohort's target_id, so reading it needs
    # that target rather than the consumer's own. cpg_flow hands every required stage's output
    # to every consumer regardless of target type, so the lookup is the only difference.
    path = inputs.as_path(
        cpg_flow.inputs.get_multicohort(),
        SelectOffDesignDefiningSites,
        key='bed',
    )
    return typing.cast('hailtop.batch.resource.Resource', cpg_utils.hail_batch.get_batch().read_input(str(path)))
