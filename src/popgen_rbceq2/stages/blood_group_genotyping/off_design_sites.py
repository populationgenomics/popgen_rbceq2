"""Which blood-group defining sites an exome cohort's capture design never targeted."""

import typing

import cpg_flow.inputs
import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import posthoc_genotype


class SelectOffDesignDefiningSites(cpg_flow.stage.MultiCohortStage):
    """Subtract the capture design's intervals from the defining sites, once for the run.

    The result is the set of defining coordinates the capture never targeted, which is what
    bounds the post-hoc fill: `FilterAndConvertGvcfsForRbceq2` may fill a hole only at a site in
    this set. See `_merge_posthoc_commands` for the rest of the rule, and `docs/PRODUCT.md` for
    why the design bounds the fill at all.

    Run once per workflow run, not per sequencing group, because the answer depends only on the
    configured design and the committed defining sites, both fixed for the whole run. The design
    is a segment of every exome output path, this one included (see
    `stage_support._release_tree`), so repointing it starts a fresh tree rather than reusing the
    previous design's BED or anything built from it.

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
        genome = cpg_utils.config.genome_build()
        prefix = stage_support.get_multicohort_output_prefix(multicohort, stage_name=self.name, category='tmp')
        return {'bed': prefix / f'off_design_defining_sites.{genome}.bed'}

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
        design_key, design_path = stage_support.exome_design_bed()

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
        j.command(
            'set -euxo pipefail\n'
            + _subtraction_commands(str(sites_bed), str(design_bed), str(j.bed), design_key, str(design_path)),
        )
        b.write_output(j.bed, str(outputs['bed']))
        return self.make_outputs(multicohort, data=outputs, jobs=[j])


def _subtraction_commands(sites_bed: str, design_bed: str, out_bed: str, design_key: str, design_path: str) -> str:
    """Shell that writes the defining sites outside the design, or fails on a design that cannot be right.

    A pure function of its arguments so the test suite can run it under real bedtools; the
    subtraction is bedtools' interval semantics and nothing of ours, and a test that spelled
    the command out again would stay green if this one changed.

    An empty *result* is legitimate: a design that targets every defining site leaves nothing
    to fill, and every hole then reaches the QC as NOCOV. Two inputs are not, and both would
    otherwise surface the same way, as every defining site off-design, which the per-sample
    check on DRAGEN's records then fails on one sample at a time, naming the design file rather
    than what was wrong with it. Both are caught here instead, once, before any of that:

    - an empty design BED;
    - a design that shares no interval with any site, which for a real exome design means the
      two BEDs name their contigs differently (`1` against `chr1`). bedtools warns about that
      on stderr and exits 0 with every site off-design.

    Args:
        sites_bed: Localised `bg_defining_sites.<genome>.bed`.
        design_bed: Localised capture design BED.
        out_bed: Where to write the off-design sites.
        design_key: The `[references]` key the design came from, for the messages.
        design_path: The path it resolved to, for the messages.

    Returns:
        Shell lines, for the job command.
    """
    return f"""
            if [ ! -s {design_bed} ]; then
                echo "ERROR: the capture design BED is empty." >&2
                echo "{stage_support.DESIGN_CONFIG_PATH} = {design_key}" >&2
                echo "resolved to {design_path}" >&2
                exit 1
            fi
            bedtools intersect -v -a {sites_bed} -b {design_bed} > {out_bed}
            n_off=$(wc -l < {out_bed} | tr -d ' ')
            n_sites=$(wc -l < {sites_bed} | tr -d ' ')
            if [ "$n_off" -eq "$n_sites" ]; then
                echo "ERROR: every defining site is outside the capture design, which no exome design leaves." >&2
                echo "The two BEDs most likely name their contigs differently:" >&2
                echo "  defining sites: $(cut -f1 {sites_bed} | sort -u | tr '\\n' ' ')" >&2
                design_contigs=$(awk '!/^(track|browser|#)/ {{print $1}}' {design_bed} | sort -u | tr '\\n' ' ')
                echo "  design:         $design_contigs" >&2
                echo "{stage_support.DESIGN_CONFIG_PATH} = {design_key}" >&2
                echo "resolved to {design_path}" >&2
                exit 1
            fi
            echo "off-design: $n_off of $n_sites defining site(s) are outside {design_key}," >&2
            echo "and so are the only sites a post-hoc record may fill" >&2
    """


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
