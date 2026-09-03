"""The capture design that bounds the fill: how it is demanded, subtracted, and enforced.

Post-hoc calls fill only defining sites the capture design never targeted. Which design that
is cannot be inferred from the gVCF, so it is configured, as a `[references]` key resolved
through `reference_path` like every other reference this pipeline reads.

The subtraction itself runs once per workflow run, in `SelectOffDesignDefiningSites`, because
its answer depends only on the configured design and the committed defining sites. Doing it
per sequencing group meant an awk pass over a vendor BED's ~230k intervals in every conversion
job, for an answer identical across the cohort.

Four things have to hold and none has a type error to catch it. An exome run that does not name
a design must fail while the graph is being built, not on the first sequencing group, because
the whole cohort is wasted either way and only one of those says why. A genome run must never
read the key at all. Repointing the key must not reuse the previous design's output. And the
conversion job must still refuse to run when DRAGEN's own records contradict the named design.

No batch runs here: each job's command is read off the mock the stage builds it on.
"""

from pathlib import Path
from unittest.mock import MagicMock

import cpg_utils.config
import pytest

from popgen_rbceq2.stages import pipeline
from popgen_rbceq2.stages.blood_group_genotyping.off_design_sites import (
    DESIGN_CONFIG_PATH,
    EXOME_DESIGN_KEY,
)
from tests.helpers import set_config
from tests.test_output_namespacing import outputs_of

pytestmark = pytest.mark.fast

TWIST_KEY = 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'
TWIST_BED = 'gs://cpg-common-main/references/exome-probesets/hg38/Twist_VCGS_Exome_Covered_Targets_hg38.bed'
CREV2_KEY = 'exome_probesets_hg38/agilent_crev2_covered_by_probes_bed'
CREV2_BED = 'gs://cpg-common-main/references/exome-probesets/hg38/S07604514_Covered.bed'
DESIGN_SECTION = 'select_off_design_defining_sites'
OFF_DESIGN_BED = 'gs://bucket-tmp/off_design_defining_sites.GRCh38.bed'


def _config(shm_tmp_path: Path, sequencing_type: str, design_key: str | None) -> None:
    """Write a config for one sequencing type, with or without the design key."""
    stage_section: dict[str, str] = {} if design_key is None else {EXOME_DESIGN_KEY: design_key}
    set_config(
        {
            'references': {
                'genome_build': 'GRCh38',
                'exome_probesets_hg38': {
                    'twist_vcgs_custom_exome_covered_targets_bed': TWIST_BED,
                    'agilent_crev2_covered_by_probes_bed': CREV2_BED,
                },
            },
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': sequencing_type,
                'driver_image': 'stub-driver:1.0',
                DESIGN_SECTION: stage_section,
            },
        },
        shm_tmp_path / f'{sequencing_type}-{(design_key or "none").rsplit("/", 1)[-1]}.toml',
    )


def _queue_subtraction(mocker, multicohort, sequencing_groups) -> MagicMock:
    """Queue the run-level subtraction over a run holding these sequencing groups."""
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    multicohort.get_sequencing_groups.return_value = list(sequencing_groups)
    pipeline.SelectOffDesignDefiningSites().queue_jobs(multicohort, MagicMock())
    return batch


def _queue_conversion(mocker, sequencing_group, multicohort) -> MagicMock:
    """Queue the conversion stage's jobs on a mock batch and return that batch."""
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    mocker.patch('cpg_flow.inputs.get_multicohort', return_value=multicohort)
    inputs = MagicMock()
    inputs.as_dict.return_value = {
        'gvcf': 'gs://bucket/SG000001.posthoc.g.vcf.gz',
        'index': 'gs://bucket/SG000001.posthoc.g.vcf.gz.tbi',
    }
    inputs.as_path.return_value = OFF_DESIGN_BED
    pipeline.FilterAndConvertGvcfsForRbceq2().queue_jobs(sequencing_group, inputs)
    return batch


# --- demanding the design ---


def test_an_exome_run_without_a_design_names_the_key_it_wants(
    mocker,
    exome_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # Failing here is the point. The alternative to a required key is a default, and every
    # default is wrong for some cohort: filling holes inside the capture would answer, with a
    # second caller, a question about that sample's DRAGEN run.
    _config(shm_tmp_path, 'exome', design_key=None)

    with pytest.raises(cpg_utils.config.ConfigError) as raised:
        _queue_subtraction(mocker, mock_multicohort, [exome_sequencing_group])

    assert DESIGN_CONFIG_PATH in str(raised.value)


def test_an_exome_run_naming_a_design_that_is_not_configured_fails(
    mocker,
    exome_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # A typo in the key resolves to nothing. reference_path raising is what makes that a
    # graph-build failure rather than a job that localises an empty file and fills nothing.
    _config(shm_tmp_path, 'exome', design_key='exome_probesets_hg38/no_such_design_bed')

    with pytest.raises(cpg_utils.config.ConfigError):
        _queue_subtraction(mocker, mock_multicohort, [exome_sequencing_group])


def test_a_run_with_no_exome_sequencing_group_produces_no_design_subtraction(
    mocker,
    mock_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # Genome-only runs have no capture edge and must not be made to carry an exome-only key.
    # Skipped rather than failed, on the same predicate the post-hoc stage uses, so the two
    # cannot disagree about whether a run needs a design.
    _config(shm_tmp_path, 'genome', design_key=None)
    mock_multicohort.get_sequencing_groups.return_value = [mock_sequencing_group]

    assert pipeline.SelectOffDesignDefiningSites().expected_outputs(mock_multicohort) is None
    assert _queue_subtraction(mocker, mock_multicohort, [mock_sequencing_group]).new_bash_job.call_count == 0


# --- doing the subtraction ---


def test_the_subtraction_runs_bedtools_over_the_configured_design(
    mocker,
    exome_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # bedtools rather than the awk this replaced: `intersect -v` is the whole rule, on a tool
    # whose half-open semantics and track-line handling match our BEDs.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    batch = _queue_subtraction(mocker, mock_multicohort, [exome_sequencing_group])

    assert TWIST_BED in [call.args[0] for call in batch.read_input.call_args_list]
    command = batch.new_bash_job.return_value.command.call_args.args[0]
    assert 'bedtools intersect -v' in command


def test_an_empty_design_bed_fails_the_subtraction(
    mocker,
    exome_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # An empty design would report every defining site as off-design. That is then caught by
    # each sample's DRAGEN-records check, one sample at a time, naming the design rather than
    # the fact that it had no intervals. Caught here instead, once, before any of that.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    command = _queue_subtraction(
        mocker,
        mock_multicohort,
        [exome_sequencing_group],
    ).new_bash_job.return_value.command.call_args.args[0]

    assert 'the capture design BED is empty' in command
    assert 'exit 1' in command


def test_repointing_the_design_key_writes_to_a_different_path(
    mocker,  # noqa: ARG001
    exome_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # The release segment above this output cannot see a config change, so without the design
    # in the path a re-run against a different design would find the previous run's BED already
    # written and reuse it. Every sample would then be gated on the wrong design, silently.
    mock_multicohort.get_sequencing_groups.return_value = [exome_sequencing_group]
    stage = pipeline.SelectOffDesignDefiningSites()

    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)
    twist = outputs_of(stage, mock_multicohort)
    _config(shm_tmp_path, 'exome', design_key=CREV2_KEY)
    crev2 = outputs_of(stage, mock_multicohort)

    assert twist['bed'] != crev2['bed']


# --- enforcing it, per sequencing group ---


def test_an_exome_merge_gates_the_fill_on_the_subtraction_it_was_handed(
    mocker,
    exome_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # The fill is the intersection: holes are narrowed to off-design sites, never the reverse.
    # The conversion job no longer localises the design BED at all, only the subtraction of it.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    batch = _queue_conversion(mocker, exome_sequencing_group, mock_multicohort)

    localised = [call.args[0] for call in batch.read_input.call_args_list]
    assert TWIST_BED not in localised
    assert OFF_DESIGN_BED in localised
    command = batch.new_bash_job.return_value.command.call_args.args[0]
    assert '> uncovered.bed' in command
    # The design subtraction itself is gone from this job: it is an input now, not a step.
    assert 'bedtools' not in command
    assert '> off_design_sites.bed' not in command


def test_an_exome_merge_fails_when_dragen_called_outside_the_configured_design(
    mocker,
    exome_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # DRAGEN emits over exactly the target BED with no padding, so a record at an off-design
    # defining site means the configured file is not the one the gVCF was called against.
    # Measured on the validation cohorts: a target-regions file where the gVCF used the probe
    # footprint leaves 309 such sites and silently gates off three quarters of the recoveries.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    command = _queue_conversion(
        mocker,
        exome_sequencing_group,
        mock_multicohort,
    ).new_bash_job.return_value.command.call_args.args[0]

    assert 'off_design_called.bed' in command
    assert DESIGN_CONFIG_PATH in command
    assert 'exit 1' in command


def test_a_genome_conversion_never_asks_for_a_design(
    mocker,
    mock_sequencing_group,
    mock_multicohort,
    shm_tmp_path: Path,
):
    # A genome gVCF has no capture edge, gets no post-hoc input, and must not read an
    # exome-only key. Its merge is a plain rename.
    _config(shm_tmp_path, 'genome', design_key=None)

    command = _queue_conversion(
        mocker,
        mock_sequencing_group,
        mock_multicohort,
    ).new_bash_job.return_value.command.call_args.args[0]

    assert 'off_design' not in command
    assert 'mv dragen.vcf.gz merged.vcf.gz' in command
