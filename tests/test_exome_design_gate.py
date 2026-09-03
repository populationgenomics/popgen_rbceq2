"""The capture design an exome run fills holes outside of, and how the stage demands it.

Post-hoc calls fill only defining sites the capture design never targeted. Which design that
is cannot be inferred from the gVCF, so it is configured, as a `[references]` key resolved
through `reference_path` like every other reference this pipeline reads.

Two things have to hold and neither has a type error to catch it. An exome run that does not
name a design must fail while the graph is being built, not on the first sequencing group's
merge, because the whole cohort is wasted either way and only one of those says why. A genome
run must never read the key at all, since it has no capture edge and no post-hoc input.

No batch runs here: the job's command is read off the mock the stage builds it on.
"""

from pathlib import Path
from unittest.mock import MagicMock

import cpg_utils.config
import pytest

from popgen_rbceq2.stages import pipeline
from popgen_rbceq2.stages.blood_group_genotyping.filter_and_convert import EXOME_DESIGN_KEY
from tests.helpers import set_config

pytestmark = pytest.mark.fast

TWIST_KEY = 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'
TWIST_BED = 'gs://cpg-common-main/references/exome-probesets/hg38/Twist_VCGS_Exome_Covered_Targets_hg38.bed'


def _config(shm_tmp_path: Path, sequencing_type: str, design_key: str | None) -> None:
    """Write a config for one sequencing type, with or without the design key."""
    stage_section: dict[str, str] = {} if design_key is None else {EXOME_DESIGN_KEY: design_key}
    set_config(
        {
            'references': {
                'genome_build': 'GRCh38',
                'exome_probesets_hg38': {'twist_vcgs_custom_exome_covered_targets_bed': TWIST_BED},
            },
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': sequencing_type,
                'driver_image': 'stub-driver:1.0',
                'filter_and_convert_gvcfs_for_rbceq2': stage_section,
            },
        },
        shm_tmp_path / f'{sequencing_type}-design.toml',
    )


def _queue(mocker, sequencing_group) -> MagicMock:
    """Queue the conversion stage's jobs on a mock batch and return that batch."""
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    inputs = MagicMock()
    inputs.as_dict.return_value = {
        'gvcf': 'gs://bucket/SG000001.posthoc.g.vcf.gz',
        'index': 'gs://bucket/SG000001.posthoc.g.vcf.gz.tbi',
    }
    pipeline.FilterAndConvertGvcfsForRbceq2().queue_jobs(sequencing_group, inputs)
    return batch


def test_an_exome_run_without_a_design_names_the_key_it_wants(mocker, exome_sequencing_group, shm_tmp_path: Path):
    # Failing here is the point. The alternative to a required key is a default, and every
    # default is wrong for some cohort: filling holes inside the capture would answer, with a
    # second caller, a question about that sample's DRAGEN run.
    _config(shm_tmp_path, 'exome', design_key=None)

    with pytest.raises(cpg_utils.config.ConfigError) as raised:
        _queue(mocker, exome_sequencing_group)

    assert EXOME_DESIGN_KEY in str(raised.value)


def test_an_exome_run_naming_a_design_that_is_not_configured_fails(mocker, exome_sequencing_group, shm_tmp_path: Path):
    # A typo in the key resolves to nothing. reference_path raising is what makes that a
    # graph-build failure rather than a job that localises an empty file and fills nothing.
    _config(shm_tmp_path, 'exome', design_key='exome_probesets_hg38/no_such_design_bed')

    with pytest.raises(cpg_utils.config.ConfigError):
        _queue(mocker, exome_sequencing_group)


def test_an_exome_merge_reads_the_configured_design_and_gates_the_fill_on_it(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # The BED has to reach the job, and the fill has to be the sites outside it that DRAGEN
    # also said nothing at — not the holes alone, which is what the stage did before.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    batch = _queue(mocker, exome_sequencing_group)

    assert TWIST_BED in [call.args[0] for call in batch.read_input.call_args_list]
    command = batch.new_bash_job.return_value.command.call_args.args[0]
    assert 'off_design_sites.bed' in command
    # The fill is the intersection: holes are narrowed to off-design sites, never the reverse.
    assert 'covered.bed off_design_sites.bed > uncovered.bed' in command


def test_an_exome_merge_fails_when_dragen_called_outside_the_configured_design(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # DRAGEN emits over exactly the target BED with no padding, so a record at an off-design
    # defining site means the configured file is not the one the gVCF was called against.
    # Measured on the validation cohorts: a target-regions file where the gVCF used the probe
    # footprint leaves 309 such sites and silently gates off three quarters of the recoveries.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    command = _queue(mocker, exome_sequencing_group).new_bash_job.return_value.command.call_args.args[0]

    assert 'off_design_called.bed' in command
    assert 'exit 1' in command


def test_a_genome_run_never_asks_for_a_design(mocker, mock_sequencing_group, shm_tmp_path: Path):
    # A genome gVCF has no capture edge, gets no post-hoc input, and must not be made to
    # carry an exome-only key. Its command is unchanged by the gate existing.
    _config(shm_tmp_path, 'genome', design_key=None)

    command = _queue(mocker, mock_sequencing_group).new_bash_job.return_value.command.call_args.args[0]

    assert 'off_design_sites.bed' not in command
    assert 'mv dragen.vcf.gz merged.vcf.gz' in command
