"""How the post-hoc caller's JVM heap is sized, and when a bad pairing fails.

`configure_job` reads this stage's container memory from its own config section, so a heap
derived from `cpu` alone can ask for more than the container was given. The JVM is then killed
for exceeding its limit, and an OOM kill says only that: nothing in the log names the heap. So
the heap comes from the same tier the container does, and a pairing that cannot fit a usable
heap fails while the graph is being built.

No batch runs here: the job's command is read off the mock the stage builds it on.
"""

from pathlib import Path
from unittest.mock import MagicMock

import cpg_utils.config
import pytest

from popgen_rbceq2 import stage_support
from popgen_rbceq2.stages import pipeline
from tests.helpers import set_config

pytestmark = pytest.mark.fast


def _posthoc_config(shm_tmp_path: Path, stage_section: dict[str, object]) -> None:
    """Write a config setting the post-hoc stage's compute keys."""
    set_config(
        {
            'references': {
                'genome_build': 'GRCh38',
                'broad': {'ref_fasta': 'gs://bucket/ref.fasta'},
            },
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'exome',
                'driver_image': 'stub-driver:1.0',
                'posthoc_genotype_off_target_sites': stage_section,
                # An exome run's output tree is keyed on its design, so naming where this
                # stage writes needs one even though the stage never reads the design itself.
                stage_support.EXOME_DESIGN_KEY: 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed',
            },
        },
        shm_tmp_path / 'posthoc-compute.toml',
    )


def _posthoc_command(mocker, sequencing_group) -> str:
    """Queue the post-hoc stage and return the command it built."""
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    pipeline.PosthocGenotypeOffTargetSites().queue_jobs(sequencing_group, MagicMock())
    return batch.new_bash_job.return_value.command.call_args.args[0]


@pytest.mark.parametrize(
    ('memory', 'cpu', 'expected_heap'),
    [('standard', 1, 2), ('standard', 2, 5), ('highmem', 2, 11)],
)
def test_the_jvm_heap_is_sized_from_the_memory_tier_the_job_was_given(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
    memory: str,
    cpu: int,
    expected_heap: int,
):
    # The heap used to come from cpu alone while configure_job took the container's memory from
    # this same section. Nothing kept the two in step, so cpu = 2 with memory = "lowmem" asked
    # the JVM for 3g inside a ~2Gb container and the job was killed with nothing naming the heap.
    _posthoc_config(shm_tmp_path, {'memory': memory, 'cpu': cpu})

    assert f'-Xmx{expected_heap}g' in _posthoc_command(mocker, exome_sequencing_group)


def test_lowmem_is_refused_before_the_job_runs(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # Hail grants lowmem about 0.9 GiB per core, not the 1 an earlier table credited it with,
    # so the heap that table sized overran the container by the JVM's own overhead. Failing
    # here beats being killed mid-run: an OOM kill on a JVM says only that the container was
    # exceeded.
    _posthoc_config(shm_tmp_path, {'memory': 'lowmem', 'cpu': 4})

    with pytest.raises(cpg_utils.config.ConfigError) as raised:
        _posthoc_command(mocker, exome_sequencing_group)

    assert 'lowmem' in str(raised.value)


def test_a_memory_value_that_is_not_a_tier_fails_rather_than_guessing_a_heap(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # Hail accepts an explicit size here, and this stage does not: there would be nothing
    # keeping -Xmx inside it, which is the whole failure being fixed.
    _posthoc_config(shm_tmp_path, {'memory': '8Gi', 'cpu': 2})

    with pytest.raises(cpg_utils.config.ConfigError):
        _posthoc_command(mocker, exome_sequencing_group)
