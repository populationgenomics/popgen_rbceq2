"""The post-hoc caller's retry on the GCS egress quota, run for real in bash.

Streaming CRAMs for two exome cohorts at once can exceed the project's GoogleEgressBandwidth
quota for a moment, and the job that draws the 429 fails while the rest of the batch carries on.
The retry exists so that job waits and succeeds instead of needing a top-up run. It must retry
only that failure: anything else has to fail on the first attempt, as it did before.

A stub stands in for HaplotypeCaller, counting its attempts in a file and failing as told. The
backoff base is 0, so no test sleeps.
"""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from popgen_rbceq2 import stage_support
from popgen_rbceq2.stages import pipeline
from popgen_rbceq2.stages.blood_group_genotyping import posthoc_genotype
from tests.helpers import set_config

pytestmark = pytest.mark.fast

QUOTA_MESSAGE = (
    'Caused by: java.io.IOException: com.google.cloud.storage.StorageException: 429 Too Many Requests\n'
    'This workload is drawing too much egress bandwidth from Cloud Storage and has exceeded the '
    'GoogleEgressBandwidth (google_egress_bandwidth) Quota in the australia-southeast1 region.'
)


def _stub(tmp_path: Path, failures: int, message: str) -> Path:
    """Write a stub that fails `failures` times printing `message`, then succeeds."""
    count = tmp_path / 'attempts'
    count.write_text('0')
    (tmp_path / 'message').write_text(message)
    stub = tmp_path / 'stub.sh'
    stub.write_text(
        f"""
n=$(( $(cat {count}) + 1 ))
echo "$n" > {count}
if [ "$n" -le {failures} ]; then
    cat {tmp_path / 'message'} >&2
    exit 3
fi
echo "stub succeeded on attempt $n"
"""
    )
    return stub


def _run(tmp_path: Path, stub: Path, retries: int = 3) -> tuple[subprocess.CompletedProcess[str], int]:
    """Run the stub under the retry, and return the result and how many attempts were made."""
    script = 'set -euo pipefail\n' + posthoc_genotype.retry_on_egress_quota(
        f'bash {stub}', retries=retries, base_seconds=0
    )
    result = subprocess.run(['bash', '-c', script], capture_output=True, text=True, check=False)  # noqa: S603, S607
    return result, int((tmp_path / 'attempts').read_text())


def test_a_quota_failure_is_retried_until_the_command_succeeds(tmp_path):
    result, attempts = _run(tmp_path, _stub(tmp_path, failures=2, message=QUOTA_MESSAGE))

    assert result.returncode == 0
    assert attempts == 3
    assert result.stderr.count('hit the GCS egress quota; retrying') == 2


def test_any_other_failure_is_not_retried(tmp_path):
    # A wrong reference or a corrupt CRAM fails the same way on every attempt, so retrying it
    # would only spend four HaplotypeCaller runs finding that out.
    result, attempts = _run(tmp_path, _stub(tmp_path, failures=1, message='A USER ERROR has occurred: bad input'))

    assert result.returncode == 3
    assert attempts == 1
    assert 'retrying' not in result.stderr


def test_the_quota_on_every_attempt_fails_after_the_last_retry_and_says_so(tmp_path):
    result, attempts = _run(tmp_path, _stub(tmp_path, failures=10, message=QUOTA_MESSAGE))

    assert result.returncode == 3
    assert attempts == posthoc_genotype.EGRESS_RETRIES + 1
    assert 'still over the GCS egress quota after 3 retries' in result.stderr


def test_the_output_of_a_failed_attempt_still_reaches_the_job_log(tmp_path):
    # The attempt's output is kept to classify the failure, and must not be swallowed doing so:
    # the job log is where anyone reads why HaplotypeCaller failed.
    result, _ = _run(tmp_path, _stub(tmp_path, failures=1, message='A USER ERROR has occurred: bad input'))

    assert 'A USER ERROR has occurred: bad input' in result.stdout


def test_the_stage_runs_haplotypecaller_under_the_retry_and_checks_the_index_after_it(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    set_config(
        {
            'references': {'genome_build': 'GRCh38', 'broad': {'ref_fasta': 'gs://bucket/ref.fasta'}},
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'exome',
                'driver_image': 'stub-driver:1.0',
                stage_support.EXOME_DESIGN_KEY: 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed',
            },
        },
        shm_tmp_path / 'posthoc-retry.toml',
    )
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    pipeline.PosthocGenotypeOffTargetSites().queue_jobs(exome_sequencing_group, MagicMock())
    command = batch.new_bash_job.return_value.command.call_args.args[0]

    loop = command.index('for attempt in $(seq 1 4)')
    caller = command.index('HaplotypeCaller')
    index_check = command.index('wrote no index')
    assert loop < caller < index_check
