"""Fixtures shared by the test suite.

Every test here builds objects and inspects them; none runs a batch or reaches a cloud, so
there is no Hail initialisation and no credentials are needed.
"""

import os
from pathlib import Path
from unittest import mock
from unittest.mock import MagicMock

import pytest
from cpg_utils.cloud import DockerImage
from google.auth import environment_vars

from tests.helpers import set_config


@pytest.fixture(scope='session', autouse=True)
def stub_image_lookups():
    """Resolve image_path(key, version) without Artifact Registry.

    Passing a version makes image_path list the cpg-common repository, which needs credentials
    CI does not have. Tests build jobs but never run a batch, so no image is ever pulled and the
    tag only has to be a plausible string.
    """

    def fake_find_image(repository: str | None, name: str, version: str) -> DockerImage:
        repo = f'images-{repository}' if repository is not None else 'images'
        uri = f'australia-southeast1-docker.pkg.dev/cpg-common/{repo}/{name}'
        return DockerImage(name, uri, f'{uri}:{version}', '0', '')

    with mock.patch('cpg_utils.config.find_image', fake_find_image):
        yield


@pytest.fixture(scope='session')
def test_storage_root(tmp_path_factory) -> Path:
    """Root directory for test scratch, on a ramdisk where there is one."""
    shm_base = Path('/dev/shm/popgen_rbceq2_tests')  # noqa: S108
    try:
        if shm_base.parent.exists():
            shm_base.mkdir(exist_ok=True, parents=True)
            return shm_base
    except (PermissionError, OSError):
        pass
    return tmp_path_factory.mktemp('popgen_rbceq2_tests')


@pytest.fixture
def shm_tmp_path(test_storage_root, tmp_path) -> Path:
    """Replacement for tmp_path that prefers /dev/shm.

    Uses the unique name from pytest's tmp_path to avoid collisions.
    """
    path = test_storage_root / tmp_path.name
    path.mkdir(exist_ok=True, parents=True)
    return path


@pytest.fixture(autouse=True)
def dummy_gcp_project():
    """A project name, so code reaching for one does not fail on its absence."""
    with mock.patch.dict(os.environ, {environment_vars.PROJECT: 'dummy-project-for-tests'}):
        yield


@pytest.fixture
def mock_cohort(mocker, shm_tmp_path: Path):
    """A cohort in a dataset with known bucket prefixes, and an active workflow to name them."""
    set_config(
        {'workflow': {'name': 'popgen_rbceq2', 'version': 'v1', 'sequencing_type': 'genome'}},
        shm_tmp_path / 'config.toml',
    )

    mock_wf = MagicMock()
    mock_wf.name = 'popgen_rbceq2'
    mock_wf.status_reporter = MagicMock()
    mocker.patch('cpg_flow.workflow.get_workflow', return_value=mock_wf)
    mocker.patch('cpg_flow.stage.get_workflow', return_value=mock_wf)

    mock_dataset = MagicMock()
    mock_dataset.prefix.side_effect = lambda category=None: Path(
        'gs://bucket-tmp' if category == 'tmp' else 'gs://bucket',
    )

    cohort = MagicMock()
    cohort.dataset = mock_dataset
    cohort.name = 'test-cohort'
    cohort.id = 'test-cohort'
    cohort.get_sequencing_groups.return_value = [MagicMock() for _ in range(3)]
    return cohort


@pytest.fixture
def mock_sequencing_group(mock_cohort):
    """A sequencing group with a gVCF, for the SequencingGroupStage outputs."""
    sg = MagicMock()
    sg.dataset = mock_cohort.dataset
    sg.id = 'SG000001'
    sg.gvcf = 'gs://bucket/SG000001.g.vcf.gz'
    return sg
