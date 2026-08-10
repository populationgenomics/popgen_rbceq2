"""Registration surface of the rbceq2 blood-group outputs.

Covers the `update_analysis_meta` callbacks that turn a written TSV into Metamist Analysis
meta, and the cohort concatenation the combine stage drives. Neither had coverage before
the QC TSV was added to them.
"""

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from click.testing import CliRunner

from popgen_rbceq2 import analysis_meta
from popgen_rbceq2.jobs.rbceq2_gather_job import MANIFEST_KEYS, concat_tsvs, main
from popgen_rbceq2.stages import pipeline
from tests.helpers import set_config

pytestmark = pytest.mark.fast

QC_TSV = 'UUID: abc123\tJK\tVEL\tFUT2\nSG000001\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)\tNA\n'


def _write_config(build: str, config_path: Path) -> None:
    """Point the config at one reference build."""
    set_config({'references': {'genome_build': build}}, config_path)


@pytest.fixture
def qc_config(shm_tmp_path: Path):
    """Config supplying the reference build that the QC meta records."""
    _write_config('GRCh38', shm_tmp_path / 'config.toml')


@pytest.mark.usefixtures('qc_config')
def test_call_qc_meta_parses_the_flags(shm_tmp_path: Path):
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text(QC_TSV)

    meta = analysis_meta.call_qc(str(qc_path))

    assert meta['blood_group_qc_flags'] == {
        'JK': 'PASS',
        'VEL': 'LOWQ:1:3774964(A>G,DP=8,GQ=45)',
        'FUT2': 'NA',
    }
    assert meta['reference_genome'] == 'GRCh38'


def test_call_qc_stage_meta_records_the_thresholds_it_gave_the_job(mocker, mock_sequencing_group, shm_tmp_path: Path):
    # The flag is only interpretable against the thresholds it was produced at, so the
    # Analysis has to carry both or a reader cannot tell what LOWQ meant on that run. They
    # are recorded by queue_jobs from the same values it hands the job — not re-read from
    # config by the meta callback — so the recorded and applied thresholds cannot diverge.
    set_config(
        {
            'references': {'genome_build': 'GRCh38'},
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'genome',
                'driver_image': 'stub-driver:1.0',
                'flag_blood_group_call_qc': {'min_depth': 12, 'min_gq': 25},
            },
        },
        shm_tmp_path / 'thresholds.toml',
    )
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    inputs = MagicMock()
    inputs.as_str.return_value = 'gs://bucket/SG000001.geno.tsv'

    output = pipeline.FlagBloodGroupCallQc().queue_jobs(mock_sequencing_group, inputs)

    assert output is not None
    assert output.meta is not None
    assert (output.meta['min_depth'], output.meta['min_gq']) == (12, 25)
    command = batch.new_bash_job.return_value.command.call_args.args[0]
    assert '--min-depth 12' in command
    assert '--min-gq 25' in command


@pytest.mark.usefixtures('qc_config')
def test_call_qc_meta_does_not_set_stage(shm_tmp_path: Path):
    # cpg-flow injects meta.stage as the Stage class name via get_job_attrs; setting it
    # here would override that with a hand-typed string.
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text(QC_TSV)

    assert 'stage' not in analysis_meta.call_qc(str(qc_path))


@pytest.mark.usefixtures('qc_config')
def test_call_qc_meta_raises_on_a_header_only_tsv(shm_tmp_path: Path):
    # Registration must fail loudly rather than record an empty flag map.
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text('UUID: abc123\tJK\tVEL\n')

    with pytest.raises(ValueError, match='Missing rows'):
        analysis_meta.call_qc(str(qc_path))


def test_call_qc_meta_tracks_the_configured_reference_build(shm_tmp_path: Path):
    # The build has to come from references.genome_build, the same source the stages pass to
    # rbceq2 and to the site extract. Read from anywhere else and the Analysis can record a
    # build the run did not use.
    _write_config('GRCh37', shm_tmp_path / 'grch37.toml')
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text(QC_TSV)

    assert analysis_meta.call_qc(str(qc_path))['reference_genome'] == 'GRCh37'


def test_blood_group_calls_meta_tracks_the_configured_reference_build(shm_tmp_path: Path):
    # Same source for the calls Analysis, so the two per-SG records cannot disagree about
    # which build produced them.
    _write_config('GRCh37', shm_tmp_path / 'grch37_sg.toml')
    geno_path = shm_tmp_path / 'SG000001.geno.tsv'
    geno_path.write_text('UUID: abc123\tJK\nSG000001\tJK*01\n')
    (shm_tmp_path / 'SG000001.pheno_alphanumeric.tsv').write_text('UUID: abc123\tJK\nSG000001\tJk(a+b-)\n')

    meta = analysis_meta.blood_group_calls(str(geno_path))

    assert meta['reference_genome'] == 'GRCh37'
    assert meta['blood_group_genotypes'] == {'JK': 'JK*01'}
    assert meta['blood_group_phenotypes'] == {'JK': 'Jk(a+b-)'}


def test_cohort_calls_meta_points_at_the_sibling_qc_tsv():
    # The Analysis is registered against the geno TSV, so the QC TSV is only discoverable
    # through the meta.
    meta = analysis_meta.cohort_calls('gs://bucket/combined.COH1.geno.tsv')

    assert meta['qc_path'] == 'gs://bucket/combined.COH1.qc.tsv'
    assert meta['pheno_numeric_path'] == 'gs://bucket/combined.COH1.pheno_numeric.tsv'


def test_concat_tsvs_keeps_one_header_and_every_sample_row():
    first = 'UUID: a\tJK\tVEL\nSG000001\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)\n'
    second = 'UUID: b\tJK\tVEL\nSG000002\tPASS\tPASS\n'

    assert concat_tsvs([first, second]).splitlines() == [
        'UUID: a\tJK\tVEL',
        'SG000001\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)',
        'SG000002\tPASS\tPASS',
    ]


def _write_per_sg(tmp_path: Path, sg: str, key: str, cell: str) -> Path:
    path = tmp_path / f'{sg}.{key}.tsv'
    path.write_text(f'UUID: {sg}\tJK\tVEL\n{sg}\tPASS\t{cell}\n')
    return path


def _manifest_for(tmp_path: Path, sgs: tuple[str, ...]) -> dict[str, list[str]]:
    """Write per-SG TSVs and return the {key: [path, ...]} manifest dict the gather job reads."""
    return {key: [str(_write_per_sg(tmp_path, sg, key, 'PASS')) for sg in sgs] for key in MANIFEST_KEYS}


def _gather_args(tmp_path: Path, manifest: dict[str, list[str]]) -> list[str]:
    args: list[str] = []
    for key in MANIFEST_KEYS:
        args += [f'--output-{key.replace("_", "-")}', str(tmp_path / f'combined.{key}.tsv')]
    manifest_path = tmp_path / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest))
    return [*args, '--manifest', str(manifest_path)]


def test_gather_writes_a_combined_qc_tsv_alongside_the_call_tsvs(shm_tmp_path: Path):
    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, _manifest_for(shm_tmp_path, ('SG000001', 'SG000002'))))

    assert result.exit_code == 0, result.output
    for key in ('geno', 'pheno_numeric', 'pheno_alphanumeric', 'qc'):
        combined = (shm_tmp_path / f'combined.{key}.tsv').read_text().splitlines()
        assert len(combined) == 3, f'{key}: expected one header plus two sample rows'


def test_gather_keeps_the_detailed_site_level_flag_in_the_cohort_qc_tsv(shm_tmp_path: Path):
    # The cohort file is the one most readers open, so it carries the full flag rather than
    # a coarse category — the category is recoverable as the text before the first colon.
    flag = 'LOWQ:1:3774964(A>G,DP=8,GQ=45)'
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    manifest['qc'].append(str(_write_per_sg(shm_tmp_path, 'SG000002', 'qc', flag)))
    for key in ('geno', 'pheno_numeric', 'pheno_alphanumeric'):
        manifest[key].append(str(_write_per_sg(shm_tmp_path, 'SG000002', key, 'PASS')))

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code == 0, result.output
    assert (shm_tmp_path / 'combined.qc.tsv').read_text().splitlines()[-1] == f'SG000002\tPASS\t{flag}'


def test_gather_raises_when_a_qc_input_is_missing(shm_tmp_path: Path):
    # A silently short cohort QC TSV would read as "every sample passed".
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    manifest['qc'].append(str(shm_tmp_path / 'SG000002.qc.tsv'))

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'blood-group qc inputs are missing or empty' in str(result.exception)


def test_gather_requires_at_least_one_qc_input(shm_tmp_path: Path):
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    manifest['qc'] = []

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'No blood-group qc inputs provided' in str(result.exception)


def test_gather_rejects_a_manifest_with_the_wrong_keys(shm_tmp_path: Path):
    # A missing key means the stage and the job disagree about the manifest shape; combining
    # the remaining three would silently drop an output.
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    del manifest['qc']

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'exactly the keys' in str(result.exception)


def test_combine_stage_writes_the_manifest_the_job_reads(mocker, mock_cohort, shm_tmp_path: Path):
    # The per-SG paths travel to the job in a manifest, not argv: at four repeated flags per
    # sequencing group the command line would outgrow ARG_MAX around ~3,000 of them.
    set_config(
        {
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'genome',
                'driver_image': 'stub-driver:1.0',
            },
        },
        shm_tmp_path / 'combine.toml',
    )
    mock_cohort.dataset.prefix.side_effect = lambda category=None: (
        shm_tmp_path if category == 'tmp' else Path('gs://bucket')
    )
    sgs = []
    for sg_id in ('SG000001', 'SG000002'):
        sg = MagicMock()
        sg.id = sg_id
        sg.gvcf = f'gs://bucket/{sg_id}.g.vcf.gz'
        sgs.append(sg)
    mock_cohort.get_sequencing_groups.return_value = sgs
    inputs = MagicMock()
    inputs.as_dict_by_target.side_effect = [
        {sg.id: {key: f'gs://bucket/{sg.id}.{key}.tsv' for key in MANIFEST_KEYS[:3]} for sg in sgs},
        {sg.id: {'qc': f'gs://bucket/{sg.id}.qc.tsv'} for sg in sgs},
    ]
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    manifest_path = (
        shm_tmp_path / 'popgen_rbceq2' / 'CombineRbceq2OutputsPerCohort' / 'test-cohort' / 'v1'
    ) / 'test-cohort.manifest.json'
    manifest_path.parent.mkdir(parents=True)  # gs:// needs no directories; this local stand-in does

    output = pipeline.CombineRbceq2OutputsPerCohort().queue_jobs(mock_cohort, inputs)

    assert output is not None
    assert json.loads(manifest_path.read_text()) == {
        key: [f'gs://bucket/SG000001.{key}.tsv', f'gs://bucket/SG000002.{key}.tsv'] for key in MANIFEST_KEYS
    }
    command = batch.new_bash_job.return_value.command.call_args.args[0]
    assert f'--manifest {manifest_path}' in command


def test_combine_stage_raises_on_a_cohort_with_no_gvcf_sequencing_groups(mock_cohort):
    # The manifest moved the per-SG paths out of argv, so build_python_command's empty-list
    # guard no longer turns "no eligible sequencing groups" into a graph-build error; the
    # stage checks it itself, naming the cause instead of a CLI flag.
    sg = MagicMock()
    sg.id = 'SG000001'
    sg.gvcf = None
    mock_cohort.get_sequencing_groups.return_value = [sg]
    inputs = MagicMock()
    inputs.as_dict_by_target.side_effect = [{}, {}]

    with pytest.raises(ValueError, match='no sequencing groups with a gVCF'):
        pipeline.CombineRbceq2OutputsPerCohort().queue_jobs(mock_cohort, inputs)
