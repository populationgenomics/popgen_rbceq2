"""Registration surface of the rbceq2 blood-group outputs.

Covers the `update_analysis_meta` callbacks that turn a written TSV into Metamist Analysis
meta, and the cohort concatenation the combine stage drives. Neither had coverage before
the QC TSV was added to them.
"""

from pathlib import Path

import pytest
from click.testing import CliRunner

from popgen_rbceq2 import analysis_meta
from popgen_rbceq2.jobs.rbceq2_gather_job import concat_tsvs, main
from tests.helpers import set_config

pytestmark = pytest.mark.fast

QC_TSV = 'UUID: abc123\tJK\tVEL\tFUT2\nSG000001\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)\tNA\n'


def _write_config(build: str, config_path: Path) -> None:
    """Point the config at one reference build, with the thresholds the QC meta records."""
    set_config(
        {
            'references': {'genome_build': build},
            'workflow': {'flag_blood_group_call_qc': {'min_depth': 12, 'min_gq': 25}},
        },
        config_path,
    )


@pytest.fixture
def qc_config(shm_tmp_path: Path):
    """Config supplying the thresholds and reference build that the QC meta records."""
    _write_config('GRCh38', shm_tmp_path / 'config.toml')


@pytest.mark.usefixtures('qc_config')
def test_call_qc_meta_parses_the_flags_and_records_the_thresholds(shm_tmp_path: Path):
    # The flag is only interpretable against the thresholds it was produced at, so the
    # Analysis has to carry both or a reader cannot tell what LOWQ meant on that run.
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text(QC_TSV)

    meta = analysis_meta.call_qc(str(qc_path))

    assert meta['blood_group_qc_flags'] == {
        'JK': 'PASS',
        'VEL': 'LOWQ:1:3774964(A>G,DP=8,GQ=45)',
        'FUT2': 'NA',
    }
    assert (meta['min_depth'], meta['min_gq']) == (12, 25)
    assert meta['reference_genome'] == 'GRCh38'


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


def _gather_args(tmp_path: Path, sgs: tuple[str, ...]) -> list[str]:
    keys = ('geno', 'pheno_numeric', 'pheno_alphanumeric', 'qc')
    args: list[str] = []
    for key in keys:
        args += [f'--output-{key.replace("_", "-")}', str(tmp_path / f'combined.{key}.tsv')]
    for key in keys:
        for sg in sgs:
            args += [f'--{key.replace("_", "-")}', str(_write_per_sg(tmp_path, sg, key, 'PASS'))]
    return args


def test_gather_writes_a_combined_qc_tsv_alongside_the_call_tsvs(shm_tmp_path: Path):
    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, ('SG000001', 'SG000002')))

    assert result.exit_code == 0, result.output
    for key in ('geno', 'pheno_numeric', 'pheno_alphanumeric', 'qc'):
        combined = (shm_tmp_path / f'combined.{key}.tsv').read_text().splitlines()
        assert len(combined) == 3, f'{key}: expected one header plus two sample rows'


def test_gather_keeps_the_detailed_site_level_flag_in_the_cohort_qc_tsv(shm_tmp_path: Path):
    # The cohort file is the one most readers open, so it carries the full flag rather than
    # a coarse category — the category is recoverable as the text before the first colon.
    flag = 'LOWQ:1:3774964(A>G,DP=8,GQ=45)'
    args = _gather_args(shm_tmp_path, ('SG000001',))
    qc_path = _write_per_sg(shm_tmp_path, 'SG000002', 'qc', flag)
    for key in ('geno', 'pheno_numeric', 'pheno_alphanumeric'):
        args += [f'--{key.replace("_", "-")}', str(_write_per_sg(shm_tmp_path, 'SG000002', key, 'PASS'))]
    args += ['--qc', str(qc_path)]

    result = CliRunner().invoke(main, args)

    assert result.exit_code == 0, result.output
    assert (shm_tmp_path / 'combined.qc.tsv').read_text().splitlines()[-1] == f'SG000002\tPASS\t{flag}'


def test_gather_raises_when_a_qc_input_is_missing(shm_tmp_path: Path):
    # A silently short cohort QC TSV would read as "every sample passed".
    args = _gather_args(shm_tmp_path, ('SG000001',))
    args += ['--qc', str(shm_tmp_path / 'SG000002.qc.tsv')]

    result = CliRunner().invoke(main, args)

    assert result.exit_code != 0
    assert 'blood-group qc inputs are missing or empty' in str(result.exception)


def test_gather_requires_at_least_one_qc_input(shm_tmp_path: Path):
    args = _gather_args(shm_tmp_path, ('SG000001',))
    qc_flag_at = args.index('--qc')
    del args[qc_flag_at : qc_flag_at + 2]

    result = CliRunner().invoke(main, args)

    assert result.exit_code != 0
    assert 'No blood-group qc inputs provided' in str(result.exception)
