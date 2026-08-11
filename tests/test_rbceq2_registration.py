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
from popgen_rbceq2.jobs.rbceq2_gather_job import (
    MANIFEST_KEYS,
    NOT_REPORTED,
    SEQUENCING_GROUPS_KEY,
    concat_tsvs,
    main,
)
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


# rbceq2 labels a row with the base name of the VCF it read, not the sample column, so a
# per-SG TSV says `<sg>.converted.vcf`. Tests feed that form so they exercise the relabelling
# every real input goes through.
def _rbceq2_tsv(uuid: str, systems: str, sg: str, cells: str) -> str:
    return f'UUID: {uuid}\t{systems}\n{sg}.converted.vcf\t{cells}\n'


def test_concat_tsvs_keys_rows_on_the_sequencing_group_not_rbceq2s_label():
    # The whole point of keying the manifest by sequencing group: rbceq2 writes an
    # intermediate file name in column 1, which no Metamist query can join to.
    first = _rbceq2_tsv('a', 'JK\tVEL', 'SG000001', 'PASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)')
    second = _rbceq2_tsv('b', 'JK\tVEL', 'SG000002', 'PASS\tPASS')

    assert concat_tsvs([('SG000001', first), ('SG000002', second)]).splitlines() == [
        'UUID: a\tJK\tVEL',
        'SG000001\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)',
        'SG000002\tPASS\tPASS',
    ]


def test_concat_tsvs_raises_when_rbceq2s_label_names_another_sequencing_group():
    # A path keyed under the wrong sequencing group, or a converted VCF renamed out from under
    # the assumption. Relabelling regardless would put one sample's calls under another's ID,
    # which is worse than the shift this whole job exists to prevent.
    tsv = _rbceq2_tsv('a', 'JK\tVEL', 'SG000002', 'PASS\tPASS')

    with pytest.raises(ValueError, match='does not name this sequencing group'):
        concat_tsvs([('SG000001', tsv)])


def test_concat_tsvs_raises_on_a_file_holding_more_than_one_sample_row():
    # One sequencing group is keyed to one file, so a second row has no ID to be keyed on.
    tsv = 'UUID: a\tJK\nSG000001.converted.vcf\tJK*01/JK*01\nSG000002.converted.vcf\tJK*01/JK*02\n'

    with pytest.raises(ValueError, match='Expected one sample row'):
        concat_tsvs([('SG000001', tsv)])


def test_concat_tsvs_aligns_cells_by_system_when_the_column_sets_differ():
    # rbceq2 emits a column only for a blood group it found a matching allele for, one frame
    # per run, so two sequencing groups can disagree on the column set. Appending the second
    # row under the first sample's header would file VEL's flag under FUT2.
    first = _rbceq2_tsv('a', 'FUT2\tJK\tVEL', 'SG000001', 'PASS\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)')
    second = _rbceq2_tsv('b', 'JK\tVEL', 'SG000002', 'PASS\tNA')

    assert concat_tsvs([('SG000001', first), ('SG000002', second)]).splitlines() == [
        'UUID: a\tFUT2\tJK\tVEL',
        'SG000001\tPASS\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)',
        f'SG000002\t{NOT_REPORTED}\tPASS\tNA',
    ]


def test_concat_tsvs_covers_a_system_only_a_later_sample_reports():
    # The union is taken across every input, not just the first, so a system the first sample
    # is missing still gets a column rather than dropping that sample's call on the floor.
    first = _rbceq2_tsv('a', 'JK', 'SG000001', 'JK*01/JK*02')
    second = _rbceq2_tsv('b', 'JK\tVEL', 'SG000002', 'JK*01/JK*01\tVEL*01/VEL*01')

    assert concat_tsvs([('SG000001', first), ('SG000002', second)]).splitlines() == [
        'UUID: a\tJK\tVEL',
        f'SG000001\tJK*01/JK*02\t{NOT_REPORTED}',
        'SG000002\tJK*01/JK*01\tVEL*01/VEL*01',
    ]


def test_concat_tsvs_sorts_the_union_rather_than_appending_unseen_systems():
    # The union has to come out in rbceq2's own alphabetical order, not first-seen order. Every
    # other case here has the first sample's columns already sorting the same as the union, so
    # an implementation that appended unseen systems to the end would pass all of them and
    # still write a cohort header disagreeing with every per-sample file it was built from.
    first = _rbceq2_tsv('a', 'VEL', 'SG000001', 'VEL*01/VEL*01')
    second = _rbceq2_tsv('b', 'ABO\tVEL', 'SG000002', 'ABO*A1.01/ABO*O.01.01\tVEL*01/VEL*02')

    assert concat_tsvs([('SG000001', first), ('SG000002', second)]).splitlines() == [
        'UUID: a\tABO\tVEL',
        f'SG000001\t{NOT_REPORTED}\tVEL*01/VEL*01',
        'SG000002\tABO*A1.01/ABO*O.01.01\tVEL*01/VEL*02',
    ]


def test_concat_tsvs_realigns_a_sample_whose_columns_are_in_a_different_order():
    # The same column set in a different order is the cell shift in its purest form: the row
    # lengths match, so nothing about the file looks wrong, and every cell lands under the
    # wrong system.
    first = _rbceq2_tsv('a', 'JK\tVEL', 'SG000001', 'JK*01/JK*02\tVEL*01/VEL*01')
    second = _rbceq2_tsv('b', 'VEL\tJK', 'SG000002', 'VEL*01/VEL*02\tJK*01/JK*01')

    assert concat_tsvs([('SG000001', first), ('SG000002', second)]).splitlines() == [
        'UUID: a\tJK\tVEL',
        'SG000001\tJK*01/JK*02\tVEL*01/VEL*01',
        'SG000002\tJK*01/JK*01\tVEL*01/VEL*02',
    ]


def test_concat_tsvs_raises_on_a_row_that_does_not_match_its_own_header():
    # Differing column sets between samples are expected; a row that disagrees with the
    # header it was written under is a corrupt file, and filling it in would invent cells.
    with pytest.raises(ValueError, match='Ragged TSV'):
        concat_tsvs([('SG000001', 'UUID: a\tJK\tVEL\nSG000001.converted.vcf\tPASS\n')])


def test_concat_tsvs_raises_on_a_repeated_system_column():
    # Matching cells to columns by name means a repeated name has no answer: keeping either
    # cell drops the other, and the row is the right length, so nothing else here would catch
    # it. The row count check would not; the union would just be one column short.
    tsv = _rbceq2_tsv('a', 'JK\tJK', 'SG000001', 'JK*01/JK*02\tJK*01/JK*01')

    with pytest.raises(ValueError, match='Repeated system column'):
        concat_tsvs([('SG000001', tsv)])


def test_concat_tsvs_raises_on_a_sample_with_no_system_columns():
    # The union fill covers samples whose column sets differ, not a sample rbceq2 reported
    # nothing for. Filling such a row would launder a broken per-sample run into all four
    # cohort TSVs as a row of NOT_REPORTED.
    first = _rbceq2_tsv('a', 'JK', 'SG000001', 'JK*01/JK*02')
    second = 'UUID: b\nSG000002.converted.vcf\n'

    with pytest.raises(ValueError, match='SG000002: the TSV has no system columns'):
        concat_tsvs([('SG000001', first), ('SG000002', second)])


def test_concat_tsvs_warning_names_the_samples_it_filled_capped_at_ten(caplog: pytest.LogCaptureFixture):
    # Whoever chases the warning needs a sequencing group to start from, not just a count,
    # but a cohort-sized list would swamp the log — the list caps at ten plus an ellipsis.
    sgs = [f'SG{i:06d}' for i in range(1, 13)]
    contents = [(sgs[0], _rbceq2_tsv('a', 'JK\tVEL', sgs[0], 'JK*01/JK*02\tVEL*01/VEL*01'))]
    contents += [(sg, _rbceq2_tsv('b', 'JK', sg, 'JK*01/JK*01')) for sg in sgs[1:]]

    with caplog.at_level('WARNING', logger='popgen_rbceq2.jobs.rbceq2_gather_job'):
        concat_tsvs(contents)

    (message,) = [r.message for r in caplog.records if NOT_REPORTED in r.message]
    assert f'VEL (11/12 samples: {", ".join(sgs[1:11])} ...)' in message


def test_concat_tsvs_is_empty_when_given_no_inputs():
    assert concat_tsvs([]) == ''


def _write_per_sg(tmp_path: Path, sg: str, key: str, cell: str, label: str | None = None) -> str:
    """Write one per-SG TSV and return its path.

    Args:
        tmp_path: Directory to write into.
        sg: The sequencing group the file belongs to.
        key: Output type, naming the file.
        cell: The VEL cell, so a test can put a flag there.
        label: What rbceq2 wrote in column 1, defaulting to the `<sg>.converted.vcf` it
            really writes. Set it to something else to stand in for a miskeyed file.
    """
    path = tmp_path / f'{sg}.{key}.tsv'
    path.write_text(f'UUID: {sg}\tJK\tVEL\n{label or f"{sg}.converted.vcf"}\tPASS\t{cell}\n')
    return str(path)


def _manifest_for(tmp_path: Path, sgs: tuple[str, ...]) -> dict:
    """Write per-SG TSVs and return the manifest dict the gather job reads."""
    return {SEQUENCING_GROUPS_KEY: list(sgs)} | {
        key: {sg: _write_per_sg(tmp_path, sg, key, 'PASS') for sg in sgs} for key in MANIFEST_KEYS
    }


def _gather_args(tmp_path: Path, manifest: dict) -> list[str]:
    args: list[str] = []
    for key in MANIFEST_KEYS:
        args += [f'--output-{key.replace("_", "-")}', str(tmp_path / f'combined.{key}.tsv')]
    manifest_path = tmp_path / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest))
    return [*args, '--manifest', str(manifest_path)]


def test_gather_writes_a_combined_qc_tsv_alongside_the_call_tsvs(shm_tmp_path: Path):
    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, _manifest_for(shm_tmp_path, ('SG000001', 'SG000002'))))

    assert result.exit_code == 0, result.output
    for key in MANIFEST_KEYS:
        combined = (shm_tmp_path / f'combined.{key}.tsv').read_text().splitlines()
        assert len(combined) == 3, f'{key}: expected one header plus two sample rows'
        assert [line.split('\t')[0] for line in combined[1:]] == ['SG000001', 'SG000002']


def test_gather_keeps_the_detailed_site_level_flag_in_the_cohort_qc_tsv(shm_tmp_path: Path):
    # The cohort file is the one most readers open, so it carries the full flag rather than
    # a coarse category — the category is recoverable as the text before the first colon.
    flag = 'LOWQ:1:3774964(A>G,DP=8,GQ=45)'
    manifest = _manifest_for(shm_tmp_path, ('SG000001', 'SG000002'))
    manifest['qc']['SG000002'] = _write_per_sg(shm_tmp_path, 'SG000002', 'qc', flag)

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code == 0, result.output
    assert (shm_tmp_path / 'combined.qc.tsv').read_text().splitlines()[-1] == f'SG000002\tPASS\t{flag}'


def test_gather_raises_when_a_qc_input_is_missing(shm_tmp_path: Path):
    # A silently short cohort QC TSV would read as "every sample passed".
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    manifest['qc']['SG000002'] = str(shm_tmp_path / 'SG000002.qc.tsv')

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'blood-group qc inputs are missing or empty' in str(result.exception)


def test_gather_requires_at_least_one_qc_input(shm_tmp_path: Path):
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    manifest['qc'] = {}

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


def test_gather_rejects_a_manifest_without_the_sequencing_groups_key(shm_tmp_path: Path):
    # Without it there is nothing to hold the four types to the same set, and a short type
    # would write a plausible file covering fewer samples than the cohort holds.
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    del manifest[SEQUENCING_GROUPS_KEY]

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'exactly the keys' in str(result.exception)


def test_gather_raises_when_a_tsv_carries_another_sequencing_groups_calls(shm_tmp_path: Path):
    # rbceq2's label is the only evidence in the file of which sample it describes, so a path
    # keyed under the wrong sequencing group would otherwise be relabelled into place.
    manifest = _manifest_for(shm_tmp_path, ('SG000001',))
    manifest['geno']['SG000001'] = _write_per_sg(
        shm_tmp_path, 'SG000001', 'geno', 'PASS', label='SG000999.converted.vcf'
    )

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'does not name this sequencing group' in str(result.exception)


def test_gather_raises_when_one_output_type_covers_fewer_samples_than_the_others(shm_tmp_path: Path):
    # The QC TSVs and the call TSVs are gathered independently, so nothing but the manifest's
    # sequencing_groups holds them to the same samples. A cohort QC file short of a sample
    # reads as though that sample passed.
    manifest = _manifest_for(shm_tmp_path, ('SG000001', 'SG000002'))
    del manifest['qc']['SG000001']  # drop its QC input, keeping its calls

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'Combined qc TSV covers the wrong sequencing groups' in str(result.exception)
    assert '1 missing (SG000001), 0 unexpected (none)' in str(result.exception)


def test_gather_raises_when_an_output_type_carries_a_sample_the_cohort_does_not_hold(shm_tmp_path: Path):
    # The other direction: a keyed path for a sequencing group outside `sequencing_groups`,
    # as a stale entry from an inactivated sample would be. Its file exists and parses, so
    # only the coverage check can name it.
    manifest = _manifest_for(shm_tmp_path, ('SG000001', 'SG000002'))
    manifest[SEQUENCING_GROUPS_KEY] = ['SG000001']

    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, manifest))

    assert result.exit_code != 0
    assert 'Combined geno TSV covers the wrong sequencing groups' in str(result.exception)
    assert '0 missing (none), 1 unexpected (SG000002)' in str(result.exception)


def test_gather_row_order_is_sorted_not_manifest_order(shm_tmp_path: Path):
    # The manifest maps' insertion order follows the stage's dict-of-targets iteration, which
    # nothing pins across runs; the combined file has to be byte-identical anyway.
    result = CliRunner().invoke(main, _gather_args(shm_tmp_path, _manifest_for(shm_tmp_path, ('SG000002', 'SG000001'))))

    assert result.exit_code == 0, result.output
    combined = (shm_tmp_path / 'combined.geno.tsv').read_text().splitlines()
    assert [line.split('\t')[0] for line in combined[1:]] == ['SG000001', 'SG000002']


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
        SEQUENCING_GROUPS_KEY: ['SG000001', 'SG000002'],
    } | {key: {sg: f'gs://bucket/{sg}.{key}.tsv' for sg in ('SG000001', 'SG000002')} for key in MANIFEST_KEYS}
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


def test_combine_stage_raises_at_graph_construction_when_an_upstream_map_is_short(mock_cohort):
    # The gather job's coverage check would also catch this, but only after every per-SG job
    # has run; the stage holds both sets when it builds the manifest, so a divergence between
    # the upstream skip predicates fails before anything is scheduled.
    sgs = []
    for sg_id in ('SG000001', 'SG000002'):
        sg = MagicMock()
        sg.id = sg_id
        sg.gvcf = f'gs://bucket/{sg_id}.g.vcf.gz'
        sgs.append(sg)
    mock_cohort.get_sequencing_groups.return_value = sgs
    inputs = MagicMock()
    inputs.as_dict_by_target.side_effect = [
        # rbceq2 TSVs for SG000001 only; QC TSVs for both.
        {'SG000001': {key: f'gs://bucket/SG000001.{key}.tsv' for key in MANIFEST_KEYS[:3]}},
        {sg.id: {'qc': f'gs://bucket/{sg.id}.qc.tsv'} for sg in sgs},
    ]

    with pytest.raises(ValueError, match=r'no upstream geno TSV for 1 sequencing group\(s\).*SG000002'):
        pipeline.CombineRbceq2OutputsPerCohort().queue_jobs(mock_cohort, inputs)
