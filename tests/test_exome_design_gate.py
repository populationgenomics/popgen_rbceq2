"""The capture design that bounds the fill: how it is demanded, resolved, and enforced.

Post-hoc calls fill only defining sites the capture design never targeted. Which design that
is cannot be inferred from the gVCF, so it is configured, as the design's key in the references
repo. The pipeline never opens the design: the defining sites outside each known design are
subtracted once, by `scripts/gen_off_design_sites.py`, and committed under `resources/`, so a
run resolves a file for its key. The committed files themselves are covered in
test_off_design_resources.

Five things have to hold and none has a type error to catch it. An exome run that does not name
a design, or names one with no committed subtraction, must fail while the graph is being built,
not on the first sequencing group, because the whole cohort is wasted either way and only one of
those says why. A genome run must never read the key at all. Repointing the key must not reuse
any output built from the previous design, which is every exome output from the conversion
onward. The conversion job must still refuse to run when DRAGEN's own records contradict the
named design. And the QC job must be handed the same off-design sites, so it can disregard a
post-hoc record at a site the merge was not allowed to fill.

No batch runs here: each job's command is read off the mock the stage builds it on.
"""

from pathlib import Path
from unittest.mock import MagicMock

import cpg_utils.config
import pytest

from popgen_rbceq2 import constants, off_design, stage_support
from popgen_rbceq2.stage_support import DESIGN_CONFIG_PATH, EXOME_DESIGN_KEY
from popgen_rbceq2.stages import pipeline
from tests.helpers import set_config
from tests.test_output_namespacing import outputs_of

pytestmark = pytest.mark.fast

TWIST_KEY = 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'
CREV2_KEY = 'exome_probesets_hg38/agilent_sureselect_clinical_research_exome_v2_covered_by_probes_bed'
# Names nothing committed under resources/.
UNKNOWN_KEY = 'exome_probesets_hg38/no_such_design_bed'


def _config(shm_tmp_path: Path, sequencing_type: str, design_key: str | None) -> None:
    """Write a config for one sequencing type, with or without the design key."""
    design: dict[str, str] = {} if design_key is None else {EXOME_DESIGN_KEY: design_key}
    set_config(
        {
            'references': {'genome_build': 'GRCh38'},
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': sequencing_type,
                'driver_image': 'stub-driver:1.0',
                **design,
            },
        },
        shm_tmp_path / f'{sequencing_type}-{(design_key or "none").rsplit("/", 1)[-1]}.toml',
    )


def _queue_conversion(mocker, sequencing_group) -> MagicMock:
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


def _queue_qc(mocker, sequencing_group) -> MagicMock:
    """Queue the QC stage's job on a mock batch and return that batch."""
    batch = MagicMock()
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=batch)
    inputs = MagicMock()
    inputs.as_str.return_value = 'gs://bucket/SG000001.some-input'
    # Localised files render under a fixed name so the command can be asserted on.
    batch.read_input.return_value = '/io/localised.bed'
    pipeline.FlagBloodGroupCallQc().queue_jobs(sequencing_group, inputs)
    return batch


def _command(batch: MagicMock) -> str:
    return batch.new_bash_job.return_value.command.call_args.args[0]


def _localised(batch: MagicMock) -> list[str]:
    return [call.args[0] for call in batch.read_input.call_args_list]


# --- demanding the design ---


def test_an_exome_run_without_a_design_fails_at_the_first_output_path(
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # Failing here is the point. The alternative to a required key is a default, and every
    # default is wrong for some cohort: filling holes inside the capture would answer, with a
    # second caller, a question about that sample's DRAGEN run. The design is in every exome
    # output path, so a run that omits it cannot even name where its first stage would write,
    # which is the earliest a graph build can fail. It holds for an exome run whose sequencing
    # groups happen to have no CRAM, where nothing would ever open the off-design sites.
    _config(shm_tmp_path, 'exome', design_key=None)
    exome_sequencing_group.cram = None

    with pytest.raises(cpg_utils.config.ConfigError) as raised:
        pipeline.FilterAndConvertGvcfsForRbceq2().expected_outputs(exome_sequencing_group)

    assert DESIGN_CONFIG_PATH in str(raised.value)


def test_an_exome_run_naming_a_design_with_no_committed_subtraction_fails(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # A typo, or a design nobody has run the generator for, resolves to no file. Failing at
    # graph build with the generator named is what stops it becoming a job that localises
    # nothing and fills nothing.
    _config(shm_tmp_path, 'exome', design_key=UNKNOWN_KEY)

    with pytest.raises(FileNotFoundError) as raised:
        _queue_conversion(mocker, exome_sequencing_group)

    assert f'{DESIGN_CONFIG_PATH} = {UNKNOWN_KEY}' in str(raised.value)
    assert off_design.GENERATOR in str(raised.value)


def test_the_cohort_says_which_exome_sequencing_groups_the_recall_is_skipped_for(
    mocker,
    exome_sequencing_group,
    mock_cohort,
    shm_tmp_path: Path,
    caplog,
):
    # A sequencing group without a CRAM is skipped, not failed, and its QC then reads NOCOV at
    # every off-design site exactly as before the recall existed. Correct per sample, but with
    # nothing said a cohort whose CRAMs were never registered would run green with the recall
    # silently off. So the cohort stage names them, once, at graph build.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)
    no_cram = MagicMock()
    no_cram.dataset = mock_cohort.dataset
    no_cram.id = 'SG000002'
    no_cram.gvcf = 'gs://bucket/SG000002.g.vcf.gz'
    no_cram.cram = None
    no_cram.sequencing_type = 'exome'
    sgs = [exome_sequencing_group, no_cram]
    mock_cohort.get_sequencing_groups.return_value = sgs
    # The combine stage writes its manifest where tmp outputs go; a local stand-in for gs://.
    mock_cohort.dataset.prefix.side_effect = lambda category=None: (
        shm_tmp_path if category == 'tmp' else Path('gs://bucket')
    )
    inputs = MagicMock()
    inputs.as_dict_by_target.side_effect = [
        {sg.id: {key: f'gs://bucket/{sg.id}.{key}.tsv' for key in constants.RBCEQ2_TSV_KEYS} for sg in sgs},
        {sg.id: {'qc': f'gs://bucket/{sg.id}.qc.tsv'} for sg in sgs},
    ]
    mocker.patch('cpg_utils.hail_batch.get_batch', return_value=MagicMock())
    stage = pipeline.CombineRbceq2OutputsPerCohort()
    # gs:// needs no directories; the local stand-in does.
    Path(str(stage_support.get_output_prefix(mock_cohort, stage.name, category='tmp'))).mkdir(parents=True)

    with caplog.at_level('INFO'):
        stage.queue_jobs(mock_cohort, inputs)

    assert 'post-hoc calling applies to 1 of 2 sequencing group(s) in cohort test-cohort' in caplog.text
    assert 'skipped for 1 exome sequencing group(s)' in caplog.text
    assert 'SG000002' in caplog.text
    assert 'SG000001' not in caplog.text.split('skipped')[1]


def test_repointing_the_design_key_moves_every_output_of_an_exome_run(
    exome_sequencing_group,
    mock_cohort,
    shm_tmp_path: Path,
):
    # cpg_flow reuses a stage whose expected outputs exist and asks nothing about how they were
    # made, and the release segment cannot see a config change. Every exome output depends on
    # the design: the merge fills the holes it names, so the converted VCF, the genotypes, the
    # QC flags and the cohort tables all do too. So the design is a segment of the whole exome
    # tree, and every stage's output moves with it.
    stages_and_targets = [
        (pipeline.PosthocGenotypeOffTargetSites(), exome_sequencing_group),
        (pipeline.FilterAndConvertGvcfsForRbceq2(), exome_sequencing_group),
        (pipeline.GenotypeBloodGroupsWithRbceq2(), exome_sequencing_group),
        (pipeline.FlagBloodGroupCallQc(), exome_sequencing_group),
        (pipeline.CombineRbceq2OutputsPerCohort(), mock_cohort),
    ]

    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)
    twist = [outputs_of(stage, target) for stage, target in stages_and_targets]
    _config(shm_tmp_path, 'exome', design_key=CREV2_KEY)
    crev2 = [outputs_of(stage, target) for stage, target in stages_and_targets]

    for (stage, _), before, after in zip(stages_and_targets, twist, crev2, strict=True):
        for key in before:
            assert str(before[key]) != str(after[key]), f'{stage.name}[{key}] did not move with the design'


# --- enforcing it, per sequencing group ---


def test_an_exome_merge_gates_the_fill_on_the_committed_off_design_sites(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # The fill is the intersection: holes are narrowed to off-design sites, never the reverse.
    # The conversion job localises the committed subtraction for its design, never the vendor
    # BED, and does no subtraction of its own.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    batch = _queue_conversion(mocker, exome_sequencing_group)

    localised = _localised(batch)
    assert off_design.resource_path(TWIST_KEY) in localised
    assert not any(name.endswith('Twist_VCGS_Exome_Covered_Targets_hg38.bed') for name in localised)
    command = _command(batch)
    assert '> uncovered.bed' in command
    assert 'bedtools' not in command
    assert '> off_design_sites.bed' not in command


def test_an_exome_merge_fails_when_dragen_called_outside_the_configured_design(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # DRAGEN emits reference blocks over exactly the target BED with no padding, so a block
    # reaching an off-design defining site means the configured file is not the one the gVCF
    # was called against. Measured on the validation cohorts: a target-regions file where the
    # gVCF used the probe footprint leaves 309 such sites and silently gates off three quarters
    # of the recoveries. The gate itself runs under real bcftools in test_posthoc_trespass; this
    # checks the stage puts it in the job.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    command = _command(_queue_conversion(mocker, exome_sequencing_group))

    assert 'off_design_in_blocks.bed' in command
    assert DESIGN_CONFIG_PATH in command
    assert 'exit 1' in command


# --- handing the same sites to the QC ---


def test_an_exome_qc_is_handed_the_fillable_sites_the_merge_wrote(
    mocker,
    exome_sequencing_group,
    shm_tmp_path: Path,
):
    # A post-hoc reference block is kept whole, so one selected for an off-design hole can be
    # the only record at an in-design hole beside it. The QC keeps that hole NOCOV only if it
    # knows which sites were fillable, and the only right answer is the list the merge used.
    #
    # That is a stage input, not the committed BED: the merge's single-copy chrX gate reads
    # this sample's own genotypes, so the committed file can name a site this sample never
    # allowed to be filled, and the QC would then trust a post-hoc record there.
    _config(shm_tmp_path, 'exome', design_key=TWIST_KEY)

    batch = _queue_qc(mocker, exome_sequencing_group)

    assert off_design.resource_path(TWIST_KEY) not in _localised(batch)
    assert '--fillable-sites gs://bucket/SG000001.some-input' in _command(batch)


def test_a_genome_qc_is_handed_no_fillable_set(
    mocker,
    mock_sequencing_group,
    shm_tmp_path: Path,
):
    # Nothing was merged, so there is no set to hand over, and the job refuses a post-hoc record
    # rather than assuming one. Reading the design key here would also break every genome run.
    _config(shm_tmp_path, 'genome', design_key=None)

    assert '--fillable-sites' not in _command(_queue_qc(mocker, mock_sequencing_group))


def test_a_genome_conversion_never_asks_for_a_design(
    mocker,
    mock_sequencing_group,
    shm_tmp_path: Path,
):
    # A genome gVCF has no capture edge, gets no post-hoc input, and must not read an
    # exome-only key. Its merge is a plain rename.
    _config(shm_tmp_path, 'genome', design_key=None)

    command = _command(_queue_conversion(mocker, mock_sequencing_group))

    assert 'off_design' not in command
    assert 'mv dragen.vcf.gz merged.vcf.gz' in command
