"""Where each stage writes, and what it calls the things it writes.

Output keys are part of the contract between stages — downstream resolves them by name through
inputs.as_path(key=...) — and the paths are what makes a re-run reuse or rebuild. Neither has a
type error to catch it, so they are asserted literally here.
"""

import re
from pathlib import Path

import cpg_utils.config
import pytest

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.stages import pipeline
from popgen_rbceq2.stages.blood_group_genotyping import posthoc_genotype
from tests.helpers import set_config

pytestmark = pytest.mark.fast

# The combined tool+release version segment the prefix helpers derive; tests build the same
# string so a tool bump moves the expectations with it.
VERSION_SEGMENT = f'rbceq2_{constants.RBCEQ2_VERSION.replace(".", "_")}_v1'


def outputs_of(stage, target) -> dict:
    """The stage's expected_outputs, asserted to be the keyed dict these stages return.

    cpg_flow types expected_outputs as a union that also holds a bare path, a string and None,
    so narrowing here is what lets a test index the result by key.
    """
    outputs = stage.expected_outputs(target)
    assert isinstance(outputs, dict)
    return outputs


def analysis_keys_of(stage) -> set[str]:
    """The keys the stage registers an Analysis against, asserted present."""
    assert stage.analysis_keys is not None
    return set(stage.analysis_keys)


def test_filter_and_convert_output_namespacing(mock_sequencing_group):
    # Downstream stages resolve both keys by name via inputs.as_path(key=...), so renaming
    # one breaks the rbceq2 run or the QC flag with no type error. The converted VCF is an
    # intermediate only the next two stages read, so it lands in tmp.
    output = outputs_of(pipeline.FilterAndConvertGvcfsForRbceq2(), mock_sequencing_group)
    prefix = Path('gs://bucket-tmp') / 'popgen_rbceq2' / VERSION_SEGMENT / 'FilterAndConvertGvcfsForRbceq2' / 'SG000001'
    assert str(output['vcf']) == str(prefix / 'SG000001.converted.vcf.gz')
    assert str(output['defining_sites']) == str(prefix / 'SG000001.defining_sites.tsv')


def test_an_exome_run_writes_under_its_capture_design(exome_sequencing_group, mock_cohort, shm_tmp_path):
    # The design sits directly under the release segment, once for the whole run, rather than
    # in each stage that reads it: every exome output from the conversion onward is built from
    # the holes the design leaves, whether or not the stage itself opens the BED. A genome run
    # never reads the key and its tree is unchanged (every other test in this file).
    design_key = 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'
    set_config(
        {
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'exome',
                stage_support.EXOME_DESIGN_KEY: design_key,
            },
        },
        shm_tmp_path / 'exome.toml',
    )
    design = 'exome_probesets_hg38_twist_vcgs_custom_exome_covered_targets_bed'

    sg_output = outputs_of(pipeline.GenotypeBloodGroupsWithRbceq2(), exome_sequencing_group)
    sg_prefix = (
        Path('gs://bucket') / 'popgen_rbceq2' / VERSION_SEGMENT / design / 'GenotypeBloodGroupsWithRbceq2' / 'SG000001'
    )
    assert str(sg_output['geno']) == str(sg_prefix / 'SG000001.geno.tsv')

    cohort_output = outputs_of(pipeline.CombineRbceq2OutputsPerCohort(), mock_cohort)
    cohort_prefix = Path('gs://bucket') / 'popgen_rbceq2' / VERSION_SEGMENT / design / 'CombineRbceq2OutputsPerCohort'
    assert str(cohort_output['geno']) == str(cohort_prefix / 'test-cohort' / 'combined.test-cohort.geno.tsv')


def test_posthoc_output_namespacing(exome_sequencing_group):
    # The conversion stage resolves the gVCF by this key, and derives the .tbi path by
    # appending to it. Both live in tmp: the supplement is an intermediate, and the record of
    # what it contributed is the POSTHOC flag in the QC TSV, not this file.
    output = outputs_of(pipeline.PosthocGenotypeOffTargetSites(), exome_sequencing_group)
    prefix = Path('gs://bucket-tmp') / 'popgen_rbceq2' / VERSION_SEGMENT / 'PosthocGenotypeOffTargetSites' / 'SG000001'
    assert str(output['gvcf']) == str(prefix / 'SG000001.posthoc.g.vcf.gz')
    assert str(output['index']) == str(prefix / 'SG000001.posthoc.g.vcf.gz.tbi')


@pytest.mark.usefixtures('mock_cohort')
def test_posthoc_registers_no_analysis():
    # The supplement is consumed through the cpg-flow graph by path, so a Metamist record
    # would have nothing reading it. Provenance reaches Metamist as a POSTHOC flag on the QC
    # TSV instead, and the caller version in that stage's Analysis meta.
    assert pipeline.PosthocGenotypeOffTargetSites().analysis_type is None


def test_posthoc_is_skipped_for_a_genome(mock_sequencing_group):
    # A genome gVCF is not called against a capture-target BED, so it has no edge for defining
    # sites to fall outside of. Recalling from its CRAM would re-derive what DRAGEN already
    # called, at the cost of a job per sequencing group.
    assert pipeline.PosthocGenotypeOffTargetSites().expected_outputs(mock_sequencing_group) is None


@pytest.mark.parametrize('missing', ['cram', 'gvcf'])
def test_posthoc_is_skipped_when_an_exome_lacks_an_input(exome_sequencing_group, missing):
    # Skipped, not failed, matching how the rest of the pipeline treats a sequencing group it
    # cannot process. Without a CRAM there is nothing to call from; without a gVCF there are no
    # primary calls to supplement, and the conversion stage produces nothing either.
    setattr(exome_sequencing_group, missing, None)
    assert pipeline.PosthocGenotypeOffTargetSites().expected_outputs(exome_sequencing_group) is None


def test_the_conversion_stage_merges_posthoc_calls_for_exactly_the_sequencing_groups_it_runs_for(
    mock_sequencing_group,
    exome_sequencing_group,
):
    # The two stages gate on one predicate. If they could disagree, the conversion stage would
    # ask cpg_flow for an output the post-hoc stage never produced, and the run would die at
    # graph-build time on a subset of a cohort.
    posthoc = pipeline.PosthocGenotypeOffTargetSites()
    assert posthoc_genotype.applies_to(exome_sequencing_group)
    assert posthoc.expected_outputs(exome_sequencing_group) is not None
    assert not posthoc_genotype.applies_to(mock_sequencing_group)
    assert posthoc.expected_outputs(mock_sequencing_group) is None


def test_genotype_output_namespacing(mock_sequencing_group):
    # rbceq2 names its own files <out>_<key>.tsv; these are where they land after write_output,
    # and analysis_meta.blood_group_calls derives the pheno path from the geno one.
    output = outputs_of(pipeline.GenotypeBloodGroupsWithRbceq2(), mock_sequencing_group)
    prefix = Path('gs://bucket') / 'popgen_rbceq2' / VERSION_SEGMENT / 'GenotypeBloodGroupsWithRbceq2' / 'SG000001'
    assert str(output['geno']) == str(prefix / 'SG000001.geno.tsv')
    assert str(output['pheno_numeric']) == str(prefix / 'SG000001.pheno_numeric.tsv')
    assert str(output['pheno_alphanumeric']) == str(prefix / 'SG000001.pheno_alphanumeric.tsv')
    # rbceq2's run log, renamed off its uuid4 name by the job. In the main prefix beside the
    # TSVs, not tmp, so it is still there when someone asks why a call was made.
    assert str(output['log']) == str(prefix / 'SG000001.log.txt')


def test_the_rbceq2_log_is_not_registered_in_metamist(mock_sequencing_group):
    # Team decision, docs/rbceq2_debug_log/SPEC.md: the log is written but not registered.
    # Adding it to analysis_keys would not just add a row — cpg-flow runs every key through the
    # same update_analysis_meta, and blood_group_calls parses the geno TSV, so it would fail.
    genotype_stage = pipeline.GenotypeBloodGroupsWithRbceq2()
    assert 'log' in outputs_of(genotype_stage, mock_sequencing_group)
    assert 'log' not in analysis_keys_of(genotype_stage)


def test_call_qc_output_namespacing(mock_sequencing_group):
    output = outputs_of(pipeline.FlagBloodGroupCallQc(), mock_sequencing_group)
    prefix = Path('gs://bucket') / 'popgen_rbceq2' / VERSION_SEGMENT / 'FlagBloodGroupCallQc' / 'SG000001'
    assert str(output['qc']) == str(prefix / 'SG000001.qc.tsv')


def test_combine_output_namespacing(mock_cohort):
    # The QC TSV is a fourth key beside rbceq2's three; the gather job resolves each by
    # name, and analysis_meta.cohort_calls derives the QC path from the geno filename, so
    # the shared `combined.<cohort>.<key>.tsv` shape is load-bearing.
    output = outputs_of(pipeline.CombineRbceq2OutputsPerCohort(), mock_cohort)
    prefix = Path('gs://bucket') / 'popgen_rbceq2' / VERSION_SEGMENT / 'CombineRbceq2OutputsPerCohort' / 'test-cohort'
    assert str(output['geno']) == str(prefix / 'combined.test-cohort.geno.tsv')
    assert str(output['qc']) == str(prefix / 'combined.test-cohort.qc.tsv')


def test_analysis_keys_name_real_output_keys(mock_cohort, mock_sequencing_group):
    # cpg-flow raises at graph-build time if analysis_keys is not a subset of the
    # expected_outputs keys, so a renamed output key breaks registration, not just a path.
    qc_stage = pipeline.FlagBloodGroupCallQc()
    assert qc_stage.analysis_type == 'blood_group_qc'
    assert analysis_keys_of(qc_stage) <= set(outputs_of(qc_stage, mock_sequencing_group))

    genotype_stage = pipeline.GenotypeBloodGroupsWithRbceq2()
    assert genotype_stage.analysis_type == 'blood_group_genotyping'
    assert analysis_keys_of(genotype_stage) <= set(outputs_of(genotype_stage, mock_sequencing_group))

    combine_stage = pipeline.CombineRbceq2OutputsPerCohort()
    assert combine_stage.analysis_type == 'blood_group_genotyping'
    assert analysis_keys_of(combine_stage) <= set(outputs_of(combine_stage, mock_cohort))


def test_sg_stages_emit_nothing_without_a_gvcf(mock_sequencing_group):
    # The whole branch starts from the gVCF, so a sequencing group without one is skipped
    # rather than failed.
    mock_sequencing_group.gvcf = None
    assert pipeline.FilterAndConvertGvcfsForRbceq2().expected_outputs(mock_sequencing_group) is None
    assert pipeline.GenotypeBloodGroupsWithRbceq2().expected_outputs(mock_sequencing_group) is None
    assert pipeline.FlagBloodGroupCallQc().expected_outputs(mock_sequencing_group) is None


def test_output_version_can_be_pinned_per_stage(mock_sequencing_group, shm_tmp_path):
    # Bumping workflow.version moves every stage's outputs; output_versions pins one, so a
    # single stage can be re-run into a fresh tree without orphaning the rest.
    set_config(
        {
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'genome',
                'output_versions': {'FlagBloodGroupCallQc': 'v2'},
            },
        },
        shm_tmp_path / 'pinned.toml',
    )
    output = outputs_of(pipeline.FlagBloodGroupCallQc(), mock_sequencing_group)
    assert str(output['qc']).endswith(
        f'rbceq2_{constants.RBCEQ2_VERSION.replace(".", "_")}_v2/FlagBloodGroupCallQc/SG000001/SG000001.qc.tsv',
    )


def _pinned(release: object) -> dict:
    """A genome config whose only oddity is what FlagBloodGroupCallQc is pinned to."""
    return {
        'name': 'popgen_rbceq2',
        'version': 'v1',
        'sequencing_type': 'genome',
        'output_versions': {'FlagBloodGroupCallQc': release},
    }


@pytest.mark.parametrize(
    ('workflow_config', 'expected_in_message'),
    [
        # No release at all. Defaulting one would name an output tree and assert itself on
        # every Metamist row under a version the run never declared.
        ({'name': 'popgen_rbceq2', 'sequencing_type': 'genome'}, 'workflow.version'),
        # A blank workflow.version, which merging a config over the defaults can produce.
        ({'name': 'popgen_rbceq2', 'version': '', 'sequencing_type': 'genome'}, 'workflow.version'),
        # A blank pin. Falling back to workflow.version here would write the pinned stage into
        # the shared tree and record the shared release — the opposite of what pinning is for.
        (_pinned('  '), 'workflow.output_versions.FlagBloodGroupCallQc'),
        # false is the one to watch: bool subclasses int, so a str/int check that does not
        # exclude it first would release the stage under 'False'.
        (_pinned(False), 'workflow.output_versions.FlagBloodGroupCallQc'),
        (_pinned(1.5), 'workflow.output_versions.FlagBloodGroupCallQc'),
        (_pinned(['v5']), 'workflow.output_versions.FlagBloodGroupCallQc'),
        # Anything that would put a separator or a space in the tree, rather than being
        # silently rewritten into a segment the meta then disagrees with.
        (_pinned('v5/rerun'), 'workflow.output_versions.FlagBloodGroupCallQc'),
        (_pinned('v5 rerun'), 'workflow.output_versions.FlagBloodGroupCallQc'),
        # output_versions written as a scalar instead of a table. cpg_utils walks the key list
        # with `in`, a substring test on a str, so this reads as "no pin" for most stage names
        # and raises a bare TypeError for one whose name is a substring of the value.
        (
            {'name': 'popgen_rbceq2', 'version': 'v1', 'sequencing_type': 'genome', 'output_versions': 'v5'},
            'workflow.output_versions',
        ),
        (
            {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'genome',
                'output_versions': 'FlagBloodGroupCallQc',
            },
            'workflow.output_versions',
        ),
    ],
    ids=[
        'no_version',
        'blank_version',
        'blank_pin',
        'false_pin',
        'float_pin',
        'list_pin',
        'pin_with_separator',
        'pin_with_space',
        'scalar_output_versions',
        'scalar_output_versions_matching_a_stage_name',
    ],
)
def test_a_release_is_never_invented(workflow_config, expected_in_message, mock_sequencing_group, shm_tmp_path):
    set_config({'workflow': workflow_config}, shm_tmp_path / 'bad_release.toml')

    with pytest.raises(cpg_utils.config.ConfigError, match=re.escape(expected_in_message)):
        outputs_of(pipeline.FlagBloodGroupCallQc(), mock_sequencing_group)


def test_a_pin_is_used_as_written_not_for_being_truthy(mock_sequencing_group, shm_tmp_path):
    # A TOML integer pin is stringified, so meta.stage_version holds the same type whatever the
    # config author typed, and the path segment it produces is the one the meta names.
    set_config(
        {
            'references': {'genome_build': 'GRCh38'},
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'genome',
                'output_versions': {'FlagBloodGroupCallQc': 2},
            },
        },
        shm_tmp_path / 'int_pin.toml',
    )
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text('UUID: abc123\tJK\nSG000001\tPASS\n')
    stage = pipeline.FlagBloodGroupCallQc()
    assert stage.update_analysis_meta is not None

    meta = stage.update_analysis_meta(str(qc_path))

    assert meta['stage_version'] == '2'
    assert f'_{meta["stage_version"]}/' in str(outputs_of(stage, mock_sequencing_group)['qc'])


def test_the_analysis_meta_records_the_release_its_outputs_went_to(mock_sequencing_group, shm_tmp_path):
    # Metamist keeps a row per release and retires none, so stage_version is how a reader tells
    # two releases of one cohort apart. It is read from the same pin the path is, hence the tie
    # to the path here: a stage with an output_versions pin must not record the unpinned
    # release.
    set_config(
        {
            'references': {'genome_build': 'GRCh38'},
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'genome',
                'output_versions': {'FlagBloodGroupCallQc': 'v2'},
            },
        },
        shm_tmp_path / 'meta_release.toml',
    )
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text('UUID: abc123\tJK\nSG000001\tPASS\n')
    stage = pipeline.FlagBloodGroupCallQc()
    assert stage.update_analysis_meta is not None

    meta = stage.update_analysis_meta(str(qc_path))

    assert meta['stage_version'] == 'v2'
    assert f'_{meta["stage_version"]}/' in str(outputs_of(stage, mock_sequencing_group)['qc'])


def test_the_analysis_meta_records_the_design_an_exome_was_called_against(
    exome_sequencing_group,
    mock_sequencing_group,
    shm_tmp_path,
):
    # The design is the third axis of an exome output tree, and the only one a reader could
    # not otherwise recover from the meta: two rows for one sequencing group called under two
    # designs are identical on every other key. Recorded as the references-repo key, not the
    # sanitised path segment, so it can be looked up; the tie to the path is asserted through
    # design_segment for the same reason the release tie is.
    design_key = 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'
    set_config(
        {
            'references': {'genome_build': 'GRCh38'},
            'workflow': {
                'name': 'popgen_rbceq2',
                'version': 'v1',
                'sequencing_type': 'exome',
                stage_support.EXOME_DESIGN_KEY: design_key,
            },
        },
        shm_tmp_path / 'exome_meta.toml',
    )
    qc_path = shm_tmp_path / 'SG000001.qc.tsv'
    qc_path.write_text('UUID: abc123\tJK\nSG000001\tPASS\n')
    stage = pipeline.FlagBloodGroupCallQc()
    assert stage.update_analysis_meta is not None

    meta = stage.update_analysis_meta(str(qc_path))

    assert meta['exome_design'] == design_key
    path = str(outputs_of(stage, exome_sequencing_group)['qc'])
    assert f'/{stage_support.design_segment(meta["exome_design"])}/' in path

    # A genome run carries the key holding None rather than dropping it, so the meta shape
    # does not depend on sequencing type and a query can filter on it either way.
    set_config(
        {
            'references': {'genome_build': 'GRCh38'},
            'workflow': {'name': 'popgen_rbceq2', 'version': 'v1', 'sequencing_type': 'genome'},
        },
        shm_tmp_path / 'genome_meta.toml',
    )
    genome_meta = stage.update_analysis_meta(str(qc_path))
    assert genome_meta['exome_design'] is None
    assert 'twist' not in str(outputs_of(stage, mock_sequencing_group)['qc'])
