"""Where each stage writes, and what it calls the things it writes.

Output keys are part of the contract between stages — downstream resolves them by name through
inputs.as_path(key=...) — and the paths are what makes a re-run reuse or rebuild. Neither has a
type error to catch it, so they are asserted literally here.
"""

from pathlib import Path

import pytest

from popgen_rbceq2 import constants
from popgen_rbceq2.stages import pipeline
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
                'output_versions': {'FlagBloodGroupCallQc': 'v2'},
            },
        },
        shm_tmp_path / 'pinned.toml',
    )
    output = outputs_of(pipeline.FlagBloodGroupCallQc(), mock_sequencing_group)
    assert str(output['qc']).endswith(
        f'rbceq2_{constants.RBCEQ2_VERSION.replace(".", "_")}_v2/FlagBloodGroupCallQc/SG000001/SG000001.qc.tsv',
    )
