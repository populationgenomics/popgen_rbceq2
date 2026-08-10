"""Declares the pipeline DAG: which stages run, and what each depends on.

Stage dependencies are declared only here. The stage classes in stages/blood_group_genotyping/
and stages/blood_group_qc/ are implementation-only. Read this file top to bottom to see the
pipeline shape.

Everything cpg_flow needs to know about a stage is declared in its ``wire`` call: what it
depends on via ``requires``, and what it records in Metamist via ``analysis_type`` /
``analysis_keys`` / ``update_analysis_meta``.

cpg_flow records each stage's name in ``Analysis.meta`` itself; ``wire`` re-asserts it after
the stage's meta function runs, so a meta function that sets ``stage`` cannot override the
framework's value with a stale hand-typed name.
"""

import cpg_flow.stage

from popgen_rbceq2 import analysis_meta, stage_support
from popgen_rbceq2.stages.blood_group_genotyping import combine, filter_and_convert, genotype
from popgen_rbceq2.stages.blood_group_qc import call_qc

# --- Blood-group genotyping ----------------------------------------------------------
# Starts from the DRAGEN gVCF on each sequencing group, so the branch has no upstream stage:
# FilterAndConvertGvcfsForRbceq2 reads sequencing_group.gvcf directly and needs no `requires`.
# A sequencing group with no gVCF is skipped, not failed — its stages emit no outputs.
#
# The implementations refer to each other by class rather than importing the wired objects
# below: cpg_flow matches stages on ``__name__``, which functools.wraps preserves, so the
# reference resolves and the import stays one-way.

FilterAndConvertGvcfsForRbceq2: cpg_flow.stage.StageDecorator = stage_support.wire(
    filter_and_convert.FilterAndConvertGvcfsForRbceq2,
)
# Per-SG calls Analysis, output = the geno TSV.
GenotypeBloodGroupsWithRbceq2: cpg_flow.stage.StageDecorator = stage_support.wire(
    genotype.GenotypeBloodGroupsWithRbceq2,
    requires=[FilterAndConvertGvcfsForRbceq2],
    analysis_type='blood_group_genotyping',
    analysis_keys=['geno'],
    update_analysis_meta=analysis_meta.blood_group_calls,
)
# Requires the conversion stage for the DP/GQ extract, and the genotyping stage for the columns
# the QC TSV has to match. Registers its own Analysis, separate from the calls above.
FlagBloodGroupCallQc: cpg_flow.stage.StageDecorator = stage_support.wire(
    call_qc.FlagBloodGroupCallQc,
    requires=[FilterAndConvertGvcfsForRbceq2, GenotypeBloodGroupsWithRbceq2],
    analysis_type='blood_group_qc',
    analysis_keys=['qc'],
    update_analysis_meta=analysis_meta.call_qc,
)
# Cohort calls Analysis, registered against the geno TSV alone. The cohort QC TSV gets no
# Analysis of its own — cpg_flow's @stage takes a single analysis_type — and is reachable only
# through the qc_path that analysis_meta.cohort_calls records in the meta.
CombineRbceq2OutputsPerCohort: cpg_flow.stage.StageDecorator = stage_support.wire(
    combine.CombineRbceq2OutputsPerCohort,
    requires=[GenotypeBloodGroupsWithRbceq2, FlagBloodGroupCallQc],
    analysis_type='blood_group_genotyping',
    analysis_keys=['geno'],
    update_analysis_meta=analysis_meta.cohort_calls,
)


# --- what the workflow runs ----------------------------------------------------------
# Handed to cpg_flow's run_workflow by the entry point. cpg_flow walks each entry's
# required_stages and pulls in what they need, so a stage reached only as a dependency does not
# have to be listed — the whole branch arrives via CombineRbceq2OutputsPerCohort.
#
# Adding a stage means adding it here or to the `requires` of something here. A stage that is
# neither is silently absent from the run; test_stage_support asserts that can't happen.
REQUESTED_STAGES: list[cpg_flow.stage.StageDecorator] = [
    CombineRbceq2OutputsPerCohort,
]
