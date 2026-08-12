"""Stage implementations, grouped by the part of the pipeline they belong to.

NO STAGES ARE DEFINED IN THIS FILE, and none are re-exported from it. Each stage is one
undecorated class in its own module under two subpackages

- `blood_group_genotyping/`
- `blood_group_qc/`

and the DAG connecting them is declared in `stages/pipeline.py` via `stage_support.wire`.

Import the module and qualify the class (`genotype.GenotypeBloodGroupsWithRbceq2`) rather than
re-exporting it. `wire` does not modify the class it is given: cpg_flow's `@stage` returns a new
function wrapping it, carrying the same `__name__`. A re-export would therefore put two
different objects of the same name in reach, only one of which can be instantiated. Qualifying
also reads as what it is — the implementation, not the DAG node.

Stage naming follows the CPG convention: PascalCase, verb + subject + preposition + object,
with only the first character of an initialism capitalised (`GenotypeBloodGroupsWithRbceq2`,
`FilterAndConvertGvcfsForRbceq2`). See
https://cpg-populationanalysis.atlassian.net/wiki/spaces/ST/pages/185597962/Pipeline+Naming+Convention+Specification

A stage's class name also derives its `[workflow.<section>]` config section, so renaming one
means renaming that section — see stage_support.config_section.
test
"""
