"""Blood-group call QC: one stage implementation per module.

A subpackage separate from `blood_group_genotyping/` because the boundary is whether we own
the logic or just call rbceq2, not position in the pipeline: FlagBloodGroupCallQc sits between
genotype and combine in the DAG, but every line of it is ours, so it lives on its own.

As in `blood_group_genotyping/`, the class is undecorated and knows nothing about its place in
the graph; the DAG is declared in stages/pipeline.py.
"""
