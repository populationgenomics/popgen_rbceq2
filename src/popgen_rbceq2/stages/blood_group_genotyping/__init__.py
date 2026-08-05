"""Blood-group genotyping with rbceq2: one stage implementation per module.

    filter_and_convert -> genotype -> call_qc -> combine

The classes here are undecorated and know nothing about each other's place in the graph; the
DAG is declared in stages/pipeline.py. Where an implementation has to name another stage — to
resolve an input by stage — it refers to the implementation class, not the wired object, so the
import stays one-way. cpg_flow matches stages on `__name__`, which is preserved by the
functools.wraps in its `@stage`, so the reference still resolves.
"""
