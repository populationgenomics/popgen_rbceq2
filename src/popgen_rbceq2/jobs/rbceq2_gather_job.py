#!/usr/bin/env python3
"""Concatenate a cohort's per-SG rbceq2 and blood-group QC TSVs into combined cohort TSVs.

The ``CombineRbceq2OutputsPerCohort`` stage resolves this cohort's per-SG TSV
paths from the workflow dependency graph and passes them in as repeated ``--geno`` / ``--pheno-numeric`` /
``--pheno-alphanumeric`` / ``--qc`` flags. For each type we concatenate the files keeping a
single header and write one combined TSV.

Every input type has the same wide one-row-per-sample layout, so one concatenation serves
all four. The QC cells are copied verbatim, keeping the detailed site-level flag in the
cohort file rather than coarsening it to a category.
"""

import logging

import click
from cpg_utils import to_path

from popgen_rbceq2.logging_setup import setup_logging

logger = logging.getLogger(__name__)


def concat_tsvs(contents: list[str]) -> str:
    """Concatenate TSV file contents, keeping only the first input's header."""
    out: list[str] = []
    for i, text in enumerate(contents):
        lines = text.splitlines()
        if not lines:
            continue
        out.extend(lines if i == 0 else lines[1:])
    return '\n'.join(out) + '\n' if out else ''


def _combine(paths: tuple[str, ...], out_path: str, label: str) -> None:
    if not paths:
        raise ValueError(f'No blood-group {label} inputs provided for {out_path}')
    contents: list[str] = []
    bad: list[str] = []
    for p in paths:
        tp = to_path(p)
        if not tp.exists():
            bad.append(f'missing: {p}')
            continue
        text = tp.read_text()
        # A valid per-SG TSV is a header row plus the single sample row.
        if len([line for line in text.splitlines() if line.strip()]) < 2:
            bad.append(f'empty/header-only: {p}')
            continue
        contents.append(text)
    if bad:
        raise ValueError(f'{len(bad)}/{len(paths)} blood-group {label} inputs are missing or empty:\n' + '\n'.join(bad))
    to_path(out_path).write_text(concat_tsvs(contents))
    logger.info(f'Wrote combined {label} TSV with {len(contents)} samples -> {out_path}')


@click.command()
@click.option('--output-geno', required=True)
@click.option('--output-pheno-numeric', required=True)
@click.option('--output-pheno-alphanumeric', required=True)
@click.option('--output-qc', required=True)
@click.option('--geno', 'geno_paths', multiple=True, help='Per-SG <sg>.geno.tsv path (repeatable)')
@click.option('--pheno-numeric', 'pheno_numeric_paths', multiple=True, help='Per-SG pheno_numeric TSV (repeatable)')
@click.option(
    '--pheno-alphanumeric',
    'pheno_alphanumeric_paths',
    multiple=True,
    help='Per-SG pheno_alphanumeric TSV (repeatable)',
)
@click.option('--qc', 'qc_paths', multiple=True, help='Per-SG <sg>.qc.tsv path (repeatable)')
def main(
    output_geno: str,
    output_pheno_numeric: str,
    output_pheno_alphanumeric: str,
    output_qc: str,
    geno_paths: tuple[str, ...],
    pheno_numeric_paths: tuple[str, ...],
    pheno_alphanumeric_paths: tuple[str, ...],
    qc_paths: tuple[str, ...],
) -> None:
    setup_logging(force=True)
    _combine(geno_paths, output_geno, 'geno')
    _combine(pheno_numeric_paths, output_pheno_numeric, 'pheno_numeric')
    _combine(pheno_alphanumeric_paths, output_pheno_alphanumeric, 'pheno_alphanumeric')
    _combine(qc_paths, output_qc, 'qc')


if __name__ == '__main__':
    main()
