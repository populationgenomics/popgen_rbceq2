#!/usr/bin/env python3
"""Concatenate a cohort's per-SG rbceq2 and blood-group QC TSVs into combined cohort TSVs.

The ``CombineRbceq2OutputsPerCohort`` stage resolves this cohort's per-SG TSV paths from the
workflow dependency graph and writes them to a JSON manifest, passed here as ``--manifest``:
one list of paths per key ``geno`` / ``pheno_numeric`` / ``pheno_alphanumeric`` / ``qc``.
(Repeated per-path flags would outgrow ARG_MAX around ~3,000 sequencing groups.) For each key
we concatenate the files keeping a single header and write one combined TSV.

Every input type has the same wide one-row-per-sample layout, so one concatenation serves
all four. The QC cells are copied verbatim, keeping the detailed site-level flag in the
cohort file rather than coarsening it to a category.
"""

import json
import logging

import click
from cpg_utils import to_path

from popgen_rbceq2.logging_setup import setup_logging

logger = logging.getLogger(__name__)

MANIFEST_KEYS = ('geno', 'pheno_numeric', 'pheno_alphanumeric', 'qc')


def read_manifest(path: str) -> dict[str, list[str]]:
    """Read the ``{key: [per-SG TSV path, ...]}`` manifest, requiring exactly the four keys.

    A missing or unexpected key means the stage and this job disagree about the manifest
    shape, so it raises rather than combining a subset.
    """
    data = json.loads(to_path(path).read_text())
    if sorted(data) != sorted(MANIFEST_KEYS):
        raise ValueError(f'Manifest {path} must have exactly the keys {sorted(MANIFEST_KEYS)}, got {sorted(data)}')
    return data


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
@click.option('--manifest', required=True, help='JSON manifest: one list of per-SG TSV paths per MANIFEST_KEYS key')
def main(
    output_geno: str,
    output_pheno_numeric: str,
    output_pheno_alphanumeric: str,
    output_qc: str,
    manifest: str,
) -> None:
    setup_logging(force=True)
    paths = read_manifest(manifest)
    _combine(tuple(paths['geno']), output_geno, 'geno')
    _combine(tuple(paths['pheno_numeric']), output_pheno_numeric, 'pheno_numeric')
    _combine(tuple(paths['pheno_alphanumeric']), output_pheno_alphanumeric, 'pheno_alphanumeric')
    _combine(tuple(paths['qc']), output_qc, 'qc')


if __name__ == '__main__':
    main()
