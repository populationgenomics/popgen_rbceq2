#!/usr/bin/env python3
"""Concatenate a cohort's per-SG rbceq2 and blood-group QC TSVs into combined cohort TSVs.

The ``CombineRbceq2OutputsPerCohort`` stage resolves this cohort's per-SG TSV paths from the
workflow dependency graph and writes them to a JSON manifest, passed here as ``--manifest``:
one ``{sequencing group: path}`` map per key ``geno`` / ``pheno_numeric`` /
``pheno_alphanumeric`` / ``qc``, plus the cohort's ``sequencing_groups``. (Repeated per-path
flags would outgrow ARG_MAX around ~3,000 sequencing groups.) For each key we concatenate the
files keeping a single header and write one combined TSV.

Every input type has the same wide one-row-per-sample layout, so one concatenation serves
all four. The QC cells are copied verbatim, keeping the detailed site-level flag in the
cohort file rather than coarsening it to a category.

Rows are keyed by sequencing group, which the manifest supplies rather than the file. RBCeq2
does not read the VCF's sample column: it labels a row with the base name of the VCF it was
given — ``Path.stem``, which strips only the final extension, so the conversion stage's
``<sg>.converted.vcf.gz`` appears as ``<sequencing group>.converted.vcf`` — an intermediate
file name, which carries the sequencing-group ID only because
``FilterAndConvertGvcfsForRbceq2`` happens to name its output that way. A cohort file keyed
on it would not join to Metamist, so column 1 is rewritten to the sequencing group the
manifest paired the file with, and rbceq2's own label is checked against it rather than
trusted.

Only the cohort files are relabelled. The per-SG TSVs keep whatever rbceq2 wrote, including
the QC TSV, which copies the cell from the geno TSV. Each of those is registered against one
sequencing group, so its own Analysis already says whose file it is and nothing reads that
cell: ``analysis_meta.parse_single_row_rbceq2_tsv`` drops column 0 by position.

Cells are matched to columns by system name, not by position. RBCeq2 emits a column only
for a blood group one of whose alleles had all its defining variants in that sample's VCF,
and writes one frame per run (``main.py`` builds each result frame with
``pd.DataFrame.from_dict`` over that run's found blood groups), so a per-sample
TSV can carry that sample's systems and no others. Two sequencing groups could therefore
disagree on the column set, and appending a row under another sample's header would shift
every cell after the first difference. Observed rbceq2 output is a constant 48 systems for
every sample, so this is a latent shift rather than one seen in a run — the fill is defensive
and `concat_tsvs` warns whenever it has to use it.
"""

import json
import logging
from collections import defaultdict
from dataclasses import dataclass

import click
from cpg_utils import to_path

from popgen_rbceq2.logging_setup import setup_logging

logger = logging.getLogger(__name__)

MANIFEST_KEYS = ('geno', 'pheno_numeric', 'pheno_alphanumeric', 'qc')

# The sequencing groups every combined TSV has to carry a row for. Its own manifest key rather
# than the union of the four maps, so a type going short of the cohort is an error rather than
# a shorter file: a cohort QC TSV missing a sample reads as though that sample passed.
SEQUENCING_GROUPS_KEY = 'sequencing_groups'

# A system some other sample in the cohort has a column for, that RBCeq2 did not report for
# this one. Distinct from the QC job's `NA`, which marks a system RBCeq2 did call but that
# has no assessable defining site: here the sample has no value of any kind.
NOT_REPORTED = 'NOT_REPORTED'


def read_manifest(path: str) -> tuple[frozenset[str], dict[str, dict[str, str]]]:
    """Read the manifest the stage wrote.

    Args:
        path: The JSON manifest.

    Returns:
        The cohort's sequencing groups, and `{key: {sequencing group: path}}` for each of
        MANIFEST_KEYS.

    Raises:
        ValueError: The manifest does not carry exactly the expected keys. That means the
            stage and this job disagree about its shape, so it raises rather than combining
            a subset.
    """
    data = json.loads(to_path(path).read_text())
    expected_keys = sorted((SEQUENCING_GROUPS_KEY, *MANIFEST_KEYS))
    if sorted(data) != expected_keys:
        raise ValueError(f'Manifest {path} must have exactly the keys {expected_keys}, got {sorted(data)}')
    return frozenset(data[SEQUENCING_GROUPS_KEY]), {key: data[key] for key in MANIFEST_KEYS}


@dataclass(frozen=True, slots=True)
class SampleRow:
    """One sample's cells, keyed by the blood-group system its column named.

    Keyed rather than positional because the column set is a property of the sample, not of
    the cohort: two sequencing groups can disagree on it, so a cell means nothing without
    the system its column was headed by.
    """

    sample_id: str
    cell_by_system: dict[str, str]


def _parse_sample_tsv(text: str) -> tuple[str, list[SampleRow]]:
    """Split one per-SG TSV into its leading header cell and its sample rows.

    Args:
        text: TSV contents, a header row plus one sample row.

    Returns:
        RBCeq2's leading header cell, and one SampleRow per data row.

    Raises:
        ValueError: The header repeats a system, or a row has a different number of cells
            from that header. Either means the file is malformed rather than just covering
            different systems: a repeated column would silently drop one of its two cells,
            and a ragged row cannot be assigned to columns at all.
    """
    rows = [line.split('\t') for line in text.splitlines() if line.strip()]
    if not rows:
        return '', []
    header, systems = rows[0], rows[0][1:]
    if len(set(systems)) != len(systems):
        raise ValueError(f'Repeated system column in the header: {systems}')
    parsed = []
    for row in rows[1:]:
        if len(row) != len(header):
            raise ValueError(f'Ragged TSV: header has {len(header)} cells but the row for {row[0]!r} has {len(row)}')
        parsed.append(SampleRow(sample_id=row[0], cell_by_system=dict(zip(systems, row[1:], strict=True))))
    return header[0], parsed


def _keyed_row(sg_id: str, text: str) -> tuple[str, SampleRow]:
    """Read one per-SG TSV and key its row on the sequencing group it belongs to.

    Args:
        sg_id: The sequencing group the manifest paired this file with.
        text: TSV contents, a header row plus one sample row.

    Returns:
        RBCeq2's leading header cell, and the sample row keyed on `sg_id`.

    Raises:
        ValueError: The file does not hold exactly one sample row, the row carries no system
            columns at all, or rbceq2's own label for that row does not name `sg_id`.
    """
    header_cell, rows = _parse_sample_tsv(text)
    if len(rows) != 1:
        raise ValueError(f'Expected one sample row for {sg_id}, got {len(rows)}')
    if not rows[0].cell_by_system:
        # The union fill covers samples whose column sets differ, not a sample rbceq2
        # reported nothing for: a row of nothing but NOT_REPORTED is a broken run, not data.
        raise ValueError(f'{sg_id}: the TSV has no system columns, only the ID column')
    # rbceq2 labels the row with the base name of the VCF it read, so this is
    # `<sequencing group>.converted.vcf`. Checking the part before the first `.` catches a
    # file paired with the wrong sequencing group. It pins only the `<sg>.` prefix: a rename
    # of the conversion output that drops the prefix fails loudly rather than relabelling a
    # cohort silently, while extension drift (`.filtered.vcf`, a lost `.gz`) still passes.
    labelled = rows[0].sample_id.split('.', 1)[0]
    if labelled != sg_id:
        raise ValueError(
            f'{sg_id}: rbceq2 labelled its row {rows[0].sample_id!r}, which does not name this '
            f'sequencing group. It labels a row with the base name of the VCF it read, which '
            f'FilterAndConvertGvcfsForRbceq2 writes as <sequencing group>.converted.vcf.gz, so '
            f'this is either the wrong file for this sequencing group or a change to that name.'
        )
    return header_cell, SampleRow(sample_id=sg_id, cell_by_system=rows[0].cell_by_system)


def concat_tsvs(contents: list[tuple[str, str]]) -> str:
    """Concatenate per-sample TSVs onto the union of their columns, keyed by sequencing group.

    Args:
        contents: `(sequencing group id, TSV contents)` pairs, each TSV a header row plus one
            sample row. The id is what column 1 of the combined file carries, replacing
            rbceq2's own label — see the module docstring.

    Returns:
        The combined TSV: the first input's leading header cell, then every system any input
        has a column for, sorted. A sample with no column for a system gets NOT_REPORTED.
        Sorting matches RBCeq2, which alphabetises its columns before writing
        (`IO.record_data.save_df`), so a cohort whose samples all report the same systems
        keeps the per-sample column order. Empty when given no inputs.

    Raises:
        ValueError: An input is malformed or is not the file for the sequencing group it was
            paired with — see `_parse_sample_tsv` and `_keyed_row`.
    """
    id_header = ''
    rows: list[SampleRow] = []
    for sg_id, text in contents:
        header_cell, row = _keyed_row(sg_id, text)
        # RBCeq2's `UUID: <uuid>` label for column 0, from whichever input came first. It names
        # one sample's run rather than the cohort's, and is kept anyway so the combined TSVs
        # have the same shape as the per-sample ones they are built from — analysis_meta reads
        # column 0 by position, so relabelling it would be a schema change for no gain.
        id_header = id_header or header_cell
        rows.append(row)
    if not rows:
        return ''

    systems = sorted({system for row in rows for system in row.cell_by_system})
    # system -> the samples RBCeq2 reported no column for, to report once at the end rather
    # than once per row: one sample carrying an extra system otherwise warns for every other
    # sample in the cohort.
    absent_from: dict[str, list[str]] = defaultdict(list)
    out = ['\t'.join([id_header, *systems])]
    for row in rows:
        cells = []
        for system in systems:
            cell = row.cell_by_system.get(system)
            if cell is None:
                absent_from[system].append(row.sample_id)
                cell = NOT_REPORTED
            cells.append(cell)
        out.append('\t'.join([row.sample_id, *cells]))
    if absent_from:
        counted = '; '.join(
            f'{system} ({len(ids)}/{len(rows)} samples: {_named(ids)})' for system, ids in sorted(absent_from.items())
        )
        logger.warning(
            f'Filled with {NOT_REPORTED} where rbceq2 reported no column — {len(absent_from)} system(s): {counted}'
        )
    return '\n'.join(out) + '\n'


def _named(ids: list[str], limit: int = 10) -> str:
    """Render a sequencing-group list for an error message, capped so a cohort cannot fill it."""
    if not ids:
        return 'none'
    return ', '.join(ids[:limit]) + (' ...' if len(ids) > limit else '')


def _combine(paths: dict[str, str], out_path: str, label: str, expected: frozenset[str]) -> None:
    """Concatenate one type's per-SG TSVs and write the combined cohort file.

    Args:
        paths: `{sequencing group: path}` for this type, from the manifest.
        out_path: Combined TSV to write.
        label: The type, for error messages and logging.
        expected: Sequencing groups this cohort must end up with a row for.

    Raises:
        ValueError: No inputs were given, an input is missing or header-only, or the combined
            file does not cover exactly `expected`.
    """
    if not paths:
        raise ValueError(f'No blood-group {label} inputs provided for {out_path}')
    contents: list[tuple[str, str]] = []
    bad: list[str] = []
    # Sorted so the combined file is byte-identical across runs: the manifest map's insertion
    # order comes from the stage's dict-of-targets iteration, which nothing pins.
    for sg_id, path in sorted(paths.items()):
        tp = to_path(path)
        if not tp.exists():
            bad.append(f'missing: {path}')
            continue
        text = tp.read_text()
        # A valid per-SG TSV is a header row plus the single sample row.
        if len([line for line in text.splitlines() if line.strip()]) < 2:
            bad.append(f'empty/header-only: {path}')
            continue
        contents.append((sg_id, text))
    if bad:
        raise ValueError(f'{len(bad)}/{len(paths)} blood-group {label} inputs are missing or empty:\n' + '\n'.join(bad))

    combined = concat_tsvs(contents)
    # Rows are keyed by sequencing group, so this cannot catch a mislabelled row — `_keyed_row`
    # does that. What it catches is a type covering the wrong set: a stage that produced
    # nothing for a sequencing group the cohort holds, or the four types disagreeing with each
    # other. A cohort QC file short of the calls it describes reads as though those samples
    # passed, which is the failure worth being loud about.
    observed = {line.split('\t', 1)[0] for line in combined.splitlines()[1:]}
    if observed != expected:
        missing, unexpected = sorted(expected - observed), sorted(observed - expected)
        raise ValueError(
            f'Combined {label} TSV covers the wrong sequencing groups: '
            f'{len(missing)} missing ({_named(missing)}), {len(unexpected)} unexpected ({_named(unexpected)}).'
        )
    to_path(out_path).write_text(combined)
    logger.info(f'Wrote combined {label} TSV with {len(observed)} samples -> {out_path}')


@click.command()
@click.option('--output-geno', required=True)
@click.option('--output-pheno-numeric', required=True)
@click.option('--output-pheno-alphanumeric', required=True)
@click.option('--output-qc', required=True)
@click.option(
    '--manifest',
    required=True,
    help=f'JSON manifest: {SEQUENCING_GROUPS_KEY}, plus one {{sequencing group: path}} map per key: '
    + ', '.join(MANIFEST_KEYS),
)
def main(
    output_geno: str,
    output_pheno_numeric: str,
    output_pheno_alphanumeric: str,
    output_qc: str,
    manifest: str,
) -> None:
    setup_logging(force=True)
    # One expected set for all four types, so it also holds them to the same sequencing groups
    # as each other, not just to a plausible-looking file of their own.
    expected, paths = read_manifest(manifest)
    _combine(paths['geno'], output_geno, 'geno', expected)
    _combine(paths['pheno_numeric'], output_pheno_numeric, 'pheno_numeric', expected)
    _combine(paths['pheno_alphanumeric'], output_pheno_alphanumeric, 'pheno_alphanumeric', expected)
    _combine(paths['qc'], output_qc, 'qc', expected)


if __name__ == '__main__':
    main()
