#!/usr/bin/env python3
"""Generate the committed blood-group site resources from RBCeq2's allele DB.

Run by hand when the pinned RBCeq2 image changes, against the `db.tsv` from the *same*
version the image installs, and commit the results:

    gen_bg_resources.py <db.tsv> GRCh38 src/popgen_rbceq2/resources

Writes three files, all from one parse (see `bg_db`):

    bg_regions.<genome>.bed         merged ±flank intervals; restricts the GVCF conversion
    bg_defining_sites.<genome>.bed  the exact defining coordinates the QC pass extracts
    bg_site_systems.<genome>.tsv    chrom/pos/ref/alt/kind/system, for per-system flagging

The regions BED must stay a strict superset of every coordinate RBCeq2 queries or
blood-group calls go silently wrong, which is why it and the defining-sites BED are
generated together and never edited by hand.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

from popgen_rbceq2.scripts import bg_db

logger = logging.getLogger(__name__)

DEFAULT_FLANK = 500_000


def write_regions_bed(rows: list[dict[str, str]], genome: str, flank: int, out_path: Path) -> int:
    """Write the merged ±flank regions BED.

    Args:
        rows: Parsed db rows.
        genome: Coordinate column to read, `GRCh37` or `GRCh38`.
        flank: Bases to extend each defining position by on either side.
        out_path: File to write.

    Returns:
        The number of intervals written.
    """
    merged = bg_db.build_intervals(rows, genome, flank)
    lines = [
        f'{chrom}\t{start}\t{end}' for chrom in sorted(merged, key=bg_db.chrom_key) for start, end in merged[chrom]
    ]
    out_path.write_text('\n'.join(lines) + '\n')
    return len(lines)


def write_defining_sites_bed(rows: list[dict[str, str]], genome: str, out_path: Path) -> int:
    """Write the defining-site coordinates as a 0-based half-open BED.

    Args:
        rows: Parsed db rows.
        genome: Coordinate column to read, `GRCh37` or `GRCh38`.
        out_path: File to write.

    Returns:
        The number of sites written.
    """
    sites = bg_db.defining_site_positions(rows, genome)
    out_path.write_text('\n'.join(f'{chrom}\t{pos - 1}\t{pos}' for chrom, pos in sites) + '\n')
    return len(sites)


def write_site_system_map(rows: list[dict[str, str]], genome: str, out_path: Path) -> int:
    """Write the site -> blood-group-system map, with a header row.

    Args:
        rows: Parsed db rows.
        genome: Coordinate column to read, `GRCh37` or `GRCh38`.
        out_path: File to write.

    Returns:
        The number of data rows written.
    """
    sites = bg_db.site_system_map(rows, genome)
    with out_path.open('w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(bg_db.SITE_SYSTEM_COLUMNS)
        for site in sites:
            writer.writerow([site.chrom, site.pos, site.ref, site.alt, site.kind, site.system])
    return len(sites)


def main(argv: list[str] | None = None) -> None:
    """Write all three blood-group site resources into the given directory.

    Args:
        argv: Command-line arguments; defaults to `sys.argv[1:]`.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('db_tsv', type=Path, help="RBCeq2's db.tsv, from the pinned image's version")
    parser.add_argument('genome', choices=['GRCh37', 'GRCh38'], help='db coordinate column to read')
    parser.add_argument('out_dir', type=Path, help='resources directory to write into')
    parser.add_argument('--flank', type=int, default=DEFAULT_FLANK, help='regions BED flank (bp)')
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s', stream=sys.stderr)
    with args.db_tsv.open(newline='') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    logger.info(f'Read {len(rows)} allele rows from {args.db_tsv}')

    n_regions = write_regions_bed(rows, args.genome, args.flank, args.out_dir / f'bg_regions.{args.genome}.bed')
    n_sites = write_defining_sites_bed(rows, args.genome, args.out_dir / f'bg_defining_sites.{args.genome}.bed')
    n_map = write_site_system_map(rows, args.genome, args.out_dir / f'bg_site_systems.{args.genome}.tsv')
    assessed = bg_db.site_system_map(rows, args.genome)
    systems = {s.system for s in assessed}
    logger.info(
        f'Wrote {n_regions} merged regions, {n_sites} defining sites, '
        f'{n_map} site-system rows across {len(systems)} blood-group systems'
    )

    # Name the systems the resources cannot speak for, so a system reported as not assessed
    # traces back to a db row rather than looking like a gap in the pipeline.
    all_sites = list(bg_db.iter_defining_sites(rows, args.genome))
    sv_sites = {s for s in all_sites if s.kind == 'sv'}
    if sv_sites:
        unassessable = sorted({s.system for s in sv_sites} - systems)
        named = f'; no assessable site remains for {", ".join(unassessable)}' if unassessable else ''
        logger.info(f'Excluded {len(sv_sites)} structural-variant site(s) that one base of DP/GQ cannot assess{named}')


if __name__ == '__main__':
    main()
