#!/usr/bin/env python3
"""Generate a synthetic DRAGEN-shaped gVCF for one sequencing group, for local checks.

This repo ships no test data and cannot: a real DRAGEN gVCF is 12-15Gb and names a real
individual. The test suite therefore builds every fixture inline, which is right for unit
tests and useless for the one question that needs a whole sample: **what do the rbceq2 TSVs
say when the sex chromosomes are encoded the way DRAGEN encodes them?**

That question is why this exists. DRAGEN calls non-PAR chrX and chrY at their real ploidy in
a male sample, writing a one-token `GT`, and three blood-group loci sit there: XK, GATA1 and
ATP11C. `--diploidise` writes the same sample the way `bcftools +fixploidy` used to leave it,
and `--mixed` writes the half-rewritten state rbceq2 2.4.4 refuses. Running all three through
the pipeline and diffing the TSVs is what the haploid-encoding section of the README reports.

Coordinates are read from the committed `bg_site_systems.<genome>.tsv` rather than written
here, so a fixture cannot describe a site the pipeline no longer ships.

    gen_synthetic_gvcf.py GRCh38 out.g.vcf
    gen_synthetic_gvcf.py GRCh38 out.diploidised.g.vcf --diploidise
    gen_synthetic_gvcf.py GRCh38 out.mixed.g.vcf --mixed

The output is a plain VCF. `bgzip` and index it before handing it to the conversion stage.
It is NOT a substitute for running a real male genome before trusting a ploidy change; it
fixes the encoding under test rather than discovering what DRAGEN actually emitted.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

from popgen_rbceq2 import stage_support
from popgen_rbceq2.scripts import bg_db

logger = logging.getLogger(__name__)

SAMPLE = 'CPGSYNTH1'

# GRCh38 non-PAR chrX. A male sample is haploid here and diploid either side of it, which is
# the distinction the fixture exists to carry. PAR1 ends at 2,781,479 and PAR2 begins at
# 155,701,383; XG and CD99 sit inside PAR1, XK, GATA1 and ATP11C between them.
NON_PAR_X = (2_781_480, 155_701_382)

# Which defining site to call, and with how many copies of the ALT, per blood-group system.
# The index selects from that system's sites in coordinate order, so a system whose site list
# grows keeps describing the same allele. Two XK entries so --mixed has a second one to
# disagree with; the rest are one apiece, enough to produce a multi-system report.
CALLS = (
    # (system, index into that system's sites, ALT copies)
    ('XK', 5, 1),  # non-PAR: hemizygous null in a male
    ('GATA1', 0, 1),  # non-PAR: hemizygous
    ('XG', 0, 1),  # PAR1: genuinely diploid, must stay two-token
    ('FY', 3, 2),  # autosomal: homozygous
    ('FY', 9, 1),  # autosomal: heterozygous
)

# The XK site --mixed writes diploid while the rest of non-PAR chrX stays haploid.
MIXED_SYSTEM, MIXED_INDEX = 'XK', 7

# Length of the reference block placed before each called site. rbceq2 never sees these -- the
# <NON_REF> drop removes them -- but the QC extract reads them, so the conversion stage has to
# survive a file that has both.
BLOCK_LEN = 20

HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowDepth,Description="Site filtered because of low read depth">
##FILTER=<ID=DRAGENSnpHardQUAL,Description="Set if true:QUAL < 10.41">
##FILTER=<ID=DRAGENIndelHardQUAL,Description="Set if true:QUAL < 7.83">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the block">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chrX,length=156040895>
"""


def is_non_par_x(chrom: str, pos: int) -> bool:
    """Whether a coordinate is one DRAGEN would call haploid in a male sample.

    Args:
        chrom: Contig name, `chr`-prefixed.
        pos: 1-based position.

    Returns:
        True for non-PAR chrX, which is where the haploid encoding appears.
    """
    return chrom == 'chrX' and NON_PAR_X[0] <= pos <= NON_PAR_X[1]


def render_gt(chrom: str, pos: int, copies: int, *, diploidise: bool) -> str:
    """Build a GT string at the ploidy of its coordinate.

    Args:
        chrom: Contig name.
        pos: 1-based position.
        copies: How many copies of the ALT allele the sample carries, 0 to 2.
        diploidise: Write a haploid site the way `bcftools +fixploidy` left it -- allele
            duplication with a phased separator -- instead of as DRAGEN called it.

    Returns:
        A one-token GT on non-PAR chrX unless `diploidise`, and a two-token GT elsewhere.
    """
    if is_non_par_x(chrom, pos):
        token = '1' if copies >= 1 else '0'
        return f'{token}|{token}' if diploidise else token
    return ('0/0', '0/1', '1/1')[copies]


def sites_by_system(genome: str) -> dict[str, list[bg_db.DefiningSite]]:
    """Read the committed site map and group its variant sites by blood-group system.

    Args:
        genome: Genome build, selecting which committed resource to read.

    Returns:
        System name to its `var` sites in coordinate order. `ref` sites are excluded: they
        define an allele by the reference base being present, so there is no ALT to call.
    """
    path = stage_support.blood_group_resource(f'bg_site_systems.{genome}.tsv')
    grouped: dict[str, list[bg_db.DefiningSite]] = {}
    with Path(path).open(newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['kind'] != 'var':
                continue
            site = bg_db.DefiningSite(
                chrom=row['chrom'],
                pos=int(row['pos']),
                ref=row['ref'],
                alt=row['alt'],
                kind='var',
                system=row['system'],
            )
            grouped.setdefault(site.system, []).append(site)
    return grouped


def build(genome: str, *, diploidise: bool = False, mixed: bool = False) -> str:
    """Render the whole gVCF.

    Args:
        genome: Genome build, selecting the committed site map to read coordinates from.
        diploidise: Write every non-PAR chrX call pre-expanded, as `+fixploidy` left it.
        mixed: Write one extra XK site diploid while the rest of non-PAR chrX stays haploid,
            reproducing a half-applied ploidy rewrite.

    Returns:
        The complete VCF text, records sorted by contig then position.

    Raises:
        KeyError: A system named in `CALLS` has no variant site in the committed map, which
            means the resources moved and this script's indices need revisiting.
    """
    grouped = sites_by_system(genome)
    calls = [*CALLS, (MIXED_SYSTEM, MIXED_INDEX, 1)] if mixed else list(CALLS)

    rows: list[tuple[str, int, str]] = []
    contigs: set[str] = set()
    for system, index, copies in calls:
        if system not in grouped:
            raise KeyError(f'{system} has no variant site in bg_site_systems.{genome}.tsv')
        site = grouped[system][index]
        contigs.add(site.chrom)
        # --mixed diploidises only the extra site, which is what makes the file disagree
        # with itself; --diploidise applies to every haploid coordinate.
        site_diploidise = diploidise or (mixed and (system, index) == (MIXED_SYSTEM, MIXED_INDEX))

        block_start = site.pos - BLOCK_LEN - 1
        block_end = block_start + BLOCK_LEN - 1
        block_gt = render_gt(site.chrom, block_start, 0, diploidise=site_diploidise)
        block_record = (
            f'{site.chrom}\t{block_start}\t.\t{site.ref[0]}\t<NON_REF>\t.\t.\t'
            f'END={block_end}\tGT:DP:GQ:MIN_DP\t{block_gt}:42:45:38'
        )
        rows.append((site.chrom, block_start, block_record))
        gt = render_gt(site.chrom, site.pos, copies, diploidise=site_diploidise)
        rows.append(
            (
                site.chrom,
                site.pos,
                f'{site.chrom}\t{site.pos}\t.\t{site.ref}\t{site.alt},<NON_REF>\t320\tPASS\t.\tGT:DP:GQ\t{gt}:50:48',
            )
        )

    rows.sort(key=lambda r: (bg_db.chrom_key(r[0]), r[1]))
    header = HEADER + f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n'
    return header + ''.join(f'{r[2]}\n' for r in rows)


def main(argv: list[str] | None = None) -> int:
    """Write one fixture gVCF.

    Args:
        argv: Command-line arguments, for testing.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('genome', help='Genome build, e.g. GRCh38')
    parser.add_argument('out', type=Path, help='Path to write the uncompressed gVCF to')
    shape = parser.add_mutually_exclusive_group()
    shape.add_argument(
        '--diploidise',
        action='store_true',
        help='Pre-expand every haploid call, as the removed `bcftools +fixploidy` step left it',
    )
    shape.add_argument(
        '--mixed',
        action='store_true',
        help='Expand one non-PAR chrX call and not the rest, the state rbceq2 2.4.4 refuses',
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    text = build(args.genome, diploidise=args.diploidise, mixed=args.mixed)
    args.out.write_text(text)
    shape_name = 'partially diploidised' if args.mixed else ('diploidised' if args.diploidise else 'native haploid')
    logger.info(f'Wrote {args.out} ({shape_name}), {text.count(chr(10)) - text.count("##") - 1} records')
    return 0


if __name__ == '__main__':
    sys.exit(main())
