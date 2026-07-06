#!/usr/bin/env python3
"""Emit one BED feature per blood-group determinant SITE from RBCeq2's allele DB.

Unlike gen_bg_bed.py (which pads +/-500 kb to build a VCF *read* filter), this emits
the determinant positions themselves - the exact sites RBCeq2 needs a call at.
Use it to find which sites your capture kit misses:

    gen_bg_sites.py db.tsv GRCh38 [flank] > sites.bed
    bedtools intersect -v -a sites.bed -b targets.bed   # sites NOT covered by targets

Coordinates are 0-based half-open BED (a 1-based DB position P -> `P-1  P`).
Column 4 lists the allele(s) that use each position. Optional `flank` (default 0)
adds a few bp either side - handy for splice-site / small-indel tolerance.
"""
import csv, sys
from collections import defaultdict


def parse_positions(cell):
    out = []
    if cell is None:
        return out
    for tok in str(cell).split(","):
        if "_" in tok:
            p = tok.split("_", 1)[0]
            if p.isdigit():
                out.append(int(p))
    return out


def norm_chrom(c):                       # chrx -> chrX, match gVCF contigs
    return "chr" + c.removeprefix("chr").upper()


def chrom_key(c):
    c = c.removeprefix("chr")
    return (int(c), "") if c.isdigit() else (99, c)


db_path, genome = sys.argv[1], sys.argv[2]
flank = int(sys.argv[3]) if len(sys.argv) > 3 else 0

# (chrom, pos) -> set of allele names that use this position
sites = defaultdict(set)
with open(db_path, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        chrom = norm_chrom(r["Chrom"])
        allele = r["Genotype"]
        for pos in parse_positions(r[genome]):
            sites[(chrom, pos)].add(allele)

# emit sorted, one 0-based half-open feature per site
for chrom, pos in sorted(sites, key=lambda cp: (chrom_key(cp[0]), cp[1])):
    start = max(0, pos - 1 - flank)      # 1-based P -> 0-based P-1, minus optional flank
    end = pos + flank
    labels = ",".join(sorted(sites[(chrom, pos)]))
    print(f"{chrom}\t{start}\t{end}\t{labels}")