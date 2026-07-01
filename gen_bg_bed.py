#!/usr/bin/env python3
"""Generate the blood-group regions BED from RBCeq2's allele DB.

Mirrors rbceq2.IO.vcf.build_intervals()/parse_positions(): this is the exact
site set RBCeq2 reads, so the output is a guaranteed superset for upstream
`bcftools view -R`. Usage: gen_bg_bed.py db.tsv GRCh38 [flank] > regions.bed
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

def build_intervals(rows, genome, flank):
    iv = defaultdict(list)
    for r in rows:
        chrom = norm_chrom(r["Chrom"])
        for pos in parse_positions(r[genome]):
            iv[chrom].append((max(0, pos - flank), pos + flank))
    merged = {}
    for chrom, lst in iv.items():
        lst.sort()
        m = []
        for s, e in lst:
            if not m or s > m[-1][1]:
                m.append([s, e])
            else:
                m[-1][1] = max(m[-1][1], e)
        merged[chrom] = m
    return merged

def chrom_key(c):
    c = c.removeprefix("chr")
    return (int(c), "") if c.isdigit() else (99, c)

db_path, genome = sys.argv[1], sys.argv[2]
flank = int(sys.argv[3]) if len(sys.argv) > 3 else 500_000
with open(db_path, newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
merged = build_intervals(rows, genome, flank)
for chrom in sorted(merged, key=chrom_key):
    for s, e in merged[chrom]:
        print(f"{chrom}\t{s}\t{e}")