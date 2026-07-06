#!/usr/bin/env python3
"""Classify RBCeq2 db.tsv allele rows by transcript region, from the HGVS `Coding`
column (c. notation). Tells you which blood-group determinants are exonic vs
intronic / promoter / UTR - i.e. which are reachable by whole-exome capture.

    classify_db_regions.py db.tsv [--by-bg] [--out classified.tsv]

Region is inferred per allele row from HGVS c. syntax:
  c.-67T>C          5'UTR / promoter   (negative position, upstream of ATG)
  c.*28G>A          3'UTR              (after the stop codon)
  c.123+5A>G        intronic           (+/-N offset from an exon boundary)
  c.123-2A>G        intronic
  c.679G>C          exonic (coding)
  (empty)           no_coding          (usually reference alleles / large SVs)

A row can list several changes (';'-separated) spanning >1 region; the single
label uses precedence intronic > promoter > 3'UTR > exonic so the "hardest for
WES" wins. Independent per-pattern counts (which may overlap) are also printed
to reproduce the raw grep tallies.
"""
import csv
import re
import sys
from collections import Counter, defaultdict

# patterns operate on the LOWERCASED Coding cell (so C.895C>T matches like c.895C>T)
RE_INTRONIC = re.compile(r"c\.[-0-9]+[+-][0-9]+")   # c.NNN+/-N offset from exon = intronic
RE_IVS = re.compile(r"ivs[-0-9]+[+-][0-9]+")        # old IVS28+1g>a style = intronic
RE_PROMOTER = re.compile(r"c\.-[0-9]")              # c.-NN = 5'UTR / promoter
RE_UTR3 = re.compile(r"c\.\*[0-9]")                 # c.*NN = 3'UTR
RE_CODING = re.compile(r"c\.[0-9]|^[0-9]+[acgt]")   # c.NN or bare cDNA pos (e.g. 3218A>G)
RE_STRUCTURAL = re.compile(r"\bg\.|nc_0000|exon|intron|gene|deletion|triplicat|del|dup")


def classify(coding):
    """Return a single region label for one Coding cell (precedence applied)."""
    c = (coding or "").strip()
    if not c:
        return "no_coding"
    cl = c.lower()
    if RE_INTRONIC.search(cl) or RE_IVS.search(cl):
        return "intronic"
    if RE_PROMOTER.search(cl):
        return "promoter_5utr"
    if RE_UTR3.search(cl):
        return "utr3"
    if RE_CODING.search(cl):
        return "exonic"
    if RE_STRUCTURAL.search(cl):
        return "structural"     # gene/exon del/dup, g. notation - large, WES-missed
    return "other"              # notation we still did not parse


def bg_of(genotype):
    return genotype.split("*")[0]


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    by_bg = "--by-bg" in args
    out_path = args[args.index("--out") + 1] if "--out" in args else None
    db_path = args[0]

    rows_out = []
    labels = Counter()
    raw = Counter()                     # independent (overlapping) pattern hits
    by_bg_counts = defaultdict(Counter)

    with open(db_path, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            coding = (r.get("Coding") or "").strip()
            label = classify(coding)
            geno = r["Genotype"]
            labels[label] += 1
            by_bg_counts[bg_of(geno)][label] += 1
            rows_out.append((geno, coding, label))

            # independent tallies (overlap allowed) - reproduces the grep counts
            if not coding:
                raw["empty_coding"] += 1
            if RE_INTRONIC.search(coding):
                raw["intronic(+/-offset)"] += 1
            if RE_PROMOTER.search(coding):
                raw["promoter_5utr(c.-N)"] += 1
            if RE_UTR3.search(coding):
                raw["utr3(c.*N)"] += 1

    total = sum(labels.values())
    print(f"# {total} allele rows\n")
    print("## Exclusive classification (precedence: intronic > promoter > 3'UTR > exonic)")
    for lab, n in labels.most_common():
        print(f"  {lab:14s} {n:5d}  ({100 * n / total:4.1f}%)")

    print("\n## Independent pattern tallies (overlapping - matches the grep counts)")
    for lab, n in raw.most_common():
        print(f"  {lab:22s} {n:5d}")

    if by_bg:
        print("\n## Non-exonic determinants by blood group")
        for bg in sorted(by_bg_counts):
            c = by_bg_counts[bg]
            detail = ", ".join(f"{k}:{v}" for k, v in c.items() if k != "exonic")
            if detail:
                print(f"  {bg:8s} {detail}")

    if out_path:
        with open(out_path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["Genotype", "Coding", "region"])
            w.writerows(rows_out)
        print(f"\nwrote per-allele table -> {out_path}")


if __name__ == "__main__":
    main()