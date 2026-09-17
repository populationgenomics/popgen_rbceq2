# Manta SV versus DRAGEN CNV in the blood-group regions: 50 OurDNA genomes

**Date:** 2026-09-17, replicated 2026-09-18 on a second 50 (§10) with a read-level check (§11) ·
**Contains no individual-level data:** aggregate observations only, no sequencing-group IDs, no per-genome tables. · **Inputs:** `gs://cpg-ourdna-main/ica/dragen_3_7_8/output/dragen_metrics/<sg>/<sg>.{sv,cnv}.vcf.gz`
for 50 sequencing groups drawn with a fixed seed from the 1,255 present (not in Metamist; bulk
download). **Regions:** the committed `bg_regions.GRCh38.bed` (50 intervals, 52.6 Mb).
**Matcher:** rbceq2 2.4.4 `SvReader`/`SvMatcher`/`select_best_per_vcf`, imported directly, with
`load_db_defs(db, svtoken_col='GRCh38')` (176 definitions, 33 RH). Sex mix by DRAGEN's own
estimate: 30 XX, 18 XY, 1 X0, 1 XYY.

The question was SPEC Q2: when the same deletion arrives from both callers, do we have to
deduplicate before rbceq2, and if so on what rule? Answer: **yes, and a size split at 10 kb is
the right rule**, but for a different reason than the 2.4.4 tied-evidence error. Details follow.

## 1. The two callers occupy different size bands

PASS records in the blood-group regions, 50 samples pooled:

| Source and type | <1 kb | 1–10 kb | 10–50 kb | 50–200 kb | ≥200 kb |
|---|---|---|---|---|---|
| CNV, PASS | 0 | 0 | 144 | 20 | 5 |
| CNV, `cnvLength` (any combo) | 0 | 770 | 0 | 0 | 0 |
| CNV, `cnvQual` only | 0 | 0 | 5 | 0 | 0 |
| Manta DEL, PASS | 4296 | 879 | 41 | 3 | 0 |
| Manta DUP, PASS | 0 | 96 | 5 | 0 | 28 |
| Manta INS, PASS | 7008 | 0 | 0 | 0 | 0 |

Per sample: CNV PASS events 1 to 7 (median 3); Manta PASS DEL/DUP 80 to 133 (median 106);
CNV non-PASS 9 to 25 (median 16).

- Below 10 kb the CNV caller emits nothing that passes; every record is `cnvLength`. Manta is
  the only source. This covers 22 of the 67 target alleles (GE, PIGG, MAM, ABO*O.16, XK*N.03,
  XK*N.05, small GYP).
- At 10 kb and above the CNV caller dominates. Manta's 28 PASS DUPs of 200 kb or more are five
  whole-arm `DUP:TANDEM` artefacts (chr3 75–130 Mb in 24 samples, chr11, chr8, chr1 with
  `MaxDepth`), not events. Its 3 DELs of 50–200 kb are real-looking but the CNV caller saw
  them too.

## 2. Where both callers report one event

PASS DEL/DUP with reciprocal overlap ≥0.5 and the same direction:

| | <1 kb | 1–10 kb | 10–50 kb | 50–200 kb | ≥200 kb |
|---|---|---|---|---|---|
| Both | 0 | 0 | 46 | 0 | 0 |
| CNV only | 0 | 0 | 98 | 20 | 5 |
| Manta only | 4296 | 973 | 2 | 3 | 28 |

The overlap is confined to 10–50 kb, and in that band Manta adds 2 events the CNV caller
missed against 98 the reverse. The 46 shared events are ten loci; 25 are one common HLA-DRB
deletion at chr6:32.49 Mb. Across them the breakpoints **never agree**: |POS delta| 224 to
3816 bp (median 515), |length delta| 28 to 8149 bp (median 3775), zero identical pairs.
Genotype agrees in 36 of 46. The CNV caller's bins are 1.1 to 2 kb wide in these regions
(0.49 to 0.94 bins per kb, median 0.86), so its breakpoints are bin-snapped and its lengths
run kb-scale longer than Manta's.

## 3. Running rbceq2's matcher on the combined records

For each sample, CNV records were rewritten as the SPEC prescribes (drop `DRAGEN:REF:`,
`SVTYPE=CNV` → `DEL`/`DUP` from ALT) and concatenated with the Manta records (BND dropped).
Five policies:

| Policy | Records | Samples raising `ambiguous_equal_best_sv_evidence` | Selected db matches |
|---|---|---|---|
| all | every record, any FILTER | 0 / 50 | 1 |
| pass | PASS only, both sources | 0 / 50 | 1 |
| split | PASS; Manta <10 kb plus all INS, CNV ≥10 kb | 0 / 50 | 1 |
| cnv_only | PASS CNV | 0 / 50 | 0 |
| sv_only | PASS Manta | 0 / 50 | 1 |

The one match is **a heterozygous GE\*01.-02.01 in one genome**, from Manta: chr2:126690214,
`SVLEN=-3609`, QUAL 506, PASS, exact length of the db token. The CNV caller saw the same
event as `DRAGEN:LOSS:chr2:126688900-126693860`, 4961 bp, CN=1, `cnvLength`, QUAL 29. This
is the 3.6 kb Gerbich exon 2 deletion (Ge:-2) and is a **real positive control** for the merge
path, which the five NA12878 replicates could not provide. No other sample carries a Manta
record within 2 kb of it.

No tie fired on real data because no two records ever had identical geometry (§2).

## 4. The tie error needs identical breakpoints; the real hazard is double-matching

Synthetic test built on the real GE record, adding a CNV-style record for the same deletion:

| Added CNV record | Result |
|---|---|
| identical breakpoints, same GT, PASS | OK, one match (Manta) |
| identical breakpoints, GT 1/1 | **raises** |
| identical breakpoints, FILTER `cnvLength` | **raises** |
| start +300 bp, same length | OK, **two matches: GE\*01.-02.01 (Manta) and GE\*01.-03.03 (CNV)** |
| start +300 bp, length +200 bp | OK, two matches as above |
| start −1 kb, end +1 kb (bin-snapped) | OK, one match (Manta); CNV record rejected by the length gate |
| two CNV records, ±400 bp, no Manta | **raises** |

So:

- The 2.4.4 error is real but fires only when two records share exact coordinates and
  differ in GT or FILTER. Manta and the CNV caller never produce that pair (§2). It could arise
  if we ever concatenated a record with a copy of itself.
- The damaging case is quieter. `select_best_per_vcf` keeps at most one db definition per
  *event*, not per locus. Two records for one deletion, offset by a few hundred bases, each pick
  their nearest db token and rbceq2 receives **two different GE alleles for one deletion**. The
  seven GE 3.6 kb alleles differ only by breakpoint, so the CNV record's bin-snapped start lands
  on a different one than Manta's exact start. The same applies to the three A4GALT deletions
  (21, 26, 33 kb, overlapping) and the GYP cluster. This is the reason to deduplicate.
- Bin-snapped CNV breakpoints are also the wrong evidence for choosing *which* sub-10 kb
  allele; Manta's are exact to within `CIPOS` (0,50 here). That independently argues for
  Manta owning the sub-10 kb band even where the CNV caller has a record.

## 5. Recommended Q2 resolution

Deduplicate in the merge stage, by size band, before rbceq2 sees the file:

1. **Manta supplies every DEL/DUP under 10 kb and every INS.** A CNV record under 10 kb is
   dropped when a reciprocal Manta DEL/DUP exists. *(Amended by §11:)* one with **no** Manta
   partner and **at least 3 bins** is kept with FILTER rewritten to PASS and flagged by the QC
   as low-resolution evidence; fewer than 3 bins is dropped.
2. **The CNV caller supplies every DEL/DUP of 10 kb or more.** Manta DEL/DUP records of 10 kb
   or more are dropped when a same-direction PASS CNV record overlaps them reciprocally ≥0.5;
   the 2 in 50 samples that have no CNV partner are kept.
3. **Drop Manta records over 200 kb** unconditionally; they are whole-arm `DUP:TANDEM`
   artefacts and the length gate would reject them anyway.
4. Never emit two records with identical CHROM/POS/END/SVTYPE.

Rule 2 costs nothing observed: no 10 kb+ target allele had a Manta-only hit in 100 samples.
Rule 1 is what makes the GE call unambiguous where both callers see it, and what keeps the
one Manta missed (§11). All four are decidable from the two files alone.

## 6. Sex chromosomes: the concern in §5.6 of the SPEC is confirmed

Large (≥100 kb) PASS CNV events on chrX occurred in two of the 50 genomes, and only in the
two whose DRAGEN ploidy estimate was neither XX nor XY:

- The genome estimated **X0** (X/autosome coverage ratio 0.51, Y/autosome 0.18) carried four
  heterozygous CN=1 deletions of 1.4 to 10.6 Mb. Between them they cover the PAR1 boundary with
  CD99 and XG, the XK locus, and the ATP11C locus.
- The genome estimated **XYY** (Y/autosome 0.93) carried one CN=3 duplication of 646 kb across
  PAR1, covering CD99 and XG.

The 18 XY genomes produced no such events, so DRAGEN's CNV caller handles a normal male
correctly. The X0 estimate describes a genome with one X and a partial Y signal (probably a male with
mosaic Y loss or a Y-coverage artefact) that the caller has treated as diploid-X, emitting
CN=1 across the chromosome as heterozygous deletions. The XYY genome has three PAR1 copies,
which is a correct dosage but reads as a CD99/XG duplication. Neither produced a false allele
here only because the length gate rejected megabase events against 11–219 kb tokens; a shorter
segment would not be rejected. The merge therefore needs the per-sample ploidy estimate, and
should either drop chrX/chrY CNV records for samples whose estimate is not XX or XY, or flag
them, rather than pass them through.

Haploid GT was seen on exactly one chrX CNV record in 50 samples; not a practical concern
for the SV/CNV files.

## 7. What the CNV VCF says about assessability, per target

For each non-RH target of 1 kb or more: in how many samples do CNV records (REF or event)
tile the whole interval, and how many bins fall inside it (median, BC scaled by overlap)?

| Target(s) | Length | Tiled by REF | Touched by an event | Gap | Bins |
|---|---|---|---|---|---|
| ATP11C\*01N.01 | 219 kb | 49 | 1 (the X0 sample) | 0 | 188 |
| GYPB\*05N.05 / \*01N / \*05N.01 / \*05N.02 / \*05N.04 | 96–224 kb | 49 | 0–1 | 1 | 82–192 |
| GYPA\*01N / \*28N.01, GYP\*01N / \*201 / \*203 | 101–124 kb | 49 | 0 | 1 | 87–106 |
| ABCC4\*01N.01 | 68 kb | 50 | 0 | 0 | 61 |
| XK\*N.01 (18 alleles, one token) | 53 kb | 49 | 1 (X0) | 0 | 47 |
| GCNT2\*01N.06 | 41 kb | 50 | 0 | 0 | 37 |
| CTL2\*01N.01 | 37 kb | 50 | 0 | 0 | 29 |
| A4GALT\*N.01 / .02 / .03 | 21–33 kb | 48–49 | 1–2 | 0 | 18–29 |
| RHAG\*01N.15, LU\*02N.06, C4A\*N.01, GYPB\*05N.03 | 19–32 kb | 49–50 | 0 | 0–1 | 16–29 |
| XG\*01N.02 / .03, CD99\*01N.01 / .02 | 11–114 kb | 48 | 2 (X0, XYY) | 0 | 9–92 |
| MAM\*01N.04 / .05, XK\*N.05, GYP\*401 / \*402 | 8–8.5 kb | 49–50 | 0–1 | 0–1 | 7 |
| PIGG\*01N.05 / .06 / .07, ABO\*O.16 | 4–6.3 kb | 50 | 0 | 0 | 3–5 |
| GE 3.6 kb alleles (7), GYP\*501 / \*504, XK\*N.03 | 3.5–3.7 kb | 49 | 1 (the GE carrier) | 0 | 3 |
| GE\*01N.01, ABCG2\*01N.27, GYP\*503 / \*301, LU\*02N.02 | 1–2.3 kb | 50 | 0 | 0 | 1–2 |
| ABCC1\*01N.01 | 21 kb | 0 | 0 | 50 | — (wrong chromosome in the db) |

Two things for the QC design:

- For every 10 kb+ target the CNV VCF answers "did the caller look here" in every sample: a
  REF record tiles it with tens to hundreds of bins, or an event touches it. The one GYP gap
  is the 117 kb duplication described in §8. This is the analogue of a gVCF reference block, with `BC` and
  `SM` as its quality.
- Below 10 kb a target has 1 to 7 bins. The CNV caller cannot resolve it (which is why
  `cnvLength` exists), so "tiled by REF" means nothing there. Assessability for GE, PIGG, MAM,
  ABO*O.16 and the small GYP/XK alleles has to come from somewhere other than the CNV VCF.

## 8. GYP locus, for the paralog question

Not a target hit, but relevant: one genome carries a PASS `<DUP>` of 117 kb at
chr4:143999799–144117036 (CN=3, 43 bins) spanning GYPB 3′ to GYPA. It matches no db
definition (the only DUP token is GYP\*505's 20 kb). Another carries a 12 kb deletion at
chr4:144.28 Mb seen by both callers, between GYPA and GYPE, also not in the db. And 36 of 50
samples carry a common 3.46 kb Manta deletion at chr4:144395293, 250 kb downstream of GYPA,
not in the db. So the CNV caller does report dosage change across the GYP paralogs at
whole-gene scale, and the GYPB deletion targets are tiled with 80 to 190 bins in 49 of 50
samples. Whether a true whole-GYPB deletion (U−) is called correctly still needs a carrier.

## 9. Reproducing

Scripts live in the session scratchpad (`overlap.py`, `matcher_exp.py`, `near_miss.py`,
`tie_test.py`, `assess.py`); they read `bcftools view -R bg_regions.GRCh38.bed` extracts of the
two VCFs and rbceq2 2.4.4 from PyPI. Sample list: 50 IDs drawn with `random.seed(20260917)`
from the sorted `dragen_metrics/` listing. No sequencing-group identifier appears in this note by
design: these are cohort observations, not individual-level data, and none of the genomes
described here is to be used as a test fixture. The scripts should be committed under `tests/` or a
`scripts/research/` directory if this becomes the merge stage's concordance check.

## 10. Replication on a disjoint second 50 (seed 20260918)

Same pipeline, 50 further genomes with no overlap with the first set. Sex mix 32 XX, 18 XY.

| Measure | First 50 | Second 50 |
|---|---|---|
| CNV PASS records under 10 kb | 0 | 0 |
| CNV PASS records 10 kb+ | 169 | 160 |
| Manta PASS DEL/DUP 10–200 kb with a CNV partner | 46 | 43 |
| Manta PASS DEL/DUP 10–200 kb with no CNV partner | 5 | 6 |
| Shared events with identical breakpoints | 0 | 0 |
| Shared events, median POS / length delta | 515 bp / 3.8 kb | 415 bp / 6.5 kb |
| Samples raising the 2.4.4 tie error, any policy | 0 | 0 |
| Manta PASS DUP ≥200 kb (whole-arm artefacts, same chr3/chr11 loci) | 28 | 28 |
| Large PASS CNV events on chrX | 5, in the X0 and XYY samples | 0 |
| Haploid GT on chrX SV/CNV records | 1 | 2 |
| 10 kb+ targets tiled by CNV records in ≥48/50 | all | all (50/50) |
| db matches under the `pass` or `split` policy | 1 (GE\*01.-02.01, Manta) | 0 |

Every structural claim in §§1–7 held. The six Manta-only 10 kb+ events in the second set are
the same recurrent loci as the first (HLA-DRB 88 kb, chr7:100.73 Mb 13 kb, a chr9:133.06 Mb
20 kb duplication 170 kb upstream of ABO in 5 samples), each with a `cnvLength` or `cnvQual`
CNV record or a REF record under it; none touches a target.

**One new observation.** Under the `all` policy one second-set genome matched GE\*01.-02.02
from a `cnvLength` CNV record at exactly the coordinates of the first set's Manta-called deletion
(`DRAGEN:LOSS:chr2:126688900-126693860`, 4961 bp, CN=1, SM 0.57, 3 bins, QUAL 21) with **no
Manta record of any kind within 20 kb**. The first set's genome had the same CNV record
(SM 0.53, QUAL 29) beside a QUAL 506 Manta call. So either Manta missed a real 3.6 kb Gerbich
deletion here, or the CNV caller produced a 3-bin false positive. Across both sets, sub-10 kb
`cnvLength` deletions of 3 to 4 bins have a reciprocal Manta deletion in 75–79% of cases, 5 or
more bins in 74–78%, 1 to 2 bins in 45%. See §11 for the read-level check.

## 11. Read-level check of the CNV-only Gerbich call: Manta missed a real deletion

Reads were streamed from three CRAMs over chr2:126,684,000–126,698,000 (the db deletion is
126,690,214–126,693,823): the genome with the CNV-only call, the genome with the Manta call,
and a genome with neither. Summarised, not tabulated per genome:

- In the CNV-only genome, MAPQ≥20 depth fell to roughly half its flanking level across the
  four 1 kb windows that coincide with the db interval and returned to baseline on both sides.
  Nine read pairs with inserts over 2 kb and two split reads spanned the interval.
- In the Manta-called genome the same shape appeared, with about four times as many
  long-insert pairs and three split reads.
- In the genome with neither call the depth was flat across the window and no long-insert
  pairs or split reads were present.

So the CNV-only call is a real heterozygous ~3.6 kb Gerbich deletion. Depth halves over exactly
the db interval; the paired-end evidence is present but thin, and Manta, which wants paired and
split support together, emitted nothing. The CNV caller found it on 3 bins and filtered it
`cnvLength`.

**Consequences for the design.**

1. Manta is not sensitive enough to own the sub-10 kb band alone. Two Gerbich carriers in 100
   genomes; Manta called one. Rule 1 in §5 has to change: a sub-10 kb `cnvLength` CNV
   deletion with **no** reciprocal Manta record must be **kept**, its FILTER rewritten to PASS
   (the per-record mechanism the SPEC's Q1 already names; not `--no_filter`). Where Manta does
   have a partner, Manta's record still wins, because its breakpoints are exact and the CNV
   caller's are bin-snapped (the CNV-only record matched GE\*01.-02.02, the Manta record
   GE\*01.-02.01, for what is presumably the same allele).
2. Kept CNV-only sub-10 kb records are lower-grade evidence. Across both sets, 3–4-bin
   `cnvLength` deletions have a Manta partner 75–79% of the time and 5+-bin ones 74–78%; the
   remainder are a mix of real events Manta missed (this one) and bin-level noise that only reads
   can separate. The QC has to mark an allele called from a CNV-only sub-10 kb record as
   provisional (a `SVLOWRES`-type flag naming the record, its `BC`, `SM` and `QUAL`), and a
   1–2-bin record (Manta partner in 45%) should not be kept at all.
3. A minimum bin count of 3 keeps the true positive here and drops the noisiest tier. It is a
   threshold set from two carriers; it should be recorded in the Analysis meta like `min_depth`
   and revisited when more carriers or a long-read truth set exist.
