# DRAGEN SV caller versus DRAGEN CNV caller in the blood-group regions: 50 OurDNA genomes

**Date:** 2026-09-17, replicated 2026-09-18 on a second 50 (§10) with a read-level check (§11), a third 50 (§12) and the review-round checks (§13) ·
**Caller naming:** the SV VCF is written by the DRAGEN 3.7.8 SV caller, which integrates and extends Manta; record IDs keep the SV caller's prefix. This note says "SV caller" or "SV", not "Manta", and claims about the caller are observations of its output, not of the SV caller's algorithm. ·
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
| SV caller DEL, PASS | 4296 | 879 | 41 | 3 | 0 |
| SV caller DUP, PASS | 0 | 96 | 5 | 0 | 28 |
| SV caller INS, PASS | 7008 | 0 | 0 | 0 | 0 |

Per sample: CNV PASS events 1 to 7 (median 3); SV PASS DEL/DUP 80 to 133 (median 106);
CNV non-PASS 9 to 25 (median 16).

- Below 10 kb the CNV caller emits nothing that passes; every record is `cnvLength`. The SV
  caller is the only source. This covers 22 of the 67 target alleles (GE, PIGG, MAM, ABO*O.16, XK*N.03,
  XK*N.05, small GYP).
- At 10 kb and above the CNV caller dominates. the SV caller's 28 PASS DUPs of 200 kb or more are five
  whole-arm `DUP:TANDEM` artefacts (chr3 75–130 Mb in 24 samples, chr11, chr8, chr1 with
  `MaxDepth`), not events. Its 3 DELs of 50–200 kb are real-looking but the CNV caller saw
  them too.

## 2. Where both callers report one event

PASS DEL/DUP with reciprocal overlap ≥0.5 and the same direction:

| | <1 kb | 1–10 kb | 10–50 kb | 50–200 kb | ≥200 kb |
|---|---|---|---|---|---|
| Both | 0 | 0 | 46 | 0 | 0 |
| CNV only | 0 | 0 | 98 | 20 | 5 |
| SV caller only | 4296 | 973 | 2 | 3 | 28 |

The overlap is confined to 10–50 kb, and in that band the SV caller adds 2 events the CNV caller
missed against 98 the reverse. The 46 shared events are ten loci; 25 are one common HLA-DRB
deletion at chr6:32.49 Mb. Across them the breakpoints **never agree**: |POS delta| 224 to
3816 bp (median 515), |length delta| 28 to 8149 bp (median 3775), zero identical pairs.
Genotype agrees in 36 of 46. The CNV caller's bins are 1.1 to 2 kb wide in these regions
(0.49 to 0.94 bins per kb, median 0.86), so its breakpoints are bin-snapped and its lengths
run kb-scale longer than the SV caller's.

## 3. Running rbceq2's matcher on the combined records

For each sample, CNV records were rewritten as the SPEC prescribes (drop `DRAGEN:REF:`,
`SVTYPE=CNV` → `DEL`/`DUP` from ALT) and concatenated with the SV records (BND dropped).
Five policies:

| Policy | Records | Samples raising `ambiguous_equal_best_sv_evidence` | Selected db matches |
|---|---|---|---|
| all | every record, any FILTER | 0 / 50 | 1 |
| pass | PASS only, both sources | 0 / 50 | 1 |
| split | PASS; SV <10 kb plus all INS, CNV ≥10 kb | 0 / 50 | 1 |
| cnv_only | PASS CNV | 0 / 50 | 0 |
| sv_only | PASS SV | 0 / 50 | 1 |

The one match is **a heterozygous GE\*01.-02.01 in one genome**, from the SV caller: chr2:126690214,
`SVLEN=-3609`, QUAL 506, PASS, exact length of the db token. The CNV caller saw the same
event as `DRAGEN:LOSS:chr2:126688900-126693860`, 4961 bp, CN=1, `cnvLength`, QUAL 29. This
is the 3.6 kb Gerbich exon 2 deletion (Ge:-2) and is a **real positive control** for the merge
path, which the five NA12878 replicates could not provide. No other sample carries an SV
record within 2 kb of it.

No tie fired on real data because no two records ever had identical geometry (§2).

## 4. The tie error needs identical breakpoints; the real hazard is double-matching

Synthetic test built on the real GE record, adding a CNV-style record for the same deletion:

| Added CNV record | Result |
|---|---|
| identical breakpoints, same GT, PASS | OK, one match (SV) |
| identical breakpoints, GT 1/1 | **raises** |
| identical breakpoints, FILTER `cnvLength` | **raises** |
| start +300 bp, same length | OK, **two matches: GE\*01.-02.01 (SV) and GE\*01.-03.03 (CNV)** |
| start +300 bp, length +200 bp | OK, two matches as above |
| start −1 kb, end +1 kb (bin-snapped) | OK, one match (SV); CNV record rejected by the length gate |
| two CNV records, ±400 bp, no SV record | **raises** |

So:

- The 2.4.4 error is real but fires only when two records share exact coordinates and
  differ in GT or FILTER. The SV caller and the CNV caller never produce that pair (§2). It could arise
  if we ever concatenated a record with a copy of itself.
- The damaging case is quieter. `select_best_per_vcf` keeps at most one db definition per
  *event*, not per locus. Two records for one deletion, offset by a few hundred bases, each pick
  their nearest db token and rbceq2 receives **two different GE alleles for one deletion**. The
  seven GE 3.6 kb alleles differ only by breakpoint, so the CNV record's bin-snapped start lands
  on a different one than the SV caller's exact start. The same applies to the three A4GALT deletions
  (21, 26, 33 kb, overlapping) and the GYP cluster. This is the reason to deduplicate.
- Bin-snapped CNV breakpoints are also the wrong evidence for choosing *which* sub-10 kb
  allele; the SV caller's are exact to within `CIPOS` (0,50 here). That independently argues for
  the SV record winning the sub-10 kb band even where the CNV caller has a record.

## 5. Recommended Q2 resolution

Superseded by SPEC §5.5 (2026-09-18), which frames this as triage rather than caller ownership:
restrict both structural VCFs to records that overlap a db SV definition or span a defining
SNV site, keep every survivor with at least 3 bins, and where two records describe one event
keep the SV caller's for its exact breakpoints. The observations that drove it: the callers' size bands
barely overlap (§1–2), the double-matching hazard (§4), a real sub-10 kb deletion only the CNV
caller reported (§11), and rbceq2's silence on deletions it cannot match to the db (SPEC §6).

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
samples carry a common 3.46 kb SV-caller deletion at chr4:144395293, 250 kb downstream of GYPA,
not in the db. So the CNV caller does report dosage change across the GYP paralogs at
whole-gene scale, and the GYPB deletion targets are tiled with 80 to 190 bins in 49 of 50
samples. Whether a true whole-GYPB deletion (U−) is called correctly still needs a carrier.

## 9. Reproducing

Scripts live in the session scratchpad (`overlap.py`, `matcher_exp.py`, `near_miss.py`,
`tie_test.py`, `assess.py`, `triage.py`, and for §13 `adjacent.py` and `adjacent_value.py`); they read `bcftools view -R bg_regions.GRCh38.bed` extracts of the
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
| SV caller PASS DEL/DUP 10–200 kb with a CNV partner | 46 | 43 |
| SV caller PASS DEL/DUP 10–200 kb with no CNV partner | 5 | 6 |
| Shared events with identical breakpoints | 0 | 0 |
| Shared events, median POS / length delta | 515 bp / 3.8 kb | 415 bp / 6.5 kb |
| Samples raising the 2.4.4 tie error, any policy | 0 | 0 |
| SV caller PASS DUP ≥200 kb (whole-arm artefacts, same chr3/chr11 loci) | 28 | 28 |
| Large PASS CNV events on chrX | 5, in the X0 and XYY samples | 0 |
| Haploid GT on chrX SV/CNV records | 1 | 2 |
| 10 kb+ targets tiled by CNV records in ≥48/50 | all | all (50/50) |
| db matches under the `pass` or `split` policy | 1 (GE\*01.-02.01, Manta) | 0 |

Every structural claim in §§1–7 held. The six SV-only 10 kb+ events in the second set are
the same recurrent loci as the first (HLA-DRB 88 kb, chr7:100.73 Mb 13 kb, a chr9:133.06 Mb
20 kb duplication 170 kb upstream of ABO in 5 samples), each with a `cnvLength` or `cnvQual`
CNV record or a REF record under it; none touches a target.

**One new observation.** Under the `all` policy one second-set genome matched GE\*01.-02.02
from a `cnvLength` CNV record at exactly the coordinates of the first set's SV-called deletion
(`DRAGEN:LOSS:chr2:126688900-126693860`, 4961 bp, CN=1, SM 0.57, 3 bins, QUAL 21) with **no
SV record of any kind within 20 kb**. The first set's genome had the same CNV record
(SM 0.53, QUAL 29) beside a QUAL 506 SV call. So either the SV caller missed a real 3.6 kb Gerbich
deletion here, or the CNV caller produced a 3-bin false positive. Across both sets, sub-10 kb
`cnvLength` deletions of 3 to 4 bins have a reciprocal SV deletion in 75–79% of cases, 5 or
more bins in 74–78%, 1 to 2 bins in 45%. See §11 for the read-level check.

## 11. Read-level check of the CNV-only Gerbich call: the SV caller missed a real deletion

Reads were streamed from three CRAMs over chr2:126,684,000–126,698,000 (the db deletion is
126,690,214–126,693,823): the genome with the CNV-only call, the genome with the SV call,
and a genome with neither. Summarised, not tabulated per genome:

- In the CNV-only genome, MAPQ≥20 depth fell to roughly half its flanking level across the
  four 1 kb windows that coincide with the db interval and returned to baseline on both sides.
  Nine read pairs with inserts over 2 kb and two split reads spanned the interval.
- In the SV-called genome the same shape appeared, with about four times as many
  long-insert pairs and three split reads.
- In the genome with neither call the depth was flat across the window and no long-insert
  pairs or split reads were present.

So the CNV-only call is a real heterozygous ~3.6 kb Gerbich deletion. Depth halves over exactly
the db interval; the paired-end evidence is present but thin, and the SV caller emitted
nothing. Manta scores paired and split support together and DRAGEN's caller extends Manta,
but why this caller stayed silent on nine long-insert pairs and two split reads is inferred
from its output, not known from its algorithm (review, 2026-09-18). The CNV caller found it
on 3 bins and filtered it `cnvLength`.

**Consequences for the design.**

1. The SV caller is not sensitive enough to own the sub-10 kb band alone. Two Gerbich carriers in 100
   genomes; the SV caller called one. Rule 1 in §5 has to change: a sub-10 kb `cnvLength` CNV
   deletion with **no** reciprocal SV record must be **kept**, its FILTER rewritten to PASS
   (the per-record mechanism the SPEC's Q1 already names; not `--no_filter`). Where the SV caller does
   have a partner, its record still wins, because its breakpoints are exact and the CNV
   caller's are bin-snapped (the CNV-only record matched GE\*01.-02.02, the SV record
   GE\*01.-02.01, for what is presumably the same allele).
2. Kept CNV-only sub-10 kb records are lower-grade evidence. Across both sets, 3–4-bin
   `cnvLength` deletions have a SV partner 75–79% of the time and 5+-bin ones 74–78%; the
   remainder are a mix of real events the SV caller missed (this one) and bin-level noise that only reads
   can separate. The QC has to mark an allele called from a CNV-only sub-10 kb record as
   provisional (a `SVLOWRES`-type flag naming the record, its `BC`, `SM` and `QUAL`), and a
   1–2-bin record (SV partner in 45%) should not be kept at all.
3. A minimum bin count of 3 keeps the true positive here and drops the noisiest tier. It is a
   threshold set from two carriers; it should be recorded in the Analysis meta like `min_depth`
   and revisited when more carriers or a long-read truth set exist.

## 12. Third disjoint 50 (seed 20260919) and the triage simulation

Same pipeline, 50 further genomes disjoint from both earlier sets. Sex mix 24 XX, 26 XY.

| Measure | First 50 | Second 50 | Third 50 |
|---|---|---|---|
| CNV PASS records under 10 kb | 0 | 0 | 0 |
| CNV PASS records 10 kb+ | 169 | 160 | 161 |
| SV caller 10–200 kb with / without a CNV partner | 46 / 5 | 43 / 6 | 42 / 1 |
| Shared events with identical breakpoints | 0 | 0 | 0 |
| Shared events, median POS / length delta | 515 bp / 3.8 kb | 415 bp / 6.5 kb | 493 bp / 6.0 kb |
| Samples raising the 2.4.4 tie error, any policy | 0 | 0 | 0 |
| SV caller PASS DUP ≥200 kb, same recurrent loci | 28 | 28 | 28 |
| Large PASS CNV events on chrX | 5 (X0, XYY) | 0 | 0 (all XX or XY) |
| Haploid GT on chrX SV/CNV records | 1 | 2 | 2 |
| 10 kb+ targets tiled by CNV records in ≥48/50 | all | all | all |
| `cnvLength` 3–4-bin deletions with SV-caller support | 75% | 79% | 76% |
| db matches, `pass` policy | 1 (GE) | 0 | 0 |

Every claim in §§1–7 and §10 held for a third time. The GYP 117 kb PASS duplication (§8)
recurred in one genome of this set too, so it is a polymorphism the CNV caller sees at the
paralog locus, still matching no db allele.

**Triage simulation (SPEC §5.5 rule 1) over all 150 genomes.** Records within `SvMatcher`
tolerance of any db SV definition: 3 in total, all Gerbich (two SV-visible, one CNV-only);
0 to 2 per genome. PASS deletions under 1 Mb spanning a defining SNV/indel site but matching
no db SV: 0 to 3 per genome, and almost all one thing: a recurrent **8.5–12.8 kb deletion at
chr19:48.69 Mb spanning every FUT2 defining site**, in 3 of 150 genomes, each reported by both
callers (SV 9.3–10.1 kb PASS, QUAL 400–797; CNV 8.5–12.8 kb, CN=1 on 6–9 bins, PASS in two and
`cnvLength` in one). Two of the three carriers share identical SV breakpoints
(chr19:48,697,201–48,706,493), consistent with one recurrent allele. (An earlier revision of this
paragraph said "5 of 150 genomes": the triage script counted PASS records, five, not carriers,
three. Corrected 2026-09-18.) It matches no db allele. A
heterozygous whole-FUT2 deletion leaves one FUT2 copy; whatever the gVCF reports at the
secretor-defining sites is then a single-copy call read as homozygous. rbceq2 does nothing with
an unmatched deletion, so this is exactly the case the `SVDEL` flag in SPEC §6 exists for. The
remaining rule-1(b) records were the X0 genome's megabase events (handled by the karyotype gate)
and one 126 Mb `MaxDepth` SV record, which is why rule 1(b) is PASS-only and capped at 1 Mb.

## 13. Review-round checks over all 150 genomes (2026-09-18)

Four questions from the first review of SPEC `popgen_rbceq2#16`, answered on the same extracts.

**Is the 3-bin Gerbich call three adjacent records or one?** One. In both carriers it is a single
record, `DRAGEN:LOSS:chr2:126688900-126693860`, `SVLEN=-4961`, `BC=3`, GT `0/1`, `CN=1`. `BC` is
the bin count inside one segment, not a count of records.

**Does either caller split one event into adjacent same-direction records?** Consecutive
DEL/DUP records on one chromosome, same direction, gap of at most 2.5 kb (about one bin):

| Caller | DEL/DUP records | Adjacent pairs | Pairs with the same `CN` | Pairs whose merged span matched a db allele that neither part matched |
|---|---|---|---|---|
| CNV | 2799 | 254 | 3 | 0 |
| SV (records ≥1 kb) | 17175 | 1 | n/a | 0 |

Widening the gap to 10 kb adds 12 CNV pairs and 1 SV pair and changes nothing else. The CNV
adjacencies sit at recurrent loci (chr6:32.4 Mb, 112 pairs; chr2:126.9 Mb, 38; chr8:74.4 Mb, 25;
chr1:158.8 Mb, 22) and 251 of 254 join records of different `CN`: a heterozygous loss beside a
homozygous one, or a `cnvLength` fragment beside a PASS segment. They are the CNV caller's
segmentation reporting a copy-number step, not one deletion cut in two. Merging them would
combine different copy states, and no merge produced a db match. So the SPEC adds no
merge-adjacent step (§7 there), while recording the reviewer's point that a target split across
fragments would fail rbceq2's length gate on each fragment alone (35% of the larger length,
50% at reciprocal overlap of 10% or more), which is the observation to re-test if a carrier ever
shows it.

**Do the FORMAT fields need harmonising before `bcftools concat`?** No. Over all 150 genomes the
SV VCF used exactly two FORMAT strings, `GT:FT:GQ:PL:PR:SR` (43,280 records) and `GT:FT:GQ:PL:PR`
(2,343), and the CNV VCF exactly one, `GT:SM:CN:BC:PE` (12,069). `GT` is first in all of them,
which is rbceq2's only requirement (`IO/vcf.py`, `_require_first_gt`). The six INFO and FORMAT
tags the two headers share (`GT`, `END`, `SVTYPE`, `SVLEN`, `CIPOS`, `CIEND`) have identical
Number and Type, and a real `bcftools concat -a` of one genome's region-restricted SV and CNV
files (sample renamed to match) ran with no header warning, sorted and indexed. Genotypes on
CNV event records: every PASS `<DUP>` carries `./1` (244 of 244), every PASS `<DEL>` `0/1` (138)
or `1/1` (108); five non-PASS records carry a haploid `1`.

**Can the sex source fail, and what agrees with it?** The `Ploidy estimation` field read `XX` in
86, `XY` in 62, `XYY` in 1 and `X0` in 1 of 150. Its ratios separate cleanly: XX genomes had
X/autosomal 0.97–1.02 and Y/autosomal 0.00; XY had 0.50–0.52 and 0.42–0.52; the XYY genome had
Y 0.93; the X0 genome Y 0.18. Review asked for a fallback because DRAGEN omits the ploidy file
on genomes with uneven coverage; the SPEC takes the somalier signals `single_sample_qc_popgen`
already registers in metamist rather than a coverage-ratio classifier, so these bands are
recorded here as a calibration and not used (§5.6 there). The CNV caller's own `SEX GENOTYPER`
line in `cnv_metrics.csv` agreed with the ploidy estimate in 49 of the 50 first-set genomes and
was blank for the X0 genome, so it is not a usable fallback either.
