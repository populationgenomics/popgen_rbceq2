# What the SNV/indel gVCF resolves today, and what a merged SNV/indel/SV/CNV VCF adds

**Date:** 2026-09-17 · **Database:** rbceq2 2.4.4 (bundled db 2.5.1), `GRCh38` column ·
**Method:** every allele's coordinate tokens classified by the largest event they name.
Script and raw output are reproducible from `db.tsv` alone (see §4).

RH (RHD, RHCE) is listed but **excluded from the gain**: the maintainer's advice on the
2.4.4 release is that DRAGEN SV/CNV does not reliably detect the RH hybrids, so `--RH`
stays off. This note is about every other system.

## 1. Allele classes

| Class | Definition | Where it can come from |
|---|---|---|
| **A small** | Every token is a SNV, a `_ref` lane site, or an indel under 50 bp | DRAGEN SNV/indel gVCF: **what we run today** |
| **B mid indel** | Largest event 50 bp to 1 kb | Sometimes the gVCF (DRAGEN calls indels up to a few hundred bp), otherwise the Manta SV VCF. Assume **not** reliably called today. |
| **C large, single** | One event of 1 kb or more: a deletion, duplication or insertion, not paired | 1 kb to 10 kb: Manta SV VCF only (DRAGEN CNV filters these as `cnvLength`). 10 kb and up: CNV VCF, and usually Manta too. **This is the gain.** |
| **D hybrid** | Two or more large events, a deletion paired with an insertion or duplication | Gene-conversion products. Long read only. |

rbceq2 reads B, C and D from the same VCF via `SvReader`, matching on `SVTYPE` and
fuzzy position/length. 2.4.4 did **not** change `SvMatcher.compatible`: `SVTYPE` must
equal the db type, and DRAGEN's `SVTYPE=CNV` still fails against `DEL`/`DUP`, so the
CNV rewrite in the SPEC (§5.4) is still required.

## 2. Per-system counts (alleles with a GRCh38 coordinate)

Only systems with at least one non-A allele are listed. Every other system, including all
35 HPA systems, ABCB6, CROM, DO, FUT1/2/3, KEL, KLF, KN, LW, SC, VEL, YT and the
rest, is **entirely class A**, so the gVCF already resolves everything the db defines.

| System | Total | A small (now) | B mid indel | C large single (gain) | D hybrid | Notes |
|---|---|---|---|---|---|---|
| ABCC1 | 1 | 0 | 0 | 1 | 0 | 21 kb del. **Uncallable on any input:** the row carries ABCC4's chr13 coordinates on chr16, an upstream db error still present in 2.4.4 (README, "Structural-variant entries"). Counted here, gained nowhere. |
| ATP11C | 1 | 0 | 0 | 1 | 0 | 219 kb whole-gene del, non-PAR chrX. **Blind today.** |
| CD99 | 3 | 1 | 0 | 2 | 0 | 11 kb and 20 kb dels, PAR1. **Effectively blind today.** |
| XG | 3 | 1 | 0 | 2 | 0 | 32 kb and 114 kb dels, PAR1 boundary. |
| XK | 67 | 44 | 3 | 20 | 0 | 19 alleles share one 53 kb token (McLeod cytogenetic dels); plus 8 kb and 3.5 kb dels. Non-PAR chrX. |
| GYPB | 40 | 33 | 1 | 6 | 0 | Whole-GYPB dels, 19 kb to 224 kb: the U- / S-s-U- alleles. Paralog locus. |
| GYPA | 26 | 24 | 0 | 2 | 0 | 101 kb and 119 kb dels. Paralog locus. |
| GYP | 18 | 2 | 4 | 8 | 4 | The MNS hybrid system. 3 of the C are the 122 kb del of GYP*201/202/203 (with SNVs); 5 are 1.8 kb to 3.6 kb explicit-sequence dels. Paralog locus. |
| GE | 20 | 10 | 2 | 8 | 0 | Seven 3.6 kb exon dels (Ge:-2 / Ge:-3, Melanesian) and one 2.3 kb del. All **under 10 kb**: Manta only. |
| A4GALT | 43 | 40 | 0 | 3 | 0 | 21 kb to 33 kb dels (P1PK null). |
| PIGG | 8 | 4 | 1 | 3 | 0 | 4 kb to 6.3 kb dels. Manta only. |
| MAM | 5 | 1 | 2 | 2 | 0 | 8.5 kb dels. Manta only. |
| LU | 38 | 36 | 0 | 2 | 0 | 1 kb and 27 kb dels. |
| ABO | 207 | 206 | 0 | 1 | 0 | ABO*O.16, 6 kb del. Manta only. |
| ABCG2 | 32 | 31 | 0 | 1 | 0 | 1.8 kb delins. Manta only. |
| ABCC4 | 4 | 3 | 0 | 1 | 0 | 68 kb whole-gene del. |
| C4A | 3 | 2 | 0 | 1 | 0 | 20 kb del, MHC segmental duplication with C4B. |
| CTL2 | 4 | 3 | 0 | 1 | 0 | 37 kb del. |
| GCNT2 | 13 | 12 | 0 | 1 | 0 | 41 kb del. |
| RHAG | 47 | 46 | 0 | 1 | 0 | 32 kb del. |
| CO | 10 | 9 | 1 | 0 | 0 | 384 bp del. |
| FY | 32 | 31 | 1 | 0 | 0 | 192 bp delins. |
| JK | 88 | 87 | 1 | 0 | 0 | 884 bp del. |
| DI | 25 | 24 | 0 | 0 | 0 | One allele's token has no coordinate (`unparsed`). |
| *RHD* | *444* | *426* | *3* | *2* | *12* | *Excluded. `--RH` off.* |
| *RHCE* | *186* | *179* | *0* | *1* | *6* | *Excluded. `--RH` off.* |
| **All non-RH (86 systems)** | **1388** | **1300** | **16** | **67** | **4** | plus 1 unparsed |

Totals across the whole db, RH included: 2018 alleles, 1905 A, 19 B, 70 C, 22 D, 2 unparsed.

System tally, non-RH: 86 systems in the db (the glossary's 48 predates the 35 HPA systems).
62 are entirely class A (35 HPA plus 27 others). 20 have at least one class-C gain allele.
3 (CO, FY, JK) have only a mid-size indel beyond class A. 1 (DI) has one allele with no
coordinate. 62 + 20 + 3 + 1 = 86.

## 3. Reading the numbers

**What we resolve now.** 1300 of 1388 non-RH alleles (93.7%) are class A. 62 of the 86
non-RH systems are fully class A. For those, the gVCF plus the existing QC is the
whole story and a merged VCF changes nothing.

**What a merged VCF gains.** 67 single-event large alleles across 20 systems (66 across 19 once ABCC1 is set aside), plus a fair
share of the 16 mid-size indels. In system terms:

- **Two systems go from blind to callable:** ATP11C and CD99, plus XG, which has one SNV
  allele but its two null alleles are dels. Today the QC reports these `NA` and rbceq2
  emits the reference phenotype for every sample. ABCC1 is also `NA` today but stays
  uncallable: its only allele has the wrong coordinates in the db.
- **Null (antigen-negative) alleles become visible in 16 systems** that are otherwise
  callable. These are the clinically interesting ones for rare-donor work: whole-GYPB
  deletion (U-), XK deletion (McLeod), A4GALT deletion (p phenotype), GE exon deletions
  (Ge:-2, Ge:-3; common in PNG and Melanesia, relevant to our cohorts), C4A deletion (Ch-),
  LU, RHAG, CTL2, GCNT2, ABCC4 nulls.
- **Not gained:** the 4 GYP hybrids and the 18 RH hybrids (class D), and anything whose
  locus mapping defeats short reads. GYPA/GYPB/GYPE and C4A/C4B are paralog pairs; a
  whole-gene deletion is a dosage loss the bin-based CNV caller can see even where
  breakpoints are unmappable, but it needs verifying on a known positive rather than
  assumed.

**Two size bands, two sources.** 22 of the 67 C alleles are under 10 kb, which DRAGEN CNV
filters as `cnvLength`; 45 are 10 kb or more. GE, PIGG, MAM, ABO*O.16, XK*N.05, ABCG2 and the small GYP dels
therefore depend on the **Manta SV VCF**, not the CNV VCF. The 10 kb and larger events
(whole-gene dels: GYPB, XK, ATP11C, A4GALT, CTL2, GCNT2, ABCC4, RHAG, LU*02N.06, CD99,
XG, C4A) can come from either. This is why the SPEC merges all three files
rather than adding just the CNV VCF.

**Exomes.** Everything above assumes genome DRAGEN outputs. The mackenzie exome runs
may have no `sv.vcf.gz`/`cnv.vcf.gz`, or a CNV VCF that needs a panel of normals; that
has to be checked before promising any exome gain.

## 4. Reproducing

```
python3 classify.py rbceq2-2.4.4/rbceq2/resources/db.tsv
```
where `classify.py` splits each `GRCh38` cell on `,`, parses `<pos>_<REF>_<ALT>` as an
indel of `|len(ALT)-len(REF)|`, `<pos>_<del|dup|ins>_<len>` as a word-form SV, and
`<pos>_ref` as a lane site; then buckets an allele by its largest event (50 bp and 1 kb
cut-offs) and calls it hybrid when it has two or more events of 1 kb or more. The system
name follows `rbceq2.core_logic.alleles.Allele.blood_group` (KLF1 becomes KLF).

## 5. What changed since the SPEC was written

- **rbceq2 2.4.4 (released 2026-09-16) handles haploid GT natively.** The SPEC's Q4 and
  our `bcftools +fixploidy` step were a workaround for a crash on haploid chrX/chrY calls.
  2.4.4 also says it keeps "chromosome-copy counts ... distinct", so `1`→`1|1` may now
  overstate dosage where 2.4.4 would have read hemizygosity correctly. Re-test with and
  without `+fixploidy` on a male sample when we bump the pin.
- **2.4.4 refuses tied SV evidence** with a named sample error
  (`SvMatcher.match/ambiguous_equal_best_sv_evidence`). The SV∩CNV overlap the SPEC's Q2
  left open (the same deletion arriving from both Manta and the CNV caller) may now fail
  the sample rather than pick one. Q2's recommendation to defer to rbceq2 no longer holds
  as written; the design PR revisits it.
- **2.4.4 keeps the selected SV event's GT, phase and FILTER together and logs the source
  event and db variant in debug output.** We already capture the debug log
  (`docs/rbceq2_debug_log/SPEC.md`); it can tell the QC which record a structural call
  rests on.
- **The QC stage now exists** (`FlagBloodGroupCallQc`), and it deliberately skips SV sites:
  a base's DP/GQ says nothing about a 21 kb deletion. See §6.
- **db notation has drifted.** RH tokens are now bare base counts (`_DEL_59419`) and the
  GE/PIGG/MAM deletions are spelled as explicit sequences, which `bg_db.py` already
  tolerates but the SPEC's §2.2 token counts no longer match.

## 6. QC for structural calls: the open design question

The current QC asks, for each defining SNV/indel site, "did the caller look here, and how
well?" using the gVCF's variant records and reference blocks. The same two questions have
different answers for the structural inputs:

| | Did the caller assess the region? | How confident is the event call? |
|---|---|---|
| **CNV VCF** | Yes, answerable: `DRAGEN:REF:` records span every assessed interval with `BC` (bin count) and `SM` (segment mean). A gene span covered by a REF record with adequate `BC` is the analogue of a reference block; a gap is the analogue of `NOCOV`. | `QUAL`, `FILTER` (`cnvLength`, `cnvQual`, `cnvCopyRatio`, `cnvBinSupportRatio`), `CN`, `SM`, `BC`. |
| **SV VCF (Manta)** | **Not answerable from the VCF.** Manta emits events only; absence is silence, the same trap as a gVCF hole. Assessability would have to come from the CRAM (depth and MAPQ0 fraction over the gene, e.g. mosdepth), which is a new input and a new job. | `QUAL`, `FILTER` (`MinQUAL`, `MaxDepth`, `MaxMQ0Frac`, `NoPairSupport`, `Ploidy`), `PR`/`SR` read counts, `CIPOS`/`CIEND`. |

So the sub-10 kb Manta-only alleles (GE, PIGG, MAM, ABO*O.16) can be *called* from the
merged VCF, but the QC cannot say whether their absence was assessed without a
CRAM-derived region metric, while the CNV-band alleles can get both answers from the CNV
VCF alone. How the QC stage should report structural calls is the design PR's question.

## 7. Cost side, briefly (inputs and compute only; the QC design is not costed here)

- Inputs: two more files per SG to resolve, and metamist analysis types for them if
  `dragen_align` does not already register them (unverified since the SPEC).
- Compute: the merge is bcftools on three small region-restricted files; negligible next
  to the existing conversion. The CRAM-derived assessability metric for Manta-band QC is
  the only new cost of note, and it reads the same ~50 regions the conversion already
  restricts to.
- Validation: no positive control among the five NA12878 replicates. Concordance across
  replicates is free; a true-positive needs a known U- or McLeod or Ge:-2 sample, or a
  spiked record.
