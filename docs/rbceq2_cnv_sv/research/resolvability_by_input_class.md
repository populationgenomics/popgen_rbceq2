# What the SNV/indel gVCF resolves today, and what a merged SNV/indel/SV/CNV VCF adds

**Date:** 2026-09-17 · **Database:** rbceq2 2.4.4 (bundled db 2.5.1), `GRCh38` column ·
**Method:** every allele's coordinate tokens classified by the largest event they name.
Script and raw output are reproducible from `db.tsv` alone (see §4).

The gVCF already resolves 1300 of the 1388 non-RH alleles (93.7%), and every allele in 62 of
the 86 non-RH systems. A merged VCF would add 67 alleles defined by one large structural
event, in 20 systems. One of them, ABCC1, has the wrong coordinates in the database, so the
real gain is 66 alleles in 19 systems. ATP11C, CD99 and XG go from blind to callable, and
null alleles become visible in 16 more systems. Alleles under 10 kb can come only from the SV
VCF, which is why the SPEC merges all three files. Hybrid alleles are not gained.

RH (RHD, RHCE) is listed but excluded from the gain. The maintainer's advice with the 2.4.4
release is that DRAGEN SV/CNV does not reliably detect the RH hybrids, so `--RH` stays off.

**Callers.** The SV VCF is written by the DRAGEN 3.7.8 SV caller, which Illumina describes as
integrating and extending Manta; its record IDs keep Manta's prefix (`MantaDEL:`,
`MantaINS:`). This note says "SV caller" and "SV VCF", not "Manta", because the caller is
DRAGEN's and claims about Manta's internals may not hold for it. The CNV VCF is written by
DRAGEN's bin-based CNV caller.

## 1. Allele classes

| Class | Definition | Where it can come from |
|---|---|---|
| **A small** | Every token is a SNV, a `_ref` lane site, or an indel under 50 bp | DRAGEN SNV/indel gVCF, which is what we run today |
| **B mid indel** | Largest event 50 bp to 1 kb | Sometimes the gVCF (DRAGEN calls indels up to a few hundred bp), otherwise the SV VCF. Assume not reliably called today. |
| **C large, single** | One event of 1 kb or more: a deletion, duplication or insertion, not paired | 1 kb to 10 kb: SV VCF only (DRAGEN CNV filters these as `cnvLength`). 10 kb and up: CNV VCF, and usually the SV VCF too. This is the gain. |
| **D hybrid** | Two or more large events, a deletion paired with an insertion or duplication | Gene-conversion products. Long read only. |

rbceq2 reads classes B, C and D from the same VCF through `SvReader`, matching on `SVTYPE`
and fuzzy position and length. 2.4.4 did not change `SvMatcher.compatible`: `SVTYPE` must
equal the db type, so DRAGEN's `SVTYPE=CNV` still fails against `DEL` and `DUP`, and the
SPEC's CNV rewrite (§5.4) is still required.

## 2. Per-system counts (alleles with a GRCh38 coordinate)

Only systems with at least one allele outside class A are listed. Every other system,
including all 35 HPA systems, ABCB6, CROM, DO, FUT1/2/3, KEL, KLF, KN, LW, SC, VEL and YT,
is entirely class A, so the gVCF already resolves everything the db defines for it.

| System | Total | A small (now) | B mid indel | C large single (gain) | D hybrid | Notes |
|---|---|---|---|---|---|---|
| ABCC1 | 1 | 0 | 0 | 1 | 0 | 21 kb del. Uncallable on any input: the row carries ABCC4's chr13 coordinates on chr16, an upstream db error still present in 2.4.4 (README, "Structural-variant entries"). Counted here, gained nowhere. |
| ATP11C | 1 | 0 | 0 | 1 | 0 | 219 kb whole-gene del, non-PAR chrX. Blind today. |
| CD99 | 3 | 1 | 0 | 2 | 0 | 11 kb and 20 kb dels, PAR1. Blind in practice today. |
| XG | 3 | 1 | 0 | 2 | 0 | 32 kb and 114 kb dels, PAR1 boundary. |
| XK | 67 | 44 | 3 | 20 | 0 | 19 alleles share one 53 kb token (McLeod cytogenetic dels); plus 8 kb and 3.5 kb dels. Non-PAR chrX. |
| GYPB | 40 | 33 | 1 | 6 | 0 | Whole-GYPB dels, 19 kb to 224 kb: the U- / S-s-U- alleles. Paralog locus. |
| GYPA | 26 | 24 | 0 | 2 | 0 | 101 kb and 119 kb dels. Paralog locus. |
| GYP | 18 | 2 | 4 | 8 | 4 | The MNS hybrid system. 3 of the C are the 122 kb del of GYP*201/202/203 (with SNVs); 5 are 1.8 kb to 3.6 kb explicit-sequence dels. Paralog locus. |
| GE | 20 | 10 | 2 | 8 | 0 | Seven 3.6 kb exon dels (Ge:-2 / Ge:-3, Melanesian) and one 2.3 kb del. All under 10 kb: SV VCF only. |
| A4GALT | 43 | 40 | 0 | 3 | 0 | 21 kb to 33 kb dels (P1PK null). |
| PIGG | 8 | 4 | 1 | 3 | 0 | 4 kb to 6.3 kb dels. SV VCF only. |
| MAM | 5 | 1 | 2 | 2 | 0 | 8.5 kb dels. SV VCF only. |
| LU | 38 | 36 | 0 | 2 | 0 | 1 kb and 27 kb dels. |
| ABO | 207 | 206 | 0 | 1 | 0 | ABO*O.16, 6 kb del. SV VCF only. |
| ABCG2 | 32 | 31 | 0 | 1 | 0 | 1.8 kb delins. SV VCF only. |
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

Of the 86 non-RH systems, 62 are entirely class A (35 HPA plus 27 others) and 20 have at
least one class C allele. Three (CO, FY, JK) have only a mid-size indel beyond class A, and
one (DI) has an allele with no coordinate. 62 + 20 + 3 + 1 = 86.

## 3. Reading the numbers

**What we resolve now.** For the 62 all-class-A systems, the gVCF plus the existing QC is the
whole story, and a merged VCF changes nothing.

**What a merged VCF gains.** The 66 callable class C alleles, plus some of the 16 mid-size
indels, depending on which the gVCF already calls. In system terms:

- **ATP11C, CD99 and XG go from blind to callable.** ATP11C has no small-variant allele, and
  CD99 and XG have one each, with their null alleles all deletions. Today the QC reports
  these `NA` and rbceq2 emits the reference phenotype for every sample. ABCC1 is also `NA`
  today, and stays uncallable because of its db coordinates.
- **Null (antigen-negative) alleles become visible in 16 systems** that are otherwise
  callable. These matter most for rare-donor work: whole-GYPB deletion (U-), XK deletion
  (McLeod), A4GALT deletion (p phenotype), the GE exon deletions (Ge:-2, Ge:-3; common in
  Papua New Guinea and Melanesia, and relevant to our cohorts), C4A deletion (Ch-), and the
  LU, RHAG, CTL2, GCNT2 and ABCC4 nulls.
- **Not gained:** the 4 GYP hybrids and the 18 RH hybrids (class D). The paralog pairs,
  GYPA/GYPB/GYPE and C4A/C4B, are also at risk. A whole-gene deletion there is a dosage loss
  the bin-based CNV caller can see even where the breakpoints are unmappable, but that needs
  verifying on a known positive, not assuming.

**Two size bands, two sources.** 22 of the 67 class C alleles are under 10 kb, which DRAGEN
CNV filters as `cnvLength`, and 45 are 10 kb or more. The small ones (GE, PIGG, MAM,
ABO*O.16, XK*N.05, ABCG2 and the small GYP deletions) therefore depend on the SV VCF. The
large ones, mostly whole-gene deletions (GYPB, XK, ATP11C, A4GALT, CTL2, GCNT2, ABCC4, RHAG,
LU*02N.06, CD99, XG, C4A), can come from either file. This is why the SPEC merges all three
files rather than adding only the CNV VCF.

**Exomes.** Everything above assumes genome DRAGEN outputs. On 2026-09-18 we checked the
buckets: every one of the 11,917 production mackenzie exomes (`cpg-mackenzie-main`, DRAGEN
3.7.8) has an `sv.vcf.gz`, a `cnv.vcf.gz`, a `ploidy_estimation_metrics.csv` and a
`wgs_coverage_metrics.csv`.

The exome CNV VCF is a different file from the genome one. It is called per capture target
against a panel of 100 normals and holds events only, with no `DRAGEN:REF:` records. It has
no `cnvLength` filter, so sub-10 kb events can PASS, and its `BC` counts capture targets.
One production exome held 1,131 CNV event records, 30 of them in the blood-group regions,
and 65 SV records, 4 in the regions.

The 10 exomes of the earlier test run (`cpg-mackenzie-test`) have the SV VCF but no CNV VCF.
Their `cnv_metrics.csv` stops after "Number of target intervals", so that run's caller
counted reads but never segmented.

So production exomes carry both structural inputs, and both size bands are in reach in
principle. The exome callers' sensitivity in the blood-group regions is unmeasured. The SPEC
(§3, §4, §5.1) puts exomes in scope on those terms, and reads whether a CNV VCF is expected
from `cnv_metrics.csv` rather than from the sequencing type.

## 4. Reproducing

```
python3 classify.py rbceq2-2.4.4/rbceq2/resources/db.tsv
```

`classify.py` splits each `GRCh38` cell on `,`. It parses `<pos>_<REF>_<ALT>` as an indel of
`|len(ALT)-len(REF)|`, `<pos>_<del|dup|ins>_<len>` as a word-form SV, and `<pos>_ref` as a
lane site. It then buckets each allele by its largest event, with cut-offs at 50 bp and
1 kb, and calls it a hybrid when it has two or more events of 1 kb or more. The system name
follows `rbceq2.core_logic.alleles.Allele.blood_group` (KLF1 becomes KLF).

## 5. What changed since the SPEC was written

- **rbceq2 2.4.4 (released 2026-09-16) handles haploid GT natively.** The SPEC's Q4 and our
  `bcftools +fixploidy` step were a workaround for a crash on haploid chrX/chrY calls. 2.4.4
  also says it keeps "chromosome-copy counts ... distinct", so expanding `1` to `1|1` may
  overstate dosage where 2.4.4 would read hemizygosity correctly. The pin bump re-tested
  this: `popgen_rbceq2#20` (merged 2026-09-22) removed `+fixploidy`, since 2.4.4 renders a
  hemizygous null as `XK*01.02/-` where the expansion had made it a homozygote. Phenotypes
  were unchanged in the synthetic check.
- **2.4.4 refuses tied SV evidence** with a named sample error
  (`SvMatcher.match/ambiguous_equal_best_sv_evidence`). The SPEC's Q2 left open what happens
  when the same deletion arrives from both callers; under 2.4.4 it may fail the sample rather
  than pick one. So Q2's recommendation to defer to rbceq2 no longer holds as written, and
  the design PR (`popgen_rbceq2#16`) revisits it.
- **2.4.4 keeps the selected SV event's GT, phase and FILTER together, and logs the source
  event and db variant in debug output.** We already capture the debug log
  (`docs/rbceq2_debug_log/SPEC.md`), so it can tell the QC which record a structural call
  rests on.
- **The QC stage now exists** (`FlagBloodGroupCallQc`), and it deliberately skips SV sites,
  since a base's DP and GQ say nothing about a 21 kb deletion. See §6.
- **The db notation has drifted.** RH tokens are now bare base counts (`_DEL_59419`), and
  the GE, PIGG and MAM deletions are spelled as explicit sequences. `bg_db.py` already
  tolerates both, but the SPEC's §2.2 token counts no longer match.

## 6. QC for structural calls: the open design question

The current QC asks two questions of each defining SNV/indel site: did the caller look here,
and how well? It answers them from the gVCF's variant records and reference blocks. For the
structural inputs, the two files give different answers:

| | Did the caller assess the region? | How confident is the event call? |
|---|---|---|
| **CNV VCF** | Yes. `DRAGEN:REF:` records span every assessed interval with `BC` (bin count) and `SM` (segment mean). A gene span covered by a REF record with adequate `BC` is the analogue of a reference block; a gap is the analogue of `NOCOV`. | `QUAL`, `FILTER` (`cnvLength`, `cnvQual`, `cnvCopyRatio`, `cnvBinSupportRatio`), `CN`, `SM`, `BC`. |
| **SV VCF (DRAGEN SV caller)** | Not from the VCF. The file holds event records only, with no reference or assessed-interval records (observed across 150 genomes and one exome), so absence is silence, the same trap as a gVCF hole. Assessability would have to come from the CRAM (depth and MAPQ0 fraction over the gene, for example from mosdepth), which is a new input and a new job. | `QUAL`, `FILTER` (`MinQUAL`, `MaxDepth`, `MaxMQ0Frac`, `NoPairSupport`, `Ploidy`, `MinGQ`, `HomRef`, `SampleFT`, as declared in the 3.7.8 header), `PR`/`SR` read counts, `CIPOS`/`CIEND`. |

So the sub-10 kb alleles that only the SV VCF carries (GE, PIGG, MAM, ABO*O.16) can be
called from the merged VCF, but without a CRAM-derived region metric the QC cannot say
whether their absence was assessed. The CNV-band alleles get both answers from the CNV VCF
alone. How the QC stage should report structural calls is for the design PR to decide.

## 7. Cost, briefly (inputs and compute only; the QC design is not costed here)

- **Inputs:** two more files per SG to resolve. `dragen_align` does not register them in
  metamist and will not before its Nextflow refactor (review, 2026-09-18), so the design PR
  derives their paths from the registered gVCF's.
- **Compute:** the merge is bcftools on three small region-restricted files, negligible next
  to the existing conversion. The only new cost of note is the CRAM-derived assessability
  metric for the SV-band QC, and it reads the same ~50 regions the conversion already
  restricts to.
- **Validation:** there is no positive control among the five NA12878 replicates. Concordance
  across replicates is free; a true positive needs a known U-, McLeod or Ge:-2 sample, or a
  spiked record.
