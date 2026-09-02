# Validation runs: post-hoc exome recall on mackenzie test exomes

End-to-end runs of `PosthocGenotypeOffTargetSites` and the merge on real data, one per capture
design. Recorded here because the design in [`SPEC.md`](./SPEC.md) rests on measurements, and
this is where they were checked against live cohorts rather than synthetic inputs.

## What was run

| | Twist | CREv2 |
|---|---|---|
| Batch | [1136502](https://batch.hail.populationgenomics.org.au/batches/1136502) | [1136522](https://batch.hail.populationgenomics.org.au/batches/1136522) |
| Cohort | COH13420 | COH13446 |
| Library | `TwistWES1VCGS1` (VCGS) | `AgilentCREv2WES` (NSW HP) |
| Samples | 10 exomes | 10 exomes |
| Jobs | 62/62 succeeded | 62/62 succeeded |
| Wall clock / cost | 8.9 min / $0.12 | 5.4 min / $0.09 |

Both on mackenzie at access-level test, image
`images-dev/popgen_rbceq2:feat-posthoc-exome-recall-0.1.0-8`, taking the DRAGEN 3.7.8 CRAM and
its matching recal gVCF from `gs://cpg-mackenzie-test/ica/dragen_3_7_8/`.

Each cohort is homogeneous in capture design, confirmed both from Metamist assay metadata and
empirically — see the design comparison below.

Two cohorts that look usable are not: COH12786 and COH12794 have sequencing groups in the main
Metamist project, so a test-level run cannot read them. A cohort spanning both projects is worse
than useless here, because the combine stage filters its expected set to sequencing groups that
have a gVCF, so the missing half is dropped with a warning rather than an error.

## The post-hoc caller is cheap

This was the main open question in the spec — streaming the CRAM over GCS NIO had only been
probed from a laptop, never timed in batch.

| | per sample |
|---|---|
| HaplotypeCaller wall clock | ~130 s (40 s of it traversal) |
| cost | $0.004 |
| regions traversed | 744 |
| reads read from the streamed CRAM | ~45,000 |

Localising an exome CRAM would have moved gigabytes per sample for the same 136kb of intervals.
Streaming stands.

## Recovery, per design

Over 1,625 assessable defining coordinates x 10 samples = 16,250 site-resolutions each:

| | Twist | CREv2 |
|---|---|---|
| holes per sample (no DRAGEN record) | 167–169 | 165–168 |
| **has a covering record, without recall** | 89.7% | 89.8% |
| **has a covering record, with recall** | **97.4%** | **97.2%** |
| and also clears DP/GQ thresholds | 83.4% → 89.6% | 80.2% → 86.7% |
| still no record from either caller | 428 (2.6%) | 450 (2.8%) |

The hole counts match the ~165 the 2026-08 coverage analysis predicted per design, from a
different set of exomes. The two designs gain almost exactly the same amount, 7.7 and 7.4
percentage points.

## The two designs miss different sites

This is the result that shows the recall is tracking real capture boundaries rather than
something incidental. Comparing the set of coordinates DRAGEN is silent at:

| Jaccard of hole sets | min | median |
|---|---|---|
| within Twist (45 pairs) | 0.99 | 0.99 |
| within CREv2 (45 pairs) | 0.97 | 0.99 |
| **across designs (100 pairs)** | **0.33** | **0.34** |

Of the coordinates silent in *every* sample of a design: 84 are shared blind spots of both
designs, 83 are Twist-only and 81 are CREv2-only. So roughly half of each design's blind spots
are unique to it, and the aggregate similarity in the table above conceals two quite different
site sets. A capture design is legible in this data, which is why validating one said nothing
about the other.

## Quality of the recovered genotypes

The question these runs had to answer is whether recovered sites are *good* calls or merely
present. They are good: a recovered site that passes is indistinguishable from a DRAGEN-called
one, on both designs.

| site class | Twist n | DP median | GQ median | CREv2 n | DP median | GQ median |
|---|---|---|---|---|---|---|
| primary, passing | 13,553 | 76 | 99 | 13,037 | 150 | 99 |
| **post-hoc, passing** | **1,010** | **73** | **99** | **1,046** | **141** | **99** |
| post-hoc, below threshold | 236 | 2 | 6 | 161 | 2 | 3 |

Recovered depth tracks each cohort's own baseline — 73 against 76 for Twist, 141 against 150
for the deeper CREv2 cohort — with identical median GQ of 99 in both. The sub-threshold sites
are genuinely marginal at median depth 2, and are correctly reported `LOWQ` with `src=` naming
the caller, not `POSTHOC`.

## Effect on the calls

Combined cohort QC tables, 10 samples x 51 systems = 510 cells each:

| category | Twist | CREv2 |
|---|---|---|
| PASS | 200 | 201 |
| POSTHOC (rests on a recovered site) | 159 | 157 |
| LOWQ | 81 | 82 |
| NOCOV | 40 | 40 |
| NA | 30 | 30 |

Around 31% of cells in both cohorts now rest on a recovered site. Before this change they would
have read `NOCOV`, meaning rbceq2 was calling those systems reference on nothing. 17 distinct
systems were recovered in at least one sample in each cohort; 15 in all ten Twist samples, 12 in
all ten CREv2 samples.

**FY is the headline**, being the original motivation. Both cohorts now carry a supported Duffy
call in all ten samples, resting on the GATA promoter sites that sit outside both captures and
previously came back empty — Twist 4x Fy(a+b-), 4x Fy(a+b+), 2x Fy(a-b+); CREv2 8x Fy(a+b+),
1x Fy(a+b-), 1x Fy(a-b+). Every one flagged `POSTHOC` so the provenance is visible.

`KLF` in the Twist cohort is a worked example of the severity order on real data: a recovered
site that passed plus two that failed, rendering `POSTHOC:...;LOWQ:...;LOWQ:...`. The provenance
flag does not mask the quality problem.

## What the runs confirmed about the design

**The sample-name check had to be a warning.** Four of ten Twist CRAMs carry a retired
sequencing-group ID in their read group while their gVCFs carry the current one, from the
upstream test-set reheadering bug. All four warned, relabelled and completed. The hard failure
the spec originally proposed would have killed 40% of a known-good cohort. No CREv2 sample
warned, so the affected inputs are a subset, not the norm — which is exactly why a hard failure
would have looked fine right up until it did not. See SPEC §10.

**`NA` is unrelated to the recall.** All 30 `NA` cells in each cohort are ABCC1, ATP11C and CD99
across ten samples. Those three have zero rows in the committed site-system map because every
one of their defining alleles is a structural variant, which one base of depth and GQ cannot
assess. The recall cannot help: the problem is not missing reads, it is that a single-base
measurement is the wrong instrument. See [`../rbceq2_cnv_sv/SPEC.md`](../rbceq2_cnv_sv/SPEC.md).

## Observations worth following up

**Zero-depth *primary* records.** 768 Twist site-resolutions had a DRAGEN record reporting
`DP=0`. Those sites are not holes — a record covers them — so the merge never offers post-hoc
data there. Most is real biology: 600 of the 768 are two samples showing ~300 of 329 RHD
coordinates at zero depth, a whole-gene RHD deletion, and 2 of 10 is in line with RhD-negative
population frequency. The residue is roughly 17 sites per sample in RHD, RHCE, C4A and C4B, some
of which may be capture edges the recall could reach if the hole rule also treated a zero-depth
primary record as silence. Worth measuring before acting on — it would change what `NOCOV` means
for the primary caller too.

**HPA remains untouched.** Per-sample logs report 51 of 85 site-map systems quality-flagged,
which looks worse than the cohort tables do; 35 of those are HPA, whose database coordinates are
GRCh37 and were never lifted. rbceq2 emits no columns for them, so they do not appear in the
cohort tables at all. Unrelated to this work.
