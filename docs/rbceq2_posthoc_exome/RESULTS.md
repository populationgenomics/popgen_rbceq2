# Validation run: post-hoc exome recall on mackenzie test exomes

First end-to-end run of `PosthocGenotypeOffTargetSites` and the merge, on real data.
Recorded here because the design in [`SPEC.md`](./SPEC.md) rests on measurements, and this is
where they were checked against a live cohort rather than against synthetic inputs.

## What was run

| | |
|---|---|
| Batch | [1136502](https://batch.hail.populationgenomics.org.au/batches/1136502) |
| Date | 2026-09-02 |
| Cohort | COH13420, `mackenzie-test-exome-VCGS-twist-20260609-dragen`, 10 exomes |
| Dataset / access level | mackenzie, test |
| Image | `images-dev/popgen_rbceq2:feat-posthoc-exome-recall-0.1.0-8` |
| Inputs | DRAGEN 3.7.8 CRAM + matching recal gVCF, both under `gs://cpg-mackenzie-test/ica/dragen_3_7_8/` |

COH13420 is the only usable cohort at test level. The two 20260528 dragen cohorts look
equivalent but their sequencing groups sit in the main Metamist project, so a test-level run
cannot read them.

## Outcome

62 of 62 jobs succeeded, no failures, 8.9 minutes wall clock, $0.12 for the cohort.

The post-hoc caller is cheap, which was the main open question in the spec — streaming the CRAM
over GCS NIO had only ever been probed from a laptop, never timed in batch:

| | per sample |
|---|---|
| HaplotypeCaller wall clock | 132 s (40 s of it traversal) |
| cost | $0.004 |
| regions traversed | 744 |
| reads read from the streamed CRAM | ~45,000 |

Localising an exome CRAM instead would have moved gigabytes per sample for the same 136kb of
intervals. Streaming stands.

## Recovery

The merge found **167 to 169 defining sites per sample with no DRAGEN record at all**, against
the ~165 the 2026-08 coverage analysis predicted from a different set of exomes. That is an
independent check that the hole-finding measures what the analysis measured.

Not every hole can be filled, and that is the point of dropping zero-depth records: a site with
no reads stays `NOCOV` rather than being dressed up as low-quality data.

Pooled over the 10 samples, across 1,625 assessable defining coordinates each:

| how the site resolved | site-resolutions | share |
|---|---|---|
| primary (DRAGEN) | 14,576 | 89.7% |
| post-hoc, passing thresholds | 1,010 | 6.2% |
| post-hoc, below thresholds | 236 | 1.5% |
| no record from either caller | 428 | 2.6% |

So the recall lifted per-sample assessable coordinates from 89.7% to 97.4%, close to the ~98%
the spec projected.

## Quality of the recovered genotypes

The question this run had to answer is whether recovered sites are *good* calls or merely
present. They are good: a recovered site that passes is statistically indistinguishable from a
DRAGEN-called one.

| site class | n | DP median | DP p25 | DP p75 | GQ median | MIN_DP median |
|---|---|---|---|---|---|---|
| primary, passing | 13,553 | 76 | 54 | 101 | 99 | 52 |
| **post-hoc, passing** | **1,010** | **73** | **49** | **93** | **99** | **36** |
| post-hoc, below threshold | 236 | 2 | 1 | 4 | 6 | 2 |

Median depth 73 against 76, and identical median GQ of 99. The lower `MIN_DP` (36 against 52)
is expected and is exactly why the flag reports it: these blocks sit at capture edges, so the
shallowest base in the band is further below the band median than it is in well-targeted
sequence. The site is still judged on `DP`, the band median.

The failing 236 are genuinely marginal — median depth 2, median GQ 6 — and are correctly
reported `LOWQ` with `src=` naming the caller, not `POSTHOC`.

## Effect on the calls

Combined cohort QC table, 10 samples x 51 systems = 510 cells:

| category | cells |
|---|---|
| PASS | 200 |
| POSTHOC (rests on a recovered site) | 159 |
| LOWQ | 81 |
| NOCOV | 40 |
| NA | 30 |

The 159 `POSTHOC` cells are systems whose defining sites had no DRAGEN record at all. Before
this change they would have read `NOCOV`, meaning rbceq2 was calling them reference on nothing.
Recovery is consistent per sample: 16 systems in nine samples, 15 in the tenth.

Fifteen systems were recovered in all ten samples — ABCB6, ABCG2, AUG, CO, DO, FY, GIL, GYP,
GYPA, GYPB, JK, KEL, KLF, RHAG, SID — plus ABO in 8 of 10 and PIGG in 1.

**FY is the headline**, being the original motivation. All ten samples now carry a supported
Duffy call resting on the GATA promoter sites that sit outside the capture and previously came
back empty: four Fy(a+b-), four Fy(a+b+), two Fy(a-b+), every one flagged `POSTHOC` so the
provenance is visible rather than implied.

`KLF` is worth reading as a worked example of the severity order on real data. Its cell holds a
recovered site that passed and two more that failed, so it reads
`POSTHOC:...;LOWQ:...;LOWQ:...` — the provenance flag does not mask the quality problem.

## Two things this run confirmed about the design

**The sample-name check had to be a warning.** Four of the ten CRAMs carry a retired
sequencing-group ID in their read group while their gVCFs carry the current one, from the
upstream test-set reheadering bug. All four warned, relabelled and completed. The hard failure
the spec originally proposed would have killed 40% of a known-good cohort. See SPEC §10.

**`NA` is unrelated to the recall.** All 30 `NA` cells are ABCC1, ATP11C and CD99 across the
ten samples. Those three have zero rows in the committed site-system map because every one of
their defining alleles is a structural variant, which one base of depth and GQ cannot assess.
The recall cannot help: the problem is not missing reads, it is that a single-base measurement
is the wrong instrument. See [`../rbceq2_cnv_sv/SPEC.md`](../rbceq2_cnv_sv/SPEC.md).

## Observations worth following up

**Zero-depth *primary* records.** 768 site-resolutions across the cohort had a DRAGEN record
reporting `DP=0`. Those sites are not holes — a record covers them — so the merge never offers
post-hoc data there. Most of that is real biology rather than a gap: 600 of the 768 are two
samples showing ~300 of 329 RHD coordinates at zero depth, which is a whole-gene RHD deletion,
and 2 of 10 is in line with RhD-negative population frequency. The residue is roughly 17 sites
per sample in RHD, RHCE, C4A and C4B, some of which may be capture edges the recall could reach
if the hole rule also treated a zero-depth primary record as silence. Worth measuring before
acting on — it would change what `NOCOV` means for the primary caller too.

**HPA remains untouched.** Per-sample logs report 51 of 85 site-map systems quality-flagged,
which looks worse than the cohort table does; 35 of those are HPA, whose database coordinates
are GRCh37 and were never lifted. rbceq2 emits no columns for them, so they do not appear in
the cohort table at all. Unrelated to this work.
