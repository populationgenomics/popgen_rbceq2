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

## The capture design bounds the fill, and the runs say which BED is the design

The merge fills only defining sites outside the cohort's capture design (README, "Merging the
post-hoc calls"). That gate needs a BED, and these runs were used to pick it and to check that
picking wrongly is detectable. For every defining coordinate in all 20 samples, whether it sits
inside the design BED and whether DRAGEN has a record there:

| design BED | in design, DRAGEN record | in design, silent | outside, DRAGEN record | outside, silent |
|---|---|---|---|---|
| Twist VCGS covered targets | 14,576 | 4 | **0** | 1,670 |
| CREv2 `S30409818_Covered` (probe footprint) | 14,593 | 7 | **0** | 1,650 |
| CREv2 `S30409818_Regions` (targets) | 14,284 | 966 | **309** | 691 |

DRAGEN's emitted footprint is exactly the Twist covered-targets file and exactly the CREv2
**Covered** file. Zero off-design records in 32,500 site-resolutions is what makes the job's
footprint check safe to fail on, and it is why the `exome_design_bed` for an Agilent design must
name `Covered` and not `Regions`. Configuring `Regions` would leave 907 of that cohort's 1,207
recoveries gated off, which no log would have shown; instead the run now stops and names the key.

The gate costs 4 Twist and 7 CREv2 site-resolutions, in C4B, RHD, RHCE and A4GALT. Those are
sites the capture targeted and DRAGEN was still silent at. They now read `NOCOV` rather than
`POSTHOC`, which is the intent: a targeted site with no DRAGEN record is a fact about that
DRAGEN run, and a second caller answering it hides the fact. Everything else in this document is
unaffected, the recovery figures above included.

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

Combined cohort QC tables, 10 samples x 51 systems = 510 cells each. Recomputed from each
sample's own defining-sites extract under the current flag scheme, in which a flag name carries
both a severity and, joined with `+`, `POSTHOC` when the post-hoc caller supplied the site:

| category | Twist | CREv2 |
|---|---|---|
| PASS | 200 | 201 |
| `POSTHOC` only, every site passing | 143 | 111 |
| a `LOWQ+POSTHOC` site | 64 | 83 |
| `POSTHOC` beside a `NOCOV` site | 30 | 41 |
| **any POSTHOC token** | **237** | **235** |
| LOWQ, no recovery involved | 20 | 19 |
| NOCOV, no recovery involved | 23 | 25 |
| NA | 30 | 30 |

An earlier version of this table reported a single `POSTHOC` row of 159 and 157. Those figures
cannot be reproduced from the stored cohort tables under any flag precedence, so they have been
replaced by the breakdown above rather than carried forward. The rows here were each recomputed
from the run artefacts.

**46% of cells rest on a recovered site** in both cohorts, 237 and 235 of 510. Before this
change they would have read `NOCOV`, meaning rbceq2 was calling those systems reference on
nothing. Of those, 143 and 111 are clean recoveries with every listed site passing; the rest
are recoveries that also carry a quality problem, which the cell now states alongside the
provenance instead of in place of it.

That last point is the whole reason the token left the severity ladder. Reading the runs' own
cohort tables, 35 Twist and 25 CREv2 cells rest on a recovered site (`src=` is present) while
where no flag would carry `POSTHOC` if severity outranked it, because every recovered site in
them is sub-threshold and would report `LOWQ` alone. That is 60 of the 472 reliant cells in the two tables.

The loss is larger over the systems the QC job assesses but rbceq2 emits no column for, which
reach the run log rather than the tables: there, 130 Twist and 104 CREv2 reliant systems were
silent, 234 of 681. Either way the mark is now on every one of them.

17 distinct systems were recovered in at least one sample in each cohort; 15 in all ten Twist
samples, 12 in all ten CREv2 samples.

**FY is the headline**, being the original motivation. Both cohorts now carry a supported Duffy
call in all ten samples, resting on the GATA promoter sites that sit outside both captures and
previously came back empty — Twist 4x Fy(a+b-), 4x Fy(a+b+), 2x Fy(a-b+); CREv2 8x Fy(a+b+),
1x Fy(a+b-), 1x Fy(a-b+). Every one carries the `POSTHOC` token so the provenance is visible.

`KLF` in the Twist cohort is a worked example on real data: a recovered site that passed plus
two that failed, rendering `POSTHOC:...;LOWQ+POSTHOC:...;LOWQ+POSTHOC:...`. Neither finding
masks the other.

## What the runs confirmed about the design

**Four Twist CRAMs disagree with their gVCF on the sample name.** They carry a retired
sequencing-group ID in their read group while their gVCFs carry the current one, from the
upstream test-set reheadering bug; no CREv2 sample does. These runs were made while the merge
downgraded that to a warning and relabelled, so all four completed, and their measurements are
sound: the post-hoc calls came from the right CRAM, only its header name was stale.

That downgrade has since been reverted, so **the Twist cohort can no longer be run as it
stands** — 4 of its 10 sequencing groups now fail the merge. The measurements above stand as a
record; reproducing them needs the CRAM headers fixed upstream. CREv2 is unaffected and
re-runnable. Sampling the test bucket shows the quirk is confined to an older block of CRAMs,
5 of the first 12, with newer additions consistent. See SPEC §10 for why accommodating it in
the pipeline was rejected.

**`NA` is unrelated to the recall.** All 30 `NA` cells in each cohort are ABCC1, ATP11C and CD99
across ten samples. Those three have zero rows in the committed site-system map because every
one of their defining alleles is a structural variant, which one base of depth and GQ cannot
assess. The recall cannot help: the problem is not missing reads, it is that a single-base
measurement is the wrong instrument. See [`../rbceq2_cnv_sv/SPEC.md`](../rbceq2_cnv_sv/SPEC.md).

## Observations worth following up

**Zero-depth *primary* records cannot be recovered.** Both callers were asked the same
question, at the same coordinates, on the same reads. 768 Twist and 823 CREv2 site-resolutions
had a DRAGEN record reporting `DP=0`, `GQ=0`, `GT=./.`. Those sites are not holes, since a
record covers them, so the merge never offers post-hoc data there, and the question was whether
the hole rule should treat such a record as silence. Measured: it should not.

| | Twist | CREv2 |
|---|---|---|
| zero-depth primary site-resolutions | 768 | 823 |
| of which whole-gene RHD deletion (2 and 3 samples, ~270-300 of 329 RHD coords each) | 602 | 790 |
| residue | 168 | 44 |
| residue where the post-hoc caller also reports `DP=0` | 165 | 44 |
| residue where it reports DP 1-9 | 3 | 0 |
| **residue that would clear DP>=10, GQ>=20 if recalled** | **0** | **0** |

The residue sits in RHD (Twist only, three samples with 30-54 zero-depth RHD coordinates, the
signature of a partial deletion or RHD-CE-D hybrid), RHCE, C4A and C4B. RhD-negative at 2 and 3
of 10 is population-plausible; C4A/C4B null alleles are common. One C4B coordinate is zero-depth
in all 20 samples on both designs, and the CRAM explains it: the reads are there, 150 across the
DRAGEN zero-depth block, and every one has mapping quality 0, because C4A and C4B are
near-identical paralogs. Both callers discard those reads, and no short-read caller will genotype
that site. So the zero-depth records are what they say, no usable reads, and changing the hole
rule would only have relabelled `LOWQ(DP=0)` as `NOCOV` at sites where the recall also has
nothing. Measured on the two validation runs above, from the `defining_sites` extracts and the
post-hoc gVCFs, with the QC stage's own `resolve_coverage`.

**HPA remains untouched.** Per-sample logs report 51 of 85 site-map systems quality-flagged,
which looks worse than the cohort tables do; 35 of those are HPA, whose database coordinates are
GRCh37 and were never lifted. rbceq2 emits no columns for them, so they do not appear in the
cohort tables at all. Unrelated to this work.
