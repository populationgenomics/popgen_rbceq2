# Design Spec: post-hoc genotyping of off-target defining sites in exomes

**Status:** IMPLEMENTED and validated on a real cohort — see
[`RESULTS.md`](./RESULTS.md) for the run and its measurements. Kept as the record of why.
**Area:** rbceq2 blood-group genotyping pipeline (`FilterAndConvertGvcfsForRbceq2` →
`GenotypeBloodGroupsWithRbceq2` → `FlagBloodGroupCallQc`), in `popgen_rbceq2`.
**Author:** Joshua Schmidt · **Reviewers:** Michael Harper (PR #14); agent-assisted reviews
2026-09-03 and 2026-09-04, findings and fixes recorded in PR #14.

Design agreed in outline 2026-08-14, from the mackenzie exome coverage analysis. This spec
maps that agreement onto the pipeline as it stands today.

---

## 1. Problem

Exome DRAGEN gVCFs are called with a capture-target BED and no `--vc-target-bed-padding`,
so gVCF emission hard-stops at the BED edges. Any blood-group defining coordinate outside
the capture footprint gets **no record at all** — not a low-quality one, none.

That silence is dangerous twice over:

- rbceq2 reads a defining site absent from its input as a *confident homozygous reference
  call* (see PRODUCT.md and `FlagBloodGroupCallQc`'s docstring).
- Our QC flags it `NOCOV`, so the call is at least visibly unsupported — but the system
  stays unassessable, and for exomes that is a lot of systems.

Measured on 20 DRAGEN 3.7.8 exomes (5 Twist VCGS, 15 Agilent CREv2):

- ~165 of 1,599 non-HPA defining coordinates are off-target per capture design.
- ~105 of those sit 1–100 bp from a target edge and have 40–110× MAPQ≥20 depth in the
  CRAMs (measured by streaming). That includes the FY GATA Duffy-null promoter sites at
  34–210×, which matter directly for the rare-blood-group use case.
- Sites more than ~250 bp off both capture designs have ~0× and are unrecoverable.

The reads exist; only the caller stopped early. Post-hoc calling of the off-target sites
from the CRAMs lifts assessable non-HPA coordinates from 91.7% to ~98%. Re-running DRAGEN
over the cohort is not affordable.

Out of scope: the 35 HPA systems (their db coordinates are GRCh37 and were never lifted —
a separate problem), and SV-defined alleles (see `docs/rbceq2_cnv_sv/SPEC.md`).

## 2. Design overview

One new per-sequencing-group stage that genotypes the defining-site neighbourhoods
directly from the CRAM, and a merge inside the existing conversion stage that keeps the
new calls **only where the DRAGEN gVCF is silent**. Downstream, rbceq2 sees one merged
VCF, and the QC flags any call resting on a post-hoc record so it can never pass silently.

```
sequencing_group.cram ──> PosthocGenotypeOffTargetSites ──> <sg>.posthoc.g.vcf.gz
                                                                  │
exome_design_bed ──> resources/bg_off_design_sites.<design>.bed (committed, once per design)
                                                                  │
sequencing_group.gvcf ──> FilterAndConvertGvcfsForRbceq2 <────────┘
                             (merge: DRAGEN wins; only off-design holes are filled;
                              supplement tagged INFO/POSTHOC)
                                    │                    │
                                converted VCF        DP/GQ extract (+ posthoc column)
                                    │                    │
                          GenotypeBloodGroupsWithRbceq2  FlagBloodGroupCallQc <── off-design sites
                                                         (POSTHOC flag / src=gatk-hc-4.6.2.0)
```

## 3. The new stage: `PosthocGenotypeOffTargetSites`

Per sequencing group. **Exome only**: `expected_outputs` returns None (skip, not fail)
when `sequencing_group.sequencing_type != 'exome'`, or the sequencing group has no CRAM,
or no gVCF (no primary calls means nothing to supplement).

- **Caller:** GATK HaplotypeCaller in GVCF mode with `--dragen-mode`, to match the DRAGEN
  3.7.8 primary calls as closely as a re-call can. No DRAGstr table initially, and no
  DRAGMAP — the reads are used as aligned, never re-aligned.
- **Streaming, not localising:** the CRAM is read from `gs://` via GATK's NIO support
  (pass the `gs://` path directly), restricted with `-L` to a committed padded
  defining-sites BED. Total interval span is ~1,600 sites × ~500 bp, well under 1 Mb, so
  the job reads megabytes of CRAM rather than localising 15–20 GB.
- **Reference:** the exact reference the CRAMs were aligned to — the DRAGEN masked
  assembly38 in the references bucket. CRAM decoding requires the aligner's reference;
  Broad's standard fasta is not it. Path comes from `references.broad.ref_fasta`, which
  already names that assembly (§9, question 1).
- **Image:** `image_path('gatk', ...)` — 4.6.2.0-2 is current in cpg-common.
- **Output:** `<sg>.posthoc.g.vcf.gz` + `.tbi`, under the tmp category like the
  conversion stage's intermediates (§9, question 3: nothing downstream reads it).

### New committed resource

`bg_defining_sites_padded.<genome>.bed`: each defining coordinate ±250 bp, merged.
Written by `scripts/gen_bg_resources.py` from the same single parse as the other three
resources, so the padded BED can never drift from the sites the QC assesses. The 250 bp
comes from the coverage analysis: beyond it there is nothing to recover. The padding is
baked into the committed file, not config — changing it means regenerating and committing,
the same rule as every other resource.

## 4. The merge, inside `FilterAndConvertGvcfsForRbceq2`

For exome sequencing groups the conversion job gains a merge step ahead of the existing
extract-and-convert. Genome runs are untouched.

**Judge silence empirically; bound the fill by the capture design.** A post-hoc record is
kept only at a coordinate that satisfies both tests: the DRAGEN gVCF has no record covering
it (reference blocks count as covering), *and* it lies outside the cohort's capture design.

The first test is per sample and reads no metadata, which is what makes the merge
self-calibrating: a design BED that lies about the real footprint still cannot cause a
DRAGEN call to be overwritten or a hole to be invented. DRAGEN wins wherever both speak.

The second test was added after the validation runs, at the author's direction, and only
ever *narrows* what the first test offers. Its argument is that "DRAGEN was silent" and
"the capture never targeted this site" are different findings that were being conflated.
A site the design did target and DRAGEN still said nothing about is a fact about that
sample's DRAGEN run; answering it with a second caller hides it. Such a site is now left
as `NOCOV`. Measured on the validation cohorts this is 4 of 1,674 Twist recoveries and 7 of
1,657 CREv2 ones, in C4B, RHD, RHCE and A4GALT — a small set, and the point is which
question each flag answers rather than the count.

Leaving it as `NOCOV` takes two steps, because a post-hoc record kept for an off-design
hole is kept whole and can reach an in-design hole beside it. The merge's `COVERED` mark
is drawn from every site the merge may not fill, the in-design holes included, so a
post-hoc *variant* reaching one is dropped. A post-hoc *reference block* reaching one is
kept for the hole it was selected for, and is then the only record at the in-design hole;
`FlagBloodGroupCallQc` reads the same off-design BED and counts a post-hoc record only at
a site in it, so the QC reports that hole `NOCOV` (§5).

The design is named by `workflow.exome_design_bed`, spelled as the design's key in the
references repo. It is deliberately not defaulted: every default is the wrong design for some
cohort, and being wrong is silent rather than fatal, so an exome run that does not name one
fails while the graph is built.

The pipeline never opens the design. Its subtraction from the defining sites depends only on
two fixed files, so it is a **committed resource**, `resources/bg_off_design_sites.<design>.
<genome>.bed`, written once per design by `scripts/gen_off_design_sites.py` in
`bedtools intersect -v`, a standard tool for exactly this operation, rather than the
hand-written awk containment loop it replaced (§10). `off_design.resource_path` resolves the
file for the configured key at graph build, and refuses a key with no committed file, or a
file whose manifest row was built from other defining sites than the shipped ones (§9.5). The
design key is also a segment of every exome output path, directly under the release (§7): the
release segment cannot see a config change, and every output from the conversion onward is
built from the holes the design leaves, so a repointed key must not reuse any of them.

**The configured design is checked against the gVCF.** DRAGEN emits reference blocks over
exactly the target BED with no padding, so a DRAGEN reference block reaching a defining
site outside the configured design means the file is not the one the gVCF was called
against, and the job fails naming the key. The check reads blocks only: a variant is
anchored inside the target but its REF can run past the edge, so a deletion at a capture
edge can cover an off-design site legitimately. That site is covered, so not filled, and
the QC reports `DEL` from DRAGEN's record. This is not hypothetical: Agilent ships both `Regions` and `Covered`
for CREv2 and they differ. Across the two validation cohorts, DRAGEN's records fall inside
Twist covered-targets and CREv2 `Covered` with **zero** off-design records in 32,500
site-resolutions, while CREv2 `Regions` leaves 309 off-design records and would gate off
907 of 1,207 recoveries — without failing, and without looking wrong in any log.

Mechanically, as implemented:

0. Compare the sample names in the two input headers and fail on a mismatch (§10).
   First in the job, reading the raw gVCF, so it depends on nothing computed and fails in
   seconds.
1. Build the DRAGEN bg-regions intermediate (unchanged first step) and index it.
2. Emit the covered span of every record reaching a defining site,
   `bcftools query -T <sites.bed> --targets-overlap 1 -f '%CHROM\t%POS0\t%END\n'`.
   Mode **1**, not the 2 the extract uses: this needs records whose *span* reaches a site,
   which is what `%END` reports and what `GvcfRecord.covers` counts as covering. Mode 2
   asks whether the *variant* overlaps and drops a deletion anchored on the site itself,
   which would make a covered site look like a hole.
3. Subtract the covered spans from the run's off-design sites to get what may be filled.
   The off-design set is the committed resource for the design and is not recomputed
   here; this subtraction is the awk one, over a few hundred spans (see §10). Separately,
   the off-design sites inside a DRAGEN *reference block*
   (`bcftools query -T <off_design.bed> --targets-overlap 1 -i 'INFO/END!="."'`) must be
   none — any fails the job, since it means the configured design is not the gVCF's. An
   off-design site under a DRAGEN variant's span is covered and unfilled, not an error.
4. Subset the post-hoc gVCF to records reaching a hole
   (`bcftools view -T <uncovered.bed> --targets-overlap 2 -e 'FORMAT/DP=0'`), set FILTER
   (`bcftools filter -s LowDepth -e 'FORMAT/DP<=1'`, which writes `PASS` on every other
   record), strip it to the fields the pipeline reads, and split multiallelics. FILTER has
   to be set because HaplotypeCaller leaves it `.` and rbceq2 uses an allele only when its
   defining variant is literally `PASS`; left alone, every recovered alternate allele was
   discarded before genotyping. `LowDepth` at DP<=1 is DRAGEN's own rule on these cohorts
   and the only one that means the same on both callers. DRAGEN's QUAL rule is not copied:
   the recalibrated gVCFs score QUAL on an ML scale (threshold 3) that is not comparable to
   HaplotypeCaller's, so post-hoc variants are used at any QUAL and their DP and GQ reach the
   QC flags, as a `PASS` DRAGEN variant's already do.
5. Drop the post-hoc *variants* that reach a base DRAGEN called (see **Boundary overlap**
   below), then tag every surviving record `INFO/POSTHOC=gatk-hc-<version>`.
6. Require the supplement and the primary gVCF to name the same sample, failing if they do
   not (see §10), then `bcftools concat -a`, which merges the two indexed files in
   coordinate order.
7. Run the existing extract and conversion over `merged.vcf.gz`.

Step 4's strip (`bcftools annotate -x`, keeping only END/GT/DP/GQ/MIN_DP) is what stops
`concat` having to reconcile two callers' definitions of tags nothing downstream reads.
Its `-e 'FORMAT/DP=0'` is what preserves `NOCOV` — see §10.

Step 6's check is what lets `concat` run at all, since it requires identical sample sets.
Relabelling the supplement would satisfy that too, and is exactly what must not happen: it
would turn a CRAM registered against the wrong sequencing group into a silent merge of
another individual's genotypes. A mismatch is fatal — see §10.

**Boundary overlap.** A record is selected for reaching a hole and is selected whole, so one
anchored in a hole can extend over a neighbouring defining site DRAGEN did call. Defining
sites are dense — most have another within 20bp — so near a capture edge this is the
ordinary case, and the two record types need opposite treatment.

A post-hoc **reference block** that straddles the edge is kept whole. It asserts nothing
rbceq2 sees, since the conversion drops every `<NON_REF>`-only record first, and dropping it
would throw away the hole it was kept for. The merged stream then carries two records at the
covered site, which is why `resolve_coverage` must prefer the DRAGEN one — see §5. Splitting
blocks on the boundary is the alternative and is not worth it.

A post-hoc **variant** that reaches a called base is dropped, and this was missed until
review. Kept, it puts two callers' alleles on one base in the file rbceq2 reads: a DRAGEN SNP
at a site, and a post-hoc deletion whose REF swallows it. The QC cannot report the conflict
either, because `resolve_coverage` prefers the primary record and so reads the base as an
ordinary `PASS`. The hole the dropped record would have filled returns to `NOCOV`, which is
the honest answer — DRAGEN wins wherever both speak, and here both spoke.

Two bcftools details make the drop work, and both were confirmed against the pinned 1.24
rather than reasoned about:

- `annotate -a <unfillable_sites.bed.gz> -m COVERED` marks on a record's whole span, not its POS,
  so a deletion anchored on a hole and reaching a called site one base away is marked. The
  annotation source is the defining sites DRAGEN covered (every site, less the holes), not the
  DRAGEN records' spans: marking on spans would also drop a post-hoc variant that merely clips
  the tail of a long reference block, reaching no called defining site and contradicting
  nothing rbceq2 reads.
- The mark and the drop run *before* `norm -m -any`, and the drop is
  `COVERED=1 && (N_ALT>1 || ALT!="<NON_REF>")`: a gVCF variant is still
  `<real ALT>,<NON_REF>` there and a block is `<NON_REF>` alone, so the alleles tell them
  apart and no INFO tag is consulted. The first version ran after the split and dropped a
  marked record with no `INFO/END`. That rested on HaplotypeCaller's habit of not writing
  `END` on a variant, and a variant that carried one survived; and switching to the ALT alone
  did not fix it, because after the split a variant's `<NON_REF>` twin inherits its REF span
  and its `END`, passes as a block, and fills the hole with an apparent hom-ref call. Before
  the split there is no twin. It also makes the logged drop count one per variant rather than
  one per allele `norm` would have split it into.

Note also that `--targets-overlap 2` in step 4 does *not* drop a deletion anchored on the
hole, the way it does in the extract. The record is still multiallelic there, and the
`<NON_REF>` allele makes bcftools match on the whole record span. That is why the drop is
needed rather than falling out of the selection mode.

### The extract gains a `posthoc` column

`%INFO/POSTHOC` is appended to the extract format **unconditionally** — genome runs and
primary records render `.`. One extract format everywhere, no per-sequencing-type
branching in the parser. `parse_extract`'s column-count check fails loudly on any extract
written before this change, which is the designed behaviour: delete stale extracts so the
conversion stage rewrites them (in practice the `workflow.version` bump in §7 makes this
moot — new tree, no stale files).

## 5. QC changes, in `FlagBloodGroupCallQc`

- **New flag prefix `POSTHOC`** for a site that resolves from a post-hoc record and
  passes the DP/GQ thresholds:
  `POSTHOC:1:159204893(T>C,src=gatk-hc-4.6.2.0,DP=42,GQ=99)`. A post-hoc site is never a
  silent `PASS` — a reviewer must be able to see the call rests on a re-call, not on the
  primary caller.
- **`src=` names the caller and version** rather than a bare `posthoc`, so a QC TSV says
  *which* caller stood in and a version bump is visible in the output rather than only in the
  code.
- **Severity order:** `NOCOV > DEL > LOWQ > PASS`, applied per site.
- **`POSTHOC` is joined to the severity, not ranked against it.** A flag name states two
  independent findings about the site: its severity (`NOCOV`, `DEL` or `LOWQ`, or absent when
  it clears both thresholds) and its provenance (`POSTHOC` when the post-hoc caller supplied
  the record). The two are joined with `+`, so a recovered site that passes is `POSTHOC` and
  one that is also sub-threshold is `LOWQ+POSTHOC`.

  Both are properties of the site, so both belong in the site's own flag. Ranking them makes
  one displace the other, and severity would win: a recovered site that was also sub-threshold
  would report only `LOWQ`, and nothing would say the call rested on a recovery. On the
  validation cohorts that would hide **234 of 681 reliant systems** — a third of the calls
  that exist only because of the recall. `rests_on_posthoc` reads the answer off a whole cell.

  A recovered site is listed even when it passes, because the antigen resting on a re-call of
  untargeted reads is itself the finding. A site that is neither poor nor recovered has no
  flag name at all and is not listed, which is what keeps a bare `PASS` cell meaning "nothing
  to report". `NOCOV` never carries `POSTHOC`: no record from either caller means no caller
  to name.

  A `POSTHOC` site counts as flagged (the cell is not `PASS`), but logs count quality-flagged
  systems, clean recoveries and reliant systems separately, so a run summary distinguishes
  "problems" from "recovered" and still says how many calls needed the recall.

- **`NOCOV` keeps its one meaning:** no record from *either* caller covers the site.
- **`resolve_coverage` prefers primary records.** Where both a DRAGEN record and a
  post-hoc record cover a site (the §4 boundary case), the DRAGEN one is chosen —
  "DRAGEN wins wherever both speak" applies to the QC's view as well as the merge.
- **A post-hoc record counts only at a fillable site.** The QC stage is handed the
  off-design BED the merge filled from, and `flags_by_system` disregards a post-hoc
  record at any site not in it. That is what keeps an in-design hole `NOCOV` when a kept
  post-hoc reference block spans it (§4). A genome run hands over no BED, and the job
  refuses an extract carrying post-hoc records without one rather than trusting them.

## 6. DAG and config

`stages/pipeline.py`:

```python
PosthocGenotypeOffTargetSites = stage_support.wire(
    posthoc_genotype.PosthocGenotypeOffTargetSites,
)  # reads sequencing_group.cram directly; no requires, no Metamist Analysis
FilterAndConvertGvcfsForRbceq2 = stage_support.wire(
    filter_and_convert.FilterAndConvertGvcfsForRbceq2,
    requires=[PosthocGenotypeOffTargetSites],
)  # also reads the committed off-design BED for the configured design (off_design.resource_path)
FlagBloodGroupCallQc = stage_support.wire(
    call_qc.FlagBloodGroupCallQc,
    requires=[FilterAndConvertGvcfsForRbceq2, GenotypeBloodGroupsWithRbceq2],
    ...
)  # the QC reads the same off-design BED the merge filled from (§5)
```

No Metamist Analysis for the new stage: its output is an intermediate consumed by the
graph, like the converted VCF (same reasoning as the debug-log decision). Provenance
reaches Metamist through the QC Analysis meta instead (§5, and record the caller version
in the QC stage's meta alongside the thresholds).

The conversion stage requests the post-hoc input only for exome sequencing groups, and for a
genome run the new stage produces nothing, so nothing is merged and the merge is a plain
rename. The genome conversion job is *not* byte-identical to today's, though: every run gains
the `INFO/POSTHOC` header line on the intermediate and a trailing `POSTHOC` column in the
extract, because `bcftools query` aborts on a tag the header does not declare. That is what
the release version bump records.

New config sections, following the class-name convention:

```toml
[workflow.posthoc_genotype_off_target_sites]
cpu = 2
memory = "standard"   # or "highmem"; lowmem is refused, the JVM heap is sized from the tier
storage = "20Gi"      # streams the CRAM; disk is for the ~3Gb reference and a tiny gVCF
```

And one `[workflow]` key, required for an exome run and never read by a genome one:

```toml
[workflow]
exome_design_bed = 'exome_probesets_hg38/<design>'   # references-repo key; selects the committed subtraction
```

## 7. Versioning

The extract format and the QC flag vocabulary change for **every** run, and exome geno
TSVs change where holes get filled, so `workflow.version` is bumped in the same PR and the
new outputs land in a fresh `rbceq2_<tool>_<release>` tree where no old extract can meet the
new parser.

`v2` covered the extract's INFO/POSTHOC column and the first QC flag vocabulary. `v3` covers
three later exome-only output changes: the capture-design gate on which holes may be filled
(§4), dropping post-hoc variants that reach a base DRAGEN called (§4, **Boundary overlap**),
and the composition of provenance into a site's flag name (§5). None moves an output path, so
without the bump an exome re-run would reuse its v2 files and none would take effect. The two
validation runs in RESULTS.md predate `v3` and wrote to the v2 tree.

`v4` covers two exome-only changes from the review of 2026-09-03: post-hoc records carry a
FILTER rbceq2 accepts (`PASS`, or `LowDepth` at DP<=1), where HaplotypeCaller left `.` and
rbceq2 therefore excluded every recovered alternate allele; and the QC disregards a post-hoc
reference block at an in-design hole (§5). The v3 re-run in RESULTS.md predates both, so its
genotype tables used no recovered allele.

An exome run's tree has one more segment than a genome run's, the configured design:
`rbceq2_<tool>_<release>/<design key>/<stage>/...`. cpg_flow reuses a stage whose expected
outputs exist and asks nothing about how they were made, and the release segment cannot see a
config change. The first draft put the design only in the output path of the subtraction,
then a stage, which let a repointed key rebuild that BED and then reuse every sample's conversion, genotypes
and QC from the old design, silently. Putting it under the release moves the whole run
instead. A genome run never reads the key and its tree is unchanged, so genome outputs stay
where v4 put them; the v4 CREv2 re-run in RESULTS.md predates the design segment and sits
directly under `rbceq2_2_4_3_v4`.

## 8. Change table

| file | change |
|---|---|
| `stages/blood_group_genotyping/posthoc_genotype.py` | new stage class |
| `stages/blood_group_genotyping/filter_and_convert.py` | merge step for exome SGs; extract gains `posthoc` column |
| `stages/blood_group_qc/call_qc.py` + `jobs/rbceq2_call_qc_job.py` | `POSTHOC` flag, `src=`, primary-record preference, separated counts, `posthoc_caller` in the Analysis meta |
| `stages/pipeline.py` | wire the new stage; add to conversion's `requires` |
| `scripts/gen_bg_resources.py` + `scripts/bg_db.py` | write `bg_defining_sites_padded.<genome>.bed` |
| `resources/` | the new committed BED: 199 intervals, 135,579 bases |
| `config/popgen_rbceq2_default_config.toml` | new stage section; `version = 'v4'` |
| `off_design.py` + `scripts/gen_off_design_sites.py` | the design subtraction as a committed resource: the generator, the manifest, and the graph-build resolver |
| `resources/` | `bg_off_design_sites.<design>.<genome>.bed` for Twist and CREv2, and `bg_off_design_sites.manifest.tsv` |
| `config/config_template.toml` | the required `exome_design_bed` key, with how to pick it |
| `constants.py` | `GATK_VERSION`, `GATK_IMAGE_TAG`, `POSTHOC_CALLER` |
| `stage_support.py` | the design key and `exome_design_bed()`; an exome run's release tree gains the design as a segment |
| README / PRODUCT.md / GLOSSARY.md | document `POSTHOC`, the stage, and the fill rule |
| `tests/test_posthoc_merge.py` | new: runs the real awk against `GvcfRecord.covers` |
| `tests/test_exome_design_gate.py` | new: the design key is required, resolved and enforced |
| `tests/test_off_design_subtraction.py` + `tests/test_off_design_resources.py` | new: bedtools' interval semantics under the generator's shell; the committed resources match their manifest and the shipped sites |
| `tests/test_posthoc_trespass.py` | new: runs the real merge shell under real bcftools |
| `tests/test_posthoc_heap.py` | new: the JVM heap tracks the configured memory tier |
| tests | QC flag logic, severity, extract parsing, exome gating, resource generation |

The reference *fasta* needed no new key — see §9.1. The capture design does: an exome run
sets `workflow.exome_design_bed` to the design's references-repo key, which selects the
committed subtraction for that design rather than being resolved to the vendor BED.

## 9. Resolved questions

1. **Reference path** — settled. `references.broad.ref_fasta` already points at
   `hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta`, the masked assembly38
   the CRAMs were aligned against, so no new config key was needed. Its `.fai` and `.dict`
   sit beside it and are localised with it.
2. **DRAGstr** — start without. `--dragen-mode` with no sample DRAGstr table means STR
   genotyping is not fully DRAGEN-equivalent. Revisit if post-hoc indel calls at
   STR-adjacent defining sites look discordant with DRAGEN's where both callers reach.
3. **Keep the post-hoc gVCF?** — no, tmp. The QC flag names the caller and the merged
   extract records what it reported, so the gVCF itself has no reader.
4. **Genome CRAM recall** — deliberately not done. A genome gVCF has no capture edge to
   stop at. Noted in PRODUCT.md as resting on an untested assumption: that a genome
   `NOCOV` site is unmappable rather than merely uncalled.

5. **The design subtraction is a committed resource, not a stage.** It was first a
   MultiCohortStage, `SelectOffDesignDefiningSites`, running `bedtools intersect -v` once per
   run. MultiCohort rather than Cohort because of how a per-sequencing-group consumer reads a
   stage's output back: it has to name the target the output was filed under, and
   `cpg_flow.inputs.get_multicohort()` always returns the run's one `MultiCohort`, whereas a
   `SequencingGroup` carries no link to its cohort (one can be in several), so a CohortStage's
   output could only be found by listing the MultiCohort's cohorts and assuming there was
   one. Both shapes were wrong for the question. The answer depends on nothing about the run,
   only on the configured design and the committed defining sites, so recomputing it per run
   bought nothing and cost a stage, a bedtools image in the pipeline, and a lookup question
   whose honest answer was "the run has one cohort". It is now
   `resources/bg_off_design_sites.<design>.<genome>.bed`, written once per design by
   `scripts/gen_off_design_sites.py` and resolved at graph build by `off_design.resource_path`,
   the way the defining sites it is a subset of already were. The generator keeps the three
   input checks the stage made (empty design, mismatched contig names, zero-length rows),
   and adds one the stage lacked: a design that targets every defining site is refused
   rather than committed as an empty BED, and the resolver refuses a manifest row recording
   zero sites, because the merge hands the file to `bcftools -T`, which aborts on an empty
   targets file. Both committed files are byte-identical to what the stage wrote for the
   validation runs.
   The cost is that derived data in the repo can go stale, which the manifest's `sites_md5`
   guards: the resolver and the test suite both refuse a subtraction of defining sites that
   are not the shipped ones, so regenerating the sites without the subtractions is a red suite
   and a failed graph build, never a wrong fill. A cohort on a design with no committed file
   fails at graph build naming the generator. `exome_design_bed` moved from the stage's
   config section to `[workflow]`, there being no stage to own it.

## 10. What testing changed about the design

None of these came out of reading documentation. Each came from running the pinned images
(`bcftools:1.24-1`, `gatk:4.6.2.0-2`) against synthetic inputs and real CPG data, and each
would have reached production:

- **`bcftools query` aborts on an undeclared INFO tag** rather than rendering `.`
  (`no such tag defined in the VCF header: INFO/POSTHOC`). Only the exome-with-holes path
  was getting the declaration via the merged file, so **every genome run would have died**
  at the extract. The intermediate now declares `POSTHOC` unconditionally.
- **HaplotypeCaller emits `DP=0,GQ=0` reference blocks across untargeted stretches** of the
  `-L` interval. Merging those would have put a record over every hole and **retired
  `NOCOV` for exomes entirely** — a site with no reads would have read `LOWQ(DP=0)`, saying
  "poor data" where the truth is "no data". Zero-depth records are dropped.
- **`bcftools annotate -x` applies its keep-list caret per category.** The first version
  used one leading `^` across INFO and FORMAT, which bcftools read as "keep these INFO
  tags, *remove* these FORMAT tags" — silently stripping GT, DP, GQ and MIN_DP from every
  supplement record while keeping a tag nothing reads. Two carets is the correct form.
- **An empty `-T` targets file is a hard error**, not an empty result
  (`Failed to read the targets`). The no-holes branch is required, not an optimisation.
- **The sample-name assertion is fatal, and some test inputs trip it.** It catches a CRAM
  registered against the wrong sequencing group, which would otherwise splice another
  individual's genotypes into these calls at exactly the sites nothing else covers. Every
  mackenzie recal gVCF carries the current sequencing-group ID, but an older block of test
  CRAMs still carries a retired one for the same individual, because an upstream test-set
  script reheadered some inputs and not others; newer additions are consistent. Sampling the
  test bucket, 5 of the first 12 CRAMs disagree with their gVCF, all in that older block.

  Relabelling the supplement would satisfy `concat` and hide this, and was the design for one
  commit. It was reverted: from inside the job a real swap and a stale header are
  indistinguishable, so downgrading the check to keep an old test cohort green weakens every
  cohort. The affected inputs are the thing to fix. Identity confirmation belongs to somalier,
  whose output sits beside these CRAMs; wiring it in is the follow-up if this ever needs to
  distinguish the two cases automatically.
- **Metamist registers several CRAMs and gVCFs per sequencing group here**, and
  `cpg_flow.metamist.get_analyses_by_sgid` keeps whichever the API returns last, unordered.
  It currently selects the DRAGEN 3.7.8 CRAM and the matching recal gVCF — the pair we want —
  but nothing pins that, so it is worth re-checking after any new analysis is registered.

Environment facts worth keeping:

- **The design subtraction moved from a hand-written awk loop to bedtools.** It first ran
  inside every conversion job as an awk containment loop over the vendor BED. `bedtools
  intersect -v` is the standard tool for exactly that operation, so it is easier to read and
  maintain, and the answer is the same for every sample, so it ran once per run in its own
  stage, and now runs once per design in `scripts/gen_off_design_sites.py`, committed (§9.5).
  Three behaviours had to match the awk for the swap to be safe, and all three were checked
  against `bedtools:2.30.0-1` on both vendor designs in use: half-open interval semantics,
  skipping `track`/`browser` lines, and ignoring columns past the third. Its output diffed
  identical to the awk's. One behaviour does not match and is refused instead: a row with end
  not greater than start, which the awk ignored and bedtools counts as covering a base. Neither
  design has one, so the generator fails on such a row rather than choose a meaning for it.
- **bedtools is not in the bcftools image** (`debian:bookworm-slim` plus bcftools binaries
  only), and its awk is **mawk**, not gawk; the remaining awk program was re-run under mawk in
  that exact image. There *is* a `cpg-common/images/bedtools` at `2.30.0-1`, which the stage
  used while the subtraction was one; the pipeline pulls no bedtools image now. The per-sample
  hole-finding stays in awk inside the conversion job (§4) because it is derived from that
  job's gVCF and the image has no bedtools.
- **GATK 4.6.2.0 rejects CRAM 3.1** (`CRAM version 3.1 is not supported`). CPG CRAMs are
  3.0, so this is a future trap, not a current one.
- **The masked reference is safe for these CRAMs.** They were aligned to unmasked hg38;
  the two builds share all 3,366 contig names, and the 786 differing checksums are all HLA
  contigs. Every blood-group defining site is on chr1-22 or chrX.
- **GATK's `gs://` NIO access works**, and is cheap: 132s and $0.004 per sample in batch,
  measured on the validation run. A laptop probe could only get as far as `Starting
  traversal`, because it streamed the 3Gb reference too, which the stage does not do.

The hole-finding rule is expressed twice — in awk in the merge, and as
`GvcfRecord.covers` in the QC job — because the merge runs in the bcftools image, which
has no Python package of ours. `tests/test_posthoc_merge.py` runs the real awk and asserts
it marks exactly the sites `covers` calls uncovered, which is the only thing tying the two
together. Moving that half to bedtools would not remove the duplication, only relocate it:
bedtools is a third implementation of containment whose semantics still have to be pinned
against `covers` by a test, and the covered spans are derived from the gVCF inside a job
whose image has no bedtools. So the awk stays for hole-finding, and the one subtraction that
is sample-independent is bedtools', run by the generator and committed (§4).

One awk program served both subtractions at first. Containment in a set of `[start, end)`
spans is the same test whether the spans come from `%POS0\t%END` on DRAGEN records or from a
vendor capture BED, so `_SITES_OUTSIDE_SPANS_AWK` was run twice rather than written twice.
That is no longer true of the design half: it is `bedtools intersect -v` in the generator,
committed once per design, and the vendor BED's extra columns and `track` header line are
bedtools' problem now. Both were re-checked against `bedtools:2.30.0-1` rather than assumed
to carry over, and its subtraction of the real Twist design diffed identical to the awk's
over all 1,625 committed sites. The awk keeps the per-sample half, where its spans always
come from `bcftools query` and so are three plain columns.
