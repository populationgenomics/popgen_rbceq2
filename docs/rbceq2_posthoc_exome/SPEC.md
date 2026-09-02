# Design Spec: post-hoc genotyping of off-target defining sites in exomes

**Status:** IMPLEMENTED. Kept as the record of why.
**Area:** rbceq2 blood-group genotyping pipeline (`FilterAndConvertGvcfsForRbceq2` →
`GenotypeBloodGroupsWithRbceq2` → `FlagBloodGroupCallQc`), in `popgen_rbceq2`.
**Author:** Joshua Schmidt · **Reviewers:** (fill in)

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
sequencing_group.gvcf ──> FilterAndConvertGvcfsForRbceq2 <────────┘
                             (merge: DRAGEN wins; supplement tagged INFO/POSTHOC)
                                    │                    │
                                converted VCF        DP/GQ extract (+ posthoc column)
                                    │                    │
                          GenotypeBloodGroupsWithRbceq2  FlagBloodGroupCallQc
                                                         (POSTHOC flag / src=posthoc)
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
  Broad's standard fasta is not it. Path comes from a `references.*` config key
  (exact key: open question 1).
- **Image:** `image_path('gatk', ...)` — 4.6.2.0-2 is current in cpg-common.
- **Output:** `<sg>.posthoc.g.vcf.gz` + `.tbi`, under the tmp category like the
  conversion stage's intermediates (open question 3 argues for keeping it).

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

**Fill empirically, not by capture metadata.** No capture BED is consulted. The rule is:
a post-hoc record is kept only at coordinates the DRAGEN gVCF has no record covering
(reference blocks count as covering); DRAGEN wins wherever both speak. This makes the
merge self-calibrating per sample — a capture BED that lies about the real footprint
cannot cause a DRAGEN call to be overwritten or a hole to be missed.

Mechanically, as implemented:

1. Build the DRAGEN bg-regions intermediate (unchanged first step) and index it.
2. Emit the covered span of every record reaching a defining site,
   `bcftools query -T <sites.bed> --targets-overlap 1 -f '%CHROM\t%POS0\t%END\n'`.
   Mode **1**, not the 2 the extract uses: this needs records whose *span* reaches a site,
   which is what `%END` reports and what `GvcfRecord.covers` counts as covering. Mode 2
   asks whether the *variant* overlaps and drops a deletion anchored on the site itself,
   which would make a covered site look like a hole.
3. Subtract those spans from the defining-sites BED to get the holes (awk; see §10).
4. Subset the post-hoc gVCF to records reaching a hole
   (`bcftools view -T <uncovered.bed> --targets-overlap 2 -e 'FORMAT/DP=0'`), strip it to
   the fields the pipeline reads, split multiallelics, and tag every kept record
   `INFO/POSTHOC=gatk-hc-<version>`.
5. Relabel the supplement to the primary gVCF's sample name (warning on a mismatch, see
   §10), then `bcftools concat -a`, which merges the two indexed files in coordinate order.
6. Run the existing extract and conversion over `merged.vcf.gz`.

Step 4's strip (`bcftools annotate -x`, keeping only END/GT/DP/GQ/MIN_DP) is what stops
`concat` having to reconcile two callers' definitions of tags nothing downstream reads.
Its `-e 'FORMAT/DP=0'` is what preserves `NOCOV` — see §10.

Step 5's relabel is what `concat` needs, since it requires identical sample sets. A
mismatch is logged rather than fatal — see §10 for why the read-group name cannot serve as
an identity check on this data.

**Boundary overlap.** A post-hoc reference block can straddle a capture edge, covering
one uncovered defining site and also one DRAGEN covers. Step 3 keeps the whole record, so
the merged stream can carry overlapping records at a covered site. That is harmless for
rbceq2 (the DRAGEN record still carries the variant it read) but the QC's
`resolve_coverage` must prefer the DRAGEN record where both cover a site — see §5. The
implementation should not try to split blocks.

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
- **Failing post-hoc sites keep their existing prefix** (`LOWQ`/`DEL`), with the same
  `src=` key in the metrics, e.g. `LOWQ:...(T>C,src=gatk-hc-4.6.2.0,DP=6,GQ=12)`. The
  failure mode stays primary in the flag; the provenance is secondary. `src=` names the
  caller and version rather than a bare `posthoc`, so a QC TSV says *which* caller stood in
  and a version bump is visible in the output rather than only in the code.
- **Severity order:** `NOCOV > DEL > LOWQ > POSTHOC > PASS`. `POSTHOC` counts as flagged
  (the per-system cell is not `PASS`), but logs count quality-flagged systems and
  post-hoc-only systems separately, so a run summary distinguishes "problems" from
  "recovered".
- **`NOCOV` keeps its one meaning:** no record from *either* caller covers the site.
- **`resolve_coverage` prefers primary records.** Where both a DRAGEN record and a
  post-hoc record cover a site (the §4 boundary case), the DRAGEN one is chosen —
  "DRAGEN wins wherever both speak" applies to the QC's view as well as the merge.

## 6. DAG and config

`stages/pipeline.py`:

```python
PosthocGenotypeOffTargetSites = stage_support.wire(
    posthoc_genotype.PosthocGenotypeOffTargetSites,
)  # reads sequencing_group.cram directly; no requires, no Metamist Analysis
FilterAndConvertGvcfsForRbceq2 = stage_support.wire(
    filter_and_convert.FilterAndConvertGvcfsForRbceq2,
    requires=[PosthocGenotypeOffTargetSites],
)
```

No Metamist Analysis for the new stage: its output is an intermediate consumed by the
graph, like the converted VCF (same reasoning as the debug-log decision). Provenance
reaches Metamist through the QC Analysis meta instead (§5, and record the caller version
in the QC stage's meta alongside the thresholds).

The conversion stage requests the post-hoc input only for exome sequencing groups; for a
genome run the new stage produces nothing and the conversion job is byte-identical to
today's.

New config section, following the class-name convention:

```toml
[workflow.posthoc_genotype_off_target_sites]
cpu = 2
memory = "standard"
storage = "20Gi"   # streams the CRAM; disk is for the ~3Gb reference and a tiny gVCF
```

## 7. Versioning

The extract format and the QC flag vocabulary change for **every** run, and exome geno
TSVs change where holes get filled. Bump `workflow.version` to `v2` in the same PR, so
the new outputs land in a fresh `rbceq2_<tool>_v2` tree and no old extract can meet the
new parser.

## 8. Change table

| file | change |
|---|---|
| `stages/blood_group_genotyping/posthoc_genotype.py` | new stage class |
| `stages/blood_group_genotyping/filter_and_convert.py` | merge step for exome SGs; extract gains `posthoc` column |
| `stages/blood_group_qc/call_qc.py` + `jobs/rbceq2_call_qc_job.py` | `POSTHOC` flag, `src=`, primary-record preference, separated counts, `posthoc_caller` in the Analysis meta |
| `stages/pipeline.py` | wire the new stage; add to conversion's `requires` |
| `scripts/gen_bg_resources.py` + `scripts/bg_db.py` | write `bg_defining_sites_padded.<genome>.bed` |
| `resources/` | the new committed BED: 199 intervals, 135,579 bases |
| `config/popgen_rbceq2_default_config.toml` | new stage section; `version = 'v2'` |
| `constants.py` | `GATK_VERSION`, `GATK_IMAGE_TAG`, `POSTHOC_CALLER` |
| README / PRODUCT.md / GLOSSARY.md | document `POSTHOC`, the stage, and the fill rule |
| `tests/test_posthoc_merge.py` | new: runs the real awk against `GvcfRecord.covers` |
| tests | QC flag logic, severity, extract parsing, exome gating, resource generation |

No new `references.*` key was needed — see §9.1.

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
- **A sample-name assertion was tried and rejected on the evidence.** The intent was to catch
  a CRAM registered against the wrong sequencing group. Checking the mackenzie DRAGEN 3.7.8
  test exomes showed the signal is not usable: every recal gVCF carries the current
  sequencing-group ID, but 4 of 10 CRAMs in the target cohort still carry a retired one for
  the same individual, because an upstream test-set script reheadered some inputs and not
  others. Asserting would have failed 40% of a known-good cohort, and a genuinely swapped
  CRAM could carry a stale-but-matching name anyway. The supplement is relabelled to the
  gVCF's name and the mismatch is logged. Identity belongs to somalier, whose output sits
  beside these CRAMs; wiring it in is the follow-up if this runs on less trusted data.
- **Metamist registers several CRAMs and gVCFs per sequencing group here**, and
  `cpg_flow.metamist.get_analyses_by_sgid` keeps whichever the API returns last, unordered.
  It currently selects the DRAGEN 3.7.8 CRAM and the matching recal gVCF — the pair we want —
  but nothing pins that, so it is worth re-checking after any new analysis is registered.

Environment facts worth keeping:

- **bedtools is not in the bcftools image** (`debian:bookworm-slim` plus bcftools binaries
  only), which settles the awk question in §4/§10 — awk stays. Its awk is **mawk**, not
  gawk; both programs were re-run under mawk in that exact image.
- **GATK 4.6.2.0 rejects CRAM 3.1** (`CRAM version 3.1 is not supported`). CPG CRAMs are
  3.0, so this is a future trap, not a current one.
- **The masked reference is safe for these CRAMs.** They were aligned to unmasked hg38;
  the two builds share all 3,366 contig names, and the 786 differing checksums are all HLA
  contigs. Every blood-group defining site is on chr1-22 or chrX.
- **GATK's `gs://` NIO access works** — a probe against a real CRAM and reference reached
  `Starting traversal`. Runtime on real data is still unverified: that probe streamed the
  3Gb reference too, which the stage does not do, and it did not finish inside seven
  minutes from a laptop. Confirm on the first real batch.

The hole-finding rule is expressed twice — in awk in the merge, and as
`GvcfRecord.covers` in the QC job — because the merge runs in the bcftools image, which
has no Python package of ours. `tests/test_posthoc_merge.py` runs the real awk and asserts
it marks exactly the sites `covers` calls uncovered, which is the only thing tying the two
together.

The hole-finding rule is expressed twice — in awk in the merge, and as
`GvcfRecord.covers` in the QC job — because the merge runs in the bcftools image, which
has no Python package of ours. `tests/test_posthoc_merge.py` runs the real awk and asserts
it marks exactly the sites `covers` calls uncovered, which is the only thing tying the two
together. If bedtools turns out to be in that image, `bedtools intersect -v` replaces the
awk and removes the duplication.
