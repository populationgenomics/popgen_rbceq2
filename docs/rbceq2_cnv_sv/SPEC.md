# Design Spec: CNV/SV support for RBCeq2 blood-group calling

**Status:** DRAFT — for review. Do not implement until the design is approved.
**Revision 2026-09-18:** Q1 and Q2 resolved from the 100-genome study; §5.6 gains a karyotype gate; new §11 QC design. See the change list at the end of §7.
**Area:** rbceq2 blood-group genotyping pipeline (`FilterAndConvertGvcfsForRbceq2` →
`GenotypeBloodGroupsWithRbceq2` → `CombineRbceq2OutputsPerCohort`), in `popgen_rbceq2`.
**Author:** Joshua Schmidt · **Reviewers:** (fill in)
**Companion docs (same directory):**
- [`dragen_three_source_merge.md`](./research/dragen_three_source_merge.md) — the visual overview + mermaid diagram (read this first for the shape).
- [`implement_cnv_rbceq2_research.md`](./research/implement_cnv_rbceq2_research.md) — the underlying RBCeq2 source analysis this spec rests on.


> **Moved from `ourdna_genomic_atlas` (2026-08).** The RBCeq2 stages now live in this repo, and
> this design was never implemented in either. The class names below are still correct, but the
> file paths are not: they predate both the stage restructure and the repo split, so
> `src/ourdna_genomic_atlas/stages.py` in §8's change table now maps to
> `src/popgen_rbceq2/stages/blood_group_genotyping/`. Content is otherwise unchanged from the
> original, including the open questions.

---

## 1. Summary

Today the pipeline feeds RBCeq2 **only SNVs/indels**, converted from each sample's DRAGEN
SNV gVCF. RBCeq2's allele database also defines a large set of blood-group alleles by
**structural variants**. These are either whole-gene or multi-exon deletions/duplications (CNVs), and
paired del+ins gene-conversion products (SVs). We are yet to define how we can + will
utilise these data, therefore currently we include none of these in rbceq2.

Our DRAGEN 3.7.8 align + genotyping pipeline produces the necessary outputs for most of the
CNV/SV variants in Rbceq2. Alongside the SNV gVCF it produces a
**Manta-derived SV VCF** and a **bin-based CNV VCF** per sample. So that is three files from DRAGEN
to use.

RBCeq2 reads a single file, so structural records need to be in the *same* `--vcf` as SNVs i.e.
we need to produce a single file.

So the work is a **per-sample preprocessing + merge**:
- convert the SNV gVCF to a plain sites VCF (drop `<NON_REF>`/reference blocks - the pipeline already does this)
- normalise the SV VCF
- **fix the CNV VCF**
- concatenate the three into one sorted/bgzipped/tabixed VCF

This new concatenated VCF will then become the input for the `GenotypeBloodGroupsWithRbceq2` stage.

The one non-trivial transform is the CNV VCF: DRAGEN writes **every** CNV record as
`SVTYPE=CNV`, but RBCeq2 requires the event `SVTYPE` to equal the database token type
(`DEL`/`DUP`). Left unfixed, **every large gene deletion silently fails to match**. This
spec pins down that rewrite and the surrounding merge.

We will not modify RBCeq2. Every change here lives in `ourdna_genomic_atlas`
and is built around RBCeq2's *observed* behaviour.

---

## 2. Background — how RBCeq2 consumes structural variants

### 2.1 One VCF, no CNV/SV flag
`find_hits` (`src/rbceq2/main.py:272-305`) reads SVs out of the same dataframe as SNVs:
`SvReader(df=vcf.df, min_size=args.min_size).events()` → `SvMatcher.match(db_defs, events)`
→ `select_best_per_vcf(...)`. There is no CLI flag and no second input; structural records
must be present in the single `--vcf`.

### 2.2 What the DB encodes
Structural alleles live in the `GRCh37`/`GRCh38` columns of `db.tsv` as word-form tokens
`<pos>_<type>_<len>` (e.g. `143914828_del_110kb`). Token-type totals in the GRCh38 column:
**132 DEL, 41 INS, 1 DUP — no `INV`, no `BND`, no literal `<CNV>`.** Practical buckets:
~40 true large CNVs (whole-gene/multi-exon del/dup), ~22 hybrid/complex SVs (paired DEL+INS,
RH & GYP), ~38 large indels (<~1 kb). Hybrids are stored as *paired* tokens, not a "hybrid"
type. (Full annotated list in the research doc §3.)

### 2.3 Matching is fuzzy in position/length, but strict on type
`SvMatcher` (`large_variants.py:159-367`) uses adaptive positional/length tolerance and
reciprocal-overlap gates, so imprecise breakpoints still match — **but
`require_same_type=True` by default**, so the DB token type must equal the event `SVTYPE`.
The symbolic-ALT fallback (`<DEL>`/`<DUP>`) in `SvReader` fires **only when `SVTYPE` is
absent** (`large_variants.py:739`). This single fact is what forces the CNV rewrite in §5.4.
RBCeq2 also pre-filters records to **±500 kb** of a DB position (`vcf.py:497`).

---

## 3. The three DRAGEN sources (evidence)

Grounded on one of the OurDNA 1000 Genomes (1KG) control replicates:

| # | Source | Path (`gs://cpg-ourdna-main/ica/dragen_3_7_8/output/`) | Caller | Ready? |
|---|--------|--------------------------------------------------------|--------|--------|
| 1 | SNV gVCF | `recal_gvcf/<sg>.hard-filtered.recal.gvcf.gz` | DRAGEN SNV | ⚠️ convert to sites VCF |
| 2 | SV VCF | `dragen_metrics/<sg>/<sg>.sv.vcf.gz` | Manta-derived | ✅ compatible |
| 3 | CNV VCF | `dragen_metrics/<sg>/<sg>.cnv.vcf.gz` | DRAGEN bin CNV | ⚠️ SVTYPE rewrite |
| — | Ploidy (sex) | `dragen_metrics/<sg>/<sg>.ploidy_estimation_metrics.csv`, `<sg>.ploidy.vcf.gz` | DRAGEN | read-only |

**1 — SNV gVCF** (called variant records — already carrying GT — interleaved with `<NON_REF>` reference blocks; the head below happens to show only ref blocks):
```text
chr1  9997   .  N  <NON_REF>  .  PASS  END=10017  ...  ./.:11,0:11:0:2:...
chr1  10018  .  C  <NON_REF>  .  PASS  END=10018  ...  0/0:17,0:17:15:17:...
```

**2 — SV VCF** (Manta; header ALTs `<DEL>`,`<INS>`,`<DUP:TANDEM>`; proper SVTYPE):
```text
chr1  789481  MantaINS:…  G  <INS>  999  PASS  END=789481;SVTYPE=INS;CIPOS=0,7;…
chr1  839442  MantaDEL:…  CACC…ACA  CT  463  PASS  END=839499;SVTYPE=DEL;SVLEN=-57;CIGAR=1M1I57D
```

**3 — CNV VCF** (every record `SVTYPE=CNV`; direction only in ALT + `CN`):
```text
chr1  817861   DRAGEN:REF:…   N  .      81  PASS       END=2650427;REFLEN=1832567                     ./.:0.98:2:1434:…
chr1  2650427  DRAGEN:LOSS:…  N  <DEL>  51  cnvLength  SVLEN=-2648;SVTYPE=CNV;END=2653075;REFLEN=2648  1/1:0.058:0:2:…
chr1  3501568  DRAGEN:GAIN:…  N  <DUP>  88  cnvLength  SVLEN=1000;SVTYPE=CNV;END=3502568;REFLEN=1000   ./1:3.18:6:1:…
```

---

## 4. Scope & non-goals

**In scope**
- New per-SG inputs: resolve the SV and CNV VCF paths from `dragen_metrics/<sg>/`.
- These inputs must be **registered in metamist** as dedicated analysis types (one for the SV VCF,
  one for the CNV VCF) before the stage can resolve them. **Ideally this registration happens
  upstream in the `dragen_align` pipeline** that produces the files; failing that, this work adds a
  small stage to define the analysis type(s) and register the existing outputs.
- A preprocessing step that (a) produces an SNV sites VCF, (b) normalises the SV VCF,
  (c) **fixes the CNV VCF** (drop REF records, rewrite `SVTYPE=CNV`→`DEL`/`DUP`), and
  (d) merges all three into one sorted/bgzipped/tabixed VCF per SG.
- Feed the merged VCF into `GenotypeBloodGroupsWithRbceq2` in place of the SNV-only VCF.
- Per-SG sex/ploidy read from the DRAGEN ploidy files (not the header) for CNV direction on
  chrX/chrY.
- A cross-replicate concordance check over the five control SGs.

**Non-goals (explicit)**
- **No changes to RBCeq2** — not `SvReader`, `SvMatcher`, `db.tsv`, or the region filter.
  If an option here appears to need a RBCeq2 change, it is rejected.
- **RH / GYP hybrid alleles are out of scope.** RHD/RHCE and GYPA/B/E are near-identical
  paralogs → MAPQ0, noisy/absent short-read SV & CNV calls; `--RH` is documented long-read
  only. We target straightforward large dels/dups and large indels; hybrids are explicitly
  not claimed for short-read DRAGEN.
- No change to the blood-group science / allele definitions.
- No re-run/backfill orchestration (tracked separately).

**NB**. I havent pursued the needs re RHD/RHCE and GYPA/B/E any further, until we gather
knowledge around how many samples and when long read data is being produced.

---

## 5. Proposed design

New per-SG stage (working name **`PreprocessDragenForRbceq2`**) inserted **between** input
resolution and `GenotypeBloodGroupsWithRbceq2`. It emits one merged VCF per SG; the existing
rbceq2 stage is repointed at it. Execution is bcftools-based (same image family as the
current filter stage). See the mermaid diagram in the companion doc.

### 5.1 Resolve the two new inputs + register in metamist
The SV and CNV VCFs are not currently referenced by the pipeline. Resolve them via the
existing `get_dragen_output_path` helper (`popgen_utils.py:32`), e.g.
`get_dragen_output_path(f'dragen_metrics/{sg}/{sg}.sv.vcf.gz')`. Both are bgzipped and
tabixed in the bucket (`.tbi` present). Missing-file handling should mirror the current
gVCF guard. These must be registered as metamist analyses under the new SV/CNV analysis types —
ideally emitted upstream by `dragen_align` (see §4), not first-registered here.

### 5.2 SNV gVCF → sites VCF
A single-sample gVCF **already contains the per-sample genotypes** (GT at variant sites,
plus `<NON_REF>` reference blocks). Producing the sites VCF RBCeq2 needs is therefore pure
preprocessing — split multiallelics, drop the `<NON_REF>` symbolic allele and the
reference-only blocks, region-restrict — **not** a genotyping step. (GATK `GenotypeGVCFs`
is *joint* genotyping across samples and has no role here.)

This is **already implemented** by `FilterAndConvertGvcfsForRbceq2` (`stages.py:410`):
`bcftools norm -m -any` + region-restrict to `resources/bg_regions.GRCh38.bed`, dropping
`<NON_REF>` (bcftools only; no reference FASTA needed). Reuse it — **with one addition, now landed in
[PR #128](https://github.com/populationgenomics/ourdna_genomic_atlas/pull/128):** a final
`| bcftools +fixploidy` in the same pipe diploid-ises DRAGEN's haploid GTs, or RBCeq2 crashes
(Q4). **Invoke it with no `-s`/`-p` arguments.** On our `chr`-prefixed hg38 the built-in ploidy
table (unprefixed `X`/`Y`) never matches, so every haploid GT is expanded by allele duplication
(`1`→`1|1`, `0`→`0|0`, `.`→`./.`) while already-diploid autosome **and PAR** calls are left
untouched (verified empirically). The parameter-free default is thus both **sufficient and
necessary**: PAR stays diploid in males *because it is already diploid, not because of a mask*,
so **no sex file and no PAR mask are needed** here — and a "biologically correct" male/non-PAR-X
= ploidy-1 config would instead *keep* the call haploid and re-crash RBCeq2. The §5.6 sex/PAR
machinery is therefore **not** used by this SNV step (it is only for CNV direction, §5.4.2). The
autosomal SNV branch is otherwise done.
> **Relocate `+fixploidy` when the merge lands.** #128 places it in the SNV pipe because today's
> input is SNV-only. Once the SV/CNV merge exists, Manta SV and DRAGEN CNV records can also carry
> haploid GT on non-PAR chrX/chrY in males, so **move** the single `+fixploidy` onto the merged VCF
> (§5.5) — one invocation covering all three sources, rather than one per source.
> ⚠️ **Interaction with PR #124:** that current filter also hard-drops genotypes below
> `DP≥20`/`GQ≥30`, which PR #124 shows manufactures false wild-type calls. Whatever PR #124
> lands as the conversion behaviour is what this spec's SNV branch inherits — this spec does
> **not** re-decide it, but the merged VCF must use the post-#124 SNV VCF.

### 5.3 SV VCF → normalised
No type change needed (Manta already writes `SVTYPE=DEL/DUP/INS/BND`). Housekeeping only:
region-restrict to the bg regions (optimisation; RBCeq2 also filters internally), ensure the
sample column name matches the merged VCF, keep `CIPOS/CIEND/SVLEN/END`, sort/bgzip/tabix.
`DUP:TANDEM` reported as `<INS>` is correct (matches DB INS/dup tokens).

### 5.4 CNV VCF → fixed (the one real transform)
Three deterministic edits:
1. **Drop `DRAGEN:REF:` records** (ALT `.`, no `SVTYPE`) — non-events.
2. **Rewrite `SVTYPE=CNV` → `DEL`/`DUP`.** Primary strategy: take the direction from the
   symbolic ALT (`<DEL>`→`DEL`, `<DUP>`→`DUP`). This is why the fix is required, not
   optional — see §2.3 (`"DEL" == "CNV"` is False; the ALT fallback only fires when
   `SVTYPE` is absent). On sex chromosomes, prefer deriving direction from `CN` relative to
   the sample's **expected ploidy** (§5.6) rather than trusting a diploid-assumed ALT.
3. **`cnvLength` (<10 kb) policy — resolved (Q1).** Sub-10 kb CNV records are not dropped as a
   class. They go through §5.5's triage like every other record: kept if database-relevant and
   ≥ 3 bins, yielding to a Manta record describing the same event, FILTER rewritten per record
   and the original preserved in `INFO/SVFILTER`. Never `--no_filter`. Evidence: Manta missed
   one of the two Gerbich 3.6 kb deletions in 100 genomes and the CNV caller found it on 3 bins
   (study §11); 1–2-bin records have Manta support in 45% of cases against 75–79% for 3–4 bins
   (study §10).
4. **Karyotype gate on chrX/chrY (§5.6).** Drop, or pass through tagged `SVSRC=CNV_KARYOTYPE`,
   every chrX/chrY CNV record from a sample whose DRAGEN ploidy estimate is not `XX` or `XY`.

`CnvRewriteStats` in the skeleton is the intended QC counter shape (ref-blocks dropped,
rewritten DEL/DUP, dropped sub-threshold, unresolved).

### 5.5 Merge: triage to the database, then one record per event (Q2 resolved)
`bcftools concat` the three normalised VCFs → sort → bgzip → tabix, into
`<sg>.rbceq2_input.vcf.gz`. All three must share the same sample column name and `chr`-prefixed
hg38 contigs (RBCeq2 strips `chr` internally).

**Why the merge has to arbitrate.** rbceq2 keeps one db definition per *record*, not per locus,
so a deletion present as both a Manta and a CNV record, offset by the CNV caller's bin snapping,
is matched to two different alleles and the sample is then read as carrying two null alleles.
The study's synthetic test, drawn:

```
one heterozygous 3.6 kb Gerbich deletion in the sample
                 |------------------------------|            truth: one event, one allele

db, seven GE alleles differing only by breakpoint (each 3.6 kb):
GE*01.-02.04  |------------------------------|
GE*01.-02.02    |------------------------------|
GE*01.-02.03          |------------------------------|
GE*01.-02.01                |------------------------------|
GE*01.-03.03                  |------------------------------|
GE*01.-03.02                          |------------------------------|
GE*01.-03.01                              |------------------------------|

records reaching rbceq2 for that one event:
Manta   POS exact (CIPOS 0-50)   |------------------------------|   -> best db token: GE*01.-02.01
CNV     POS bin-snapped +300 bp     |------------------------------| -> best db token: GE*01.-03.03

select_best_per_vcf: "best db definition for THIS RECORD"      x 2 records
                                                                 = GE*01.-02.01 / GE*01.-03.03
                                                                 = compound heterozygote, Ge:-2,-3
one record (either caller) -> one allele -> GE*01.-02.01 / GE*01 heterozygote, correct
```

The tie error would need the two records to have identical coordinates; they never do. 2.4.4's `ambiguous_equal_best_sv_evidence` error does not catch this: it fires only on
identical coordinates with conflicting GT or FILTER, which the two callers never produce
together (0 of 89 shared events in 100 genomes had identical breakpoints). That is a gap in
rbceq2 worth raising upstream, and until it is closed the merge must guarantee one record per
event. Empirical basis in `research/sv_cnv_overlap_50_samples.md`.

**Triage, not caller ownership.** The first draft of this resolution gave Manta the sub-10 kb
band and the CNV caller the rest. The study's read-level check (§11 there) showed Manta missing
a real 3.6 kb Gerbich deletion that the CNV caller found on 3 bins, so neither caller can own a
band. Instead:

1. **Restrict both structural VCFs to database-relevant records.** Keep a Manta or CNV record
   only if it (a) overlaps a db SV definition within rbceq2's own positional and length
   tolerance (`SvMatcher` defaults; use the same code), or (b) is a **PASS deletion under 1 Mb**
   spanning any defining SNV/indel site in `bg_site_systems.<genome>.tsv`. Everything else is
   discarded here, which also removes the whole-arm Manta `DUP:TANDEM` artefacts and, for (b),
   the megabase karyotype events (§5.6 handles those) and a 126 Mb `MaxDepth` Manta record seen
   once over AUG and RHAG.
2. **Keep every surviving record from either caller**, whatever its FILTER, with two edits:
   a `cnvLength` CNV record is dropped if it has fewer than **3 bins** (`BC < 3`; Manta
   corroborates 45% of 1–2-bin deletions against 75–79% of 3–4-bin ones), and a surviving
   non-PASS record has its FILTER rewritten to PASS and the original recorded in
   `INFO/SVFILTER`. `--no_filter` is never used; it is global.
3. **Where two surviving records describe one event** (same direction, reciprocal overlap
   ≥ 0.5), **keep the Manta record.** Its breakpoints are exact to `CIPOS` (0–50 bp) where the
   CNV caller's are bin-snapped by 0.2–4.6 kb, and the sub-10 kb targets (seven GE alleles,
   three A4GALT, the GYP cluster) differ from each other only by breakpoint. Record the dropped
   partner's ID in `INFO/SVPARTNER` so the QC can report that both callers agreed.
4. **Tag every kept record** with `INFO/SVSRC` (`MANTA`, `CNV`, or `CNV_LOWRES` for a CNV-only
   record under 10 kb) so the QC (§11) can grade the allele it produces.
5. **Assert no two surviving records share CHROM, POS, END and SVTYPE.** Fail the sample.

In 150 genomes rule 1(a) kept 0 to 2 records per genome (3 in total: the two Gerbich records
and the Gerbich CNV-only record) and rule 1(b) another 0 to 3, so the merged file carries a
handful of structural records, not hundreds. Rule 3 fired only on the Gerbich case, and the only
allele actually called was the Gerbich deletion. Rule 1(b) is where a new observation sits: a
recurrent 9–13 kb deletion at chr19:48.69 Mb spanning a FUT2 defining site, in 5 of 150 genomes,
seen by both callers, matching no db allele. rbceq2 ignores it (§11); the QC must not.

> **Haploid GT beyond SNVs.** The SPEC previously proposed moving `+fixploidy` to the merged VCF.
> Two facts have changed: haploid GT appeared on 1–2 chrX CNV records per 50 genomes, and 2.4.4
> handles haploid GT natively (Q4). Decide the `+fixploidy` question once, on the pin bump, for
> SNV and structural records together; do not add a second invocation here.

### 5.6 Sex / ploidy
For X-linked systems (XK/Kx, XG, CD99), expected copy number depends on the sample's sex
**and the region** (PAR vs non-PAR), and RBCeq2 infers zygosity from GT without knowing it. Read per-SG sex from
`ploidy_estimation_metrics.csv` / `.ploidy.vcf.gz` (X/Y coverage ratios). **Do not** use the
`##referenceSexKaryotype` header — it is a reference/config constant reading `XXYY` for
every sample (verified 102/102 across the cohort). Feed resolved ploidy into the CNV
direction logic (§5.4.2) — that is its **only** consumer. The SNV GT diploid-isation (§5.2, Q4)
does **not** need it: parameter-free `bcftools +fixploidy` handles haploid GT with no sex or PAR
input (PR #128).

> **Expected CN is region × sex — not a blanket "chrX = 1 in males" (from `db.tsv`).** Of the
> seven chrX structural alleles, four are in **PAR1** — CD99\*01N.01/02 (~2.71 Mb) and
> XG\*01N.02/03 (~2.78 Mb, on the PAR1 boundary) — which is **diploid in males** (baseline CN=2).
> Only **XK** (37.7 Mb) and **ATP11C** (139.7 Mb) are genuinely hemizygous (CN=1 in males). A
> naive haploid-X rule would expect CN=1 in PAR and miscall the normal CD99/XG state as a
> deletion. So both the CNV direction (§5.4.2) and the Q4 fix-up (§5.2) must be **PAR-mask
> aware**: diploid in PAR1/PAR2, hemizygous only in non-PAR X/Y for males (chrY: 1 in males,
> 0 in females).

> **Karyotype gate (added 2026-09-18).** In 100 genomes the 36 XY samples produced no large chrX
> CNV event, so DRAGEN's caller handles a normal male. The one sample it estimated **X0** carried
> PASS heterozygous CN=1 deletions of 1.4–10.6 Mb across XK, CD99, XG and ATP11C, and the one
> **XYY** genome a PASS CN=3 duplication across PAR1 (CD99, XG). Both escaped a false allele only
> because the length gate rejected megabase events against 11–219 kb tokens; a shorter segment
> would not be rejected. So the merge reads `Ploidy estimation` from
> `ploidy_estimation_metrics.csv` and, for any value other than `XX` or `XY`, drops or tags the
> sample's chrX/chrY CNV records (§5.4.4) and the QC reports the X-linked systems as
> `KARYOTYPE:<estimate>` rather than assessing them.

> ⚠️ **Haploid GT on sex chromosomes (Q4 — resolved in [PR #128](https://github.com/populationgenomics/ourdna_genomic_atlas/pull/128)).**
> DRAGEN emits *true haploid* `GT` on non-PAR chrX / chrY for male samples (e.g. `1`, not the
> pseudo-diploid `1/1` GATK writes). RBCeq2 assumes diploid GT and **crashes** on haploid input
> (`assert len(GT) == 3`, `data_procesing.py:878`). The SNV branch (§5.2) diploid-ises with a
> parameter-free `bcftools +fixploidy` in the conversion pipe (`1`→`1|1`, `0`→`0|0`); already-
> diploid PAR/autosome calls are untouched, so **no** sex or PAR mask is applied — see Q4.

### 5.7 Wire into the caller stage
`GenotypeBloodGroupsWithRbceq2` (`stages.py:492`) changes only its input: `--vcf` now points
at `<sg>.rbceq2_input.vcf.gz`. **Do not** use `--no_filter` to admit sub-10 kb CNVs — it is
global (would also let through non-PASS SNVs/SVs); if Q1 resolves to keep them, rewrite those
CNV records' `FILTER`→`PASS` in preprocessing instead (§5.4.3). Leave
`--phased` off (never required for detection) and `--RH` off (out of scope, §4).

---

## 6. Touch points

| File / symbol | Change |
|---|---|
| `src/ourdna_genomic_atlas/stages.py` | New `PreprocessDragenForRbceq2` SG stage; repoint `GenotypeBloodGroupsWithRbceq2.required_stages` + `--vcf` at its output |
| `src/ourdna_genomic_atlas/jobs/` | New job module implementing the skeleton (`fix_dragen_cnv_vcf`, `normalize_sv`, `merge_variant_vcfs`, `validate_for_rbceq2`) |
| `src/ourdna_genomic_atlas/popgen_utils.py:32` | Reuse `get_dragen_output_path` for `dragen_metrics/<sg>/<sg>.{sv,cnv}.vcf.gz` and ploidy |
| `config/ourdna_genomic_atlas_default_config.toml` | New `[workflow.preprocess_dragen_for_rbceq2]` block: `keep_subthreshold_cnv`, `cnv_svtype_from`, `min_size`, resources |
| `resources/bg_regions.GRCh38.bed` | Reused unchanged for region-restrict |
| rbceq2 image `2.4.1-1` | Possibly bump — research referenced 2.4.2 (version TBC) |

---

## 7. Open questions for reviewers

- **Q1 — `cnvLength` (<10 kb) policy → RESOLVED 2026-09-18, see §5.4.3.** Original text kept for the record: Drop sub-10 kb CNVs and source small events from the
  SV VCF (RBCeq2 default PASS-only), or keep them? **`--no_filter` is not the mechanism** —
  it is *global* and would also admit non-PASS SNVs and SVs (flagged in review); to keep
  sub-10 kb CNVs, rewrite just those records' `FILTER`→`PASS` in preprocessing (§5.4.3).
  **From `db.tsv`:** only 6 DEL alleles fall in the 1–10 kb band, and 5 are RH/GYP (paralog
  loci out of scope, §4); the **one in-scope** sub-10 kb structural allele is **XK\*N.05, an
  8 kb XK deletion** (non-paralog locus, squarely in Manta's range → the SV VCF supplies it).
  (Recommendation: SV VCF for <10 kb, keep CNV PASS-only; the `cnvLength`→`PASS` rewrite is
  only needed if XK\*N.05 proves weak in the SV VCF — check in concordance, §9.)
- **Q2 — SV∩CNV overlap → RESOLVED 2026-09-18, reversed.** Triage both files to database-relevant
  records and keep one record per event, Manta's where both callers describe it (§5.5). The original recommendation (defer to `select_best_per_vcf`) was
  wrong on two counts: that function keeps one db definition per record, so a duplicated event
  yields two alleles, and 2.4.4 raises on exact-coordinate ties. Empirical basis in
  `research/sv_cnv_overlap_50_samples.md` §§2–5, replicated in §10.
- **Q3 — CNV direction source.** Trust the `<DEL>`/`<DUP>` ALT everywhere, or derive from
  `CN` vs expected ploidy on autosomes too? (Recommendation: ALT on autosomes, `CN` vs
  **region × sex** ploidy on chrX/chrY — endorsed in review. NB PAR-aware: CD99 + XG sit in
  PAR1 → diploid in males; only XK/ATP11C are hemizygous — §5.6.)
- **Q4 — Haploid GT on sex chromosomes → RESOLVED (implemented in [PR #128](https://github.com/populationgenomics/ourdna_genomic_atlas/pull/128)).**
  DRAGEN writes haploid `GT` (`1`) on non-PAR chrX/chrY for male samples; RBCeq2 assumes
  diploid GT and **hard-crashes** on haploid input. Confirmed in source (v2.4.2): the zygosity
  determiner `get_ref` does `assert len(GT) == 3` (`core_logic/data_procesing.py:878`) and runs
  for every defining variant via `make_variant_pool` (`data_procesing.py:454`), so a haploid
  `"1"` (len 1) raises `AssertionError`. Corroborating diploid assumptions: `remove_home_ref` /
  `get_variants` drop only `"0/0"`, not haploid `"0"` (`IO/vcf.py:120,261`), and
  `split_vcf_to_dfs` asserts a `/`|`|` separator at index 1 (`IO/vcf.py:299`). **Fix (PR #128):**
  a final `| bcftools +fixploidy` in the conversion pipe, **with no `-s`/`-p` arguments**. On
  `chr`-prefixed hg38 the built-in (unprefixed) ploidy table never matches, so every haploid GT
  is expanded to diploid (`1`→`1|1`, `0`→`0|0`, `.`→`./.`) while diploid autosome/PAR calls stay
  untouched. This is **necessary, not just convenient**: a sex-aware config (male non-PAR X =
  ploidy 1) leaves the call haploid and re-crashes RBCeq2 (both behaviours verified empirically).
  So the §5.6 sex/PAR machinery is **not** used here — only for CNV direction. **Caveat
  (unchanged):** `1`→`1|1` reads as HOM (dosage 2, `core_logic/alleles.py:284`), overstating a
  truly hemizygous call — fine for detection, but relevant to zygosity-dependent filters;
  RBCeq2's native HEM status is deletion-derived, not from GT. Missed by the XX 1KG fixtures
  (§9). **Residual:** #128 covers SNV records only; when the merge stage (§5.5) lands, re-apply
  `+fixploidy` to the merged VCF.

**Change list, 2026-09-18 revision** (design PR; the factual refresh of 2026-09-17 was PR #15):
- §5.4.3 Q1 resolved: `cnvLength` records are triaged like any other, not dropped as a class; 3-bin floor; per-record FILTER rewrite with the original kept.
- §5.4.4 and §5.6: karyotype gate on chrX/chrY CNV records for non-XX/XY samples.
- §5.5 Q2 resolved and reversed: triage both structural VCFs to database-relevant records, keep one record per event with Manta's breakpoints where both callers agree, 3-bin floor for CNV records, assert uniqueness. Raise the per-record-not-per-locus matching upstream.
- §5.5: no second `+fixploidy`; decide once on the 2.4.4 pin bump.
- §9: synthetic fixtures shaped like the study's observations; 100-genome concordance scripts. No cohort genome is a fixture.
- §11 new: QC for structural calls, with `SVNOCOV`, `SVLOWRES`, `SVUNASSESSED`, `SVDEL`, `KARYOTYPE` flags and thresholds in Analysis meta; rbceq2 ignores a deletion it cannot match, so the QC must report deletions over defining SNV sites itself.

**To raise upstream with the RBCeq2 maintainers** (feature requests, not things this design works
around beyond reporting):

1. **Matching is per record, not per locus.** Two records for one event yield two alleles (§5.5
   figure). Our merge guarantees one record per event; a per-locus selection in
   `select_best_per_vcf`, or a warning when two records match overlapping tokens of one system,
   would make that unnecessary.
2. **Unmatched deletions are ignored.** A deletion that matches no db token never enters the
   variant pool, so the hemizygosity adjustment for SNVs inside it never runs and a single-copy
   genotype is read as homozygous. Our QC reports it (`SVDEL`, §11); rbceq2 could apply the same
   adjustment it already applies to matched deletions.
3. **Undefined null alleles.** A deletion removing coding sequence of a blood-group gene is a null
   whether or not the db names it. The study found a recurrent 9–13 kb deletion over a FUT2
   defining site in 5 of 150 genomes, matching no token. Two asks: consider it as a candidate
   allele for the db, and consider a generic "coding deletion in system X, unnamed" call. We do
   not build gene-model resources for this here; it stays an observation to bear in mind, and a
   reason the `SVDEL` flag names the system.

---

## 8. Risks & considerations

- **Short-read paralog loci.** Beyond the RH/GYP hybrids already excluded, any blood-group
  gene with a close paralog risks mismapped/absent CNV calls. Concordance across the five
  replicates a necessary check/analysis to do for QC.
- **Bin-CNV false positives.** DRAGEN CNV is bin-based; small events are noisy (hence the
  `cnvLength` filter). Keeping sub-10 kb (Q1) raises sensitivity but risks spurious dosage calls.
- **Haploid sex-chromosome GT (confirmed — Q4).** DRAGEN's haploid `GT` on non-PAR chrX/chrY
  (male samples) makes RBCeq2 **crash** (`assert len(GT) == 3`, `data_procesing.py:878`); it
  assumes diploid GT. Mitigated by diploid-ising in preprocessing (§5.2); residual: `1`→`1/1`
  reads as HOM, overstating dosage for truly hemizygous calls.
- **Do we need positive controls?** The five replicate control samples showed *no* blood-group
  CNVs at GYP/XK/RHD in the earlier data look — good for a concordance/regression baseline,
  but does **not** exercise a true-positive large deletion. We may need a synthetic or known-
  positive sample to prove a real gene deletion is called (see §9).
- **Merge correctness.** Contig/sample-name mismatches or unsorted concat will silently break
  RBCeq2 region fetch; `validate_for_rbceq2` must assert sorted, single-sample, `chr`-prefixed,
  indexed.

## 9. Testing plan

Based on a small extract from one of OurDNA's five 1000 Genomes (1KG) control replicates — independent sequencings of `NA12878`, a public 1000 Genomes control. Note `NA12878` is **XX**, so these fixtures cannot exercise the haploid male chrX/chrY GT path (Q4) — that needs a separate male sample.

> **Fixtures, not live paths.** Unit tests run against **small checked-in extracts** of one
> replicate's DRAGEN outputs — a handful of records per source (CNV: a `DRAGEN:REF:` block plus
> `<DEL>`/`<DUP>` events; SV: a DEL/INS; gVCF: a few called sites plus one `<NON_REF>` block; plus a
> chrX record for the ploidy path), committed as fixtures under `tests/`. **CI cannot read
> `gs://cpg-ourdna-main/…`**, so no unit test may resolve a live GCS path. The five-replicate
> concordance run (below) *does* need the real files, so it is a **manual/analysis check**, not CI.

- **Unit:** `fix_dragen_cnv_vcf` on the replicate's CNV **extract** — assert REF records dropped,
  every surviving `SVTYPE` ∈ {DEL,DUP}, counts match `CnvRewriteStats`; `<DEL>`→DEL / `<DUP>`→DUP.
- **Unit:** `merge_variant_vcfs` output is sorted, single-sample, tabix-indexed, `chr`-prefixed.
- **Integration (concordance — manual, needs GCS):** run the full preprocess+rbceq2 on all five
  replicates; structural calls must be
  **identical across all five** (same individual). Any divergence is a bug or a QC signal
  (cf. PR #124's SNV divergence table).
- **Sub-10 kb path, synthetic fixtures shaped like the study's observations (no cohort data in
  the repo):** (a) a Manta DEL record of exactly 3609 bp at the GE\*01.-02.01 coordinates plus a
  4961 bp `cnvLength` CN=1 CNV record over it → the merge keeps the Manta record and records the CNV
  partner, and rbceq2 calls one GE allele, never two; (b) the CNV record alone, 3 bins → kept, PASS-rewritten, tagged,
  and the QC flags the GE call `SVLOWRES`; (c) the same with 2 bins → dropped, GE unassessed.
  The study (§§3, 10, 11) observed both (a) and (b) in real genomes; the fixtures are constructed
  records, not extracts.
- **Karyotype gate, synthetic:** a CNV record set with megabase CN=1 deletions across XK, CD99/XG
  and ATP11C paired with a ploidy metrics file reading `X0`, and a PAR1 CN=3 duplication paired
  with `XYY` → no X-linked allele is called and the QC reports `KARYOTYPE:<estimate>`.
- **Concordance across 100 genomes:** the study's five scripts, once committed, are the check
  that a code change to the merge does not alter the per-band record counts or the single GE
  hit under the `pass` policy. They read live GCS paths, so this is a manual check, not CI.
- **True-positive (DB-derived targets), still wanted:** a sample with a known 10 kb+ deletion —
  assert the allele is called. Concrete targets: **XK\*N.05** (8 kb, non-PAR chrX → must come via
  the **SV VCF**; exercises the <10 kb path, Q1) and a larger non-PAR del — **XK\*N.01** (53 kb) or
  **ATP11C\*01N.01** (219 kb) — from the **CNV VCF** (exercises the hemizygous ploidy path). Needed
  because the control shows no positive structural hit (§8).
- **PAR-vs-hemizygous ploidy (needs a male sample):** confirm the region-aware expected CN (§5.6)
  — CD99 / XG (PAR1, diploid in males) are **not** miscalled as deletions, while XK / ATP11C
  (non-PAR) are handled hemizygously. The XX 1KG replicates cannot provide this check.

## 10. Expected outcome

RBCeq2 gains access to the ~40 large-CNV and ~38 large-indel blood-group alleles it currently
cannot call, using data DRAGEN already produces, with **no change to RBCeq2** and one new,
well-tested preprocessing stage. Hybrid RH/GYP alleles remain explicitly out of reach on
short-read data.

## 11. QC for structural calls (new, 2026-09-18)

`FlagBloodGroupCallQc` asks, per defining SNV site, whether the caller looked and how well, using
gVCF records and reference blocks. Structural alleles need the same two answers from different
evidence, and the two sources answer differently:

| Question | CNV VCF (10 kb+ targets) | gVCF depth (sub-10 kb targets) |
|---|---|---|
| Did the caller assess the region? | Yes: `DRAGEN:REF:` records tile every 10 kb+ target in 98–100% of genomes with 16–195 bins (study §7). A gap, or a `cnvQual` event, is the analogue of `NOCOV`. | Not from the structural VCFs: Manta emits events only, and a sub-10 kb target holds 1–7 bins. But the **gVCF** already tiles the interval with reference blocks carrying `DP`/`MIN_DP`, and the conversion stage reads it once. Mean DP over the target against its flanks is the same signal the study's read-level check used (§11 there): about half for a heterozygous deletion, near zero for homozygous, flat for none. |
| How good is the call? | `QUAL`, `CN`, `SM`, `BC`, FILTER. | Manta: `QUAL`, `PR`/`SR`, `CIPOS`. CNV-only: `BC`, `SM`, `QUAL`, and the fact that Manta saw nothing. |

Design:

1. **The site-system map gains interval rows.** `bg_db.py` stops dropping `kind == 'sv'`; each
   SV definition becomes a `(chrom, start, end, system, allele)` row. The QC job reads the
   merged VCF's structural records (or rbceq2's debug log, which in 2.4.4 names the source
   record for each SV match) rather than the DP/GQ extract for these rows.
2. **New flags**, joined with `+` to provenance like today's:
   - `SVNOCOV:<system>(<allele>,gap=<bp>)`: a 10 kb+ target not tiled by CNV records.
   - `SVLOWRES:<system>(<allele>,src=CNV,BC=<n>,SM=<x>,QUAL=<q>)`: an allele called from a
     §5.4.3 CNV-only sub-10 kb record. Provisional by construction.
   - `SVDEPTH:<system>(<allele>,ratio=<x>,DP=<n>,flank=<n>)`: gVCF depth over a sub-10 kb target
     is below `sv_depth_ratio` of its flanking depth and no kept record explains it. The
     dosage-drop signal without a breakpoint call; what the missed Gerbich deletion would have
     raised had Manta and the CNV caller both been silent.
   - `SVUNASSESSED:<system>`: a sub-10 kb target with no gVCF record over the interval at all
     (an exome hole, or an unmapped region), so neither a call nor its absence can be judged.
   - `KARYOTYPE:<estimate>`: X-linked systems in a non-XX/XY sample (§5.6).
   - `SVDEL:<system>(<site>,del=<chrom:pos-end>,src=<caller>,GT=<gt>)`: a kept deletion, matched
     to a db allele or not, spans a defining SNV/indel site of the system. **rbceq2 does nothing
     with an unmatched deletion**: a structural record enters its variant pool only under the db
     token it matched, and the zygosity adjustment that turns a homozygous call inside a
     deletion into hemizygous (`modify_variant_pool_if_large_indel`) only sees pool entries. So a
     gVCF `A/A` under a heterozygous deletion that matches no db SV is reported homozygous, the
     same silent-wrong-call class as absent-means-reference. The QC has to say so, as it already
     does for a small deletion that removed the base (`DEL`).
   - A structural call that passes everything is not listed, so `PASS` keeps meaning "nothing to
     report".
3. **The defining-sites extract grows by the sub-10 kb target intervals.** `gen_bg_resources.py`
   emits them as a second BED from the same db parse; `FilterAndConvertGvcfsForRbceq2` extracts
   `DP`/`MIN_DP`/`END` over them and over a flank each side (say 5 kb) in the same pass it already
   makes for the SNV sites. No new input, no new job.
4. **Thresholds in Analysis meta**, as `min_depth`/`min_gq` are today: `sv_min_bins = 3`,
   `sv_recip_overlap = 0.5`, `sv_lowres_max_bp = 10000` (the size below which a CNV-only record
   is graded `CNV_LOWRES`), `sv_depth_ratio = 0.7`, `sv_depth_flank_bp = 5000`.
5. **What gVCF depth does not give**, accepted for the first cut: the MAPQ0 fraction and
   discordant-pair evidence a CRAM read would. DRAGEN's DP already excludes reads failing its
   mapping filters, so poor mappability appears as low depth rather than as its own signal, which
   a flag can live with. Breakpoint confirmation is the callers' job, not the QC's.
