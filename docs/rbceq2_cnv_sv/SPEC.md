# Design spec: structural-variant support for RBCeq2 blood-group calling

**Status:** draft, awaiting design approval; do not implement.
**Revision:** 25 September 2026, after the first review round on `popgen_rbceq2#15` and `#16` and the `+fixploidy` removal in `#20`. **Author:** Joshua Schmidt. **Reviewers:** Alexander Stuckey.
**Reader:** a pipeline engineer or scientist on the team who has not worked on the rbceq2 stages.
**Decision asked of reviewers:** approve the triage rule and karyotype gate in §5 and the thresholds in §8, or say which to change and why.
**Area:** rbceq2 blood-group genotyping pipeline (`FilterAndConvertGvcfsForRbceq2` to `GenotypeBloodGroupsWithRbceq2` to `CombineRbceq2OutputsPerCohort`), in `popgen_rbceq2`. rbceq2 is pinned at 2.4.4 (`constants.py`, `popgen_rbceq2#17`), and the database figures below are from that release.
**Pull-request references:** the stages were ported from `ourdna_genomic_atlas` in August 2026 and both repos number PRs from 1, so every PR below is written `ourdna_genomic_atlas#n` or `popgen_rbceq2#n`.
**Callers:** the SV VCF is written by the DRAGEN 3.7.8 SV caller, which Illumina describes as integrating and extending Manta; its record IDs keep Manta's prefix (`MantaDEL:`). This spec says "SV caller" and "SV VCF", not "Manta", so that nothing here rests on Manta's internals. The CNV VCF is written by DRAGEN's bin-based CNV caller.
**Companion docs (same directory):**
- [`dragen_three_source_merge.md`](./research/dragen_three_source_merge.md): the visual overview and mermaid diagram; read this first for the shape.
- [`implement_cnv_rbceq2_research.md`](./research/implement_cnv_rbceq2_research.md): the July 2026 RBCeq2 source analysis this spec rests on, read against v2.4.2; where 2.4.4 moved a fact, this spec says so.
- [`resolvability_by_input_class.md`](./research/resolvability_by_input_class.md): per-system estimate of what the gVCF resolves today and what the merge adds, from the 2.4.4 database.
- [`sv_cnv_overlap_50_samples.md`](./research/sv_cnv_overlap_50_samples.md): the 150-genome SV-caller-versus-CNV-caller study the triage rule and karyotype gate rest on.

---

## 1. Summary

RBCeq2's 2.4.4 database defines 67 non-RH alleles across 20 systems by a single large structural event. Our pipeline calls none of them today, because it feeds RBCeq2 only the SNV gVCF (`research/resolvability_by_input_class.md`). DRAGEN already produces the two structural files needed, an SV VCF and a bin-based CNV VCF, but RBCeq2 reads a single VCF. The design is therefore a per-sequencing-group (SG) preprocessing stage that merges the three files into one.

The merge does three things beyond concatenation. It rewrites DRAGEN's `SVTYPE=CNV` to `DEL` or `DUP`, without which no large deletion matches a database allele. It triages both structural files to records that could match a database allele or that span a defining SNV site. It keeps one record per event, because RBCeq2 matches per record and a deletion seen by both callers would otherwise yield two alleles. It drops or tags chrX and chrY CNV records for samples whose DRAGEN karyotype estimate is not XX or XY.

A new QC design (§6) reports whether the callers assessed each structural target and how well, and reports deletions over defining SNV sites, which RBCeq2 itself ignores.

The evidence is a 150-genome OurDNA study (`research/sv_cnv_overlap_50_samples.md`). Across those genomes the design calls one structural allele, a Gerbich deletion. The value is correctness on rare alleles and honest QC, not call volume; §10 states the expected outcome and success criteria in those terms.

RBCeq2 is not modified. RH and GYP hybrid alleles stay out of scope on short-read data, and `--RH` stays off, on the maintainer's advice that DRAGEN SV/CNV does not reliably detect the RH hybrids.

Exome SGs are in scope. Production exomes carry an SV VCF and a panel-of-normals CNV VCF (§3), so the merge treats them as it treats genomes; §5.4 and §6 handle how the exome CNV file differs. Where a run produced no CNV VCF, as the 10-exome test run did, the merge is gVCF plus SV VCF and every 10 kb+ target is reported unassessed.

---

## 2. Background: how RBCeq2 consumes structural variants

### 2.1 One VCF, no CNV/SV flag

`find_hits` (`src/rbceq2/main.py:272-305`) reads SVs out of the same dataframe as SNVs: `SvReader(df=vcf.df, min_size=args.min_size).events()`, then `SvMatcher.match(db_defs, events)`, then `select_best_per_vcf(...)`. There is no CLI flag and no second input; structural records must be present in the single `--vcf`.

### 2.2 What the DB encodes

Structural alleles live in the `GRCh37`/`GRCh38` columns of `db.tsv` as word-form tokens `<pos>_<type>_<len>` (for example `143914828_del_110kb`), or as explicit sequences `<pos>_<REFseq>_<ALTseq>` that `parse_db_token` reads as a DEL or INS of `|len(ALT)-len(REF)|`. In the 2.4.4 GRCh38 column there are 74 DEL, 22 INS and 1 DUP word-form tokens, plus 33 DEL and 8 INS sequence-form tokens of 50 bp or more; no `INV`, `BND` or literal `<CNV>`. The RH tokens are bare base counts (`_DEL_59419`) and the GE, PIGG and MAM deletions are spelled as sequences, so the v2.4.2 counts in `implement_cnv_rbceq2_research.md` §3 no longer match. Per allele across the non-RH database: 67 single large events of 1 kb or more, 16 mid indels of 50 bp to 1 kb, 4 GYP hybrids (`resolvability_by_input_class.md`). Hybrids are stored as paired tokens, not a single hybrid type.

### 2.3 Matching is fuzzy in position and length, but strict on type

`SvMatcher` (`large_variants.py:159-367`) uses adaptive positional and length tolerance and reciprocal-overlap gates, so imprecise breakpoints still match, but `require_same_type=True` by default, so the DB token type must equal the event `SVTYPE`. The symbolic-ALT fallback (`<DEL>`/`<DUP>`) in `SvReader` fires only when `SVTYPE` is absent (`large_variants.py:739` in 2.4.2, `:755-767` in 2.4.4; `SvMatcher.compatible` is unchanged between them). This single fact is what forces the CNV rewrite in §5.4. RBCeq2 also pre-filters records to within 500 kb of a DB position (`vcf.py:497`).

---

## 3. Inputs

Grounded on one of the OurDNA 1000 Genomes (1KG) control replicates:

| # | Source | Path (`gs://cpg-ourdna-main/ica/dragen_3_7_8/output/`) | Caller | Ready? |
|---|--------|--------------------------------------------------------|--------|--------|
| 1 | SNV gVCF | `recal_gvcf/<sg>.hard-filtered.recal.gvcf.gz` | DRAGEN SNV | Convert to sites VCF |
| 2 | SV VCF | `dragen_metrics/<sg>/<sg>.sv.vcf.gz` | DRAGEN SV caller (extends Manta) | Compatible |
| 3 | CNV VCF | `dragen_metrics/<sg>/<sg>.cnv.vcf.gz` | DRAGEN bin CNV | Rewrite SVTYPE |
| - | Ploidy (sex) | `sg.meta['qc']` in metamist, else `dragen_metrics/<sg>/<sg>.ploidy_estimation_metrics.csv` | DRAGEN via `single_sample_qc_popgen` | Read only (§5.6) |

What each sequencing type has, checked in the buckets on 2026-09-18 (150 OurDNA genomes; the 11,917 mackenzie exomes in `cpg-mackenzie-main`; the 10 exomes of the DRAGEN 3.7.8 test run in `cpg-mackenzie-test`):

| File | Genome | Exome, production (`-main`) | Exome, test run (`-test`) |
|---|---|---|---|
| SNV gVCF | yes | yes | yes |
| SV VCF | yes | yes, 11,917 of 11,917 (one exome: 65 records, 4 in the blood-group regions) | yes, 10 of 10 (one exome: 103 records, 9 in the regions) |
| CNV VCF | yes | yes, 11,917 of 11,917 | no, 0 of 10: `cnv_metrics.csv` stops after "Number of target intervals", so the caller counted reads but never segmented |
| `ploidy_estimation_metrics.csv` | yes | yes | yes |
| `ploidy.vcf.gz` | yes | no | no |

The exome CNV VCF is not the genome file with fewer records. One production exome's header and records show four differences:

- It is called per capture target (`--cnv-target-bed`) against a panel of 100 normals (`--cnv-normals-list`; `Number of normal samples,100` in `cnv_metrics.csv`; self-normalisation off), with HSLM segmentation.
- Every record is an event. There are no `DRAGEN:REF:` records: 1,131 records, 30 of them in the blood-group regions, all events.
- FILTER values are PASS, `cnvQual`, `cnvBinSupportRatio` and `cnvCopyRatio`. There is no `cnvLength`, so sub-10 kb events pass: 129 PASS under 1 kb and 121 PASS of 1 to 10 kb in that exome.
- `BC` counts capture targets, not 1–2 kb bins.

FORMAT is the same `GT:SM:CN:BC:PE`, and duplications carry `./1` as in genomes. §5.4 and §6 say what follows from each difference.

**1, SNV gVCF** (called variant records, already carrying GT, interleaved with `<NON_REF>` reference blocks; the head below happens to show only ref blocks):

```text
chr1  9997   .  N  <NON_REF>  .  PASS  END=10017  ...  ./.:11,0:11:0:2:...
chr1  10018  .  C  <NON_REF>  .  PASS  END=10018  ...  0/0:17,0:17:15:17:...
```

**2, SV VCF** (DRAGEN SV caller; record IDs keep Manta's prefix; header ALTs `<DEL>`,`<INS>`,`<DUP:TANDEM>`; proper SVTYPE; FORMAT `GT:FT:GQ:PL:PR:SR`):

```text
chr1  789481  MantaINS:…  G  <INS>  999  PASS  END=789481;SVTYPE=INS;CIPOS=0,7;…
chr1  839442  MantaDEL:…  CACC…ACA  CT  463  PASS  END=839499;SVTYPE=DEL;SVLEN=-57;CIGAR=1M1I57D
```

**3, CNV VCF** (every record `SVTYPE=CNV`; direction only in ALT and `CN`; FORMAT `GT:SM:CN:BC:PE`; every PASS duplication in the study carries GT `./1`, every PASS deletion `0/1` or `1/1`):

```text
chr1  817861   DRAGEN:REF:…   N  .      81  PASS       END=2650427;REFLEN=1832567                     ./.:0.98:2:1434:…
chr1  2650427  DRAGEN:LOSS:…  N  <DEL>  51  cnvLength  SVLEN=-2648;SVTYPE=CNV;END=2653075;REFLEN=2648  1/1:0.058:0:2:…
chr1  3501568  DRAGEN:GAIN:…  N  <DUP>  88  cnvLength  SVLEN=1000;SVTYPE=CNV;END=3502568;REFLEN=1000   ./1:3.18:6:1:…
```

---

## 4. Goals and non-goals

The goals of this design are to:

- give RBCeq2 access to the structural blood-group alleles it cannot currently call, using the SV and CNV VCFs DRAGEN already produces alongside the SNV gVCF;
- merge the three sources into one record per structural event, so RBCeq2's per-record matching does not read a single deletion as two alleles; and
- add QC that reports whether the callers could assess each structural target, and flags deletions RBCeq2 itself ignores.

In scope

- New per-SG inputs: the SV and CNV VCFs, `cnv_metrics.csv` and the ploidy file, derived from the registered gVCF's DRAGEN output prefix and recorded in the merged VCF's Analysis meta. No new metamist analysis types and no cpg_flow extension (§5.1; Q6 decided).
- A preprocessing step that (a) produces an SNV sites VCF, (b) normalises the SV VCF, (c) fixes the CNV VCF (drop REF records, rewrite `SVTYPE=CNV` to `DEL`/`DUP`), and (d) merges all three into one sorted, bgzipped, tabixed VCF per SG.
- Exome SGs, through the same three branches. Production exomes carry a panel-of-normals CNV VCF (§3); where a run produced none, step (c) is skipped and the 10 kb+ targets are reported unassessed (§5.1, §6).
- Feed the merged VCF into `GenotypeBloodGroupsWithRbceq2` in place of the SNV-only VCF.
- Per-SG sex and ploidy from the QC meta that `single_sample_qc_popgen` already writes to metamist (§5.6). It feeds CNV direction on chrX/chrY, the karyotype gate (§5.4, §5.6), and the `KARYOTYPE` and `SEXCHECK` QC flags (§6).
- A cross-replicate concordance check over the OurDNA 1KG control SGs (§10).

Non-goals (explicit)

- No changes to RBCeq2: not `SvReader`, `SvMatcher`, `db.tsv`, or the region filter. If an option here appears to need a RBCeq2 change, it is rejected.
- RH and GYP hybrid alleles are out of scope on short-read data. RHD/RHCE and GYPA/B/E are near-identical paralogs, so mapping quality drops to zero and short-read SV and CNV calls are noisy or absent. `--RH` is documented as long-read only, and the maintainer's comparison of DRAGEN SV/CNV against matched long reads for the 2.4.4 release found deletions useful but hybrids not reliably detected. Revisit when the volume and timing of long-read data are known.
- No CNV calling of our own. Where a DRAGEN run produced no CNV VCF (the exome test run, §3) this design adds no caller; that SG's 10 kb+ targets stay unassessed and the QC says so (§6).
- No change to the blood-group science or allele definitions.
- No re-run or backfill orchestration (tracked separately).

---

## 5. Design

A new per-SG stage, working name `PreprocessDragenForRbceq2`, runs between `FilterAndConvertGvcfsForRbceq2` and `GenotypeBloodGroupsWithRbceq2` for every SG and emits one merged VCF. The rbceq2 stage is repointed at that VCF. For an SG whose run produced no CNV VCF (§5.1), the stage skips §5.4 and §5.5 triages the SV records alone. Execution is bcftools-based, in the same image family as the current filter stage. The mermaid diagram in `dragen_three_source_merge.md` shows the shape.

### 5.1 Resolve the new inputs

This step resolves, per SG, the SV and CNV VCF paths, `cnv_metrics.csv`, and the ploidy file that §5.6 falls back to.

The SNV gVCF comes from `sequencing_group.gvcf`, which cpg_flow populates from the gVCF analysis in metamist. Review (2026-09-18) ruled out extending cpg_flow for the other files, and confirmed that `dragen_align` will not register them before its Nextflow refactor. So the stage derives them from the gVCF's path. `<prefix>/recal_gvcf/<sg>.hard-filtered.recal.gvcf.gz` gives `<prefix>/dragen_metrics/<sg>/<sg>.sv.vcf.gz`, `<sg>.cnv.vcf.gz`, `<sg>.cnv_metrics.csv` and `<sg>.ploidy_estimation_metrics.csv`, a layout verified in both buckets (§3).

The SV VCF and `cnv_metrics.csv` are expected for every SG, and a missing one fails the SG, as the gVCF guard does today. Both VCFs need their `.tbi`.

Whether a CNV VCF is expected is read from `cnv_metrics.csv`, which every run writes. It is not a property of the sequencing type (production exomes have one, the exome test run does not, §3), and it is not a config default.

- A run that segmented carries `Number of segments` (and `Number of normal samples` when a panel was used), and must have a CNV VCF.
- A run whose metrics stop after `Number of target intervals` never called, and must not have one.
- Any other combination fails the SG.

The ploidy file is read only when `sg.meta['qc']` carries no ploidy estimate, and it may be absent for a genome with uneven coverage (§5.6). An SG with no estimate in meta, no ploidy file and no somalier signals in meta fails.

The derived paths are recorded in the merged VCF's Analysis meta, as the thresholds are (§8). No new metamist analysis types are created for the raw DRAGEN files (Q6 decided, §13).

### 5.2 SNV gVCF to sites VCF

The SNV branch is the `vcf` output of the existing conversion stage, `FilterAndConvertGvcfsForRbceq2` (`stages/blood_group_genotyping/filter_and_convert.py`), taken unchanged. This spec adds nothing to that stage and re-decides nothing about it; this section describes it for the reader.

The stage reads the gVCF once, restricted to `resources/bg_regions.<genome>.bed`; the restriction is unconditional, so the gVCF `.tbi` is required. `bcftools norm -m -any` splits multiallelics. From that intermediate the stage writes the defining-sites extract the QC stage reads, and the rbceq2 input. The rbceq2 input drops every `<NON_REF>` record (the reference blocks and the split-off symbolic twin of each variant) and trims now-unused ALT alleles. No reference FASTA and no genotyping are involved: a single-sample gVCF already carries GT at every variant site.

Genotypes are never rewritten. rbceq2 2.4.4 reads DRAGEN's one-token non-PAR chrX genotypes natively, and `popgen_rbceq2#20` removed the `bcftools +fixploidy` step that expanded them to diploid for 2.4.3 (§13, the 2026-09-21 entry). On an exome the post-hoc fill is kept out of single-copy chrX, so the file carries one ploidy per region.

The stage filters nothing on `FORMAT/DP` or `FORMAT/GQ`. rbceq2 reads a defining site absent from its input as confident homozygous reference, so dropping a borderline genotype would manufacture a wild-type call rather than a no-call. Since `ourdna_genomic_atlas#136`, rbceq2's own PASS-only rule handles DRAGEN's hard-filter names, and `FlagBloodGroupCallQc` reports DP and GQ per system instead (`ourdna_genomic_atlas#138`, `#140`). The merged VCF inherits this posture, and §6 extends the same QC to structural calls.

For an exome SG the stage also merges in post-hoc calls at defining sites outside the capture design (`popgen_rbceq2#14`). That output is the exome SNV branch, taken unchanged in the same way.

### 5.3 SV VCF to normalised

The DRAGEN SV caller already writes `SVTYPE=DEL/DUP/INS/BND`, so no type change is needed and this step is housekeeping:

- drop `BND` records, which pair by `MATEID` and match no db token;
- restrict to the bg regions (an optimisation; RBCeq2 also filters internally);
- make the sample column name match the merged VCF;
- keep `CIPOS`, `CIEND`, `SVLEN` and `END`, and FORMAT as written (§5.5 says why nothing needs harmonising);
- sort, bgzip and tabix.

`DUP:TANDEM` reported as `<INS>` is correct; it matches db INS and dup tokens.

### 5.4 CNV VCF to fixed (the one real transform)

This step makes four deterministic edits to the DRAGEN CNV VCF so RBCeq2 can match its records. It depends on the sex and ploidy read in §5.6, and is skipped for an SG whose run produced no CNV VCF (§5.1).

1. Drop `DRAGEN:REF:` records (ALT `.`, no `SVTYPE`); these are non-events.
2. Rewrite `SVTYPE=CNV` to `DEL` or `DUP`. On autosomes take the direction from the symbolic ALT (`<DEL>` to `DEL`, `<DUP>` to `DUP`); on sex chromosomes derive it from `CN` against the sample's expected ploidy (§5.6).
3. Keep sub-10 kb `cnvLength` records rather than dropping them as a class. §5.5 triages them like every other record, keeping those with at least 3 bins and rewriting FILTER per record, with the original preserved in `INFO/SVFILTER` (§13, Q1 resolved).
4. Apply the karyotype gate: drop, or pass through tagged `SVSRC=CNV_KARYOTYPE`, every chrX/chrY CNV record from a sample whose DRAGEN ploidy estimate is not `XX` or `XY` (open decision, Q5).

Edit 2 is what makes the CNV VCF usable at all. `SvMatcher` runs with `require_same_type=True`, so `"DEL" == "CNV"` is false, and the ALT fallback fires only when `SVTYPE` is absent (§2.3).

Edit 3 rests on two findings. The SV caller missed one of the two Gerbich 3.6 kb deletions in 150 genomes, which the CNV caller found on 3 bins (`sv_cnv_overlap_50_samples.md` §11). And 1–2-bin records have SV-caller support in 45% of cases, against 75–79% for 3–4 bins (same note, §10). The 3-bin Gerbich call is one record with `BC=3`, not three adjacent records (same note, §13); adjacent CNV records are copy-number steps, and no merge-adjacent step is added (§7).

An exome CNV VCF (§3) changes two things. It has no `DRAGEN:REF:` records, so edit 1 drops nothing and §6 cannot read assessability from the file. And it has no `cnvLength` filter, and its `BC` counts capture targets rather than 1–2 kb bins. Edit 3's `sv_min_bins` floor applies to that `BC` as written; whether three targets is the right floor is untested (§11). The other edits apply unchanged.

The job counts what it did per SG: reference records dropped, records rewritten to DEL and DUP, sub-threshold records dropped, and records left unresolved. `CnvRewriteStats` is the proposed name for that record; no code for it exists yet in either repo.

### 5.5 Merge: triage to the database, then one record per event

The merge concatenates the normalised VCFs, then triages the structural records so RBCeq2 sees one record per event. The QC in §6 grades what survives, and the thresholds live in §8.

`bcftools concat` the three normalised VCFs (two for an SG without a CNV VCF), then sort, bgzip and tabix into `<sg>.rbceq2_input.vcf.gz`. All inputs must share the sample column name and `chr`-prefixed hg38 contigs, since RBCeq2 strips `chr` internally.

FORMAT needs no harmonisation across the sources (review asked, 2026-09-18). VCF FORMAT is per record, and rbceq2 requires only that `GT` is the first key of every row (`IO/vcf.py`, `_require_first_gt`). All three sources meet that: gVCF-derived rows carry `GT:AD:DP:GQ:...`, SV rows `GT:FT:GQ:PL:PR:SR` (or without `SR`), and CNV rows `GT:SM:CN:BC:PE`, with no other string in the 150 study genomes. `SvReader` takes `SVTYPE`, `SVLEN`, `END`, `CIPOS` and `CIEND` from INFO. The tags the SV and CNV headers share have identical Number and Type in DRAGEN 3.7.8, and a real concat ran without a header warning (`sv_cnv_overlap_50_samples.md` §13). `validate_for_rbceq2` asserts GT-first on every row, alongside sorted, single-sample, `chr`-prefixed and indexed.

One genotype fact to check at implementation: every PASS CNV duplication in the study carries GT `./1` (244 of 244; deletions carry `0/1` or `1/1`). rbceq2 accepts the string, but how it scores a half-called GT on a matched DUP record is unverified. It matters for one token, GYP\*505.

The merge has to arbitrate because rbceq2 keeps one db definition per record, not per locus. A deletion present as both an SV and a CNV record, offset by the CNV caller's bin snapping, is matched to two different alleles, and the sample reads as carrying two null alleles:

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
SV      POS exact (CIPOS 0–50)   |------------------------------|   -> best db token: GE*01.-02.01
CNV     POS bin-snapped +300 bp     |------------------------------| -> best db token: GE*01.-03.03

select_best_per_vcf: "best db definition for THIS RECORD"      x 2 records
                                                                 = GE*01.-02.01 / GE*01.-03.03
                                                                 = compound heterozygote, Ge:-2,-3
one record (either caller) -> one allele -> GE*01.-02.01 / GE*01 heterozygote, correct
```

2.4.4's `ambiguous_equal_best_sv_evidence` error does not catch this. It fires only on identical coordinates with conflicting GT or FILTER, and the two callers never produce identical coordinates: 0 of 131 shared events in 150 genomes (`sv_cnv_overlap_50_samples.md` §2, §12). That is a gap in rbceq2 to raise upstream (§12); until it is closed, the merge must guarantee one record per event.

Triage, not caller ownership, decides what survives. The SV caller missed a real 3.6 kb Gerbich deletion that the CNV caller found on 3 bins, so neither caller can own a size band (`sv_cnv_overlap_50_samples.md` §11). For an SG without a CNV VCF the rules run over the SV records alone.

#### Triage rules

| Rule | What is kept or done | Why |
|---|---|---|
| 1. Restrict both structural VCFs to db-relevant records | Keep an SV or CNV record only if it (a) overlaps a db SV definition within `SvMatcher`'s own positional and length tolerance (its defaults; use the same code), or (b) is a PASS deletion under `sv_del_max_bp` (1 Mb) spanning a defining SNV/indel site in `bg_site_systems.<genome>.tsv`. Discard everything else. | Removes whole-arm `DUP:TANDEM` artefacts and megabase karyotype events, which §5.6 handles instead. |
| 2. Keep every surviving record from either caller, whatever its FILTER | A `cnvLength` CNV record is dropped if `BC < sv_min_bins`; a surviving non-PASS record has its FILTER rewritten to PASS with the original recorded in `INFO/SVFILTER`. | Corroboration from the SV caller is far more likely at higher bin counts, so a low-bin CNV-only record is dropped rather than trusted; `--no_filter` is never used, since it is global. |
| 3. Where two surviving records describe one event, keep the SV record | Triggered by reciprocal overlap of at least `sv_recip_overlap` in the same direction; the dropped partner's ID is recorded in `INFO/SVPARTNER`. | The SV caller's breakpoints are exact; the CNV caller's are bin-snapped, and several sub-10 kb targets differ from each other only by breakpoint. |
| 4. Tag every kept record with `INFO/SVSRC` | `SV`, `CNV`, or `CNV_LOWRES` for a CNV-only record under 10 kb. | Lets the QC (§6) grade the allele it produces. |
| 5. Assert no two surviving records share CHROM, POS, END and SVTYPE | Fail the sample if they do. | Guarantees the one-record-per-event property the merge exists to provide. |

Rule 3 is the reciprocal-overlap collapse a `bedtools merge` would do (review suggested it). The job implements it with `SvMatcher`'s own tolerance code rather than bedtools, so the overlap test is the one rbceq2 applies, and the kept record can carry the partner's ID and the source tag. Review also noted that a copy-loss and a deletion record could be merged into one synthetic record. The allele is the same either way, and keeping the SV record with `SVPARTNER` preserves both callers' evidence for the QC without inventing coordinates. Rule 3 covers both directions, but in 150 genomes it fired only on deletions; the one DUP token in the database (GYP\*505) is the only case where a duplication pair could matter.

#### What the study showed

- Rule 1(a) kept 0 to 2 records per genome, 3 in total across 150 genomes: the two Gerbich records and the Gerbich CNV-only record.
- Rule 1(b) kept another 0 to 3 records per genome. They include a recurrent 8.5–12.8 kb deletion at chr19:48.69 Mb spanning every FUT2 defining site in 3 of 150 genomes, seen by both callers (5 PASS records plus one `cnvLength`). It matches no db allele, so rbceq2 ignores it, but the QC must not (§6).
- Rule 1 also removed a 126 Mb `MaxDepth` SV record seen once over AUG and RHAG.
- Rule 2's 3-bin floor (`sv_min_bins`) reflects that the SV caller corroborates 45% of 1–2-bin deletions against 75–79% of 3–4-bin ones.
- Rule 3 fired only on the Gerbich case. The SV caller's breakpoints were exact to `CIPOS` (0–50 bp) where the CNV caller's were bin-snapped by 0.2–4.6 kb, and the sub-10 kb targets (seven GE alleles, three A4GALT, the GYP cluster) differ from each other only by breakpoint.
- The only allele actually called across the study was the Gerbich deletion.

The merge must not add a `bcftools +fixploidy`, `+setGT` or any other ploidy rewrite. Since `popgen_rbceq2#20` no step in the pipe rewrites a genotype, and rbceq2 refuses a file that claims two ploidies in one region (§13, the 2026-09-21 entry). The SNV records already carry DRAGEN's true haploid GT on non-PAR chrX, and the study saw haploid GT on a few chrX CNV records too. So a haploid structural record there is consistent, and a diploid one is the contradiction `validate_for_rbceq2` should catch.

### 5.6 Sex and ploidy

This step derives each sample's sex and ploidy from DRAGEN's own estimate, for CNV direction and the karyotype gate (§5.4) and for the `KARYOTYPE` and `SEXCHECK` QC flags (§6). For the X-linked systems (XK/Kx, XG, CD99), expected copy number depends on the sample's sex and on whether the region is PAR; RBCeq2 infers zygosity from GT without knowing either.

#### Sources, in order

1. **`sg.meta['qc']`, written by `single_sample_qc_popgen`.** Its `RegisterQcMetricsToMetamist` stage registers `ploidy_estimation`, `norm_x_coverage` and `norm_y_coverage` per SG. These are DRAGEN's `Ploidy estimation`, `X median / Autosomal median` and `Y median / Autosomal median` from `ploidy_estimation_metrics.csv`, as MultiQC reads them. This is the primary source because it is already in metamist and is the only one that names X0 and XYY. When the meta is absent (a cohort not yet through that pipeline), the stage reads the ploidy file directly; the file was present for all 150 study genomes and every exome checked (§3). If the file is also absent, source 2 applies.
2. **The somalier signals in the same meta:** `f_stat_raw`, `x_het_rate`, `y_calls` and `y_n`. They are computed from genotype sketches, not from the ploidy file, so they survive the uneven-coverage genomes for which DRAGEN omits that file (review, 2026-09-18: three Garvan DSP samples in `tenk10k-phase2`, none from AGRF). `karyotype_from_signals` in `ourdna_genomic_atlas` (`ImputeSex`) already turns them and the DRAGEN estimate into a karyotype, with a loss-of-Y rescue and an `ambiguous` gate. This design reuses that function rather than deriving a second one, and treats `ambiguous` as a non-XX/XY estimate. Where the function lives so both repos can import it is an implementation decision (§9). An SG with none of sources 1 and 2 fails (§5.1).
3. **Reported sex, as a cross-check only.** `single_sample_qc_popgen` already compares DRAGEN's ploidy estimate with the participant's reported sex, and records a failure in `sg.meta['qc']['qc_checks_failed']`. The `SEXCHECK` flag (§6) reads that recorded outcome for the X-linked systems and recomputes nothing. Reported sex never drives the gate, because it cannot express X0 or XYY.

Not used, with the reason for each in §7: the chrX/chrY averages in `wgs_coverage_metrics.csv` (review suggested them as the fallback), the `##referenceSexKaryotype` header, `.ploidy.vcf.gz`, and `SEX GENOTYPER` in `cnv_metrics.csv`. The calibration the coverage ratios would have had is in `sv_cnv_overlap_50_samples.md` §13.

#### Expected copy number by region and karyotype

Expected copy number depends on region and sex together, not on a blanket 'chrX = 1 in males' rule. Of the seven chrX structural alleles, four sit in PAR1 and are diploid in males; only XK and ATP11C are hemizygous.

| Region | XX | XY |
|---|---|---|
| PAR1/PAR2 (CD99, XG) | 2 | 2 |
| Non-PAR chrX (XK, ATP11C) | 2 | 1 |
| chrY | 0 | 1 |

- PAR1: CD99\*01N.01/02 (about 2.71 Mb) and XG\*01N.02/03 (about 2.78 Mb, on the PAR1 boundary).
- Non-PAR chrX: XK (37.7 Mb) and ATP11C (139.7 Mb).

A haploid-X rule would expect CN=1 in PAR1 and miscall the normal CD99/XG state in males as a deletion. So CNV direction (§5.4) reads the table above: diploid in PAR1/PAR2, hemizygous in non-PAR chrX for males, and chrY 1 in males and 0 in females.

#### Karyotype gate

The gate reads the sample's estimate from source 1, or from source 2 when no ploidy estimate exists. For any value other than `XX` or `XY`, `ambiguous` included, it drops or tags the sample's chrX/chrY CNV records (§5.4, open decision Q5). The QC then reports the X-linked systems as `KARYOTYPE:<estimate>` rather than assessing them.

The gate exists because of what the study saw in 150 genomes:

- The 62 XY samples produced no large chrX CNV event, so DRAGEN's caller handles a normal male.
- The one genome estimated X0 carried PASS heterozygous CN=1 deletions of 1.4–10.6 Mb across XK, CD99, XG and ATP11C.
- The one XYY genome carried a PASS CN=3 duplication across PAR1 (CD99, XG).

Both escaped a false allele only because the length gate rejected megabase events against 11–219 kb tokens; a shorter segment would have matched.

### 5.7 Wire into the caller stage

`GenotypeBloodGroupsWithRbceq2` (`stages/blood_group_genotyping/genotype.py`) changes only its input: `--vcf` now points at `<sg>.rbceq2_input.vcf.gz`.

`--no_filter` is not used to admit sub-10 kb CNVs. It is global and would also let through non-PASS SNVs and SVs; instead, the kept sub-10 kb CNV records have their FILTER rewritten to PASS in preprocessing (§5.4, §5.5). `--phased` stays off, since detection never needs it, and `--RH` stays off (§4).

---

## 6. QC for structural calls

This section extends `FlagBloodGroupCallQc` to structural calls with seven new flags. For each defining SNV site the QC already asks whether the caller looked and how well. Structural alleles need the same two answers, from different evidence:

| Question | 10 kb+ targets | Sub-10 kb targets |
|---|---|---|
| Did the caller assess the region? | Yes, from the CNV VCF's `DRAGEN:REF:` records, which tile every 10 kb+ target in 98–100% of genomes with 16–195 bins. A gap, or a `cnvQual` event, is the analogue of `NOCOV`. | Not from the structural VCFs: the SV VCF holds events only, and such a target holds 1–7 CNV bins. Inferred instead from the gVCF reference blocks over the interval, whose `DP`/`MIN_DP` the conversion stage already reads. |
| How good is the call? | `QUAL`, `CN`, segment mean (`SM`), bin count (`BC`) and FILTER. | The SV caller's `QUAL`, `PR`/`SR` and `CIPOS`; for a CNV-only call, `BC`, `SM`, `QUAL`, and the fact that the SV caller saw nothing. Mean gVCF depth over the target against its flanks is the signal the study's read-level check used: about half for a heterozygous deletion, near zero for homozygous, flat for none. |

Design:

1. The site-system map gains interval rows. `bg_db.py` stops dropping `kind == 'sv'`, and each SV definition becomes a `(chrom, start, end, system, allele)` row. The QC job reads structural records from the merged VCF or from rbceq2's debug log, which in 2.4.4 names the source record for each SV match; which is the source of record is open (Q7).
2. New flags, joined with `+` to provenance like today's:
   - `SVNOCOV:<system>(<allele>,gap=<bp>)`: a 10 kb+ target not tiled by CNV records.
   - `SVLOWRES:<system>(<allele>,src=CNV,BC=<n>,SM=<x>,QUAL=<q>)`: an allele called from a §5.4 CNV-only sub-10 kb record. Provisional by construction.
   - `SVDEPTH:<system>(<allele>,ratio=<x>,DP=<n>,flank=<n>)`: gVCF depth over a sub-10 kb target falls below `sv_depth_ratio` of its flanking depth with no kept record to explain it. This is the dosage-drop signal without a breakpoint call.
   - `SVUNASSESSED:<system>`: a sub-10 kb target with no gVCF record over the interval at all, an exome hole or an unmapped region, so neither a call nor its absence can be judged.
   - `KARYOTYPE:<estimate>`: X-linked systems in a non-XX/XY sample (§5.6), `ambiguous` from the somalier fallback included.
   - `SEXCHECK:<reported>/<estimated>`: X-linked systems of a sample for which `single_sample_qc_popgen` recorded a ploidy-versus-reported-sex failure in `sg.meta['qc']['qc_checks_failed']` (§5.6).
   - `SVDEL:<system>(<site>,del=<chrom:pos-end>,src=<caller>,GT=<gt>)`: a kept deletion, matched to a db allele or not, spans a defining SNV/indel site of the system.

   A structural call that passes everything is not listed, so `PASS` keeps meaning 'nothing to report'.

   For an exome the CNV VCF has no reference records (§3), so assessment is read from the capture design instead. A 10 kb+ target with capture targets inside it counts as assessed by the panel-of-normals caller, and one with none is `SVUNASSESSED:<system>`; the design-gate resource is the one the post-hoc exome recall already reads (`popgen_rbceq2#14`). `SVNOCOV` keeps its genome meaning, a gap in the CNV caller's tiling, and is never emitted for an exome. An SG whose run produced no CNV VCF has every 10 kb+ target `SVUNASSESSED`. Sub-10 kb targets follow the gVCF depth path where the interval lies in the capture design, and are `SVUNASSESSED` otherwise.

3. The defining-sites extract grows by the sub-10 kb target intervals. `gen_bg_resources.py` emits them as a second BED from the same db parse, and `FilterAndConvertGvcfsForRbceq2` extracts `DP`/`MIN_DP`/`END` over them and a flank each side, in the pass it already makes for the SNV sites.
4. Accepted for the first cut: gVCF depth does not give the mapping-quality-zero fraction or the discordant-pair evidence a CRAM read would. DRAGEN's DP already excludes reads failing its mapping filters, so poor mappability shows as low depth rather than as its own signal, which a flag can live with. Breakpoint confirmation is the callers' job, not the QC's.

#### Why `SVDEL` exists

rbceq2 does nothing with an unmatched deletion. A structural record enters its variant pool only under the db token it matched, and the adjustment that turns a homozygous call inside a deletion into hemizygous (`modify_variant_pool_if_large_indel`) sees only pool entries. So a gVCF `A/A` under a heterozygous deletion that matches no db SV is reported homozygous, the same silent-wrong-call class as absent-means-reference. The QC has to say so, as it already does with `DEL` for a small deletion that removed the base.

---

## 7. Alternatives rejected

| Alternative | Why rejected | Where the evidence is |
|---|---|---|
| Size-band caller ownership: the SV caller owns records under 10 kb, the CNV caller owns everything above | The CNV-only Gerbich deletion shows the SV caller missed a real sub-10 kb event that the CNV caller found on 3 bins | `sv_cnv_overlap_50_samples.md` §11 |
| Merge adjacent same-direction records into one event before triage (review suggestion, 2026-09-18) | Of 254 adjacent CNV pairs within 2.5 kb in 150 genomes, only 3 shared a `CN`, and no merged span matched a db allele that neither part matched; the SV VCF had 1 adjacent pair in 17,175 records. The adjacencies are copy-number steps, not one event in pieces. Revisit if a carrier ever shows a target split across records, since each fragment would fail rbceq2's length gate alone | `sv_cnv_overlap_50_samples.md` §13 |
| Keep exome SGs out of scope | Every production exome has both structural VCFs, with records in the blood-group regions | §3, `resolvability_by_input_class.md` §3 |
| Decide whether a CNV VCF is expected from the sequencing type, or from a config flag | The two exome runs differ, so the type is the wrong key; a config flag would be a default standing in for a fact `cnv_metrics.csv` already records | §5.1 |
| Metamist reported sex as the karyotype source | Cannot express X0 or XYY, the two estimates the gate exists for; kept as a cross-check | §5.6 |
| `SEX GENOTYPER` in `cnv_metrics.csv` as the sex fallback | Agreed with the ploidy estimate in 49 of 50 genomes but was blank for the X0 genome, so it fails on exactly the sample that needs it | `sv_cnv_overlap_50_samples.md` §13 |
| `wgs_coverage_metrics.csv` chrX/chrY averages as the sex fallback (review suggestion) | Registered nowhere (`single_sample_qc_popgen` takes only mean, median and percent-over-20x coverage from the file). A classifier on them would duplicate `karyotype_from_signals` with a weaker signal, and capture design shifts the ratios: an XY exome read Y 0.35, against 0.42–0.52 for XY genomes | §5.6, `sv_cnv_overlap_50_samples.md` §13 |
| Recomputing the reported-sex check for `SEXCHECK` | `single_sample_qc_popgen` already records the ploidy-versus-reported-sex outcome in `qc_checks_failed` | §5.6, §6 |
| New metamist analysis types for the raw SV and CNV VCFs | `dragen_align` cannot register them before its Nextflow refactor, cpg_flow is not to be extended, and the paths derive from the gVCF's | §5.1 |
| `--no_filter` to admit sub-10 kb CNVs | It is global, so it would also admit non-PASS SNVs and SVs | §5.4, §5.7 |
| The `##referenceSexKaryotype` header as the sex source | It is a constant `XXYY` in 102 of 102 samples, not a per-sample estimate | §5.6 |
| `.ploidy.vcf.gz` as the sex source | Absent for exomes, and its `##estimatedSexKaryotype` header repeats the ploidy estimate already in the QC meta | §3, §5.6 |
| Deferring dedup to `select_best_per_vcf` | It keeps one db definition per record, not per locus, so two records for one deletion yield two alleles | §5.5, `sv_cnv_overlap_50_samples.md` §4 |
| The SV VCF alone | The study's `sv_only` policy found the same single allele, but the CNV-only Gerbich deletion is a real event the SV caller missed | `sv_cnv_overlap_50_samples.md` §3, §11 |
| Relocating `+fixploidy` onto the merged VCF | Would fabricate homozygotes for every structural and SNV record on non-PAR chrX; rbceq2 2.4.4 reads the haploid GT natively | §13, the 2026-09-21 entry |
| Checked-in extracts of a cohort genome's DRAGEN outputs as unit fixtures | No individual-level data in the repo; the study's observations are reproduced as constructed records instead | §10 |

---

## 8. Thresholds

| Parameter | Value | What it gates | Where it came from |
|---|---|---|---|
| `sv_min_bins` | 3 | Minimum CNV bin count to keep a sub-10 kb `cnvLength` record (§5.4, §5.5) | `sv_cnv_overlap_50_samples.md` §10–11: 3–4-bin records have SV-caller support 75–79% of the time, 1–2-bin only 45% |
| `sv_recip_overlap` | 0.5 | Reciprocal overlap required to treat two records as the same event (§5.5, rule 3) | `sv_cnv_overlap_50_samples.md` §2 |
| `sv_lowres_max_bp` | 10000 | Size below which a CNV-only record is graded `CNV_LOWRES` and flagged `SVLOWRES` | `sv_cnv_overlap_50_samples.md` §1, §11 |
| `sv_depth_ratio` | 0.7 | Depth ratio below which gVCF depth over a sub-10 kb target triggers `SVDEPTH` | §6 |
| `sv_depth_flank_bp` | 5000 | Flank size either side of a target used to compute the depth ratio | §6 |
| `sv_del_max_bp` (name proposed) | 1000000 | Upper size cap on a PASS deletion admitted by rule 1(b) for spanning a defining SNV/indel site | §5.5, rule 1(b); caps out the X0 sample's megabase events and the 126 Mb SV `MaxDepth` record |

Each threshold is recorded in the Analysis meta of the stage that applies it, as `min_depth`/`min_gq` are today; §9 says which config section holds which. The sex fallback (§5.6) carries no thresholds of its own: `karyotype_from_signals` owns them.

---

## 9. Touch points

| File / symbol | Change |
|---|---|
| `src/popgen_rbceq2/stages/blood_group_genotyping/filter_and_convert.py` (`FilterAndConvertGvcfsForRbceq2`) | Unchanged. Its `vcf` output is the SNV branch (§5.2); the merge stage lists it as a required stage |
| `src/popgen_rbceq2/stages/blood_group_genotyping/` | New `PreprocessDragenForRbceq2` per-SG stage for every SG, with the CNV branch conditional on `cnv_metrics.csv` (§5.1); repoint `GenotypeBloodGroupsWithRbceq2`'s required stages and `--vcf` (`genotype.py`) at its output |
| `src/popgen_rbceq2/jobs/` | New job module. Proposed function names, no code yet: `resolve_dragen_inputs`, `read_sample_sex`, `fix_dragen_cnv_vcf`, `normalize_sv`, `merge_variant_vcfs`, `validate_for_rbceq2` |
| SV/CNV input resolution | Paths derived from `sequencing_group.gvcf`'s DRAGEN output prefix (§5.1), CNV expectation read from `cnv_metrics.csv`, paths recorded in the merged VCF's Analysis meta (Q6 decided) |
| Sex and ploidy | Read from `sequencing_group.meta['qc']` as `single_sample_qc_popgen` writes it: `ploidy_estimation`, `norm_x_coverage`, `norm_y_coverage`, `f_stat_raw`, `x_het_rate`, `y_calls`, `y_n`, `qc_checks_failed`; the ploidy file only when the meta is absent (§5.6) |
| `karyotype_from_signals` (today in `ourdna_genomic_atlas`, `jobs/sample_qc/impute_sex_job.py`) | Reused for the somalier fallback. Import across repos or move to a shared package: decide at implementation, do not re-derive |
| `src/popgen_rbceq2/config/popgen_rbceq2_default_config.toml`, `[workflow.preprocess_dragen_for_rbceq2]` (new) | `cnv_svtype_from`, `min_size`, resources, and the merge thresholds from §8 (`sv_min_bins`, `sv_recip_overlap`, `sv_lowres_max_bp`, `sv_del_max_bp`) |
| `src/popgen_rbceq2/resources/bg_regions.GRCh38.bed` | Reused unchanged for region-restrict |
| `src/popgen_rbceq2/resources/bg_site_systems.GRCh38.tsv` | Gains SV definition rows (§6) |
| `src/popgen_rbceq2/scripts/bg_db.py` (`SiteKind`) | Stops excluding `kind == 'sv'` rows when building the site-system map |
| `src/popgen_rbceq2/scripts/gen_bg_resources.py` | Emits the sub-10 kb target BED (§6) from the same db parse |
| `src/popgen_rbceq2/stages/blood_group_qc/call_qc.py` (`FlagBloodGroupCallQc`), `src/popgen_rbceq2/jobs/rbceq2_call_qc_job.py` | Today raises `ValueError` if the site-system map carries any `kind == 'sv'` row; gains the §6 QC design once that map has SV rows |
| `[workflow.flag_blood_group_call_qc]` in the same config file | Holds `min_depth = 10`, `min_gq = 20` today; gains the QC thresholds from §8: `sv_depth_ratio`, `sv_depth_flank_bp` |
| `src/popgen_rbceq2/constants.py` (`RBCEQ2_VERSION`, `RBCEQ2_IMAGE_TAG`, `NON_PAR_X`) | Pinned at `2.4.4` / `2.4.4-1` (`popgen_rbceq2#17`); `+fixploidy` removed in `#20`, which also added the single-copy chrX bounds the post-hoc gate reads |

---

## 10. Testing plan, expected outcome and success criteria

Constructed-record tests run in CI, and concordance checks run by hand against real files in GCS. Two further tests need samples the team does not yet have.

#### What runs in CI

Unit tests use constructed records committed under `tests/`, shaped like the records in §3 and the study's observations; no extract of any cohort genome enters the repo. CI cannot read `gs://cpg-ourdna-main/…`, so no unit test may resolve a live GCS path.

| Test | Input (constructed) | Asserts |
|---|---|---|
| `fix_dragen_cnv_vcf` | A CNV VCF holding a `DRAGEN:REF:` block and `<DEL>`/`<DUP>` events | REF records dropped; every surviving `SVTYPE` in {DEL,DUP}; counts match `CnvRewriteStats` |
| `merge_variant_vcfs` | An SNV VCF, SV VCF and CNV VCF, with a chrX record for the ploidy path, each with its own FORMAT string as in §3 | Output is sorted, single-sample, tabix-indexed, `chr`-prefixed; header holds the union of the three; every row has `GT` first |
| `validate_for_rbceq2` | A merged VCF with one row whose FORMAT is `DP:GT` | Fails naming the row |
| No-CNV path | An SNV VCF and SV VCF, a `cnv_metrics.csv` that stops after the target-interval count, no CNV VCF | Merge succeeds; every 10 kb+ target is `SVUNASSESSED`; the same inputs plus a CNV VCF, or segment metrics without a CNV VCF, fail the SG |
| Exome CNV path | An exome-shaped CNV VCF: events only, no `DRAGEN:REF:` records, PASS and `cnvQual`, `BC` in target counts, with `Number of segments` in the metrics | Edit 1 drops nothing; triage keeps the db-relevant records; the QC reads assessability from the capture design and emits no `SVNOCOV` |
| Sex sources | (a) meta with `ploidy_estimation`; (b) meta without it plus a ploidy file; (c) neither, with somalier signals shaped like XX, XY and loss-of-Y; (d) none of the three; (e) meta whose `qc_checks_failed` records the ploidy-versus-reported-sex check | (a) estimate read from meta, file untouched; (b) read from the file; (c) `XX`, `XY`, `ambiguous` via `karyotype_from_signals`, the last gated; (d) the SG fails; (e) `SEXCHECK` on the X-linked systems only |
| Sub-10 kb path (a) | An SV DEL record of exactly 3609 bp at the GE\*01.-02.01 coordinates plus a 4961 bp `cnvLength` CN=1 CNV record over it | The merge keeps the SV record and records the CNV partner; rbceq2 calls one GE allele, never two |
| Sub-10 kb path (b) | The CNV record alone, 3 bins | Kept, PASS-rewritten, tagged; the QC flags the GE call `SVLOWRES` |
| Sub-10 kb path (c) | The CNV record alone, 2 bins | Dropped; GE unassessed |
| Karyotype gate | Megabase CN=1 deletions across XK, CD99/XG and ATP11C with a ploidy metrics file reading `X0`; a PAR1 CN=3 duplication with `XYY` | No X-linked allele is called; the QC reports `KARYOTYPE:<estimate>` |

The study observed sub-10 kb cases (a) and (b) in real genomes. The fixtures are constructed records shaped like them, not extracts.

#### What runs by hand against GCS

The OurDNA 1000 Genomes (1KG) control replicates are independent sequencings of one public control individual (`NA12878`, XX): five SGs when this spec was first drafted, nine by the time of `ourdna_genomic_atlas#136`.

- Replicate concordance: run the full preprocess-plus-rbceq2 on every replicate. Structural calls must be identical across all of them, since they are one individual. Any divergence is a bug or a QC signal (compare the SNV divergence table in `ourdna_genomic_atlas#124`).
- 150-genome concordance: the study's scripts, once committed, check that a code change to the merge does not alter the per-band record counts or the single GE hit under the `pass` policy.

#### Still wanted

Both need a sample the replicates cannot supply: being XX, they exercise neither the male chrX/chrY paths nor a positive structural hit.

- A true positive: a carrier of a known structural allele from the database, asserting the allele is called. Concrete targets: XK\*N.05 (8 kb, non-PAR chrX), which must come from the SV VCF and exercises the sub-10 kb path; and XK\*N.01 (53 kb) or ATP11C\*01N.01 (219 kb), larger non-PAR deletions from the CNV VCF that exercise the hemizygous ploidy path.
- PAR-versus-hemizygous ploidy, on a male genome: confirm the region-aware expected CN (§5.6), so that CD99/XG (PAR1, diploid in males) are not miscalled as deletions while XK/ATP11C (non-PAR) are handled as hemizygous.

#### Expected outcome

RBCeq2 gains the ability to call the 66 callable single-event structural alleles in 19 non-RH systems, and some of the 16 mid-size indels, that it cannot call today (`resolvability_by_input_class.md`). It does so from data DRAGEN already produces, with no change to RBCeq2 and one new preprocessing stage. Hybrid RH/GYP alleles remain out of reach on short-read data.

These alleles are rare: across the 150 study genomes the design calls one, the Gerbich deletion. The measurable gain is therefore correctness on the carriers that do occur, plus QC that reports unassessed targets and the deletions over defining sites that RBCeq2 ignores.

#### Success criteria

- The constructed-record tests above pass.
- The 150-genome per-band counts and the single GE hit reproduce under the `pass` policy.
- No X-linked allele is called in the X0 or XYY fixtures.
- The merged VCF passes `validate_for_rbceq2`.
- An SG without a CNV VCF runs end to end on the gVCF and SV VCF alone, and its QC names every 10 kb+ target `SVUNASSESSED`.
- An exome with a CNV VCF runs through §5.4, and its QC emits no `SVNOCOV`.

---

## 11. Risks

- Short-read paralog loci. Beyond the excluded RH/GYP hybrids, any blood-group gene with a close paralog risks mismapped or absent CNV calls. Concordance across the 1KG replicates is the necessary check.
- DRAGEN version. Our calls are DRAGEN 3.7.8. Review (Slack, 2026-09-18) relays Genomics England's internal view that DRAGEN 4.x SV/CNV calls outperform Canvas and Manta and are worth using, while 3.x calls are not. The maintainer's evidence that DRAGEN SV/CNV detects blood-group deletions may have come from 4.x, and the author is asking which version (Q8). If it was 4.x, this design stands for 3.7.8 on the study's own concordance and read-level checks alone, and a true-positive carrier (§10) matters more.
- DRAGEN's ploidy estimator may be off in production. The v5 comparison (`docs/release_v5_comparison.md` §1) found production runs made with `--enable-ploidy-estimator false`, so their chrX copy number follows the sex the ICA pipeline was given. If `Ploidy estimation` in those runs echoes that input rather than measuring it, source 1 of §5.6 can never read X0 or XYY, and the karyotype gate passes the chromosome-scale chrX events it exists to stop. The 150-genome study read OurDNA runs, where the estimate did name X0 and XYY. Check this before implementing §5.6: in a production run with the estimator off, does `ploidy_estimation_metrics.csv` exist, and does it disagree with the supplied sex for any sample? If it only echoes the input, the somalier signals (source 2) become the primary source.
- Sub-10 kb CNV-only records. A kept CNV-only sub-10 kb record is weaker evidence than an SV-confirmed one, so the QC grades it `SVLOWRES`, and a 1–2-bin record is dropped.
- Half-called duplication GT. Every PASS CNV duplication carries GT `./1` (§5.5). rbceq2 accepts it, but its scoring of a half-called GT on a matched DUP is unverified; only GYP\*505 is exposed.
- Exome sensitivity and thresholds. On capture data the SV caller sees breakpoints only where reads reach them, so a 3.6 kb exon deletion with intronic breakpoints outside the capture may be missed even though the exon is covered. The exome CNV caller is also a different instrument from the genome one (§3). `sv_min_bins` was set from 1–2 kb genome bins, while exome `BC` counts capture targets. And the panel-of-normals caller's false-positive rate in the blood-group regions is unmeasured: one exome had 30 event records in the regions, against 1 to 7 PASS events per genome. Until an exome study like `sv_cnv_overlap_50_samples.md` is run, the exome path makes structural alleles possible, not promised, and `SVLOWRES` and `SVDEPTH` are the flags that qualify an exome call.
- Positive controls. The 1KG replicates showed no blood-group CNVs at GYP, XK or RHD in the earlier data look, which suits a concordance baseline but tests no positive. The study's Gerbich carrier is a real positive for one sub-10 kb allele; a 10 kb or larger true positive is still wanted to exercise the CNV VCF path (§10).
- Merge correctness. A contig or sample-name mismatch, or an unsorted concat, silently breaks RBCeq2's region fetch, which is why `validate_for_rbceq2` asserts sorted, single-sample, `chr`-prefixed and indexed.
- A ploidy rewrite reintroduced on the merged VCF. Since `popgen_rbceq2#20` the file rbceq2 reads carries DRAGEN's true haploid GT and a hemizygous call renders as `XK*01.02/-`. A `+fixploidy` or `+setGT` on the concat would turn it back into a fabricated homozygote, and nothing fails or logs: the VCF stays well formed and the TSV stays plausible (§13, the 2026-09-21 entry).

---

## 12. Upstream requests

Feature requests to raise with the RBCeq2 maintainers. This design works around the first and reports the other two; none blocks it.

1. Matching is per record, not per locus. Two records for one event yield two alleles (§5.5 figure). Our merge guarantees one record per event; a per-locus selection in `select_best_per_vcf`, or a warning when two records match overlapping tokens of one system, would make that unnecessary.
2. Unmatched deletions are ignored. A deletion that matches no db token never enters the variant pool, so the hemizygosity adjustment for SNVs inside it never runs and a single-copy genotype is read as homozygous. Our QC reports it (`SVDEL`, §6); rbceq2 could apply the same adjustment it already applies to matched deletions.
3. Undefined null alleles. A deletion removing coding sequence of a blood-group gene is a null whether or not the db names it. The study found a recurrent 8.5–12.8 kb whole-FUT2 deletion over every FUT2 defining site in 3 of 150 genomes, matching no token. Two asks: consider it as a candidate allele for the db, and consider a generic 'coding deletion in system X, unnamed' call. We do not build gene-model resources for this here; it stays an observation to bear in mind, and a reason the `SVDEL` flag names the system.

---

## 13. Decision log

Newest first; the Q3 entry is undated. Each entry gives the decision, the alternatives rejected and the consequences; a later entry that reverses one says so.

**2026-09-21, `+fixploidy` removed at the pin bump (`popgen_rbceq2#20`), closing the 2026-09-18 deferral.** Decision: the SNV conversion pipe no longer rewrites any genotype. rbceq2 2.4.4 reads DRAGEN's one-token non-PAR chrX GT natively and renders a hemizygous null as `XK*01.02/-`. The `1` to `1|1` expansion had rendered it `XK*01.02/XK*01.02`, indistinguishable from a female homozygote; phenotype TSVs were unchanged in the synthetic check. rbceq2 derives one chromosome-copy count per blood group and refuses a file whose records disagree (`Undetermined`), while a record claiming fewer copies than its neighbours passes silently and can flip the phenotype. So the exome post-hoc fill is kept out of single-copy chrX (`constants.NON_PAR_X`), at the cost of 10 off-design XK sites on Twist and 2 on Agilent CREv2, which reach the QC as `NOCOV`. Alternatives rejected: relocating `+fixploidy` onto the merged VCF (the superseded 2026-09-18 position), which would re-fabricate homozygotes for every structural and SNV record on non-PAR chrX; rewriting the post-hoc caller's genotypes to match DRAGEN, which needs a ploidy inference and turns a wrong inference into a confident wrong call. Consequences for this design: the concat in §5.5 adds no ploidy rewrite; `validate_for_rbceq2` should treat a two-token GT on a non-PAR chrX structural record beside one-token SNV records as the contradiction rbceq2 will refuse; §5.2, §5.5, §9 and §11 updated. Open: the gate is sex-blind, so female exomes also lose those XK sites; a gate based on DRAGEN's encoding was raised in the `#20` review.

**2026-09-18, first review round (`popgen_rbceq2#15`, `#16`, and two Slack threads).** Six decisions from Alexander Stuckey's comments:

1. *Naming.* The SV VCF's caller is the DRAGEN SV caller, which integrates and extends Manta. Every doc in this directory now says "SV caller" or "SV VCF", keeping literal record IDs, and the `INFO/SVSRC` value `MANTA` becomes `SV`. Rejected: keeping "Manta" as shorthand, which lets claims about Manta's internals stand in for claims about DRAGEN's caller.
2. *Exomes in scope.* Every production exome has an SV VCF and a panel-of-normals CNV VCF (11,917 of 11,917 in `cpg-mackenzie-main`); the 10 exomes of the test run have the SV VCF only (§3). The merge runs for every SG, reads whether a CNV VCF is expected from `cnv_metrics.csv` (§5.1), handles the exome file's lack of reference records and of `cnvLength` (§5.4), and reads exome assessability from the capture design (§6). Rejected: keeping exomes out and correcting only the wording, which leaves structural alleles uncalled on data that carries them; keying the CNV expectation on sequencing type, which the two exome runs contradict.
3. *Q6 decided: no new metamist analysis types.* cpg_flow is not extended, and `dragen_align` cannot register the files before its Nextflow refactor. So the stage derives the paths from the registered gVCF's prefix, fails the SG on a missing expected file, and records the paths in the merged VCF's Analysis meta (§5.1). Rejected: a small registration stage in this repo, which would register another pipeline's outputs under types nobody else reads.
4. *FORMAT needs no harmonisation.* Verified on the study files and a real concat (§5.5); `validate_for_rbceq2` gains a GT-first assertion.
5. *No merge-adjacent step.* Measured over 150 genomes, adjacent CNV records are copy-number steps, not fragments, and the 3-bin Gerbich call is one record (§7; `sv_cnv_overlap_50_samples.md` §13). Rejected: an optional pre-triage merge with a gap threshold, which would have merged records of different `CN` in 251 of 254 observed cases.
6. *Sex from the QC meta already in metamist.* Sex is read from what `single_sample_qc_popgen` registers in `sg.meta['qc']`, with the ploidy file and then the somalier signals as fallbacks, and `SEXCHECK` echoing the recorded reported-sex check (§5.6). An earlier draft of this round proposed a coverage-ratio fallback from `wgs_coverage_metrics.csv` with four thresholds; the author pointed out the QC pipeline's somalier signals, and it was withdrawn the same day. The rejected sources and why are in §7.

Also recorded in this round: the `+fixploidy` question review raised at §5.2 was the one already deferred to the pin bump (the 2026-09-18 deferral entry below, since closed by the 2026-09-21 entry). The bedtools-merge suggestion is rule 3 (§5.5). Two reviewer observations stand as written: the callers' breakpoints never agree, and the 2.4.4 tie error will not fire on real pairs. The DRAGEN-version question is a new risk (§11) and open question (Q8).

**2026-09-18, Q2 reversed to triage.** Decision: deduplicate the merged structural records by triage, not by caller ownership. Both VCFs are restricted to db-relevant records, survivors are kept, and the SV record is preferred on overlap. Rejected: giving the SV caller the sub-10 kb band and the CNV caller everything else; deferring entirely to `select_best_per_vcf`. Consequences: §5.5 rules 1–5; the per-record, not per-locus, matching gap raised upstream (§12).

**2026-09-18, Q1 resolved.** Decision: sub-10 kb CNV records are triaged like any other structural record, not dropped as a class. They are kept if database-relevant or spanning a defining site and at least 3 bins, with FILTER rewritten per record and the original kept in `INFO/SVFILTER`. Rejected: dropping all sub-10 kb CNV records and sourcing small events from the SV VCF only; using `--no_filter` to admit them. Consequences: the 3-bin floor (`sv_min_bins`, §8) and the per-record FILTER rewrite in §5.4; the QC grades a kept CNV-only sub-10 kb allele `SVLOWRES`/`CNV_LOWRES`.

**2026-09-18, karyotype gate added.** Decision: chrX/chrY CNV records from a sample whose DRAGEN ploidy estimate is not XX or XY are dropped or tagged `SVSRC=CNV_KARYOTYPE` (§5.4, §5.6; the choice is still open, Q5). Rejected: passing all CNV records through untagged and relying on the length gate alone. Consequences: the X0 and XYY synthetic fixtures in §10; the `KARYOTYPE` QC flag in §6.

**2026-09-18, QC design added.** Decision: add the §6 QC design, six flags for structural calls (a seventh, `SEXCHECK`, came with the review round above), reading structural records from the merged VCF or rbceq2's debug log (source still open, Q7), plus a second BED for sub-10 kb target depth. Rejected: none tested; this is a first design. Consequences: `bg_db.py` stops excluding `kind == 'sv'` rows; `gen_bg_resources.py` gains a second BED; thresholds recorded in §8.

**2026-09-18, `+fixploidy` deferred to the pin bump. Superseded by the 2026-09-21 entry.** Decision: keep the single `+fixploidy` from `ourdna_genomic_atlas#128`, add no second one when the merge lands, and decide at the 2.4.4 pin bump, for SNV and structural records together, whether it moves onto the merged VCF or goes. 2.4.4 was described as "focused on supporting haploid encoding in VCF", so the crash it worked around might be gone. This entry superseded an earlier position, moving the invocation onto the merged VCF once the merge existed.

**2026-09-17, §5.2 restated against the conversion stage as it is.** Decision: the SNV branch is the conversion stage's `vcf` output as it exists, described rather than changed. The July 2026 passage it replaced had been overtaken four times in `ourdna_genomic_atlas`: `#124` showed the then DP 20 / GQ 30 hard filter manufacturing false wild-type calls; `#128` appended `+fixploidy`; `#136` removed the filter, taking the replicate cohort's discordant systems from eight to zero; `#138` and `#140` made the stage a single region-restricted pass with a defining-sites extract and added `FlagBloodGroupCallQc`. The stages then moved here (`popgen_rbceq2#5`) and gained the exome post-hoc merge (`popgen_rbceq2#14`). Consequences: §5.2 rewritten; exomes declared out of scope, since reversed by the 2026-09-18 review round; every PR reference repo-qualified.

**2026-09-17, factual refresh (`popgen_rbceq2#15`), and this branch rebased onto it.** Decision: refresh the database figures for 2.4.4 (§2.2, §2.3); add `resolvability_by_input_class.md`; record the maintainer's advice against `--RH` and for a single combined VCF (§1, §4); rewrite the touch points for this repo's paths and the QC stage that did not exist in July; note that 2.4.4's native haploid support makes `+fixploidy` a re-test item. `popgen_rbceq2#16` was first opened from `main` in parallel, and its restructure of this file dropped those changes; it now sits on top of `#15` with the refresh folded back in. Rejected: merging both branches to `main` independently, which would have conflicted in this file. Consequences: one SPEC, with the figures from `#15` and the design from `#16`.

**2026-08, moved from `ourdna_genomic_atlas`.** The RBCeq2 stages moved into this repo (`popgen_rbceq2`); this design was never implemented in either repo. The class names in this spec are unchanged, but the old file paths are not: `src/ourdna_genomic_atlas/stages.py` now maps to `src/popgen_rbceq2/stages/blood_group_genotyping/` (§9).

**2026-07-19, Q4 resolved in `ourdna_genomic_atlas#128`, at rbceq2 2.4.2. Superseded by the 2026-09-21 entry.** DRAGEN writes haploid `GT` (`1`) on non-PAR chrX/chrY for male samples, and rbceq2 2.4.2 crashed on it: `get_ref` asserts `len(GT) == 3` (`core_logic/data_procesing.py:878`) for every defining variant, and `IO/vcf.py` assumed a diploid separator in three more places (`:120`, `:261`, `:299`). Fix: a final parameter-free `bcftools +fixploidy`, whose unprefixed built-in ploidy table never matches `chr`-prefixed hg38, so every haploid GT became diploid (`1` to `1|1`) while PAR, already diploid from DRAGEN, stayed untouched; verified empirically. Rejected: a sex-aware config setting male non-PAR X to ploidy 1, which left the call haploid and re-crashed rbceq2. Caveat recorded then: `1|1` reads as HOM (`core_logic/alleles.py:284`), overstating a hemizygous call. The XX 1KG fixtures missed the crash (§10).

**Q3 endorsed in review, resolved.** Decision: CNV direction comes from the ALT symbol on autosomes, and from `CN` against region-by-sex ploidy on chrX/chrY, PAR-aware: CD99 and XG sit in PAR1 and are diploid in males; only XK and ATP11C are hemizygous. Rejected: trusting the ALT symbol everywhere, including the sex chromosomes. Consequences: §5.4 and §5.6.

---

## 14. Open questions for reviewers

- **Q5, the karyotype gate (§5.4, §5.6).** Drop chrX/chrY CNV records for non-XX/XY samples, or pass them through tagged `SVSRC=CNV_KARYOTYPE`? No reason to prefer either has been recorded; the QC reports `KARYOTYPE:<estimate>` for the sample under both options, so the choice only affects what the merged VCF carries.
- **Q7, the QC's source for structural records.** The merged VCF's own records, or rbceq2's debug log, which in 2.4.4 names the source record for each SV match? Recommendation, not a decision: the merged VCF's records are named first and read as the primary source, with the debug log as a secondary cross-check (§6).
- **Q8, which DRAGEN version the maintainer's evidence came from.** Genomics England's experience (relayed in review) is that DRAGEN 4.x SV/CNV is worth using and 3.x is not; we run 3.7.8 (§11). The author has asked the maintainer. If the answer is 4.x, the design's evidence for 3.7.8 is the 150-genome study and its read-level check alone, and the true-positive carrier in §10 moves up the list.

Q6 (metamist registration of the SV and CNV VCFs) was decided in the review round of 2026-09-18 (§13): no new analysis types; paths derive from the registered gVCF and are recorded in the merged VCF's Analysis meta (§5.1).
