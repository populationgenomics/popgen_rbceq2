# Design spec: structural-variant support for RBCeq2 blood-group calling

**Status:** draft, awaiting design approval; do not implement.
**Revision:** 18 September 2026, after the first review round on `popgen_rbceq2#15` and `#16`. **Author:** Joshua Schmidt. **Reviewers:** Alexander Stuckey.
**Reader:** everyone!
**Decision asked of reviewers:** approve the triage rule and karyotype gate in §5 and the thresholds in §8, or say which to change and why.
**Area:** rbceq2 blood-group genotyping pipeline (`FilterAndConvertGvcfsForRbceq2` to `GenotypeBloodGroupsWithRbceq2` to `CombineRbceq2OutputsPerCohort`), in `popgen_rbceq2`. rbceq2 is pinned at 2.4.3 (`constants.py`); database figures below are from the 2.4.4 release of 2026-09-16, which the pin bump (a separate PR) will adopt.
**Pull-request references:** the stages were ported from `ourdna_genomic_atlas` in August 2026 and both repos number PRs from 1, so every PR below is written `ourdna_genomic_atlas#n` or `popgen_rbceq2#n`.
**Callers:** the SV VCF is written by the DRAGEN 3.7.8 SV caller, which Illumina describes as integrating and extending Manta; its record IDs keep Manta's prefix (`MantaDEL:`). This spec says "SV caller" and "SV VCF", not "Manta", so that nothing here rests on Manta's internals. The CNV VCF is written by DRAGEN's bin-based CNV caller.
**Companion docs (same directory):**
- [`dragen_three_source_merge.md`](./research/dragen_three_source_merge.md): the visual overview and mermaid diagram, read this first for the shape.
- [`implement_cnv_rbceq2_research.md`](./research/implement_cnv_rbceq2_research.md): the underlying RBCeq2 source analysis this spec rests on, read against v2.4.2; where 2.4.4 moved a fact, this spec says so.
- [`resolvability_by_input_class.md`](./research/resolvability_by_input_class.md): per-system estimate of what the gVCF resolves today and what the merge adds, from the 2.4.4 database.
- [`sv_cnv_overlap_50_samples.md`](./research/sv_cnv_overlap_50_samples.md): the 150-genome SV-caller-versus-CNV-caller study the triage rule and karyotype gate rest on.

---

## 1. Summary

RBCeq2's 2.4.4 database defines 67 non-RH alleles across 20 systems by a single large structural event, and today our pipeline does not call any of them, because it feeds RBCeq2 only the SNV gVCF (`research/resolvability_by_input_class.md`). DRAGEN already produces the two structural files needed, an SV VCF from its Manta-derived SV caller and a bin-based CNV VCF, but RBCeq2 reads a single VCF.

The design is therefore a per-sequencing-group (SG) preprocessing stage that merges the three files into one. The maintainer's advice with the 2.4.4 release also suggests to leave `--RH` off because DRAGEN SV/CNV does not reliably detect the RH hybrids.

The merge does three things beyond concatenation. It rewrites DRAGEN's `SVTYPE=CNV` to `DEL` or `DUP`, without which no large deletion matches a database allele. It triages both structural files to records that could match a database allele or that span a defining SNV site, and keeps one record per event, because RBCeq2 matches per record and a deletion seen by both callers would otherwise yield two alleles. It drops or tags chrX and chrY CNV records for samples whose DRAGEN karyotype estimate is not XX or XY.

A new QC design (§6) reports whether the callers assessed each structural target and how well, and reports deletions over defining SNV sites, which RBCeq2 itself ignores.

The evidence is a 150-genome OurDNA study (`research/sv_cnv_overlap_50_samples.md`). Across those genomes the design calls one structural allele, a Gerbich deletion. The value is correctness on rare alleles and honest QC, not call volume; §10 states the expected outcome and success criteria in those terms.

RBCeq2 is not modified. RH and GYP hybrid alleles stay out of scope on short-read data. Exome sequencing groups are in scope. The production exomes carry an SV VCF and a panel-of-normals CNV VCF (§3), so the merge treats them as it treats genomes, with the exome CNV file's differences (no reference records, no `cnvLength` filter, `BC` counting capture targets) handled in §5.4 and §6. Where a run produced no CNV VCF, as the 10-exome test run did, the merge is gVCF plus SV VCF and every 10 kb+ target is reported unassessed.

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
| - | Ploidy (sex, primary) | `dragen_metrics/<sg>/<sg>.ploidy_estimation_metrics.csv` | DRAGEN | Read only |
| - | Coverage (sex, fallback) | `dragen_metrics/<sg>/<sg>.wgs_coverage_metrics.csv` | DRAGEN | Read only |

What each sequencing type has, checked in the buckets on 2026-09-18 (150 OurDNA genomes; the 11,917 mackenzie exomes in `cpg-mackenzie-main`; the 10 exomes of the DRAGEN 3.7.8 test run in `cpg-mackenzie-test`):

| File | Genome | Exome, production (`-main`) | Exome, test run (`-test`) |
|---|---|---|---|
| SNV gVCF | yes | yes | yes |
| SV VCF | yes | yes, 11,917 of 11,917 (one exome: 65 records, 4 in the blood-group regions) | yes, 10 of 10 (one exome: 103 records, 9 in the regions) |
| CNV VCF | yes | yes, 11,917 of 11,917 | no, 0 of 10: `cnv_metrics.csv` stops after "Number of target intervals", so the caller counted reads but never segmented |
| `ploidy_estimation_metrics.csv` | yes | yes | yes |
| `wgs_coverage_metrics.csv` | yes | yes | yes |
| `ploidy.vcf.gz` | yes | no | no |

The exome CNV VCF is not the genome file with fewer records. From one production exome's header and records: it is called per capture target (`--cnv-target-bed`) against a panel of 100 normals (`--cnv-normals-list`, `Number of normal samples,100` in `cnv_metrics.csv`, self-normalisation off) with HSLM segmentation. Every record is an event: no `DRAGEN:REF:` records at all (1,131 records, 30 of them in the blood-group regions, all events). FILTER values are PASS, `cnvQual`, `cnvBinSupportRatio` and `cnvCopyRatio`; there is no `cnvLength`, so sub-10 kb events pass (129 PASS under 1 kb, 121 PASS of 1 to 10 kb in that exome). `BC` counts capture targets, not 1–2 kb bins. FORMAT is the same `GT:SM:CN:BC:PE`, and duplications carry `./1` as in genomes. §5.4 and §6 say what follows from each difference.

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

- New per-SG inputs: the SV and CNV VCFs and the two sex-metric files, derived from the registered gVCF's DRAGEN output prefix and recorded in the merged VCF's Analysis meta. No new metamist analysis types and no cpg_flow extension (§5.1; Q6 decided).
- A preprocessing step that (a) produces an SNV sites VCF, (b) normalises the SV VCF, (c) fixes the CNV VCF (drop REF records, rewrite `SVTYPE=CNV` to `DEL`/`DUP`), and (d) merges all three into one sorted, bgzipped, tabixed VCF per SG.
- Exome SGs, through the same three branches. Production exomes carry a panel-of-normals CNV VCF (§3); where a run produced none, step (c) is skipped and the 10 kb+ targets are reported unassessed (§5.1, §6).
- Feed the merged VCF into `GenotypeBloodGroupsWithRbceq2` in place of the SNV-only VCF.
- Per-SG sex and ploidy from DRAGEN's ploidy estimate, with its coverage metrics as fallback and metamist's reported sex as a cross-check (§5.6), for CNV direction on chrX/chrY, the karyotype gate (§5.4, §5.6), and the `KARYOTYPE` and `SEXCHECK` QC flags (§6).
- A cross-replicate concordance check over the OurDNA 1KG control SGs (§10).

Non-goals (explicit)

- No changes to RBCeq2: not `SvReader`, `SvMatcher`, `db.tsv`, or the region filter. If an option here appears to need a RBCeq2 change, it is rejected.
- RH and GYP hybrid alleles are out of scope on short-read data. RHD/RHCE and GYPA/B/E are near-identical paralogs, so mapping quality drops to zero and short-read SV and CNV calls are noisy or absent; `--RH` is documented as long-read only, and the maintainer's own comparison of DRAGEN SV/CNV against matched long reads for the 2.4.4 release found deletions useful but hybrids not reliably detected. This is revisited when the volume and timing of long-read data are known.
- No CNV calling of our own. Where a DRAGEN run produced no CNV VCF (the exome test run, §3) this design adds no caller; that SG's 10 kb+ targets stay unassessed and the QC says so (§6).
- No change to the blood-group science or allele definitions.
- No re-run or backfill orchestration (tracked separately).

---

## 5. Design

New per-SG stage (working name `PreprocessDragenForRbceq2`) inserted between `FilterAndConvertGvcfsForRbceq2` and `GenotypeBloodGroupsWithRbceq2`, for every SG (§4). It emits one merged VCF per SG.

The existing rbceq2 stage is repointed at it for every SG. For an SG whose run produced no CNV VCF (§5.1) the stage skips §5.4 and §5.5 triages the SV records alone. Execution is bcftools-based, the same image family as the current filter stage. See the mermaid diagram in `dragen_three_source_merge.md`.

### 5.1 Resolve the new inputs

This step resolves the SV and CNV VCF paths per SG, and the two sex-metric files §5.6 reads.

The SNV gVCF is resolved from `sequencing_group.gvcf`, which cpg_flow populates from the gVCF analysis in metamist. That attribute is cpg_flow's, and review (2026-09-18) ruled out extending cpg_flow for the other files and confirmed that `dragen_align` will not register them until its Nextflow refactor. So the stage derives them from the gVCF's path: `<prefix>/recal_gvcf/<sg>.hard-filtered.recal.gvcf.gz` gives `<prefix>/dragen_metrics/<sg>/<sg>.sv.vcf.gz`, `<sg>.cnv.vcf.gz`, `<sg>.ploidy_estimation_metrics.csv` and `<sg>.wgs_coverage_metrics.csv`, the layout verified in both buckets (§3).

The SV VCF, `cnv_metrics.csv` and at least one sex-metric file are expected for every SG; a missing one fails the SG, as the gVCF guard does today. Whether a CNV VCF is expected is not a property of the sequencing type (production exomes have one, the exome test run does not, §3) and is not a config default; it is read from `cnv_metrics.csv`, which every run writes. A run that segmented and called carries `Number of segments` (and `Number of normal samples` when a panel was used) and must have a CNV VCF; a run whose metrics stop after `Number of target intervals` never called and must not. A CNV VCF present without segment metrics, or absent with them, fails the SG. The ploidy file may be absent for a genome with uneven coverage (§5.6), which is the one case the coverage fallback covers; both sex-metric files missing fails the SG. Both VCFs need their `.tbi`.

The derived paths are recorded in the merged VCF's Analysis meta, as the thresholds are (§8). No new metamist analysis types are created for the raw DRAGEN files (Q6 decided, §13).

### 5.2 SNV gVCF to sites VCF

This stage converts the SNV gVCF into a sites VCF; the merge (§5.5) takes its output as one of three inputs, unchanged.

The SNV branch is the `vcf` output of the existing conversion stage, `FilterAndConvertGvcfsForRbceq2` (`stages/blood_group_genotyping/filter_and_convert.py`). This spec adds nothing to that stage and re-decides nothing about it.

The stage reads the gVCF once, restricted to `resources/bg_regions.<genome>.bed` (the restriction is unconditional, so the gVCF `.tbi` is required). `bcftools norm -m -any` then splits multiallelics.

From that intermediate the stage writes two things: the defining-sites extract the QC stage reads, and the rbceq2 input. The rbceq2 input drops every `<NON_REF>` record (the reference blocks and the split-off symbolic twin of each variant) and trims now-unused ALT alleles.

It then ends in a parameter-free `bcftools +fixploidy`, which expands DRAGEN's haploid non-PAR chrX/chrY genotypes to diploid (`1` to `1|1`), because the pinned rbceq2 (2.4.3) asserts a three-character GT and crashes otherwise. 2.4.4 reads haploid GT natively, so whether this invocation stays is decided at the pin bump, for SNV and structural records together (§13, the 2026-09-18 entry). The full rationale, including why this parameter-free default is sufficient and necessary at 2.4.3 and why a sex-aware config re-crashes rbceq2, is in §13 (Q4) and the §7 table.

No reference FASTA is involved, and no genotyping happens here. A single-sample gVCF already carries GT at every variant site, and GATK `GenotypeGVCFs` has no role here.

The stage filters nothing on `FORMAT/DP` or `FORMAT/GQ`. rbceq2 reads a defining site that is absent from its input as confident homozygous reference, so dropping a borderline genotype manufactures a wild-type call rather than a no-call.

Since `ourdna_genomic_atlas#136`, borderline genotypes stay in. rbceq2's own PASS-only rule handles DRAGEN's hard-filter names, and `FlagBloodGroupCallQc` reports DP and GQ per system instead (`ourdna_genomic_atlas#138`, `#140`). The merged VCF inherits exactly this posture; §6 extends the same QC to structural calls.

For an exome SG the same stage also merges in post-hoc calls at defining sites outside the capture design (`popgen_rbceq2#14`). That output is the SNV branch for exomes exactly as the plain conversion output is for genomes; the merge stage takes it unchanged.

### 5.3 SV VCF to normalised

This step normalises the SV VCF so it can be concatenated with the other two branches; no type change is needed.

The DRAGEN SV caller already writes `SVTYPE=DEL/DUP/INS/BND`, so this is housekeeping only: drop `BND` records (they pair by `MATEID` and match no db token), region-restrict to the bg regions (an optimisation, RBCeq2 also filters internally), ensure the sample column name matches the merged VCF, and keep `CIPOS`/`CIEND`/`SVLEN`/`END`. FORMAT stays as written (`GT:FT:GQ:PL:PR:SR`, or without `SR`); §5.5 says why nothing needs harmonising. Sort, bgzip and tabix the result.

`DUP:TANDEM` reported as `<INS>` is correct; it matches DB INS/dup tokens.

### 5.4 CNV VCF to fixed (the one real transform)

This step fixes the DRAGEN CNV VCF so RBCeq2 can match its records; it depends on the sex and ploidy read in §5.6. It makes four deterministic edits and counts what it did per SG. It is skipped for an SG whose run produced no CNV VCF (§5.1).

For an exome CNV VCF (§3) two facts differ from the genome file. There are no `DRAGEN:REF:` records, so edit 1 drops nothing and §6 cannot read assessability from this file. And there is no `cnvLength` filter, because calling is per capture target against a panel of normals and `BC` counts targets rather than 1–2 kb bins; edit 3's `sv_min_bins` floor applies to that `BC` as written, and whether three targets is the right floor for exomes is untested (§11). Everything else, the `SVTYPE` rewrite, the sex-chromosome direction and the karyotype gate, applies unchanged.

1. Drop `DRAGEN:REF:` records (ALT `.`, no `SVTYPE`); these are non-events.
2. Rewrite `SVTYPE=CNV` to `DEL` or `DUP`. On autosomes take the direction from the symbolic ALT (`<DEL>` to `DEL`, `<DUP>` to `DUP`); on sex chromosomes derive direction from `CN` relative to the sample's expected ploidy (§5.6) instead.
3. Keep sub-10 kb `cnvLength` records rather than dropping them as a class; triage them like every other record in §5.5, keeping those with at least 3 bins, rewriting FILTER per record with the original preserved in `INFO/SVFILTER` (§13, Q1 resolved).
4. Apply the karyotype gate: drop, or pass through tagged `SVSRC=CNV_KARYOTYPE`, every chrX/chrY CNV record from a sample whose DRAGEN ploidy estimate is not `XX` or `XY` (open decision, Q5).

Edit 2 is required, not optional: `SvMatcher` runs with `require_same_type=True`, so the DB token type must equal the event `SVTYPE`, `"DEL" == "CNV"` is false, and the ALT fallback only fires when `SVTYPE` is absent (§2.3).

The evidence for edit 3: the SV caller missed one of the two Gerbich 3.6 kb deletions in 150 genomes that the CNV caller found on 3 bins (`sv_cnv_overlap_50_samples.md` §11); 1–2-bin records have SV-caller support in 45% of cases against 75–79% for 3–4 bins (same note, §10). `--no_filter` is never used for this.

The 3-bin Gerbich call is one CNV record with `BC=3`, not three adjacent records (same note, §13). Adjacent same-direction CNV records do occur, but they are copy-number steps at recurrent loci rather than one event split in pieces, and no merge-adjacent step is added (§7).

The job should count what it did, per SG: ref-blocks dropped, rewritten DEL and DUP, dropped sub-threshold, unresolved. `CnvRewriteStats` is the proposed name for that record; no code for it exists yet, in this repo or the old one.

### 5.5 Merge: triage to the database, then one record per event

The merge concatenates the three normalised VCFs into one, then triages the structural records so RBCeq2 sees one record per event. The QC in §6 grades what survives, and the thresholds it uses live in §8.

`bcftools concat` the three normalised VCFs (two for an SG without a CNV VCF), sort, bgzip and tabix into `<sg>.rbceq2_input.vcf.gz`. All must share the same sample column name and `chr`-prefixed hg38 contigs, since RBCeq2 strips `chr` internally.

FORMAT needs no harmonisation across the sources (review asked, 2026-09-18). VCF FORMAT is per record: gVCF-derived rows carry `GT:AD:DP:GQ:...`, SV rows `GT:FT:GQ:PL:PR:SR` (or without `SR`), CNV rows `GT:SM:CN:BC:PE`; across the 150 study genomes no other FORMAT string appeared in either structural file. rbceq2 requires only that `GT` is the first key of every retained row (`IO/vcf.py`, `_require_first_gt`), which holds for all three, and `SvReader` takes `SVTYPE`, `SVLEN`, `END`, `CIPOS` and `CIEND` from INFO. `bcftools concat` writes the union of the input headers; the six INFO and FORMAT tags the SV and CNV headers share (`GT`, `END`, `SVTYPE`, `SVLEN`, `CIPOS`, `CIEND`) have identical Number and Type in DRAGEN 3.7.8, and a real concat of one genome's region-restricted SV and CNV files ran without a header warning and sorted and indexed cleanly. `validate_for_rbceq2` asserts `GT`-first on every row alongside sorted, single-sample, `chr`-prefixed and indexed.

One FORMAT fact to check at implementation, not before: every PASS CNV duplication in the study carries GT `./1` (244 of 244; deletions carry `0/1` or `1/1`). rbceq2 accepts the string, but how it scores a half-called GT on a matched DUP record is unverified. It matters for exactly one token, GYP\*505.

The merge has to arbitrate because rbceq2 keeps one db definition per record, not per locus. A deletion present as both an SV and a CNV record, offset by the CNV caller's bin snapping, is matched to two different alleles, and the sample is then read as carrying two null alleles.

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

The tie error would need the two records to have identical coordinates; they never do. 2.4.4's `ambiguous_equal_best_sv_evidence` error does not catch this: it fires only on identical coordinates with conflicting GT or FILTER, which the two callers never produce together (0 of 131 shared events in 150 genomes had identical breakpoints).

That is a gap in rbceq2 worth raising upstream (§12), and until it is closed the merge must guarantee one record per event. Empirical basis in `sv_cnv_overlap_50_samples.md`.

Triage, not caller ownership, decides what survives. The study's read-level check showed the SV caller missing a real 3.6 kb Gerbich deletion that the CNV caller found on 3 bins, so neither caller can own a band by size alone (`sv_cnv_overlap_50_samples.md` §11). For an SG without a CNV VCF the rules run over the SV records alone.

#### Triage rules

| Rule | What is kept or done | Why |
|---|---|---|
| 1. Restrict both structural VCFs to db-relevant records | Keep an SV or CNV record only if it (a) overlaps a db SV definition within `SvMatcher`'s own positional and length tolerance (its defaults; use the same code), or (b) is a PASS deletion under `sv_del_max_bp` (1 Mb) spanning a defining SNV/indel site in `bg_site_systems.<genome>.tsv`. Discard everything else. | Removes whole-arm `DUP:TANDEM` artefacts and megabase karyotype events, which §5.6 handles instead. |
| 2. Keep every surviving record from either caller, whatever its FILTER | A `cnvLength` CNV record is dropped if `BC < sv_min_bins`; a surviving non-PASS record has its FILTER rewritten to PASS with the original recorded in `INFO/SVFILTER`. | Corroboration from the SV caller is far more likely at higher bin counts, so a low-bin CNV-only record is dropped rather than trusted; `--no_filter` is never used, since it is global. |
| 3. Where two surviving records describe one event, keep the SV record | Triggered by reciprocal overlap of at least `sv_recip_overlap` in the same direction; the dropped partner's ID is recorded in `INFO/SVPARTNER`. | The SV caller's breakpoints are exact; the CNV caller's are bin-snapped, and several sub-10 kb targets differ from each other only by breakpoint. |
| 4. Tag every kept record with `INFO/SVSRC` | `SV`, `CNV`, or `CNV_LOWRES` for a CNV-only record under 10 kb. | Lets the QC (§6) grade the allele it produces. |
| 5. Assert no two surviving records share CHROM, POS, END and SVTYPE | Fail the sample if they do. | Guarantees the one-record-per-event property the merge exists to provide. |

Rule 3 is the reciprocal-overlap collapse a `bedtools merge` would do (review suggested it), implemented in the job with `SvMatcher`'s own tolerance code rather than bedtools, so the kept record can carry the partner's ID and the source tag, and so the overlap test is the one rbceq2 itself applies. Review also noted that a copy-loss and a deletion record could be merged into one synthetic record rather than one kept; the outcome is the same allele either way, and keeping the SV record with `SVPARTNER` preserves both callers' evidence for the QC without inventing coordinates. Rule 3 is written for both directions, but in 150 genomes it fired only on deletions, and the one DUP token in the database (GYP\*505) is the only case where a duplication pair could matter.

#### What the study showed

- Rule 1(a) kept 0 to 2 records per genome, 3 in total across 150 genomes: the two Gerbich records and the Gerbich CNV-only record.
- Rule 1(b) kept another 0 to 3 records per genome, including a recurrent 9–13 kb deletion at chr19:48.69 Mb spanning a FUT2 defining site in 5 of 150 genomes, matching no db allele; rbceq2 ignores it, but the QC must not (§6).
- Rule 1 also removed a 126 Mb `MaxDepth` SV record seen once over AUG and RHAG.
- Rule 2's 3-bin floor (`sv_min_bins`) reflects that the SV caller corroborates 45% of 1–2-bin deletions against 75–79% of 3–4-bin ones.
- Rule 3 fired only on the Gerbich case: the SV caller's breakpoints were exact to `CIPOS` (0–50 bp) where the CNV caller's were bin-snapped by 0.2–4.6 kb, and the sub-10 kb targets (seven GE alleles, three A4GALT, the GYP cluster) differ from each other only by breakpoint.
- The only allele actually called across the study was the Gerbich deletion.

Haploid GT can also appear on non-PAR chrX/chrY SV and CNV records once the merge exists, at a similarly low rate to the one seen on chrX CNV records in the study. The merge does not add a second `bcftools +fixploidy` invocation for this.

That question is decided once, on the 2.4.4 pin bump, for SNV and structural records together, because 2.4.4 claims native haploid support that may make the invocation redundant everywhere (§13, Q4 and the 2026-09-18 entry).

### 5.6 Sex and ploidy

This step derives each sample's sex and ploidy from DRAGEN's own estimate, for use by CNV direction (§5.4), the karyotype gate (§5.4), and the `KARYOTYPE` and `SEXCHECK` QC flags (§6).

For X-linked systems (XK/Kx, XG, CD99), expected copy number depends on the sample's sex and on region, PAR versus non-PAR. RBCeq2 infers zygosity from GT without knowing either.

#### Sources, in order

1. **`Ploidy estimation` in `ploidy_estimation_metrics.csv`**, DRAGEN's karyotype call from X and Y median coverage over autosomal median coverage. Primary, because it is the only source that names X0 and XYY. Present for all 150 study genomes and all 10 test exomes (§3).
2. **`wgs_coverage_metrics.csv`**, the fallback when the ploidy file is absent. Review (2026-09-18) reports DRAGEN omitting the ploidy file for genomes with uneven coverage: three Garvan DSP samples in `tenk10k-phase2`, none from AGRF. The coverage file is written regardless and carries `Average chr X coverage over genome`, `Average chr Y coverage over genome` and `Average autosomal coverage over genome`, so the same two ratios can be formed. They are classified with the §8 bands: `XX` when X/autosomal is at least `sex_xx_x_min` and Y/autosomal is below `sex_xx_y_max`; `XY` when both ratios fall in `sex_xy_x_range` and `sex_xy_y_range`; anything else `UNDETERMINED`, which the karyotype gate treats as a non-XX/XY estimate. The bands are set from the 150 genome ploidy files, where XX genomes had X/autosomal 0.97–1.02 and Y/autosomal 0.00, XY had 0.50–0.52 and 0.42–0.52, the XYY genome had Y 0.93 and the X0 genome Y 0.18. Two exomes read: an `XX` with X 0.88 and Y 0.01, and an `XY` with X 0.50 and Y 0.35, the latter on the lower edge of `sex_xy_y_range`. Capture design shifts the Y ratio, so the fallback will return `UNDETERMINED` for some XY exomes; that gates them conservatively rather than wrongly, and the bands are revisited if the ploidy file ever goes missing on an exome.
3. **Metamist reported sex**, `sequencing_group.pedigree.sex`, which cpg_flow populates from the participant's `reportedSex` (`cpg_flow/inputs.py`). A cross-check, never the gate's input, because reported sex cannot express X0 or XYY. Disagreement with the estimate from source 1 or 2 raises `SEXCHECK:<reported>/<estimated>` on the X-linked systems (§6) and nothing else changes.

Not used: the `##referenceSexKaryotype` header, a constant reading `XXYY` for every sample, verified 102/102 (§7); `.ploidy.vcf.gz`, absent for exomes and whose `##estimatedSexKaryotype` header repeats source 1; and `SEX GENOTYPER` in `cnv_metrics.csv`, which agreed with source 1 in 49 of 50 genomes and was blank for the X0 genome, so it fails on the sample that needs it. The SNV GT diploid-isation in the conversion stage (§5.2) is not a consumer; it needs no sex or PAR input (§13, Q4).

Expected CN is region times sex, not a blanket 'chrX = 1 in males' rule. Of the seven chrX structural alleles, four sit in PAR1 and are diploid in males; only XK and ATP11C are genuinely hemizygous (CN=1 in males).

#### Expected copy number by region and karyotype

| Region | XX | XY |
|---|---|---|
| PAR1/PAR2 (CD99, XG) | 2 | 2 |
| Non-PAR chrX (XK, ATP11C) | 2 | 1 |
| chrY | 0 | 1 |

- PAR1: CD99\*01N.01/02 (about 2.71 Mb), XG\*01N.02/03 (about 2.78 Mb, on the PAR1 boundary). Diploid in males (baseline CN=2).
- Non-PAR chrX: XK (37.7 Mb), ATP11C (139.7 Mb). Hemizygous in males (CN=1).

A naive haploid-X rule would expect CN=1 in PAR and miscall the normal CD99/XG state as a deletion. So the CNV direction (§5.4) must be PAR-mask aware: diploid in PAR1/PAR2, hemizygous only in non-PAR X/Y for males (chrY: 1 in males, 0 in females).

#### Karyotype gate

The gate reads the sample's estimate from source 1, or source 2 when the ploidy file is absent. For any value other than `XX` or `XY`, `UNDETERMINED` included, it drops or tags the sample's chrX/chrY CNV records (§5.4, open decision Q5), and the QC reports the X-linked systems as `KARYOTYPE:<estimate>` rather than assessing them.

What the study observed, across the 150 genomes:

- The 62 XY samples produced no large chrX CNV event, so DRAGEN's caller handles a normal male.
- The one genome estimated X0 carried PASS heterozygous CN=1 deletions of 1.4–10.6 Mb across XK, CD99, XG and ATP11C.
- The one XYY genome carried a PASS CN=3 duplication across PAR1 (CD99, XG).

Both escaped a false allele only because the length gate rejected megabase events against 11–219 kb tokens; a shorter segment would not be rejected.

### 5.7 Wire into the caller stage

This step repoints `GenotypeBloodGroupsWithRbceq2` at the merged VCF; it depends on the FILTER rewrite for sub-10 kb CNVs in §5.4.

`GenotypeBloodGroupsWithRbceq2` (`stages/blood_group_genotyping/genotype.py`) changes only its input: `--vcf` now points at `<sg>.rbceq2_input.vcf.gz`.

Do not use `--no_filter` to admit sub-10 kb CNVs; it is global and would also let through non-PASS SNVs and SVs. Sub-10 kb CNV records that are kept have their FILTER rewritten to PASS in preprocessing instead (§5.4).

Leave `--phased` off, never required for detection, and `--RH` off, out of scope (§4).

---

## 6. QC for structural calls

This section extends `FlagBloodGroupCallQc` to structural calls, emitting seven new flags. For each defining SNV site the QC already asks whether the caller looked and how well; structural alleles need the same two answers from different evidence.

| Question | CNV VCF (10 kb+ targets) | gVCF depth (sub-10 kb targets) |
|---|---|---|
| Did the caller assess the region? | Yes, near-complete tiling | Not directly; inferred from depth |
| How good is the call? | `QUAL`, `CN`, `SM`, `BC`, FILTER | `QUAL`, `PR`/`SR`, `CIPOS`, or bin evidence |

#### Did the caller assess the region?

`DRAGEN:REF:` records tile every 10 kb+ target in 98–100% of genomes, with 16–195 bins. A gap, or a `cnvQual` event, is the analogue of `NOCOV`.

The structural VCFs give nothing directly for sub-10 kb targets: the SV VCF holds events only, and such a target holds 1–7 bins. But the gVCF already tiles the interval with reference blocks carrying `DP`/`MIN_DP`, read once by the conversion stage.

#### How good is the call?

For 10 kb+ targets, quality comes from `QUAL`, `CN`, segment mean (`SM`), bin count (`BC`) and FILTER.

For sub-10 kb targets, the SV caller gives `QUAL`, `PR`/`SR` and `CIPOS`; a CNV-only call gives `BC`, `SM`, `QUAL`, and the fact that the SV caller saw nothing. Mean depth over the target against its flanks is the same signal the study's read-level check used: about half for a heterozygous deletion, near zero for homozygous, flat for none.

Design:

1. The site-system map gains interval rows. `bg_db.py` stops dropping `kind == 'sv'`; each SV definition becomes a `(chrom, start, end, system, allele)` row. The QC job reads structural records from either the merged VCF or rbceq2's debug log, which in 2.4.4 names the source record for each SV match; which of the two is the source of record is not yet decided (Q7).
2. New flags, joined with `+` to provenance like today's:
   - `SVNOCOV:<system>(<allele>,gap=<bp>)`: a 10 kb+ target not tiled by CNV records.
   - `SVLOWRES:<system>(<allele>,src=CNV,BC=<n>,SM=<x>,QUAL=<q>)`: an allele called from a §5.4 CNV-only sub-10 kb record. Provisional by construction.
   - `SVDEPTH:<system>(<allele>,ratio=<x>,DP=<n>,flank=<n>)`: gVCF depth over a sub-10 kb target falls below `sv_depth_ratio` of its flanking depth with no kept record to explain it. This is the dosage-drop signal without a breakpoint call.
   - `SVUNASSESSED:<system>`: a sub-10 kb target with no gVCF record over the interval at all, an exome hole or an unmapped region, so neither a call nor its absence can be judged.
   - `KARYOTYPE:<estimate>`: X-linked systems in a non-XX/XY sample (§5.6), `UNDETERMINED` from the coverage fallback included.
   - `SEXCHECK:<reported>/<estimated>`: X-linked systems where metamist's reported sex disagrees with DRAGEN's estimate (§5.6).
   - `SVDEL:<system>(<site>,del=<chrom:pos-end>,src=<caller>,GT=<gt>)`: a kept deletion, matched to a db allele or not, spans a defining SNV/indel site of the system.

   A structural call that passes everything is not listed, so `PASS` keeps meaning 'nothing to report'.

   For an exome the CNV VCF has no reference records (§3), so "did the caller assess the region" is read from the capture design instead: a 10 kb+ target with capture targets inside it counts as assessed by the panel-of-normals caller, and one with none is `SVUNASSESSED:<system>`. That is the design-gate resource the post-hoc exome recall already reads (`popgen_rbceq2#14`). `SVNOCOV` keeps its genome meaning, a gap in the CNV caller's tiling, and is never emitted for an exome. An SG whose run produced no CNV VCF has every 10 kb+ target `SVUNASSESSED`. Sub-10 kb targets follow the gVCF depth path where the interval lies in the capture design and are `SVUNASSESSED` otherwise.

3. The defining-sites extract grows by the sub-10 kb target intervals. `gen_bg_resources.py` emits them as a second BED from the same db parse; `FilterAndConvertGvcfsForRbceq2` extracts `DP`/`MIN_DP`/`END` over them and over a flank each side in the same pass it already makes for the SNV sites.
4. What gVCF depth does not give, accepted for the first cut: the mapping-quality-zero fraction and discordant-pair evidence a CRAM read would. DRAGEN's DP already excludes reads failing its mapping filters, so poor mappability appears as low depth rather than as its own signal, which a flag can live with. Breakpoint confirmation is the callers' job, not the QC's.

#### Why `SVDEL` exists

`SVDEL` exists because rbceq2 does nothing with an unmatched deletion. A structural record enters its variant pool only under the db token it matched, and the zygosity adjustment that turns a homozygous call inside a deletion into hemizygous (`modify_variant_pool_if_large_indel`) only sees pool entries.

So a gVCF `A/A` under a heterozygous deletion that matches no db SV is reported homozygous, the same silent-wrong-call class as absent-means-reference. The QC has to say so, as it already does for a small deletion that removed the base (`DEL`).

---

## 7. Alternatives rejected

| Alternative | Why rejected | Where the evidence is |
|---|---|---|
| Size-band caller ownership: the SV caller owns records under 10 kb, the CNV caller owns everything above | The CNV-only Gerbich deletion shows the SV caller missed a real sub-10 kb event that the CNV caller found on 3 bins | `sv_cnv_overlap_50_samples.md` §11 |
| Merge adjacent same-direction records into one event before triage (review suggestion, 2026-09-18) | Measured over the 150 genomes: 254 adjacent CNV pairs within 2.5 kb, only 3 with the same `CN`, and no merged span matched a db allele that neither part matched; 1 adjacent SV pair in 17,175 records. The adjacencies are copy-number steps at recurrent loci, not one event in pieces. The 3-bin Gerbich call is a single record. Revisit if a carrier ever shows a target split across records, where each fragment would fail rbceq2's length gate alone | `sv_cnv_overlap_50_samples.md` §13 |
| Keep exome SGs out of scope | Every production exome has an SV VCF and a panel-of-normals CNV VCF, both with records in the blood-group regions; only the 10-exome test run lacks the CNV VCF | §3, `resolvability_by_input_class.md` §3 |
| Decide whether a CNV VCF is expected from the sequencing type, or from a config flag | Production exomes have one and the exome test run does not, so the type is the wrong key; a config flag would be a default standing in for a fact the run already records in `cnv_metrics.csv` | §5.1 |
| Metamist reported sex as the karyotype source | Cannot express X0 or XYY, the two estimates the gate exists for; kept as a cross-check | §5.6 |
| `SEX GENOTYPER` in `cnv_metrics.csv` as the sex fallback | Blank for the X0 genome, so it fails on exactly the sample that needs it; `wgs_coverage_metrics.csv` is written regardless | §5.6 |
| New metamist analysis types for the raw SV and CNV VCFs | `dragen_align` cannot register them before its Nextflow refactor, cpg_flow is not to be extended, and the paths derive from the registered gVCF's prefix | §5.1 |
| `--no_filter` to admit sub-10 kb CNVs | It is global: it would also admit non-PASS SNVs and non-PASS SVs, not only the CNV records it was meant for | §5.4, §5.7 |
| The `##referenceSexKaryotype` header as the sex source | It is a constant `XXYY` in 102 of 102 samples, not a per-sample estimate | §5.6 |
| Deferring dedup to `select_best_per_vcf` | It keeps one db definition per record, not per locus, so two records for one deletion yield two alleles | §5.5, `sv_cnv_overlap_50_samples.md` §4 |
| The SV VCF alone, no CNV VCF | The study's `sv_only` policy found the same single allele, but the CNV-only Gerbich deletion shows a real event the SV caller misses | `sv_cnv_overlap_50_samples.md` §3, §11 |
| A sex-aware `+fixploidy` config (male non-PAR X set to ploidy 1) | Leaves the call haploid and the pinned RBCeq2 crashes; verified empirically | §13 (Q4) |
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
| `sex_xx_x_min` | 0.8 | Coverage-fallback classifier (§5.6): minimum X/autosomal ratio for `XX` | 150 ploidy files: XX genomes 0.97–1.02, XY 0.50–0.52; one exome 0.88 |
| `sex_xx_y_max` | 0.1 | Coverage-fallback classifier: maximum Y/autosomal ratio for `XX` | 150 ploidy files: XX genomes 0.00; the X0 genome 0.18 |
| `sex_xy_x_range` | 0.4–0.6 | Coverage-fallback classifier: X/autosomal band for `XY` | 150 ploidy files: XY genomes 0.50–0.52 |
| `sex_xy_y_range` | 0.35–0.65 | Coverage-fallback classifier: Y/autosomal band for `XY`; the XYY genome (0.93) and the X0 genome (0.18) fall outside it | 150 ploidy files: XY genomes 0.42–0.52 |

Each threshold is recorded in the Analysis meta of the stage that applies it, as `min_depth`/`min_gq` are today; §9 says which config section holds which.

---

## 9. Touch points

| File / symbol | Change |
|---|---|
| `src/popgen_rbceq2/stages/blood_group_genotyping/filter_and_convert.py` (`FilterAndConvertGvcfsForRbceq2`) | Unchanged. Its `vcf` output is the SNV branch (§5.2); the merge stage lists it as a required stage |
| `src/popgen_rbceq2/stages/blood_group_genotyping/` | New `PreprocessDragenForRbceq2` per-SG stage for every SG, with the CNV branch conditional on sequencing type (§5.1); repoint `GenotypeBloodGroupsWithRbceq2`'s required stages and `--vcf` (`genotype.py`) at its output |
| `src/popgen_rbceq2/jobs/` | New job module. Proposed function names, no code yet: `resolve_dragen_inputs`, `read_sample_sex`, `fix_dragen_cnv_vcf`, `normalize_sv`, `merge_variant_vcfs`, `validate_for_rbceq2` |
| SV/CNV input resolution | Paths derived from `sequencing_group.gvcf`'s DRAGEN output prefix (§5.1), expected-file check by `sequencing_group.sequencing_type`, paths recorded in the merged VCF's Analysis meta. Reported sex from `sequencing_group.pedigree.sex` (Q6 decided) |
| `src/popgen_rbceq2/config/popgen_rbceq2_default_config.toml`, `[workflow.preprocess_dragen_for_rbceq2]` (new) | `cnv_svtype_from`, `min_size`, resources, the merge thresholds from §8 (`sv_min_bins`, `sv_recip_overlap`, `sv_lowres_max_bp`, `sv_del_max_bp`) and the sex-fallback bands (`sex_xx_x_min`, `sex_xx_y_max`, `sex_xy_x_range`, `sex_xy_y_range`) |
| `src/popgen_rbceq2/resources/bg_regions.GRCh38.bed` | Reused unchanged for region-restrict |
| `src/popgen_rbceq2/resources/bg_site_systems.GRCh38.tsv` | Gains SV definition rows (§6) |
| `src/popgen_rbceq2/scripts/bg_db.py` (`SiteKind`) | Stops excluding `kind == 'sv'` rows when building the site-system map |
| `src/popgen_rbceq2/scripts/gen_bg_resources.py` | Emits the sub-10 kb target BED (§6) from the same db parse |
| `src/popgen_rbceq2/stages/blood_group_qc/call_qc.py` (`FlagBloodGroupCallQc`), `src/popgen_rbceq2/jobs/rbceq2_call_qc_job.py` | Today raises `ValueError` if the site-system map carries any `kind == 'sv'` row; gains the §6 QC design once that map has SV rows |
| `[workflow.flag_blood_group_call_qc]` in the same config file | Holds `min_depth = 10`, `min_gq = 20` today; gains the QC thresholds from §8: `sv_depth_ratio`, `sv_depth_flank_bp` |
| `src/popgen_rbceq2/constants.py` (`RBCEQ2_VERSION`, `RBCEQ2_IMAGE_TAG`) | Pinned at `2.4.3` / `2.4.3-1`. The bump to 2.4.4 is its own PR: image, regenerated resources, and the `+fixploidy` re-test (§13, Q4) |

---

## 10. Testing plan, expected outcome and success criteria

Checks run in two places. Constructed-record tests run in CI; concordance checks run by hand against real files in GCS. Two further tests need samples the team does not yet have.

#### What runs in CI

Unit tests use constructed records committed under `tests/`, shaped like the records in §3 and the study's observations. No extract of any cohort genome enters the repo. CI cannot read `gs://cpg-ourdna-main/…`, so no unit test may resolve a live GCS path.

| Test | Input (constructed) | Asserts |
|---|---|---|
| `fix_dragen_cnv_vcf` | A CNV VCF holding a `DRAGEN:REF:` block and `<DEL>`/`<DUP>` events | REF records dropped; every surviving `SVTYPE` in {DEL,DUP}; counts match `CnvRewriteStats` |
| `merge_variant_vcfs` | An SNV VCF, SV VCF and CNV VCF, with a chrX record for the ploidy path, each with its own FORMAT string as in §3 | Output is sorted, single-sample, tabix-indexed, `chr`-prefixed; header holds the union of the three; every row has `GT` first |
| `validate_for_rbceq2` | A merged VCF with one row whose FORMAT is `DP:GT` | Fails naming the row |
| No-CNV path | An SNV VCF and SV VCF, a `cnv_metrics.csv` that stops after the target-interval count, no CNV VCF | Merge succeeds; every 10 kb+ target is `SVUNASSESSED`; the same inputs plus a CNV VCF, or segment metrics without a CNV VCF, fail the SG |
| Exome CNV path | An exome-shaped CNV VCF: events only, no `DRAGEN:REF:` records, PASS and `cnvQual`, `BC` in target counts, with `Number of segments` in the metrics | Edit 1 drops nothing; triage keeps the db-relevant records; the QC reads assessability from the capture design and emits no `SVNOCOV` |
| Sex sources | No ploidy file and a coverage metrics file with ratios shaped like an XX, an XY and an X0 genome; a ploidy file reading `XY` with reported sex female | `XX`, `XY`, `UNDETERMINED`; `SEXCHECK:female/XY` on the X-linked systems |
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

Both need a sample the replicates cannot supply. Being XX, they exercise neither the haploid male chrX/chrY paths nor a positive structural hit.

- A true positive from a database-derived target: a sample with a known 10 kb+ deletion, asserting the allele is called. Concrete targets: XK\*N.05 (8 kb, non-PAR chrX, must come via the SV VCF, exercises the sub-10 kb path), and a larger non-PAR deletion from the CNV VCF, XK\*N.01 (53 kb) or ATP11C\*01N.01 (219 kb), which exercises the hemizygous ploidy path.
- PAR-versus-hemizygous ploidy, needing a male genome: confirm the region-aware expected CN (§5.6). CD99/XG (PAR1, diploid in males) are not miscalled as deletions, while XK/ATP11C (non-PAR) are handled hemizygously.

#### Expected outcome

RBCeq2 gains the ability to call the roughly 40 large-CNV and 38 large-indel blood-group alleles it cannot call today, from data DRAGEN already produces, with no change to RBCeq2 and one new preprocessing stage. Hybrid RH/GYP alleles remain out of reach on short-read data.

These alleles are rare: across the 150 study genomes the design calls one, the Gerbich deletion. The measurable gain is therefore correctness on the carriers that do occur, plus QC that reports unassessed targets and the deletions over defining sites that RBCeq2 ignores.

#### Success criteria

- The constructed-record tests above pass.
- The 150-genome per-band counts and the single GE hit reproduce under the `pass` policy.
- No X-linked allele is called in the X0 or XYY fixtures.
- The merged VCF passes `validate_for_rbceq2`.
- An SG without a CNV VCF runs end to end on the gVCF and SV VCF alone and its QC names every 10 kb+ target `SVUNASSESSED`; an exome with a CNV VCF runs through §5.4 and its QC emits no `SVNOCOV`.

---

## 11. Risks

- Short-read paralog loci. Beyond the RH/GYP hybrids already excluded, any blood-group gene with a close paralog risks mismapped or absent CNV calls. Concordance across the 1KG replicates is a necessary check for QC.
- DRAGEN version. Our calls are DRAGEN 3.7.8. Review (Slack, 2026-09-18) relays Genomics England's internal view that DRAGEN 4.x SV/CNV calls outperform Canvas and Manta and are worth using, while 3.x calls are not. The maintainer's evidence that DRAGEN SV/CNV detects blood-group deletions may have come from 4.x; the author is asking which version was used (Q8). If it was 4.x, the study's own concordance and read-level checks are the evidence this design stands on for 3.7.8, and a true-positive carrier (§10) becomes more important, not less.
- Sub-10 kb CNV-only records. A kept CNV-only sub-10 kb record is lower-grade evidence than an SV-confirmed one; the QC grades it `SVLOWRES` rather than treating it as equivalent, and a 1–2-bin record is dropped rather than kept.
- Half-called duplication GT. Every PASS CNV duplication carries GT `./1` (§5.5). rbceq2 accepts it, but its scoring of a half-called GT on a matched DUP is unverified; only GYP\*505 is exposed.
- Exome sensitivity and thresholds. The SV caller on capture data sees breakpoints only where reads reach them; a 3.6 kb exon deletion whose breakpoints fall in intronic sequence outside the capture may be missed even though the exon itself is covered. The exome CNV caller is a different instrument from the genome one (§3): `sv_min_bins` was set from 1–2 kb genome bins, exome `BC` counts capture targets, and the panel-of-normals caller's false-positive rate in the blood-group regions (30 event records in one exome's regions, against 1 to 7 PASS events per genome) is unmeasured. The exome path therefore adds structural alleles as a possibility, not a promise, until an exome study like `sv_cnv_overlap_50_samples.md` is run; `SVLOWRES` and `SVDEPTH` are the flags that keep an exome call honest meanwhile.
- Positive controls. The 1KG replicate control samples showed no blood-group CNVs at GYP/XK/RHD in the earlier data look, good for a concordance/regression baseline, but the study's Gerbich carrier partly answers this for one 3.6 kb sub-10 kb allele. A 10 kb or larger true positive is still wanted to exercise the CNV-VCF path directly (§10).
- Merge correctness. Contig or sample-name mismatches, or unsorted concat, will silently break RBCeq2's region fetch; `validate_for_rbceq2` must assert sorted, single-sample, `chr`-prefixed, indexed.
- Residual haploid-GT dosage. Diploid-ising a true haploid call (`1` to `1|1`) reads as HOM (dosage 2, `core_logic/alleles.py:284`), overstating dosage for a truly hemizygous call. Fine for detection, but relevant to any future zygosity-dependent filter, and 2.4.4 may read hemizygosity correctly from the haploid GT that `+fixploidy` currently erases (§13, Q4).

---

## 12. Upstream requests

To raise upstream with the RBCeq2 maintainers, feature requests rather than things this design works around beyond reporting:

1. Matching is per record, not per locus. Two records for one event yield two alleles (§5.5 figure). Our merge guarantees one record per event; a per-locus selection in `select_best_per_vcf`, or a warning when two records match overlapping tokens of one system, would make that unnecessary.
2. Unmatched deletions are ignored. A deletion that matches no db token never enters the variant pool, so the hemizygosity adjustment for SNVs inside it never runs and a single-copy genotype is read as homozygous. Our QC reports it (`SVDEL`, §6); rbceq2 could apply the same adjustment it already applies to matched deletions.
3. Undefined null alleles. A deletion removing coding sequence of a blood-group gene is a null whether or not the db names it. The study found a recurrent 9–13 kb deletion over a FUT2 defining site in 5 of 150 genomes, matching no token. Two asks: consider it as a candidate allele for the db, and consider a generic 'coding deletion in system X, unnamed' call. We do not build gene-model resources for this here; it stays an observation to bear in mind, and a reason the `SVDEL` flag names the system.

---

## 13. Decision log

**2026-09-18, first review round (`popgen_rbceq2#15`, `#16`, and the two Slack threads).** Six decisions from Alexander Stuckey's comments:

1. *Naming.* The SV VCF's caller is the DRAGEN SV caller, which integrates and extends Manta; every doc in this directory now says "SV caller" or "SV VCF", keeping literal record IDs. `INFO/SVSRC` value `MANTA` becomes `SV`. Alternatives rejected: keeping "Manta" as shorthand, which lets claims about Manta's internals stand in for claims about DRAGEN's caller.
2. *Exomes in scope.* Every production exome has an SV VCF and a panel-of-normals CNV VCF (11,917 of 11,917 in `cpg-mackenzie-main`); the 10 exomes of the test run have the SV VCF and no CNV VCF (§3). The merge runs for every SG; whether a CNV VCF is expected is read from `cnv_metrics.csv` (§5.1); §5.4 handles the exome file's lack of reference records and of `cnvLength`; the QC reads exome assessability from the capture design (§6). Alternatives rejected: keeping exomes out and correcting only the wording, which leaves structural alleles uncalled on data that carries the calls; keying the CNV expectation on sequencing type, which the two exome runs contradict.
3. *Q6 decided: no new metamist analysis types.* cpg_flow is not extended and `dragen_align` cannot register the files before its Nextflow refactor, so the stage derives the paths from the registered gVCF's prefix, checks expected presence by sequencing type, fails on deviation, and records the paths in the merged VCF's Analysis meta (§5.1). Alternatives rejected: a small registration stage in this repo, which would register another pipeline's outputs under types nobody else reads.
4. *FORMAT needs no harmonisation.* Verified on the study files and a real concat (§5.5); `validate_for_rbceq2` gains a `GT`-first assertion. Alternatives rejected: none needed.
5. *No merge-adjacent step.* Measured over 150 genomes: adjacent CNV records are copy-number steps, not fragments, and the 3-bin Gerbich call is one record (§7, `sv_cnv_overlap_50_samples.md` §13). Alternatives rejected: an optional pre-triage merge with a gap threshold, which would have merged records of different `CN` in 251 of 254 observed cases.
6. *Sex from three sources in order.* DRAGEN's ploidy estimate, then `wgs_coverage_metrics.csv` when the ploidy file is absent (review reports this on uneven-coverage genomes), with metamist's reported sex as a cross-check raising `SEXCHECK` (§5.6, §8). Alternatives rejected: reported sex as the gate's input (cannot express X0/XYY); `SEX GENOTYPER` in `cnv_metrics.csv` (blank on the X0 genome); `.ploidy.vcf.gz` (absent for exomes).

Also recorded: the `+fixploidy` question raised again at §5.2 is the one already deferred to the pin bump (the 2026-09-18 entry below); §5.2 now points there. The bedtools-merge suggestion is rule 3 (§5.5). Two reviewer observations stand as written: the callers' breakpoints never agree, and the 2.4.4 tie error will not fire on real pairs. The DRAGEN-version question is a new risk (§11) and open question (Q8).

**2026-09-18, Q2 reversed to triage.** Decision: deduplicate the merged structural records by triage, restricting both VCFs to db-relevant records, keeping survivors, preferring the SV record on overlap, rather than by caller ownership. Alternatives rejected: giving the SV caller the sub-10 kb band and the CNV caller everything else; deferring entirely to `select_best_per_vcf`. Consequences: §5.5 rules 1–5 as currently designed; raised the per-record-not-per-locus matching gap upstream (§12).

**2026-09-18, Q1 resolved.** Decision: sub-10 kb CNV records are triaged like any other structural record, not dropped as a class; kept if database-relevant or spanning a defining site and at least 3 bins, FILTER rewritten per record with the original kept in `INFO/SVFILTER`. Alternatives rejected: dropping all sub-10 kb CNV records and sourcing small events from the SV VCF only; using `--no_filter` to admit them. Consequences: the 3-bin floor (`sv_min_bins`, §8) and the per-record FILTER rewrite in §5.4; the QC grades a kept CNV-only sub-10 kb allele `SVLOWRES`/`CNV_LOWRES`.

**2026-09-18, karyotype gate added.** Decision: chrX/chrY CNV records from a sample whose DRAGEN ploidy estimate is not XX or XY are dropped or tagged `SVSRC=CNV_KARYOTYPE` (§5.4, §5.6; choice still open, Q5). Alternatives rejected: passing all CNV records through untagged and relying on the length gate alone. Consequences: the X0 and XYY synthetic fixtures in §10; the `KARYOTYPE` QC flag in §6.

**2026-09-18, QC design added.** Decision: add the §6 QC design, six flags for structural calls (a seventh, `SEXCHECK`, came with the review round above), reading structural records from the merged VCF or rbceq2's debug log (source still open, Q7), plus a second BED for sub-10 kb target depth. Alternatives rejected: none tested; this is a first design. Consequences: `bg_db.py` stops excluding `kind == 'sv'` rows; `gen_bg_resources.py` gains a second BED; thresholds recorded in §8.

**2026-09-18, `+fixploidy` deferred to the pin bump (supersedes an earlier position).** Decision: keep the single `bcftools +fixploidy` invocation in the SNV conversion pipe from `ourdna_genomic_atlas#128`, and do not add a second invocation when the merge lands; whether it moves onto the merged VCF, or goes altogether, is decided once, on the 2.4.4 pin bump, for SNV and structural records together. The 2.4.4 release is described as "focused on supporting haploid encoding in VCF", motivated by DRAGEN and array data, and as keeping "chromosome-copy counts ... distinct", so the crash the invocation works around may be gone, and `1` to `1|1` may then overstate dosage where 2.4.4 would read hemizygosity correctly. The bump PR tests a male genome with and without `+fixploidy` before deciding. Alternatives rejected (superseded position): move the single `+fixploidy` invocation onto the merged VCF once the SV/CNV merge exists, on the reasoning that SV and CNV records can also carry haploid GT on non-PAR chrX/chrY in males. Consequences: §5.2 and §5.5 state one position; the pin-bump decision is out of scope for this PR.

**2026-09-17, §5.2 restated against the conversion stage as it is.** Decision: the SNV branch is the conversion stage's `vcf` output as it exists, and the spec describes that output rather than proposing changes to it. The passage it replaces was written in July 2026 against the stage of that time and had been overtaken four times: `ourdna_genomic_atlas#124` (design only, 2026-07-16) showed that the then hard filter at DP 20 and GQ 30 manufactured false wild-type calls across five replicate SGs of one control; `ourdna_genomic_atlas#128` (2026-07-19) appended `+fixploidy`; `ourdna_genomic_atlas#136` (2026-07-29) removed the DP/GQ filter, taking the replicate cohort's discordant systems from eight to zero; `ourdna_genomic_atlas#138` and `#140` (2026-08-03) rewrote the stage as a single pass with a mandatory region restrict and a defining-sites extract, added `FlagBloodGroupCallQc`, and registered its output. The stages were then ported here (`popgen_rbceq2#5`, 2026-08-10) and the stage gained the exome post-hoc merge (`popgen_rbceq2#14`, 2026-09-07). Alternatives rejected: none; the old text described a filter that no longer exists and an addition that had already landed. Consequences: §5.2 rewritten; exome SGs declared out of scope (§4, §5), reversed by the review round of 2026-09-18 above; §9's conversion-stage row reads "unchanged"; every PR reference in the document is repo-qualified.

**2026-09-17, factual refresh (`popgen_rbceq2#15`), and this branch rebased onto it.** Decision: refreshed the database figures for 2.4.4 (§2.2, §2.3), added `resolvability_by_input_class.md`, recorded the maintainer's advice against `--RH` and for a single combined VCF (§1, §4), rewrote the touch points for this repo's paths and the QC stage that did not exist in July, and noted that 2.4.4's native haploid support makes `+fixploidy` a re-test item. `popgen_rbceq2#16` was first opened from `main` in parallel and its restructure of this file dropped those changes; it now sits on top of `#15` with the refresh folded back into the restructured body. Alternatives rejected: merging both branches to `main` independently, which would have conflicted in this file. Consequences: one SPEC, the figures from `#15`, the design from `#16`.

**Q4 resolved in `ourdna_genomic_atlas#128` (2026-07-19), at rbceq2 2.4.2.** DRAGEN writes haploid `GT` (`1`) on non-PAR chrX/chrY for male samples; rbceq2 2.4.2 assumes diploid GT and hard-crashes on haploid input. Confirmed in source: the zygosity determiner `get_ref` does `assert len(GT) == 3` (`core_logic/data_procesing.py:878`) and runs for every defining variant via `make_variant_pool` (`data_procesing.py:454`), so a haploid `"1"` (length 1) raises `AssertionError`. Corroborating diploid assumptions: `remove_home_ref` and `get_variants` drop only `"0/0"`, not haploid `"0"` (`IO/vcf.py:120,261`), and `split_vcf_to_dfs` asserts a `/` or `|` separator at index 1 (`IO/vcf.py:299`). Fix: a final `bcftools +fixploidy` in the conversion pipe, with no `-s`/`-p` arguments. On `chr`-prefixed hg38 the built-in, unprefixed, ploidy table never matches, so every haploid GT is expanded to diploid (`1` to `1|1`, `0` to `0|0`, `.` to `./.`) while diploid autosome and PAR calls stay untouched, verified empirically. The parameter-free default is both sufficient and necessary: PAR stays diploid in males because DRAGEN already writes it diploid, not because of a mask, so no sex file and no PAR mask are needed. Alternatives rejected: a sex-aware config (male non-PAR X set to ploidy 1) leaves the call haploid and re-crashes rbceq2, verified empirically; the §5.6 sex/PAR machinery is therefore not used for this fix, only for CNV direction. Consequences: necessary and sufficient for the SNV branch at the pinned 2.4.3. Caveat: `1` to `1|1` reads as HOM (dosage 2, `core_logic/alleles.py:284`), overstating a truly hemizygous call; fine for detection, but relevant to any zygosity-dependent filter, since rbceq2's native HEM status is deletion-derived, not from GT. Missed by the XX 1KG fixtures (§10). Re-opened for the 2.4.4 bump by the 2026-09-18 entry above.

**Q3 endorsed in review, resolved.** Decision: CNV direction is taken from the ALT symbol on autosomes, and from `CN` versus region-by-sex ploidy on chrX/chrY, PAR-aware: CD99 and XG sit in PAR1 and are diploid in males, only XK and ATP11C are hemizygous. Alternatives rejected: trusting the ALT symbol everywhere, including sex chromosomes. Consequences: implemented in §5.4 and §5.6.

**2026-08, moved from `ourdna_genomic_atlas`.** Decision: the RBCeq2 stages moved into this repo (`popgen_rbceq2`); this design was never implemented in either repo. Alternatives rejected: none; this was a repository restructure. Consequences: the class names in this spec are still correct, but old file paths are not. `src/ourdna_genomic_atlas/stages.py`, referenced in this spec's earlier touch-points table, now maps to `src/popgen_rbceq2/stages/blood_group_genotyping/` (§9).

---

## 14. Open questions for reviewers

- **Q5, the karyotype gate (§5.4, §5.6).** Drop chrX/chrY CNV records for non-XX/XY samples, or pass them through tagged `SVSRC=CNV_KARYOTYPE`? No reason to prefer either has been recorded; the QC reports `KARYOTYPE:<estimate>` for the sample under both options, so the choice only affects what the merged VCF carries.
- **Q7, the QC's source for structural records.** The merged VCF's own records, or rbceq2's debug log, which in 2.4.4 names the source record for each SV match? Recommendation, not a decision: the merged VCF's records are named first and read as the primary source, with the debug log as a secondary cross-check (§6).
- **Q8, which DRAGEN version the maintainer's evidence came from.** Genomics England's experience (relayed in review) is that DRAGEN 4.x SV/CNV is worth using and 3.x is not; we run 3.7.8 (§11). The author has asked the maintainer. If the answer is 4.x, the design's evidence for 3.7.8 is the 150-genome study and its read-level check alone, and the true-positive carrier in §10 moves up the list.

Q6 (metamist registration of the SV and CNV VCFs) was decided in the review round of 2026-09-18 (§13): no new analysis types; paths derive from the registered gVCF and are recorded in the merged VCF's Analysis meta (§5.1).
