# Design spec: structural-variant support for RBCeq2 blood-group calling

**Status:** draft, awaiting design approval; do not implement.
**Revision:** 17 September 2026. **Author:** Joshua Schmidt. **Reviewers:** (fill in).
**Reader:** a pipeline engineer who has not worked on the rbceq2 stages.
**Decision asked of reviewers:** approve the triage rule and karyotype gate in §5 and the thresholds in §8, or say which to change and why.
**Area:** rbceq2 blood-group genotyping pipeline (`FilterAndConvertGvcfsForRbceq2` to `GenotypeBloodGroupsWithRbceq2` to `CombineRbceq2OutputsPerCohort`), in `popgen_rbceq2`. rbceq2 is pinned at 2.4.3 (`constants.py`); database figures below are from the 2.4.4 release of 2026-09-16, which the pin bump (a separate PR) will adopt.
**Pull-request references:** the stages were ported from `ourdna_genomic_atlas` in August 2026 and both repos number PRs from 1, so every PR below is written `ourdna_genomic_atlas#n` or `popgen_rbceq2#n`.
**Companion docs (same directory):**
- [`dragen_three_source_merge.md`](./research/dragen_three_source_merge.md): the visual overview and mermaid diagram, read this first for the shape.
- [`implement_cnv_rbceq2_research.md`](./research/implement_cnv_rbceq2_research.md): the underlying RBCeq2 source analysis this spec rests on, read against v2.4.2; where 2.4.4 moved a fact, this spec says so.
- [`resolvability_by_input_class.md`](./research/resolvability_by_input_class.md): per-system estimate of what the gVCF resolves today and what the merge adds, from the 2.4.4 database.
- [`sv_cnv_overlap_50_samples.md`](./research/sv_cnv_overlap_50_samples.md): the 150-genome Manta-versus-DRAGEN-CNV study the triage rule and karyotype gate rest on.

---

## 1. Summary

RBCeq2's 2.4.4 database defines 67 non-RH alleles across 20 systems by a single large structural event, and today the pipeline can call none of them, because it feeds RBCeq2 only the SNV gVCF (`research/resolvability_by_input_class.md`). DRAGEN already produces the two structural files needed, a Manta SV VCF and a bin-based CNV VCF, but RBCeq2 reads a single VCF. The design is therefore a per-sequencing-group (SG) preprocessing stage that merges the three files into one. The maintainer's advice with the 2.4.4 release says the same: combine the VCFs so every variant is in one file, and leave `--RH` off because DRAGEN SV/CNV does not reliably detect the RH hybrids.

The merge does three things beyond concatenation. It rewrites DRAGEN's `SVTYPE=CNV` to `DEL` or `DUP`, without which no large deletion matches a database allele. It triages both structural files to records that could match a database allele or that span a defining SNV site, and keeps one record per event, because RBCeq2 matches per record and a deletion seen by both callers would otherwise yield two alleles. It drops or tags chrX and chrY CNV records for samples whose DRAGEN karyotype estimate is not XX or XY.

A new QC design (§6) reports whether the callers assessed each structural target and how well, and reports deletions over defining SNV sites, which RBCeq2 itself ignores.

The evidence is a 150-genome OurDNA study (`research/sv_cnv_overlap_50_samples.md`). Across those genomes the design calls one structural allele, a Gerbich deletion. The value is correctness on rare alleles and honest QC, not call volume; §10 states the expected outcome and success criteria in those terms.

RBCeq2 is not modified. RH and GYP hybrid alleles stay out of scope on short-read data, and so do exome sequencing groups, which have no DRAGEN structural calls to merge.

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
| 2 | SV VCF | `dragen_metrics/<sg>/<sg>.sv.vcf.gz` | Manta-derived | Compatible |
| 3 | CNV VCF | `dragen_metrics/<sg>/<sg>.cnv.vcf.gz` | DRAGEN bin CNV | Rewrite SVTYPE |
| - | Ploidy (sex) | `dragen_metrics/<sg>/<sg>.ploidy_estimation_metrics.csv`, `<sg>.ploidy.vcf.gz` | DRAGEN | Read only |

**1, SNV gVCF** (called variant records, already carrying GT, interleaved with `<NON_REF>` reference blocks; the head below happens to show only ref blocks):

```text
chr1  9997   .  N  <NON_REF>  .  PASS  END=10017  ...  ./.:11,0:11:0:2:...
chr1  10018  .  C  <NON_REF>  .  PASS  END=10018  ...  0/0:17,0:17:15:17:...
```

**2, SV VCF** (Manta; header ALTs `<DEL>`,`<INS>`,`<DUP:TANDEM>`; proper SVTYPE):

```text
chr1  789481  MantaINS:…  G  <INS>  999  PASS  END=789481;SVTYPE=INS;CIPOS=0,7;…
chr1  839442  MantaDEL:…  CACC…ACA  CT  463  PASS  END=839499;SVTYPE=DEL;SVLEN=-57;CIGAR=1M1I57D
```

**3, CNV VCF** (every record `SVTYPE=CNV`; direction only in ALT and `CN`):

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

- New per-SG inputs: resolve the SV and CNV VCF paths from `dragen_metrics/<sg>/`.
- These inputs must be registered in metamist as dedicated analysis types, one for the SV VCF and one for the CNV VCF, before the stage can resolve them (Q6, open questions below).
- A preprocessing step that (a) produces an SNV sites VCF, (b) normalises the SV VCF, (c) fixes the CNV VCF (drop REF records, rewrite `SVTYPE=CNV` to `DEL`/`DUP`), and (d) merges all three into one sorted, bgzipped, tabixed VCF per SG.
- Feed the merged VCF into `GenotypeBloodGroupsWithRbceq2` in place of the SNV-only VCF.
- Per-SG sex and ploidy, read from the DRAGEN ploidy files rather than the header, for CNV direction on chrX/chrY, the karyotype gate (§5.4, §5.6), and the `KARYOTYPE` QC flag (§6).
- A cross-replicate concordance check over the OurDNA 1KG control SGs (§10).

Non-goals (explicit)

- No changes to RBCeq2: not `SvReader`, `SvMatcher`, `db.tsv`, or the region filter. If an option here appears to need a RBCeq2 change, it is rejected.
- RH and GYP hybrid alleles are out of scope on short-read data. RHD/RHCE and GYPA/B/E are near-identical paralogs, so mapping quality drops to zero and short-read SV and CNV calls are noisy or absent; `--RH` is documented as long-read only, and the maintainer's own comparison of DRAGEN SV/CNV against matched long reads for the 2.4.4 release found deletions useful but hybrids not reliably detected. This is revisited when the volume and timing of long-read data are known.
- Exome sequencing groups are out of scope. The merge stage runs for genome SGs only; an exome SG keeps today's path, the conversion stage's output (with its post-hoc records, `popgen_rbceq2#14`) straight into `GenotypeBloodGroupsWithRbceq2`.
- No change to the blood-group science or allele definitions.
- No re-run or backfill orchestration (tracked separately).

---

## 5. Design

New per-SG stage (working name `PreprocessDragenForRbceq2`) inserted between `FilterAndConvertGvcfsForRbceq2` and `GenotypeBloodGroupsWithRbceq2`, for genome SGs only (§4). It emits one merged VCF per SG.

The existing rbceq2 stage is repointed at it for those SGs and keeps reading the conversion output for exomes. Execution is bcftools-based, the same image family as the current filter stage. See the mermaid diagram in `dragen_three_source_merge.md`.

### 5.1 Resolve the two new inputs and register in metamist

This step resolves the SV and CNV VCF paths per SG, the two new inputs the merge needs alongside the SNV gVCF.

The SNV gVCF is already resolved from the sequencing group's registered gVCF (`sequencing_group.gvcf`). The SV and CNV VCFs need an equivalent metamist-backed resolution; the mechanism is left open (Q6). Both files are bgzipped and tabixed in the bucket (`.tbi` present), and missing-file handling should mirror the current gVCF guard.

They must be registered as metamist analyses under new SV and CNV analysis types. Ideally `dragen_align` emits them upstream; failing that, this work adds a small stage to define the analysis types and register the existing outputs.

### 5.2 SNV gVCF to sites VCF

This stage converts the SNV gVCF into a sites VCF; the merge (§5.5) takes its output as one of three inputs, unchanged.

The SNV branch is the `vcf` output of the existing conversion stage, `FilterAndConvertGvcfsForRbceq2` (`stages/blood_group_genotyping/filter_and_convert.py`). This spec adds nothing to that stage and re-decides nothing about it.

The stage reads the gVCF once, restricted to `resources/bg_regions.<genome>.bed` (the restriction is unconditional, so the gVCF `.tbi` is required). `bcftools norm -m -any` then splits multiallelics.

From that intermediate the stage writes two things: the defining-sites extract the QC stage reads, and the rbceq2 input. The rbceq2 input drops every `<NON_REF>` record (the reference blocks and the split-off symbolic twin of each variant) and trims now-unused ALT alleles.

It then ends in a parameter-free `bcftools +fixploidy`, which expands DRAGEN's haploid non-PAR chrX/chrY genotypes to diploid (`1` to `1|1`), because the pinned rbceq2 asserts a three-character GT and crashes otherwise. The full rationale, including why this parameter-free default is sufficient and necessary and why a sex-aware config re-crashes rbceq2, is in §13 (Q4) and the §7 table.

No reference FASTA is involved, and no genotyping happens here. A single-sample gVCF already carries GT at every variant site, and GATK `GenotypeGVCFs` has no role here.

The stage filters nothing on `FORMAT/DP` or `FORMAT/GQ`. rbceq2 reads a defining site that is absent from its input as confident homozygous reference, so dropping a borderline genotype manufactures a wild-type call rather than a no-call.

Since `ourdna_genomic_atlas#136`, borderline genotypes stay in. rbceq2's own PASS-only rule handles DRAGEN's hard-filter names, and `FlagBloodGroupCallQc` reports DP and GQ per system instead (`ourdna_genomic_atlas#138`, `#140`). The merged VCF inherits exactly this posture; §6 extends the same QC to structural calls.

For an exome SG the same stage also merges in post-hoc calls at defining sites outside the capture design (`popgen_rbceq2#14`). That path is unaffected by this design, because exome SGs do not enter the merge stage (§4).

### 5.3 SV VCF to normalised

This step normalises the Manta SV VCF so it can be concatenated with the other two branches; no type change is needed.

Manta already writes `SVTYPE=DEL/DUP/INS/BND`, so this is housekeeping only: region-restrict to the bg regions (an optimisation, RBCeq2 also filters internally), ensure the sample column name matches the merged VCF, and keep `CIPOS`/`CIEND`/`SVLEN`/`END`. Sort, bgzip and tabix the result.

`DUP:TANDEM` reported as `<INS>` is correct; it matches DB INS/dup tokens.

### 5.4 CNV VCF to fixed (the one real transform)

This step fixes the DRAGEN CNV VCF so RBCeq2 can match its records; it depends on the sex and ploidy read in §5.6. It makes four deterministic edits and counts what it did per SG.

1. Drop `DRAGEN:REF:` records (ALT `.`, no `SVTYPE`); these are non-events.
2. Rewrite `SVTYPE=CNV` to `DEL` or `DUP`. On autosomes take the direction from the symbolic ALT (`<DEL>` to `DEL`, `<DUP>` to `DUP`); on sex chromosomes derive direction from `CN` relative to the sample's expected ploidy (§5.6) instead.
3. Keep sub-10 kb `cnvLength` records rather than dropping them as a class; triage them like every other record in §5.5, keeping those with at least 3 bins, rewriting FILTER per record with the original preserved in `INFO/SVFILTER` (§13, Q1 resolved).
4. Apply the karyotype gate: drop, or pass through tagged `SVSRC=CNV_KARYOTYPE`, every chrX/chrY CNV record from a sample whose DRAGEN ploidy estimate is not `XX` or `XY` (open decision, Q5).

Edit 2 is required, not optional: `SvMatcher` runs with `require_same_type=True`, so the DB token type must equal the event `SVTYPE`, `"DEL" == "CNV"` is false, and the ALT fallback only fires when `SVTYPE` is absent (§2.3).

The evidence for edit 3: Manta missed one of the two Gerbich 3.6 kb deletions in 150 genomes that the CNV caller found on 3 bins (`sv_cnv_overlap_50_samples.md` §11); 1–2-bin records have Manta support in 45% of cases against 75–79% for 3–4 bins (same note, §10). `--no_filter` is never used for this.

The job should count what it did, per SG: ref-blocks dropped, rewritten DEL and DUP, dropped sub-threshold, unresolved. `CnvRewriteStats` is the proposed name for that record; no code for it exists yet, in this repo or the old one.

### 5.5 Merge: triage to the database, then one record per event

The merge concatenates the three normalised VCFs into one, then triages the structural records so RBCeq2 sees one record per event. The QC in §6 grades what survives, and the thresholds it uses live in §8.

`bcftools concat` the three normalised VCFs, sort, bgzip and tabix into `<sg>.rbceq2_input.vcf.gz`. All three must share the same sample column name and `chr`-prefixed hg38 contigs, since RBCeq2 strips `chr` internally.

The merge has to arbitrate because rbceq2 keeps one db definition per record, not per locus. A deletion present as both a Manta and a CNV record, offset by the CNV caller's bin snapping, is matched to two different alleles, and the sample is then read as carrying two null alleles.

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
Manta   POS exact (CIPOS 0–50)   |------------------------------|   -> best db token: GE*01.-02.01
CNV     POS bin-snapped +300 bp     |------------------------------| -> best db token: GE*01.-03.03

select_best_per_vcf: "best db definition for THIS RECORD"      x 2 records
                                                                 = GE*01.-02.01 / GE*01.-03.03
                                                                 = compound heterozygote, Ge:-2,-3
one record (either caller) -> one allele -> GE*01.-02.01 / GE*01 heterozygote, correct
```

The tie error would need the two records to have identical coordinates; they never do. 2.4.4's `ambiguous_equal_best_sv_evidence` error does not catch this: it fires only on identical coordinates with conflicting GT or FILTER, which the two callers never produce together (0 of 131 shared events in 150 genomes had identical breakpoints).

That is a gap in rbceq2 worth raising upstream (§12), and until it is closed the merge must guarantee one record per event. Empirical basis in `sv_cnv_overlap_50_samples.md`.

Triage, not caller ownership, decides what survives. The study's read-level check showed Manta missing a real 3.6 kb Gerbich deletion that the CNV caller found on 3 bins, so neither caller can own a band by size alone (`sv_cnv_overlap_50_samples.md` §11).

#### Triage rules

| Rule | What is kept or done | Why |
|---|---|---|
| 1. Restrict both structural VCFs to db-relevant records | Keep a Manta or CNV record only if it (a) overlaps a db SV definition within `SvMatcher`'s own positional and length tolerance (its defaults; use the same code), or (b) is a PASS deletion under `sv_del_max_bp` (1 Mb) spanning a defining SNV/indel site in `bg_site_systems.<genome>.tsv`. Discard everything else. | Removes whole-arm `DUP:TANDEM` artefacts and megabase karyotype events, which §5.6 handles instead. |
| 2. Keep every surviving record from either caller, whatever its FILTER | A `cnvLength` CNV record is dropped if `BC < sv_min_bins`; a surviving non-PASS record has its FILTER rewritten to PASS with the original recorded in `INFO/SVFILTER`. | Corroboration from Manta is far more likely at higher bin counts, so a low-bin CNV-only record is dropped rather than trusted; `--no_filter` is never used, since it is global. |
| 3. Where two surviving records describe one event, keep the Manta record | Triggered by reciprocal overlap of at least `sv_recip_overlap` in the same direction; the dropped partner's ID is recorded in `INFO/SVPARTNER`. | Manta's breakpoints are exact; the CNV caller's are bin-snapped, and several sub-10 kb targets differ from each other only by breakpoint. |
| 4. Tag every kept record with `INFO/SVSRC` | `MANTA`, `CNV`, or `CNV_LOWRES` for a CNV-only record under 10 kb. | Lets the QC (§6) grade the allele it produces. |
| 5. Assert no two surviving records share CHROM, POS, END and SVTYPE | Fail the sample if they do. | Guarantees the one-record-per-event property the merge exists to provide. |

#### What the study showed

- Rule 1(a) kept 0 to 2 records per genome, 3 in total across 150 genomes: the two Gerbich records and the Gerbich CNV-only record.
- Rule 1(b) kept another 0 to 3 records per genome, including a recurrent 9–13 kb deletion at chr19:48.69 Mb spanning a FUT2 defining site in 5 of 150 genomes, matching no db allele; rbceq2 ignores it, but the QC must not (§6).
- Rule 1 also removed a 126 Mb `MaxDepth` Manta record seen once over AUG and RHAG.
- Rule 2's 3-bin floor (`sv_min_bins`) reflects that Manta corroborates 45% of 1–2-bin deletions against 75–79% of 3–4-bin ones.
- Rule 3 fired only on the Gerbich case: Manta's breakpoints were exact to `CIPOS` (0–50 bp) where the CNV caller's were bin-snapped by 0.2–4.6 kb, and the sub-10 kb targets (seven GE alleles, three A4GALT, the GYP cluster) differ from each other only by breakpoint.
- The only allele actually called across the study was the Gerbich deletion.

Haploid GT can also appear on non-PAR chrX/chrY Manta and CNV records once the merge exists, at a similarly low rate to the one seen on chrX CNV records in the study. The merge does not add a second `bcftools +fixploidy` invocation for this.

That question is decided once, on the 2.4.4 pin bump, for SNV and structural records together, because 2.4.4 claims native haploid support that may make the invocation redundant everywhere (§13, Q4 and the 2026-09-18 entry).

### 5.6 Sex and ploidy

This step derives each sample's sex and ploidy from DRAGEN's own estimate, for use by CNV direction (§5.4), the karyotype gate (§5.4), and the `KARYOTYPE` QC flag (§6).

For X-linked systems (XK/Kx, XG, CD99), expected copy number depends on the sample's sex and on region, PAR versus non-PAR. RBCeq2 infers zygosity from GT without knowing either.

Read per-SG sex from `ploidy_estimation_metrics.csv` and `.ploidy.vcf.gz` (X/Y coverage ratios), not the `##referenceSexKaryotype` header, a constant reading `XXYY` for every sample, verified 102/102 (§7). The SNV GT diploid-isation in the conversion stage (§5.2) is not a fourth consumer; it needs no sex or PAR input (§13, Q4).

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

The gate reads `Ploidy estimation` from `ploidy_estimation_metrics.csv` per sample. For any value other than `XX` or `XY`, it drops or tags the sample's chrX/chrY CNV records (§5.4, open decision Q5), and the QC reports the X-linked systems as `KARYOTYPE:<estimate>` rather than assessing them.

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

This section extends `FlagBloodGroupCallQc` to structural calls, emitting six new flags. For each defining SNV site the QC already asks whether the caller looked and how well; structural alleles need the same two answers from different evidence.

| Question | CNV VCF (10 kb+ targets) | gVCF depth (sub-10 kb targets) |
|---|---|---|
| Did the caller assess the region? | Yes, near-complete tiling | Not directly; inferred from depth |
| How good is the call? | `QUAL`, `CN`, `SM`, `BC`, FILTER | `QUAL`, `PR`/`SR`, `CIPOS`, or bin evidence |

#### Did the caller assess the region?

`DRAGEN:REF:` records tile every 10 kb+ target in 98–100% of genomes, with 16–195 bins. A gap, or a `cnvQual` event, is the analogue of `NOCOV`.

The structural VCFs give nothing directly for sub-10 kb targets: Manta emits events only, and such a target holds 1–7 bins. But the gVCF already tiles the interval with reference blocks carrying `DP`/`MIN_DP`, read once by the conversion stage.

#### How good is the call?

For 10 kb+ targets, quality comes from `QUAL`, `CN`, segment mean (`SM`), bin count (`BC`) and FILTER.

For sub-10 kb targets, Manta gives `QUAL`, `PR`/`SR` and `CIPOS`; a CNV-only call gives `BC`, `SM`, `QUAL`, and the fact that Manta saw nothing. Mean depth over the target against its flanks is the same signal the study's read-level check used: about half for a heterozygous deletion, near zero for homozygous, flat for none.

Design:

1. The site-system map gains interval rows. `bg_db.py` stops dropping `kind == 'sv'`; each SV definition becomes a `(chrom, start, end, system, allele)` row. The QC job reads structural records from either the merged VCF or rbceq2's debug log, which in 2.4.4 names the source record for each SV match; which of the two is the source of record is not yet decided (Q7).
2. New flags, joined with `+` to provenance like today's:
   - `SVNOCOV:<system>(<allele>,gap=<bp>)`: a 10 kb+ target not tiled by CNV records.
   - `SVLOWRES:<system>(<allele>,src=CNV,BC=<n>,SM=<x>,QUAL=<q>)`: an allele called from a §5.4 CNV-only sub-10 kb record. Provisional by construction.
   - `SVDEPTH:<system>(<allele>,ratio=<x>,DP=<n>,flank=<n>)`: gVCF depth over a sub-10 kb target falls below `sv_depth_ratio` of its flanking depth with no kept record to explain it. This is the dosage-drop signal without a breakpoint call.
   - `SVUNASSESSED:<system>`: a sub-10 kb target with no gVCF record over the interval at all, an exome hole or an unmapped region, so neither a call nor its absence can be judged.
   - `KARYOTYPE:<estimate>`: X-linked systems in a non-XX/XY sample (§5.6).
   - `SVDEL:<system>(<site>,del=<chrom:pos-end>,src=<caller>,GT=<gt>)`: a kept deletion, matched to a db allele or not, spans a defining SNV/indel site of the system.

   A structural call that passes everything is not listed, so `PASS` keeps meaning 'nothing to report'.

3. The defining-sites extract grows by the sub-10 kb target intervals. `gen_bg_resources.py` emits them as a second BED from the same db parse; `FilterAndConvertGvcfsForRbceq2` extracts `DP`/`MIN_DP`/`END` over them and over a flank each side in the same pass it already makes for the SNV sites.
4. What gVCF depth does not give, accepted for the first cut: the mapping-quality-zero fraction and discordant-pair evidence a CRAM read would. DRAGEN's DP already excludes reads failing its mapping filters, so poor mappability appears as low depth rather than as its own signal, which a flag can live with. Breakpoint confirmation is the callers' job, not the QC's.

#### Why `SVDEL` exists

`SVDEL` exists because rbceq2 does nothing with an unmatched deletion. A structural record enters its variant pool only under the db token it matched, and the zygosity adjustment that turns a homozygous call inside a deletion into hemizygous (`modify_variant_pool_if_large_indel`) only sees pool entries.

So a gVCF `A/A` under a heterozygous deletion that matches no db SV is reported homozygous, the same silent-wrong-call class as absent-means-reference. The QC has to say so, as it already does for a small deletion that removed the base (`DEL`).

---

## 7. Alternatives rejected

| Alternative | Why rejected | Where the evidence is |
|---|---|---|
| Size-band caller ownership: Manta owns records under 10 kb, the CNV caller owns everything above | The CNV-only Gerbich deletion shows Manta missed a real sub-10 kb event that the CNV caller found on 3 bins | `sv_cnv_overlap_50_samples.md` §11 |
| `--no_filter` to admit sub-10 kb CNVs | It is global: it would also admit non-PASS SNVs and non-PASS SVs, not only the CNV records it was meant for | §5.4, §5.7 |
| The `##referenceSexKaryotype` header as the sex source | It is a constant `XXYY` in 102 of 102 samples, not a per-sample estimate | §5.6 |
| Deferring dedup to `select_best_per_vcf` | It keeps one db definition per record, not per locus, so two records for one deletion yield two alleles | §5.5, `sv_cnv_overlap_50_samples.md` §4 |
| Manta alone, no CNV VCF | The study's `sv_only` policy found the same single allele, but the CNV-only Gerbich deletion shows a real event Manta misses | `sv_cnv_overlap_50_samples.md` §3, §11 |
| A sex-aware `+fixploidy` config (male non-PAR X set to ploidy 1) | Leaves the call haploid and the pinned RBCeq2 crashes; verified empirically | §13 (Q4) |
| Checked-in extracts of a cohort genome's DRAGEN outputs as unit fixtures | No individual-level data in the repo; the study's observations are reproduced as constructed records instead | §10 |

---

## 8. Thresholds

| Parameter | Value | What it gates | Where it came from |
|---|---|---|---|
| `sv_min_bins` | 3 | Minimum CNV bin count to keep a sub-10 kb `cnvLength` record (§5.4, §5.5) | `sv_cnv_overlap_50_samples.md` §10–11: 3–4-bin records have Manta support 75–79% of the time, 1–2-bin only 45% |
| `sv_recip_overlap` | 0.5 | Reciprocal overlap required to treat two records as the same event (§5.5, rule 3) | `sv_cnv_overlap_50_samples.md` §2 |
| `sv_lowres_max_bp` | 10000 | Size below which a CNV-only record is graded `CNV_LOWRES` and flagged `SVLOWRES` | `sv_cnv_overlap_50_samples.md` §1, §11 |
| `sv_depth_ratio` | 0.7 | Depth ratio below which gVCF depth over a sub-10 kb target triggers `SVDEPTH` | §6 |
| `sv_depth_flank_bp` | 5000 | Flank size either side of a target used to compute the depth ratio | §6 |
| `sv_del_max_bp` (name proposed) | 1000000 | Upper size cap on a PASS deletion admitted by rule 1(b) for spanning a defining SNV/indel site | §5.5, rule 1(b); caps out the X0 sample's megabase events and the 126 Mb Manta `MaxDepth` record |

Each threshold is recorded in the Analysis meta of the stage that applies it, as `min_depth`/`min_gq` are today; §9 says which config section holds which.

---

## 9. Touch points

| File / symbol | Change |
|---|---|
| `src/popgen_rbceq2/stages/blood_group_genotyping/filter_and_convert.py` (`FilterAndConvertGvcfsForRbceq2`) | Unchanged. Its `vcf` output is the SNV branch (§5.2); the merge stage lists it as a required stage |
| `src/popgen_rbceq2/stages/blood_group_genotyping/` | New `PreprocessDragenForRbceq2` per-SG stage, genome SGs only; repoint `GenotypeBloodGroupsWithRbceq2`'s required stages and `--vcf` (`genotype.py`) at its output for those SGs |
| `src/popgen_rbceq2/jobs/` | New job module. Proposed function names, no code yet: `fix_dragen_cnv_vcf`, `normalize_sv`, `merge_variant_vcfs`, `validate_for_rbceq2` |
| SV/CNV input resolution | The SNV gVCF is resolved from `sequencing_group.gvcf`; the SV and CNV VCFs need an equivalent metamist-backed resolution. Mechanism left open (Q6); whether `dragen_align` registers them today is unverified |
| `src/popgen_rbceq2/config/popgen_rbceq2_default_config.toml`, `[workflow.preprocess_dragen_for_rbceq2]` (new) | `cnv_svtype_from`, `min_size`, resources, and the merge thresholds from §8: `sv_min_bins`, `sv_recip_overlap`, `sv_lowres_max_bp`, `sv_del_max_bp` |
| `src/popgen_rbceq2/resources/bg_regions.GRCh38.bed` | Reused unchanged for region-restrict |
| `src/popgen_rbceq2/resources/bg_site_systems.GRCh38.tsv` | Gains SV definition rows (§6) |
| `src/popgen_rbceq2/scripts/bg_db.py` (`SiteKind`) | Stops excluding `kind == 'sv'` rows when building the site-system map |
| `src/popgen_rbceq2/scripts/gen_bg_resources.py` | Emits the sub-10 kb target BED (§6) from the same db parse |
| `src/popgen_rbceq2/stages/blood_group_qc/call_qc.py` (`FlagBloodGroupCallQc`), `src/popgen_rbceq2/jobs/rbceq2_call_qc_job.py` | Today raises `ValueError` if the site-system map carries any `kind == 'sv'` row; gains the §6 QC design once that map has SV rows |
| `[workflow.flag_blood_group_call_qc]` in the same config file | Holds `min_depth = 10`, `min_gq = 20` today; gains the QC thresholds from §8: `sv_depth_ratio`, `sv_depth_flank_bp` |
| `src/popgen_rbceq2/constants.py` (`RBCEQ2_VERSION`, `RBCEQ2_IMAGE_TAG`) | Pinned at `2.4.3` / `2.4.3-1`. The bump to 2.4.4 is its own PR: image, regenerated resources, and the `+fixploidy` re-test (§13, Q4) |

---

## 10. Testing plan, expected outcome and success criteria

Two kinds of check. Unit tests run in CI on constructed records committed under `tests/`, shaped like the records in §3 and the study's observations: no extract of any cohort genome enters the repo, and CI cannot read `gs://cpg-ourdna-main/…`, so no unit test may resolve a live GCS path. Concordance checks run by hand against the real files of the OurDNA 1000 Genomes (1KG) control replicates, independent sequencings of one public control individual (`NA12878`, XX): five SGs when this spec was first drafted, nine by the time of `ourdna_genomic_atlas#136`. Being XX, the replicates cannot exercise the haploid male chrX/chrY paths; that needs a male genome.

- Unit: `fix_dragen_cnv_vcf` on a constructed CNV VCF holding a `DRAGEN:REF:` block and `<DEL>`/`<DUP>` events, asserting REF records dropped, every surviving `SVTYPE` in {DEL,DUP}, counts match `CnvRewriteStats`.
- Unit: `merge_variant_vcfs` on a constructed SNV VCF, SV VCF and CNV VCF with a chrX record for the ploidy path; output is sorted, single-sample, tabix-indexed, `chr`-prefixed.
- Integration (concordance, manual, needs GCS): run the full preprocess-plus-rbceq2 on every replicate; structural calls must be identical across all of them (same individual). Any divergence is a bug or a QC signal (compare the SNV divergence table in `ourdna_genomic_atlas#124`).
- Sub-10 kb path, synthetic fixtures shaped like the study's observations, no cohort data in the repo: (a) a Manta DEL record of exactly 3609 bp at the GE\*01.-02.01 coordinates plus a 4961 bp `cnvLength` CN=1 CNV record over it, where the merge keeps the Manta record and records the CNV partner, and rbceq2 calls one GE allele, never two; (b) the CNV record alone, 3 bins, kept, PASS-rewritten, tagged, and the QC flags the GE call `SVLOWRES`; (c) the same with 2 bins, dropped, GE unassessed. The study observed both (a) and (b) in real genomes; the fixtures are constructed records, not extracts.
- Karyotype gate, synthetic: a CNV record set with megabase CN=1 deletions across XK, CD99/XG and ATP11C paired with a ploidy metrics file reading `X0`, and a PAR1 CN=3 duplication paired with `XYY`; no X-linked allele is called and the QC reports `KARYOTYPE:<estimate>`.
- Concordance across 150 genomes: the study's scripts, once committed, are the check that a code change to the merge does not alter the per-band record counts or the single GE hit under the `pass` policy. They read live GCS paths, so this is a manual check, not CI.
- True-positive (DB-derived targets), still wanted: a sample with a known 10 kb+ deletion, asserting the allele is called. Concrete targets: XK\*N.05 (8 kb, non-PAR chrX, must come via the SV VCF, exercises the sub-10 kb path) and a larger non-PAR deletion, XK\*N.01 (53 kb) or ATP11C\*01N.01 (219 kb), from the CNV VCF (exercises the hemizygous ploidy path). Needed because the control shows no positive structural hit.
- PAR-versus-hemizygous ploidy (needs a male sample): confirm the region-aware expected CN (§5.6): CD99/XG (PAR1, diploid in males) are not miscalled as deletions, while XK/ATP11C (non-PAR) are handled hemizygously. The XX 1KG replicates cannot provide this check.

Expected outcome

RBCeq2 gains the ability to call the roughly 40 large-CNV and 38 large-indel blood-group alleles it cannot call today, from data DRAGEN already produces, with no change to RBCeq2 and one new preprocessing stage. These alleles are rare: across the 150 study genomes the design calls one, the Gerbich deletion. The measurable gain is therefore correctness on the carriers that do occur, plus QC that reports unassessed targets and the deletions over defining sites that RBCeq2 ignores. Hybrid RH/GYP alleles remain out of reach on short-read data.

Success criteria

- The synthetic fixtures above pass.
- The 150-genome per-band counts and the single GE hit reproduce under the `pass` policy.
- No X-linked allele is called in the X0 or XYY synthetic fixtures.
- The merged VCF passes `validate_for_rbceq2`.

---

## 11. Risks

- Short-read paralog loci. Beyond the RH/GYP hybrids already excluded, any blood-group gene with a close paralog risks mismapped or absent CNV calls. Concordance across the 1KG replicates is a necessary check for QC.
- Sub-10 kb CNV-only records. A kept CNV-only sub-10 kb record is lower-grade evidence than a Manta-confirmed one; the QC grades it `SVLOWRES` rather than treating it as equivalent, and a 1–2-bin record is dropped rather than kept.
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

**2026-09-18, Q2 reversed to triage.** Decision: deduplicate the merged structural records by triage, restricting both VCFs to db-relevant records, keeping survivors, preferring Manta on overlap, rather than by caller ownership. Alternatives rejected: giving Manta the sub-10 kb band and the CNV caller everything else; deferring entirely to `select_best_per_vcf`. Consequences: §5.5 rules 1–5 as currently designed; raised the per-record-not-per-locus matching gap upstream (§12).

**2026-09-18, Q1 resolved.** Decision: sub-10 kb CNV records are triaged like any other structural record, not dropped as a class; kept if database-relevant or spanning a defining site and at least 3 bins, FILTER rewritten per record with the original kept in `INFO/SVFILTER`. Alternatives rejected: dropping all sub-10 kb CNV records and sourcing small events from the SV VCF only; using `--no_filter` to admit them. Consequences: the 3-bin floor (`sv_min_bins`, §8) and the per-record FILTER rewrite in §5.4; the QC grades a kept CNV-only sub-10 kb allele `SVLOWRES`/`CNV_LOWRES`.

**2026-09-18, karyotype gate added.** Decision: chrX/chrY CNV records from a sample whose DRAGEN ploidy estimate is not XX or XY are dropped or tagged `SVSRC=CNV_KARYOTYPE` (§5.4, §5.6; choice still open, Q5). Alternatives rejected: passing all CNV records through untagged and relying on the length gate alone. Consequences: the X0 and XYY synthetic fixtures in §10; the `KARYOTYPE` QC flag in §6.

**2026-09-18, QC design added.** Decision: add the §6 QC design, six flags for structural calls, reading structural records from the merged VCF or rbceq2's debug log (source still open, Q7), plus a second BED for sub-10 kb target depth. Alternatives rejected: none tested; this is a first design. Consequences: `bg_db.py` stops excluding `kind == 'sv'` rows; `gen_bg_resources.py` gains a second BED; thresholds recorded in §8.

**2026-09-18, `+fixploidy` deferred to the pin bump (supersedes an earlier position).** Decision: keep the single `bcftools +fixploidy` invocation in the SNV conversion pipe from `ourdna_genomic_atlas#128`, and do not add a second invocation when the merge lands; whether it moves onto the merged VCF, or goes altogether, is decided once, on the 2.4.4 pin bump, for SNV and structural records together. The 2.4.4 release is described as "focused on supporting haploid encoding in VCF", motivated by DRAGEN and array data, and as keeping "chromosome-copy counts ... distinct", so the crash the invocation works around may be gone, and `1` to `1|1` may then overstate dosage where 2.4.4 would read hemizygosity correctly. The bump PR tests a male genome with and without `+fixploidy` before deciding. Alternatives rejected (superseded position): move the single `+fixploidy` invocation onto the merged VCF once the SV/CNV merge exists, on the reasoning that Manta and CNV records can also carry haploid GT on non-PAR chrX/chrY in males. Consequences: §5.2 and §5.5 state one position; the pin-bump decision is out of scope for this PR.

**2026-09-17, §5.2 restated against the conversion stage as it is.** Decision: the SNV branch is the conversion stage's `vcf` output as it exists, and the spec describes that output rather than proposing changes to it. The passage it replaces was written in July 2026 against the stage of that time and had been overtaken four times: `ourdna_genomic_atlas#124` (design only, 2026-07-16) showed that the then hard filter at DP 20 and GQ 30 manufactured false wild-type calls across five replicate SGs of one control; `ourdna_genomic_atlas#128` (2026-07-19) appended `+fixploidy`; `ourdna_genomic_atlas#136` (2026-07-29) removed the DP/GQ filter, taking the replicate cohort's discordant systems from eight to zero; `ourdna_genomic_atlas#138` and `#140` (2026-08-03) rewrote the stage as a single pass with a mandatory region restrict and a defining-sites extract, added `FlagBloodGroupCallQc`, and registered its output. The stages were then ported here (`popgen_rbceq2#5`, 2026-08-10) and the stage gained the exome post-hoc merge (`popgen_rbceq2#14`, 2026-09-07). Alternatives rejected: none; the old text described a filter that no longer exists and an addition that had already landed. Consequences: §5.2 rewritten; exome SGs declared out of scope (§4, §5); §9's conversion-stage row reads "unchanged"; every PR reference in the document is repo-qualified.

**2026-09-17, factual refresh (`popgen_rbceq2#15`), and this branch rebased onto it.** Decision: refreshed the database figures for 2.4.4 (§2.2, §2.3), added `resolvability_by_input_class.md`, recorded the maintainer's advice against `--RH` and for a single combined VCF (§1, §4), rewrote the touch points for this repo's paths and the QC stage that did not exist in July, and noted that 2.4.4's native haploid support makes `+fixploidy` a re-test item. `popgen_rbceq2#16` was first opened from `main` in parallel and its restructure of this file dropped those changes; it now sits on top of `#15` with the refresh folded back into the restructured body. Alternatives rejected: merging both branches to `main` independently, which would have conflicted in this file. Consequences: one SPEC, the figures from `#15`, the design from `#16`.

**Q4 resolved in `ourdna_genomic_atlas#128` (2026-07-19), at rbceq2 2.4.2.** DRAGEN writes haploid `GT` (`1`) on non-PAR chrX/chrY for male samples; rbceq2 2.4.2 assumes diploid GT and hard-crashes on haploid input. Confirmed in source: the zygosity determiner `get_ref` does `assert len(GT) == 3` (`core_logic/data_procesing.py:878`) and runs for every defining variant via `make_variant_pool` (`data_procesing.py:454`), so a haploid `"1"` (length 1) raises `AssertionError`. Corroborating diploid assumptions: `remove_home_ref` and `get_variants` drop only `"0/0"`, not haploid `"0"` (`IO/vcf.py:120,261`), and `split_vcf_to_dfs` asserts a `/` or `|` separator at index 1 (`IO/vcf.py:299`). Fix: a final `bcftools +fixploidy` in the conversion pipe, with no `-s`/`-p` arguments. On `chr`-prefixed hg38 the built-in, unprefixed, ploidy table never matches, so every haploid GT is expanded to diploid (`1` to `1|1`, `0` to `0|0`, `.` to `./.`) while diploid autosome and PAR calls stay untouched, verified empirically. The parameter-free default is both sufficient and necessary: PAR stays diploid in males because DRAGEN already writes it diploid, not because of a mask, so no sex file and no PAR mask are needed. Alternatives rejected: a sex-aware config (male non-PAR X set to ploidy 1) leaves the call haploid and re-crashes rbceq2, verified empirically; the §5.6 sex/PAR machinery is therefore not used for this fix, only for CNV direction. Consequences: necessary and sufficient for the SNV branch at the pinned 2.4.3. Caveat: `1` to `1|1` reads as HOM (dosage 2, `core_logic/alleles.py:284`), overstating a truly hemizygous call; fine for detection, but relevant to any zygosity-dependent filter, since rbceq2's native HEM status is deletion-derived, not from GT. Missed by the XX 1KG fixtures (§10). Re-opened for the 2.4.4 bump by the 2026-09-18 entry above.

**Q3 endorsed in review, resolved.** Decision: CNV direction is taken from the ALT symbol on autosomes, and from `CN` versus region-by-sex ploidy on chrX/chrY, PAR-aware: CD99 and XG sit in PAR1 and are diploid in males, only XK and ATP11C are hemizygous. Alternatives rejected: trusting the ALT symbol everywhere, including sex chromosomes. Consequences: implemented in §5.4 and §5.6.

**2026-08, moved from `ourdna_genomic_atlas`.** Decision: the RBCeq2 stages moved into this repo (`popgen_rbceq2`); this design was never implemented in either repo. Alternatives rejected: none; this was a repository restructure. Consequences: the class names in this spec are still correct, but old file paths are not. `src/ourdna_genomic_atlas/stages.py`, referenced in this spec's earlier touch-points table, now maps to `src/popgen_rbceq2/stages/blood_group_genotyping/` (§9).

---

## 14. Open questions for reviewers

- **Q5, the karyotype gate (§5.4, §5.6).** Drop chrX/chrY CNV records for non-XX/XY samples, or pass them through tagged `SVSRC=CNV_KARYOTYPE`? No reason to prefer either has been recorded; the QC reports `KARYOTYPE:<estimate>` for the sample under both options, so the choice only affects what the merged VCF carries.
- **Q6, metamist registration of the SV and CNV VCFs.** Upstream in `dragen_align`, or a small stage in this repo? Recommendation, not a decision: the author's stated preference is upstream in `dragen_align`, with a small stage here only as a fallback (§4, §5.1, §9).
- **Q7, the QC's source for structural records.** The merged VCF's own records, or rbceq2's debug log, which in 2.4.4 names the source record for each SV match? Recommendation, not a decision: the merged VCF's records are named first and read as the primary source, with the debug log as a secondary cross-check (§6).
