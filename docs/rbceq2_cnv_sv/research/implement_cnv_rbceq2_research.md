# Implementing CNV/SV support with RBCeq2 — research notes

**Date:** 2026-07-14
**Repo:** `RBCeq2` (v2.4.2), blood-group allele inference
**Goal:** Understand how CNV/SV data is consumed by RBCeq2 and assess compatibility
with our ICA **DRAGEN 3.7.8** outputs (separate SNV gVCF + SV VCF + CNV VCF per sample).

---

## 1. How CNV/SV data is provided to RBCeq2

**There is no separate CNV/SV input and no dedicated CLI flag.** Structural records
go **in the same VCF** passed to `--vcf`, encoded with standard VCF SV fields. RBCeq2
reads them alongside SNVs/indels and fuzzy-matches them to SV alleles in its database.

### Pipeline wiring
`src/rbceq2/main.py:291-305` (`find_hits`):

```python
reader  = SvReader(df=vcf.df, min_size=args.min_size)   # read SVs out of the VCF
events  = list(reader.events())
db_defs = load_db_defs(db.df)                            # read SV tokens out of db.tsv
matcher = SvMatcher()
matches = matcher.match(db_defs, events)                 # fuzzy match
best    = select_best_per_vcf(matches, tie_tol=1e-9)
```

### VCF side — what `SvReader` expects (`src/rbceq2/core_logic/large_variants.py:678-788`)
Caller-agnostic; keys off standard fields:
- `SVTYPE` in INFO → `DEL`, `DUP`, `INS`, `INV`, `CNV`, `BND`
- `SVLEN` (signed), `END` — optional
- `CIPOS` / `CIEND` — optional, **widen matching tolerance**
- Symbolic ALTs `<DEL>` etc. parsed **only if `SVTYPE` is absent** (line 739)
- Explicit-sequence large indels auto-detected when `abs(len(ALT)-len(REF)) >= min_size`
- `BND` paired via `MATEID`
- `chr` prefix stripped automatically; zygosity read from FORMAT/sample columns

Docstring lists supported callers: long-read (Sniffles2, SVIM, CuteSV, NanoVar, pbsv),
short-read (Manta, LUMPY, DELLY, GRIDSS, CNVnator, GATK-gCNV, DRAGEN).

Per-caller encoders live in `src/rbceq2/IO/encoders.py` (`CnvnatorEncoder`,
`GatkGcnvEncoder`, `DragenEncoder`). **Important:** these only build the display
`variant` string — they do **not** change the `SVTYPE` used for matching. `CN` copy
number is used only as a caller signature, never as dosage; every event collapses to a
`chrom:pos_svtype_size` token.

### Matching is fuzzy (`SvMatcher`, lines 159-367)
- `require_same_type=True` by default → DB token type must equal event `SVTYPE`
  (DEL/DUP/INS/INV/CNV/INDEL; INDEL is the only cross-compatible bucket).
- Adaptive positional tolerance (~50–75% of size, cap 50 kb) and a hard length gate
  (~35–50% of size, cap 100 kb), plus reciprocal ≥10% interval overlap.
- Defaults: `tol_pos=25_000`, `tol_len=10_000`, `tol_ratio=0.25`.
- So imprecise breakpoints still match — but **type must agree** (see DRAGEN CNV gotcha §5).

### Region pre-filter
`read_vcf` keeps only records within **±500 kb** of a DB variant position
(`src/rbceq2/IO/vcf.py:497`, `build_intervals`, `flank=500_000`).

### Relevant CLI options
- `--min_size` (default `10`): minimum indel/SV size for fuzzy matching.
- `--RH`: RHD/RHCE results — **"WARNING! Long read only!"**
- `--phased`: use phase info (default off).
- `--no_filter`: use all variants, not just `FILTER=PASS`.

---

## 2. Database representation of structural alleles

DB file: `src/rbceq2/resources/db.tsv` (2,029 rows, tab-separated). Columns:
`Chrom, Genotype, Genotype_alt, Coding, Protein, Phenotype, Phenotype_change,
Phenotype_alt, Phenotype_alt_change, Reference_genotype, Antithetical,
Weight_of_genotype, GRCh37, GRCh38, Transcript, Sub_type, Lane, Note`.

The matchable variant lives in the **`GRCh37` / `GRCh38`** columns as comma-separated
tokens. Three encodings (parsed by `parse_db_token`, `large_variants.py:91-129`):

1. **SNV / small indel:** `<pos>_<REF>_<ALT>` (e.g. `144920596_G_A`), or `<pos>_ref`
   (anchor "reference base here").
2. **Word-form SV:** `<pos>_<type>_<len>`, type ∈ {del,dup,ins,inv,cnv}, len unitized
   (`kb`/`bp`/`mb`). **Key CNV/SV encoding.**
3. **Sequence-form large indel:** `<pos>_<REFseq>_<ALTseq>` with literal DNA
   (used where the whole junction is spelled out, e.g. GYP*401/402).

Notation is **inconsistent between loci** but the parser handles both:
- **GYP (MNS):** lowercase, unit-suffixed → `144916340_del_110kb`, `145041698_dup_20kb`
- **RH:** UPPERCASE, bare bp → `25599038_DEL_59419`, `25696957_INS_38403`

**Hybrids are encoded as *paired* tokens, not a "hybrid" type:**
- RH hybrid = RHD `DEL` + RHCE-derived `INS`/dup
  (e.g. `RHD*01N.03` → `25611061_DEL_37390,25696957_INS_38403`)
- GYP hybrid = `del` + discriminating SNVs
  (e.g. `GYP*201.02` → `144918724_A_G,144918730_del_122kb`)

Example DB rows (Genotype | Coding | GRCh38 token):
```
GYPB*05N.01  Whole GYPB deletion            143914828_del_110kb
GYPB*05N.05  Del GYPB and GYPE (unofficial) 143832512_del_224kb
GYPA*28N.01  Del GYPA exons 2-7 GYPB exon 1 144019250_del_101kb
XK*N.01      Deletion of XK gene            37680879_del_53kb
RHD*01N.03   c.149_1227del;c.149_1227dup    25284570_DEL_37390,25370466_INS_38403
```

> Caveat: the token `svtype` is an **encoding artifact**, not the biology — a coding
> `dup` is stored as an `INS` token (e.g. GE "Triplicated Exon 3" `c.107_190dup` →
> `INS~84bp`). Classify by the `Coding` column, not the token type.
>
> The DB is still being refined — several rows carry author TODOs
> (CTL2 "Discuss ways to handle SVs", `rm,...` GYP rows, some RHD entries).

---

## 3. Annotated list of DB structural variants (CNV vs SV)

**All entries are SVs; CNV is the dosage-changing subset (DEL/DUP).** Token-type totals
in the GRCh38 column: **132 DEL, 41 INS, 1 DUP** — no `INV`, no `BND`, no literal `<CNV>`.
Practical buckets:
- **CNV** — dosage change ≥~1 kb (gene/exon loss/gain); needs a CNV/SV caller.
- **Hybrid/complex SV** — paired DEL+INS gene-conversion products (RH & GYP); needs long-read.
- **Large indel** — sub-kb (~10 bp–~1 kb); often from a normal indel caller.

### True large CNVs (whole-gene / multi-exon del & dup)
| Allele | Token(s) |
|---|---|
| `ATP11C*01N.01` | DEL ~219 kb (whole ATP11C) |
| `GYPB*05N.05` | DEL ~224 kb (GYPB+GYPE, unofficial) |
| `GYP*01N` | DEL ~124 kb (GYPA ex2-7 + GYPB ex1-5) |
| `GYPB*01N` / `GYPA*01N` | DEL ~121 kb / ~119 kb |
| `GYPB*05N.01` / `.02` | DEL ~110 kb / ~103 kb (whole GYPB) |
| `GYPA*28N.01` | DEL ~101 kb |
| `GYPB*05N.04` | DEL ~96 kb (GYPB ex2-5 + GYPE ex1) |
| `XG*01N.03` / `.02` | DEL ~114 kb / ~32 kb |
| `ABCC4*01N.01` | DEL ~68 kb (gene del) |
| `XK*N.01` + 20× `XK*N.01.00x/.10x` | DEL ~53 kb (all identical token) |
| `GCNT2*01N.06` | DEL ~41 kb |
| `CTL2*01N.01` | DEL ~37 kb |
| `A4GALT*N.01/.02/.03` | DEL ~33/21/26 kb |
| `RHAG*01N.15` | DEL ~32 kb |
| `LU*02N.06` | DEL ~27 kb |
| `ABCC1*01N.01` | DEL ~21 kb |
| `CD99*01N.02` / `.01` | DEL ~20 kb / ~11 kb |
| `GYPB*05N.03` | DEL ~19 kb |
| `RHD*01W.150` | INS ~12 kb (`c.327_487-4164dup`, dup→INS) |
| `MAM*01N.04/.05` | DEL ~8.5 kb |
| `XK*N.05` | DEL ~8 kb |
| `PIGG*01N.05/.06/.07` | DEL ~4.0/6.0/6.3 kb |
| `GE*01.-02.*` / `GE*01.-03.*` | DEL ~3.6 kb (7 alleles) |
| `XK*N.03` | DEL ~3.5 kb |
| `GE*01N.01` | DEL ~2.3 kb |
| `ABCG2*01N.27` | DEL ~1.8 kb |
| `LU*02N.02` | DEL ~1.06 kb |

### Complex / hybrid SVs (paired DEL + INS — RH & GYP)
| Allele | Tokens |
|---|---|
| `RHD*01N.02` | DEL ~49 kb + INS ~50 kb |
| `RHD*01N.03` | DEL ~37 kb + INS ~38 kb |
| `RHD*01N.04` / `RHCE*02N.08` | DEL ~31 kb + INS ~32 kb |
| `RHD*01N.42` | DEL ~25 kb + INS ~23 kb + INS ~226 bp + DEL ~148 bp |
| `RHD*01N.05` | INS ~23 kb + DEL ~22 kb |
| `RHD*01EL.44` / `RHD*03N.02` / `RHCE*01.29` | DEL ~21 kb + INS ~21.7 kb |
| `RHD*01N.43` / `RHCE*03.02` | DEL ~18.2 kb + INS ~18.3 kb |
| `RHD*01N.06/.07` / `RHD*03N.01` / `RHCE*01.34` | INS ~6.4 kb + DEL ~5.8 kb |
| `RHD*01EL.23` | INS ~5.2 kb + DEL ~5.2 kb |
| `RHCE*02.08.02` | DEL ~26.6 kb + INS ~28.6 kb |
| `RHCE*02N.07` | DEL ~27.4 kb + INS ~26.4 kb |
| `RHCE*01.44` | INS ~1.9 kb + DEL ~167 bp |
| `GYP*505` | DEL ~18 kb + **DUP ~20 kb** (only literal `dup` token in DB) |
| `GYP*201.01/.02`, `GYP*203.01` | DEL ~122 kb (+ SNVs) — A-B hybrids |
| `GYP*401` / `GYP*402` | DEL ~9/8 kb + INS ~2.9 kb (junction spelled out) |
| `GYP*301.01/.02` | DEL ~20 bp + INS ~1.79 kb |
| `GYP*501` | INS 25 bp + INS 28 bp + DEL ~3.6 kb |
| `GYP*503` | INS ~19 bp + DEL ~1.75 kb |
| `GYP*504` | INS 73 bp ×2 + DEL ~3.6 kb |
| `GYP*101.01/.02/.03` | INS ~39 bp + DEL ~2.8 kb (+SNV) |

### Large indels (<~1 kb — usually caught by an indel caller)
ABO*O.16 (DEL 725bp), XK*N.06 (391), CO*N.01 (384), XK*N.04 (263), XK*N.02 (245),
RHD*01N.67 (148), GYPB*04N.05/.06 (~97bp), GE*01.06.01/.02 (INS 84bp = dup),
PIGG*01N.08 (69), RHD*08N.01 (INS 37), FUT1*01N.29 (33), A4GALT*01N.19/*0XN.32 (26),
A4GALT*02N.25 (INS 23), KLF1*BGM51 (INS 23 = dup), RHD*01N.37 (23), RHD*01N.44 (22),
KLF*BGM71 (19), A4GALT*0XN.06 (17), FUT1*01N.33 (17), VEL*01N.01 (17),
A4GALT*01N.35 (INS 16 = dup), RHD*01N.28 (16), JK*01N.12 (15), RHAG*01N.17 (15),
FY*01N.02 (14), XK*N.14 (14), RHD*01N.30/.47 (13), ABCB6*01N.21 (12), RHD*01N.41 (12),
XK*N.38 (INS 12), ABCB6*01N.04 (INS 11 = dup), RHD*01N.36 (10), LW*05N.01 (10),
JK*02N.22 (10), GYPB*04N.04 (10), ABCB6*01N.07 (INS 9).

**Counts:** ~40 true large CNVs · ~22 hybrid/complex SVs · ~38 large indels.
Notes: the 21 `XK*N.01.*` cytogenetic dels share ONE `37540132_del_53kb` token
(indistinguishable by breakpoint); `GYP*505` holds the only literal `dup` token.

---

## 4. Do these need to be phased?

**No — phasing is opt-in (`--phased`, default off) and never required for detection/matching.**
It only affects *cis/trans* allele pairing in compound cases (hybrids, del+SNV combos).

Key behaviour (`src/rbceq2/core_logic/data_procesing.py:680-813`,
`modify_phase_of_large_indel`):
- SV callers rarely phase; RBCeq2 does **not** need the SV record phased. It **infers**
  the deletion's phase from overlapping **phased SNVs**: a het SNV called *inside* a
  hemizygous deletion must sit on the retained haplotype, so the deletion is the
  opposite phase (`flip_phase`, line 809). `_ref` anchors inferred likewise.
- Conditions to fire: `--phased` on; phased SNVs (`|` + `PS`) overlap the deletion;
  **all variants share one phase set** (`if len(all_phase_sets) != 1: return bg`,
  line 784); deletion currently unphased (`/`).
- In practice needs **long-read phased SNVs (HP/PS)** across the locus, single phase set.

Bottom line: leave `--phased` off and everything still calls; turn it on (with phased
long-read SNVs) to resolve cis/trans for hybrids and del+SNV alleles.

---

## 5. DRAGEN 3.7.8 (ICA) compatibility

Grounding files inspected (a single de-identified test sample, `SAMPLE_A`; real
paths/IDs redacted — see internal notes):
- SNV: `<bucket>/ica/dragen_3_7_8/output/recal_gvcf/SAMPLE_A.hard-filtered.recal.gvcf.gz`
- CNV: `<bucket>/ica/dragen_3_7_8/output/dragen_metrics/SAMPLE_A/SAMPLE_A.cnv.vcf.gz`
- SV (sibling): `.../SAMPLE_A.sv.vcf.gz`

DRAGEN version: `SW: 13.021.604.3.7.8f`. Genome: **hg38** (`chr`-prefixed).
Run flags of note: `--enable-sv true --enable-cnv true`, `--vc-enable-vcf-output false`
(so only a gVCF exists for SNVs).

> **Do NOT use `##referenceSexKaryotype` for per-sample sex.** It is a reference/config
> constant — it reports the karyotype the DRAGEN reference (`hg38_alt_masked_graph_v2`,
> which carries both chrX and chrY) was built to model, not the individual. It reads
> `XXYY` for **every** sample: verified `XXYY` on 102/102 gVCFs sampled across the whole
> OurDNA cohort range. Per-sample sex comes from
> `ploidy_estimation_metrics.csv` / `.ploidy.vcf.gz` (X/Y median-coverage ratios). The
> test sample used here is in fact **XX** (X median cov ≈ autosomal, Y median cov = 0).

### 5a. SNV file = gVCF (needs conversion to a plain sites VCF)
```
##fileformat=VCFv4.2
##ALT=<ID=NON_REF,...>
##source=DRAGEN_SNV
chr1  9997   .  N  <NON_REF>  .  PASS  END=10017  GT:AD:DP:GQ:MIN_DP:PL:SPL:ICNT  ./.:11,0:11:0:2:...
chr1  10018  .  C  <NON_REF>  .  PASS  END=10018  ...  0/0:17,0:17:15:17:...
```
Called variant records (already carrying GT) interleaved with `<NON_REF>` reference blocks
(the head above shows only ref blocks). **Must convert gVCF → sites VCF** — split
multiallelics, drop the `<NON_REF>` symbolic allele and reference-only blocks, region-restrict.
This is *preprocessing, not genotyping*: the per-sample genotypes already exist in the gVCF.
bcftools alone suffices (as the current `FilterAndConvertGvcfsForRbceq2` stage does); no
reference FASTA and no GATK `GenotypeGVCFs` (that is *joint* genotyping across samples).

### 5b. SV VCF (Manta-derived) — directly compatible ✅
```
##source=DRAGEN 13.021.604.3.7.8f
##ALT=<ID=DEL,...> <ID=INS,...> <ID=DUP:TANDEM,...>
SVTYPE distribution:  4983 DEL | 5774 INS | 46 DUP | 2640 BND
chr1  839442  MantaDEL:...  CACC...ACA  CT      463  PASS  END=839499;SVTYPE=DEL;SVLEN=-57;CIGAR=1M1I57D
chr1  998743  MantaDUP:TANDEM:...  T   <INS>   369  PASS  END=998743;SVTYPE=INS;SVLEN=52;DUPSVLEN=42;...
```
Proper `SVTYPE=DEL/DUP/INS/BND`, CIPOS/CIEND, MATEID. DUP:TANDEM reported as `<INS>`
(Manta convention) — fine, matches DB INS/dup tokens. **This is exactly what SvReader
was built for.**

### 5c. CNV VCF — needs preprocessing ⚠️
```
##ALT=<ID=CNV,..> <ID=DEL,..> <ID=DUP,..>
##INFO=<ID=REFLEN,..> <ID=SVLEN,..> <ID=SVTYPE,..> <ID=END,..> <ID=CIPOS,..> <ID=CIEND,..>
##FILTER=<ID=cnvLength,Description="CNV with length below 10000">
##FORMAT=<ID=GT> <ID=SM> <ID=CN> <ID=BC> <ID=PE>
SVTYPE distribution:  727 CNV   (ALL records)      FILTER: 975 PASS | 560 cnvLength | ...
chr1  2650427   DRAGEN:LOSS:...  N  <DEL>  51  cnvLength  SVLEN=-2648;SVTYPE=CNV;END=2653075;REFLEN=2648   1/1:0.058:0:2:...
chr1  13224579  DRAGEN:GAIN:...  N  <DUP>  47  PASS       SVLEN=17455;SVTYPE=CNV;END=13242034;REFLEN=17455 ./1:1.66:3:7:...
chr1  817861    DRAGEN:REF:...   N  .      81  PASS       END=2650427;REFLEN=1832567                       ./.:0.98:2:1434:...
```
Issues:
1. **Every record is `SVTYPE=CNV`** — direction only in the `<DEL>`/`<DUP>` ALT and
   `CN`/`SM`. RBCeq2's `SvMatcher.compatible()` requires DB-type == event-type; the DB
   has `DEL`/`DUP` tokens, so `"DEL" == "CNV"` → **False → large gene deletions never
   match.** The `<DEL>`/`<DUP>` fallback in `SvReader` fires only when `SVTYPE` is
   *absent* (line 739), which it isn't here. **The `DragenEncoder` does NOT rescue this**
   (it only builds the display string; matching uses raw `SVTYPE`).
   - Fix A (best): rewrite `SVTYPE=CNV`→`DEL`/`DUP` from ALT (or `CN<2`→DEL,`CN>2`→DUP).
   - Fix B: strip the `SVTYPE` INFO field so RBCeq2 reads the `<DEL>`/`<DUP>` ALT.
2. **`cnvLength` filter** drops CNVs <10 kb (612/727 here). RBCeq2 keeps only
   `FILTER=PASS` unless `--no_filter`. Source 2–50 kb events from the **SV VCF** instead,
   or run `--no_filter`.
3. `DRAGEN:REF:` records (ALT=`.`, no SVTYPE) are correctly skipped by `SvReader`.

### 5d. Test sample had no blood-group CNVs
SAMPLE_A shows no calls across GYP (chr4 ~144.8–145.0 Mb) and only small sub-2 kb
`cnvLength`-flagged noise near XK/RHD — a normal individual at these loci. Format
analysis (above) is the deliverable, not a positive hit.

---

## 6. Blockers & recommended workflow

**Blockers**
1. Three separate files → RBCeq2 takes **one VCF per sample** (or a folder of
   per-sample single-sample VCFs). Must merge.
2. SNV gVCF must be converted to a plain sites VCF (preprocessing, not genotyping).
3. CNV `SVTYPE=CNV` must be rewritten to `DEL`/`DUP`.
4. `cnvLength` filter vs default PASS-only behaviour.
5. The SV and CNV VCFs are not yet tracked in **metamist** — a pipeline stage cannot resolve
   them until they are registered under dedicated analysis types (one each), ideally upstream in
   the `dragen_align` pipeline that produces them.

**Per-sample pipeline**
1. **Convert** the gVCF → `SAMPLE.snv.vcf.gz` (drop `<NON_REF>`/ref blocks, split
   multiallelics; genotypes already present — no genotyping step).
2. **Fix CNV VCF:** drop `DRAGEN:REF:` records; rewrite `SVTYPE=CNV`→`DEL`/`DUP` per ALT;
   decide whether to keep sub-10 kb (`--no_filter`).
3. **Concatenate** SNV + SV + CNV → one sorted, bgzipped, tabixed VCF.
   (SV and CNV overlap in ~2–50 kb — dedup or let `select_best_per_vcf` choose.)
4. **Run:** `rbceq2 --vcf SAMPLE.merged.vcf.gz --out SAMPLE --reference_genome GRCh38`
   (add `--no_filter` if keeping sub-10 kb CNVs; omit `--RH` unless accepting it's
   unvalidated on short reads).

**Caveats**
- **RH / GYP hybrids unreliable on short-read DRAGEN.** RHD/RHCE and GYPA/B/E are
  near-identical paralogs → high MAPQ0, noisy/absent short-read SV & CNV calls. `--RH`
  is documented **long-read only**. Straightforward large deletions and small indels are
  fine; treat hybrid alleles as out of reach for this pipeline.
- **Sex-chromosome dosage:** for X-linked systems (XK/Kx, XG, CD99), the sample's sex
  karyotype affects expected copy number, and RBCeq2 infers zygosity from GT without
  knowing it. Determine per-sample sex from `ploidy_estimation_metrics.csv` /
  `.ploidy.vcf.gz` (X/Y coverage ratios) — **not** from the `##referenceSexKaryotype`
  header, which is a constant `XXYY` for every sample (see §5). Feed the resolved
  ploidy into the CNV direction logic rather than assuming diploid on X/Y.

---

## 7. Key source references
- `src/rbceq2/core_logic/large_variants.py` — `SvDef`, `parse_db_token`, `SvReader`,
  `SvMatcher`, `load_db_defs`, `select_best_per_vcf`
- `src/rbceq2/IO/encoders.py` — `DragenEncoder` (l.993), `MantaEncoder`,
  `CnvnatorEncoder`, `GatkGcnvEncoder`, `encode_sv_standard`, `normalize_svtype`
- `src/rbceq2/IO/vcf.py:497` — `build_intervals` (±500 kb region filter)
- `src/rbceq2/core_logic/data_procesing.py:680-813` — large-indel phase inference
- `src/rbceq2/main.py:272-305` — `find_hits` wiring; CLI args ~l.148-159
- `src/rbceq2/resources/db.tsv` — allele DB (SV tokens in GRCh37/GRCh38 columns)

## 8. Open next step (not yet done)
Write + test a reusable preprocessing script (bcftools + `SVTYPE=CNV`→DEL/DUP fixer +
merge) end-to-end on a de-identified test sample.
