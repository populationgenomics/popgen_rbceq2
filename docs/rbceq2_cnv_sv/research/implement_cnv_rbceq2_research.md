# Implementing CNV/SV support with RBCeq2: research notes

**Date:** 2026-07-14
**Repo:** `RBCeq2` (v2.4.2), blood-group allele inference
**Goal:** understand how RBCeq2 consumes CNV/SV data, and assess compatibility with our ICA
DRAGEN 3.7.8 outputs (a separate SNV gVCF, SV VCF and CNV VCF per sample).
**Status:** a July 2026 snapshot, read against rbceq2 2.4.2. Line numbers and counts are from
that release. Three parts are superseded: the allele counts in §3 by
`resolvability_by_input_class.md` (2.4.4 database); the open decisions in §5c and §6
(`cnvLength` policy, SV/CNV dedup, metamist registration) by the SPEC; and the sex source in
§5 and §6 by SPEC §5.6.

---

## 1. How CNV/SV data is provided to RBCeq2

There is no separate CNV/SV input and no dedicated CLI flag. Structural records go in the
same VCF passed to `--vcf`, encoded with standard VCF SV fields. RBCeq2 reads them alongside
SNVs and indels and fuzzy-matches them to SV alleles in its database.

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

### VCF side: what `SvReader` expects (`src/rbceq2/core_logic/large_variants.py:678-788`)

It is caller-agnostic and keys off standard fields:

- `SVTYPE` in INFO: `DEL`, `DUP`, `INS`, `INV`, `CNV` or `BND`
- `SVLEN` (signed) and `END`, both optional
- `CIPOS` / `CIEND`, optional, which widen the matching tolerance
- symbolic ALTs such as `<DEL>`, parsed only if `SVTYPE` is absent (line 739)
- explicit-sequence large indels, detected when `abs(len(ALT)-len(REF)) >= min_size`
- `BND`, paired via `MATEID`
- the `chr` prefix, stripped automatically; zygosity is read from the FORMAT/sample columns

The docstring lists supported callers: long-read (Sniffles2, SVIM, CuteSV, NanoVar, pbsv) and
short-read (Manta, LUMPY, DELLY, GRIDSS, CNVnator, GATK-gCNV, DRAGEN).

Per-caller encoders live in `src/rbceq2/IO/encoders.py` (`CnvnatorEncoder`, `GatkGcnvEncoder`,
`DragenEncoder`). They only build the display `variant` string and do not change the `SVTYPE`
used for matching. `CN` is used only as a caller signature, never as dosage, and every event
collapses to a `chrom:pos_svtype_size` token.

### Matching is fuzzy (`SvMatcher`, lines 159-367)

- `require_same_type=True` by default, so the db token type must equal the event `SVTYPE`
  (DEL/DUP/INS/INV/CNV/INDEL; INDEL is the only cross-compatible bucket).
- Positional tolerance is adaptive (about 50–75% of size, capped at 50 kb), with a hard
  length gate (about 35–50% of size, capped at 100 kb) and at least 10% reciprocal interval
  overlap.
- Defaults: `tol_pos=25_000`, `tol_len=10_000`, `tol_ratio=0.25`.
- So imprecise breakpoints still match, but the type must agree (see the DRAGEN CNV problem in §5c).

### Region pre-filter

`read_vcf` keeps only records within ±500 kb of a db variant position
(`src/rbceq2/IO/vcf.py:497`, `build_intervals`, `flank=500_000`).

### Relevant CLI options

- `--min_size` (default `10`): minimum indel/SV size for fuzzy matching.
- `--RH`: RHD/RHCE results, marked "WARNING! Long read only!".
- `--phased`: use phase information (default off).
- `--no_filter`: use all variants, not just `FILTER=PASS`.

---

## 2. Database representation of structural alleles

DB file: `src/rbceq2/resources/db.tsv` (2,029 rows, tab-separated). Columns:
`Chrom, Genotype, Genotype_alt, Coding, Protein, Phenotype, Phenotype_change,
Phenotype_alt, Phenotype_alt_change, Reference_genotype, Antithetical,
Weight_of_genotype, GRCh37, GRCh38, Transcript, Sub_type, Lane, Note`.

The matchable variant lives in the `GRCh37` / `GRCh38` columns as comma-separated tokens, in
three encodings parsed by `parse_db_token` (`large_variants.py:91-129`):

1. **SNV or small indel:** `<pos>_<REF>_<ALT>` (for example `144920596_G_A`), or `<pos>_ref`,
   an anchor meaning "reference base here".
2. **Word-form SV:** `<pos>_<type>_<len>`, with type one of del, dup, ins, inv or cnv and a
   length with a unit (`kb`, `bp`, `mb`). This is the main CNV/SV encoding.
3. **Sequence-form large indel:** `<pos>_<REFseq>_<ALTseq>` with literal DNA, used where the
   whole junction is spelled out (for example GYP*401/402).

Notation differs between loci, but the parser handles both forms:

- GYP (MNS): lower case with a unit, for example `144916340_del_110kb` and `145041698_dup_20kb`.
- RH: upper case with bare base counts, for example `25599038_DEL_59419` and `25696957_INS_38403`.

Hybrids are encoded as paired tokens, not as a hybrid type:

- an RH hybrid is an RHD `DEL` plus an RHCE-derived `INS` or dup (for example `RHD*01N.03` is
  `25611061_DEL_37390,25696957_INS_38403`);
- a GYP hybrid is a `del` plus discriminating SNVs (for example `GYP*201.02` is
  `144918724_A_G,144918730_del_122kb`).

Example db rows (Genotype | Coding | GRCh38 token):

```
GYPB*05N.01  Whole GYPB deletion            143914828_del_110kb
GYPB*05N.05  Del GYPB and GYPE (unofficial) 143832512_del_224kb
GYPA*28N.01  Del GYPA exons 2-7 GYPB exon 1 144019250_del_101kb
XK*N.01      Deletion of XK gene            37680879_del_53kb
RHD*01N.03   c.149_1227del;c.149_1227dup    25284570_DEL_37390,25370466_INS_38403
```

Two caveats. The token's `svtype` is an encoding artefact, not the biology: a coding `dup` is
stored as an `INS` token (for example GE "Triplicated Exon 3", `c.107_190dup`, is `INS~84bp`),
so classify by the `Coding` column, not the token type. And the db is still being refined:
several rows carry author TODOs (CTL2 "Discuss ways to handle SVs", the `rm,...` GYP rows,
some RHD entries).

---

## 3. Annotated list of db structural variants (CNV versus SV)

These counts are from the 2.4.2 database; `resolvability_by_input_class.md` has the 2.4.4
figures.

All the entries are SVs; CNVs are the dosage-changing subset (DEL and DUP). Token-type totals
in the GRCh38 column are 132 DEL, 41 INS and 1 DUP, with no `INV`, no `BND` and no literal
`<CNV>`. In practice they fall into three buckets:

- **CNV:** a dosage change of about 1 kb or more (gene or exon loss or gain), which needs a
  CNV or SV caller.
- **Hybrid or complex SV:** paired DEL and INS gene-conversion products (RH and GYP), which
  need long reads.
- **Large indel:** under about 1 kb (about 10 bp to 1 kb), often called by an ordinary indel
  caller.

### True large CNVs (whole-gene or multi-exon deletions and duplications)

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
| `RHD*01W.150` | INS ~12 kb (`c.327_487-4164dup`, a dup stored as INS) |
| `MAM*01N.04/.05` | DEL ~8.5 kb |
| `XK*N.05` | DEL ~8 kb |
| `PIGG*01N.05/.06/.07` | DEL ~4.0/6.0/6.3 kb |
| `GE*01.-02.*` / `GE*01.-03.*` | DEL ~3.6 kb (7 alleles) |
| `XK*N.03` | DEL ~3.5 kb |
| `GE*01N.01` | DEL ~2.3 kb |
| `ABCG2*01N.27` | DEL ~1.8 kb |
| `LU*02N.02` | DEL ~1.06 kb |

### Complex or hybrid SVs (paired DEL and INS: RH and GYP)

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
| `GYP*505` | DEL ~18 kb + DUP ~20 kb (the only literal `dup` token in the db) |
| `GYP*201.01/.02`, `GYP*203.01` | DEL ~122 kb (+ SNVs); A-B hybrids |
| `GYP*401` / `GYP*402` | DEL ~9/8 kb + INS ~2.9 kb (junction spelled out) |
| `GYP*301.01/.02` | DEL ~20 bp + INS ~1.79 kb |
| `GYP*501` | INS 25 bp + INS 28 bp + DEL ~3.6 kb |
| `GYP*503` | INS ~19 bp + DEL ~1.75 kb |
| `GYP*504` | INS 73 bp ×2 + DEL ~3.6 kb |
| `GYP*101.01/.02/.03` | INS ~39 bp + DEL ~2.8 kb (+SNV) |

### Large indels (under about 1 kb, usually caught by an indel caller)

ABO*O.16 (DEL 725bp), XK*N.06 (391), CO*N.01 (384), XK*N.04 (263), XK*N.02 (245),
RHD*01N.67 (148), GYPB*04N.05/.06 (~97bp), GE*01.06.01/.02 (INS 84bp = dup),
PIGG*01N.08 (69), RHD*08N.01 (INS 37), FUT1*01N.29 (33), A4GALT*01N.19/*0XN.32 (26),
A4GALT*02N.25 (INS 23), KLF1*BGM51 (INS 23 = dup), RHD*01N.37 (23), RHD*01N.44 (22),
KLF*BGM71 (19), A4GALT*0XN.06 (17), FUT1*01N.33 (17), VEL*01N.01 (17),
A4GALT*01N.35 (INS 16 = dup), RHD*01N.28 (16), JK*01N.12 (15), RHAG*01N.17 (15),
FY*01N.02 (14), XK*N.14 (14), RHD*01N.30/.47 (13), ABCB6*01N.21 (12), RHD*01N.41 (12),
XK*N.38 (INS 12), ABCB6*01N.04 (INS 11 = dup), RHD*01N.36 (10), LW*05N.01 (10),
JK*02N.22 (10), GYPB*04N.04 (10), ABCB6*01N.07 (INS 9).

**Counts:** about 40 true large CNVs, 22 hybrid or complex SVs and 38 large indels. The 21
`XK*N.01.*` cytogenetic deletions share one `37540132_del_53kb` token, so they cannot be
told apart by breakpoint, and `GYP*505` holds the only literal `dup` token.

---

## 4. Do these need to be phased?

No. Phasing is opt-in (`--phased`, default off) and never required for detection or
matching. It only affects cis/trans allele pairing in compound cases (hybrids, and deletion
plus SNV combinations).

The behaviour is in `modify_phase_of_large_indel` (`src/rbceq2/core_logic/data_procesing.py:680-813`):

- SV callers rarely phase, and RBCeq2 does not need the SV record phased. It infers the
  deletion's phase from overlapping phased SNVs: a het SNV called inside a hemizygous
  deletion must sit on the retained haplotype, so the deletion takes the opposite phase
  (`flip_phase`, line 809). `_ref` anchors are inferred the same way.
- It fires only when `--phased` is on, phased SNVs (`|` plus `PS`) overlap the deletion, all
  variants share one phase set (`if len(all_phase_sets) != 1: return bg`, line 784), and the
  deletion is currently unphased (`/`).
- In practice it needs long-read phased SNVs (HP/PS) across the locus, in a single phase set.

So leave `--phased` off and everything still calls; turn it on, with phased long-read SNVs, to
resolve cis/trans for hybrids and deletion-plus-SNV alleles.

---

## 5. DRAGEN 3.7.8 (ICA) compatibility

Files inspected, from one test sample:

- SNV: `<bucket>/ica/dragen_3_7_8/output/recal_gvcf/<sg>.hard-filtered.recal.gvcf.gz`
- CNV: `<bucket>/ica/dragen_3_7_8/output/dragen_metrics/<sg>/<sg>.cnv.vcf.gz`
- SV (sibling): `.../<sg>.sv.vcf.gz`

DRAGEN version: `SW: 13.021.604.3.7.8f`. Genome: hg38 (`chr`-prefixed). Run flags of note:
`--enable-sv true --enable-cnv true` and `--vc-enable-vcf-output false`, so the only SNV
output is a gVCF.

Do not use `##referenceSexKaryotype` for per-sample sex. It is a reference constant: it
reports the karyotype the DRAGEN reference (`hg38_alt_masked_graph_v2`, which carries both
chrX and chrY) was built to model, not the individual. It reads `XXYY` for every sample,
verified on 102 of 102 gVCFs sampled across the OurDNA cohort range. Per-sample sex came, at
this point, from `ploidy_estimation_metrics.csv` or `.ploidy.vcf.gz` (X/Y median-coverage
ratios); SPEC §5.6 now reads it from the QC meta in metamist.

### 5a. The SNV file is a gVCF, and needs converting to a plain sites VCF

```
##fileformat=VCFv4.2
##ALT=<ID=NON_REF,...>
##source=DRAGEN_SNV
chr1  9997   .  N  <NON_REF>  .  PASS  END=10017  GT:AD:DP:GQ:MIN_DP:PL:SPL:ICNT  ./.:11,0:11:0:2:...
chr1  10018  .  C  <NON_REF>  .  PASS  END=10018  ...  0/0:17,0:17:15:17:...
```

Called variant records, which already carry GT, are interleaved with `<NON_REF>` reference
blocks; the head above shows only reference blocks. The gVCF must be converted to a sites
VCF: split multiallelics, drop the `<NON_REF>` symbolic allele and the reference-only blocks,
and restrict to the regions. This is preprocessing, not genotyping, since the per-sample
genotypes already exist in the gVCF. bcftools alone suffices, as the current
`FilterAndConvertGvcfsForRbceq2` stage shows; no reference FASTA is needed, and GATK
`GenotypeGVCFs`, which genotypes jointly across samples, has no role.

### 5b. SV VCF (DRAGEN SV caller, which extends Manta): directly compatible

```
##source=DRAGEN 13.021.604.3.7.8f
##ALT=<ID=DEL,...> <ID=INS,...> <ID=DUP:TANDEM,...>
SVTYPE distribution:  4983 DEL | 5774 INS | 46 DUP | 2640 BND
chr1  839442  MantaDEL:...  CACC...ACA  CT      463  PASS  END=839499;SVTYPE=DEL;SVLEN=-57;CIGAR=1M1I57D
chr1  998743  MantaDUP:TANDEM:...  T   <INS>   369  PASS  END=998743;SVTYPE=INS;SVLEN=52;DUPSVLEN=42;...
```

Records carry a proper `SVTYPE=DEL/DUP/INS/BND`, `CIPOS`/`CIEND` and `MATEID`. `DUP:TANDEM` is
reported as `<INS>`, a Manta convention the DRAGEN caller keeps along with the `Manta` prefix
on record IDs, and it matches the db's INS and dup tokens. This is the input `SvReader` was
built for.

### 5c. CNV VCF: needs preprocessing

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

1. **Every record is `SVTYPE=CNV`.** Direction lives only in the `<DEL>`/`<DUP>` ALT and in
   `CN`/`SM`. RBCeq2's `SvMatcher.compatible()` requires the db type to equal the event type,
   and the db has `DEL` and `DUP` tokens, so `"DEL" == "CNV"` is false and large gene
   deletions never match. The ALT fallback in `SvReader` fires only when `SVTYPE` is absent
   (line 739), which it is not here. `DragenEncoder` does not rescue this: it only builds the
   display string, and matching uses the raw `SVTYPE`.
   - Fix A (preferred): rewrite `SVTYPE=CNV` to `DEL` or `DUP` from the ALT (or `CN < 2` to
     DEL, `CN > 2` to DUP).
   - Fix B: strip the `SVTYPE` INFO field, so RBCeq2 reads the `<DEL>`/`<DUP>` ALT.
2. **The `cnvLength` filter** drops CNVs under 10 kb (612 of 727 here), and RBCeq2 keeps only
   `FILTER=PASS` unless run with `--no_filter`. The options at this point were to source
   2–50 kb events from the SV VCF, or to run `--no_filter`. (Superseded: SPEC §5.4 triages
   these records and rewrites FILTER per record; `--no_filter` is never used.)
3. `DRAGEN:REF:` records (ALT `.`, no `SVTYPE`) are correctly skipped by `SvReader`.

### 5d. The test sample had no blood-group CNVs

The sample showed no calls across GYP (chr4, about 144.8–145.0 Mb) and only small sub-2 kb
`cnvLength` noise near XK and RHD, as expected for a normal individual at these loci. The
format analysis above is the result, not a positive hit.

---

## 6. Blockers and recommended workflow (July 2026)

The SPEC resolves each of these; this is the July view.

**Blockers**

1. There are three separate files, and RBCeq2 takes one VCF per sample (or a folder of
   single-sample VCFs), so they must be merged.
2. The SNV gVCF must be converted to a plain sites VCF (preprocessing, not genotyping).
3. The CNV VCF's `SVTYPE=CNV` must be rewritten to `DEL`/`DUP`.
4. The `cnvLength` filter conflicts with RBCeq2's default PASS-only behaviour.
5. The SV and CNV VCFs are not tracked in metamist. (Superseded: the SPEC derives their paths
   from the registered gVCF and adds no analysis types, §5.1.)

**Per-sample pipeline**

1. **Convert** the gVCF to `SAMPLE.snv.vcf.gz`: drop `<NON_REF>` and reference blocks and
   split multiallelics. The genotypes are already present, so there is no genotyping step.
2. **Fix the CNV VCF:** drop `DRAGEN:REF:` records, rewrite `SVTYPE=CNV` to `DEL`/`DUP` from
   the ALT, and decide whether to keep sub-10 kb records.
3. **Concatenate** SNV, SV and CNV into one sorted, bgzipped, tabixed VCF. The SV and CNV VCFs
   overlap at about 2–50 kb, so either deduplicate or let `select_best_per_vcf` choose.
   (Superseded: `select_best_per_vcf` cannot deduplicate, and the SPEC's triage does, §5.5.)
4. **Run** `rbceq2 --vcf SAMPLE.merged.vcf.gz --out SAMPLE --reference_genome GRCh38`,
   adding `--no_filter` if keeping sub-10 kb CNVs, and omitting `--RH`, which is unvalidated
   on short reads.

**Caveats**

- **RH and GYP hybrids are unreliable on short-read DRAGEN.** RHD/RHCE and GYPA/B/E are
  near-identical paralogs, so reads map with high MAPQ0 and short-read SV and CNV calls are
  noisy or absent. `--RH` is documented as long-read only. Large deletions and small indels
  are in reach; hybrid alleles are not, for this pipeline.
- **Sex-chromosome dosage.** For the X-linked systems (XK/Kx, XG, CD99), the sample's
  karyotype affects expected copy number, and RBCeq2 infers zygosity from GT without knowing
  it. The resolved ploidy should feed the CNV direction logic rather than a diploid
  assumption on chrX/chrY, and never come from the `##referenceSexKaryotype` header (§5).

---

## 7. Key source references

- `src/rbceq2/core_logic/large_variants.py`: `SvDef`, `parse_db_token`, `SvReader`,
  `SvMatcher`, `load_db_defs`, `select_best_per_vcf`
- `src/rbceq2/IO/encoders.py`: `DragenEncoder` (l.993), `MantaEncoder`, `CnvnatorEncoder`,
  `GatkGcnvEncoder`, `encode_sv_standard`, `normalize_svtype`
- `src/rbceq2/IO/vcf.py:497`: `build_intervals` (±500 kb region filter)
- `src/rbceq2/core_logic/data_procesing.py:680-813`: large-indel phase inference
- `src/rbceq2/main.py:272-305`: `find_hits` wiring; CLI args at about l.148-159
- `src/rbceq2/resources/db.tsv`: the allele db (SV tokens in the GRCh37/GRCh38 columns)

## 8. Next step recorded in July

Write and test a reusable preprocessing script (bcftools, an `SVTYPE=CNV` to DEL/DUP fixer,
and a merge) end to end on a test sample. The SPEC now carries this as a design.
