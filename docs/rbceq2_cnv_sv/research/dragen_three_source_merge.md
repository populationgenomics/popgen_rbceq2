# DRAGEN to RBCeq2: three data sources per sample, and how to combine them

**Date:** 2026-07-16; diagram revised 2026-09-18 to match SPEC §5.5 and the first review round. The source tables are unchanged from July apart from the caller name.
**Reader:** anyone who wants the shape of the design before reading the SPEC.

ICA DRAGEN 3.7.8 (`SW: 13.021.604.3.7.8f`, hg38, `chr`-prefixed) emits three variant files per sample. RBCeq2 (v2.4.x) reads SNVs, indels and structural variants from one VCF per sample, with no separate CNV/SV input or CLI flag. So the three DRAGEN files must be transformed and merged into one sorted, bgzipped, tabixed VCF before RBCeq2 runs.

This note shows each source as it really is, with records copied from one of the OurDNA 1KG control replicates: independent sequencings of one 1000 Genomes control sample, used to check that merged calls agree across replicates. It then shows what has to change and how the pieces combine. [`implement_cnv_rbceq2_research.md`](./implement_cnv_rbceq2_research.md) has the RBCeq2 source analysis this rests on.

The SV VCF is written by the DRAGEN 3.7.8 SV caller, which integrates and extends Manta; record IDs keep Manta's prefix. This note says "SV caller", not "Manta".

---

## The three sources (real layout, `gs://cpg-ourdna-main/ica/dragen_3_7_8/output/`)

| # | Source | Path (per SG `<sg>`) | Caller | RBCeq2-ready? |
|---|--------|----------------------|--------|---------------|
| 1 | SNV gVCF | `recal_gvcf/<sg>.hard-filtered.recal.gvcf.gz` | DRAGEN SNV | Needs converting to a sites VCF, which the existing stage does |
| 2 | SV VCF | `dragen_metrics/<sg>/<sg>.sv.vcf.gz` | DRAGEN SV caller (extends Manta) | Directly compatible |
| 3 | CNV VCF | `dragen_metrics/<sg>/<sg>.cnv.vcf.gz` | DRAGEN bin CNV (genomes); DRAGEN panel-of-normals target CNV (exomes) | Needs an `SVTYPE` rewrite; absent when the run did not call CNVs |

Sex is read from `sg.meta['qc']`, which `single_sample_qc_popgen` already populates (SPEC §5.6):

- `ploidy_estimation`, `norm_x_coverage` and `norm_y_coverage`: DRAGEN's karyotype estimate and X and Y over autosomal median coverage, from `ploidy_estimation_metrics.csv`. The file itself is read only when the meta is absent.
- The somalier signals `f_stat_raw`, `x_het_rate`, `y_calls` and `y_n`: the fallback when DRAGEN omitted the ploidy file, as it does for genomes with uneven coverage. The atlas repo's `karyotype_from_signals` turns them into a karyotype.
- `qc_checks_failed`, which records the ploidy-versus-reported-sex check.

Do not read sex from the `##referenceSexKaryotype` header: it is a reference constant that reads `XXYY` for every sample.

**Exomes** (DRAGEN 3.7.8, checked 2026-09-18). Every exome has the gVCF, the SV VCF and the ploidy file. The 11,917 production exomes also have a CNV VCF, called per capture target against a panel of 100 normals: events only, no `DRAGEN:REF:` records, no `cnvLength` filter, and `BC` counting targets. The 10 exomes of the test run have no CNV VCF; their `cnv_metrics.csv` stops after the target-interval count. So whether a CNV VCF is expected is read from `cnv_metrics.csv`, not from the sequencing type (SPEC §3, §5.1).

**Paths.** The `gs://…` paths above are shown for grounding. The stage derives them from the registered gVCF's DRAGEN output prefix, fails an SG missing an expected file, and records the paths in the merged VCF's Analysis meta. There are no new metamist analysis types: `dragen_align` cannot register the files before its Nextflow refactor, and cpg_flow is not to be extended (SPEC §5.1).

---

## 1. SNV gVCF: convert to a plain sites VCF

Real records:

```text
chr1  9997   .  N  <NON_REF>  .  PASS  END=10017  GT:AD:DP:GQ:MIN_DP:PL:SPL:ICNT  ./.:11,0:11:0:2:...
chr1  10018  .  C  <NON_REF>  .  PASS  END=10018  ...                            0/0:17,0:17:15:17:...
```

**Problem:** it is a gVCF. Called variant records, which already carry the per-sample
genotype (for example `0/1`), are interleaved with `<NON_REF>` reference blocks; the head
above happens to show only reference blocks. RBCeq2 fetches blood-group sites by coordinate
and expects a plain VCF, and the reference blocks and the symbolic `<NON_REF>` allele break it.

**Modification:** convert the gVCF to a sites VCF: split multiallelics, drop `<NON_REF>` and
the reference blocks, and restrict to the blood-group regions. This is preprocessing, not
genotyping. The genotypes are already in the gVCF, and GATK `GenotypeGVCFs`, which genotypes
jointly across samples, has no role. It needs bcftools only, and no reference FASTA.

This is already implemented in the `FilterAndConvertGvcfsForRbceq2` stage
(`src/popgen_rbceq2/stages/blood_group_genotyping/filter_and_convert.py`), with
`bcftools norm -m -any` and a region restrict to `resources/bg_regions.GRCh38.bed`. The SNV
path is done; the gap is sources 2 and 3.

---

## 2. SV VCF (DRAGEN SV caller): directly compatible

Real records:

```text
chr1  789481  MantaINS:...  G  <INS>  999  PASS  END=789481;SVTYPE=INS;CIPOS=0,7;...
chr1  839442  MantaDEL:...  CACC…ACA  CT  463  PASS  END=839499;SVTYPE=DEL;SVLEN=-57;CIGAR=1M1I57D
chr1  934064  MantaDEL:...  AGGG…    A   412  PASS  END=934904;SVTYPE=DEL;SVLEN=-840;CIPOS=0,27;HOMLEN=27
```

Header ALTs are `<DEL>`, `<INS>` and `<DUP:TANDEM>`. Records carry a proper
`SVTYPE=DEL/DUP/INS/BND`, plus `SVLEN`, `END`, `CIPOS`/`CIEND` and `MATEID`. `DUP:TANDEM` is
reported as `<INS>`, a Manta convention the DRAGEN caller keeps, which matches the db's INS
and dup tokens. FORMAT is `GT:FT:GQ:PL:PR:SR` (or without `SR`) on every record.

**Modification:** none for matching; this is the input RBCeq2's `SvReader` was built for.
Only housekeeping: drop `BND`, restrict to the blood-group regions (a size optimisation;
RBCeq2 also does this internally), make the sample column name match the merged VCF, and
sort, bgzip and tabix.

---

## 3. CNV VCF: needs an SVTYPE rewrite

Real records:

```text
chr1  817861    DRAGEN:REF:chr1:817861-2650427   N  .      81  PASS       END=2650427;REFLEN=1832567                        ./.:0.98:2:1434:...
chr1  2650427   DRAGEN:LOSS:chr1:2650428-2653075 N  <DEL>  51  cnvLength  SVLEN=-2648;SVTYPE=CNV;END=2653075;REFLEN=2648    1/1:0.058:0:2:...
chr1  3501568   DRAGEN:GAIN:chr1:3501569-3502568 N  <DUP>  88  cnvLength  SVLEN=1000;SVTYPE=CNV;END=3502568;REFLEN=1000     ./1:3.18:6:1:...
```

Header ALTs `<CNV>`/`<DEL>`/`<DUP>`; FORMAT is `GT:SM:CN:BC:PE` on every record, with `CN`
the estimated copy number and `BC` the bin count. Every PASS `<DUP>` in 150 genomes carried
GT `./1`; every PASS `<DEL>` carried `0/1` or `1/1`.

**Three problems, three fixes:**

1. **Every record is `SVTYPE=CNV`.** Direction lives only in the `<DEL>`/`<DUP>` ALT and the
   `CN`/`SM` values. RBCeq2's `SvMatcher.compatible()` requires the db token type to equal the
   event `SVTYPE`, and the db has `DEL` and `DUP` tokens, so `"DEL" == "CNV"` is false and
   large gene deletions never match. The ALT fallback in `SvReader` fires only when `SVTYPE`
   is absent (`large_variants.py:739` in 2.4.2), and `DragenEncoder` does not rescue this: it
   only builds the display string.
   - **Fix:** rewrite `SVTYPE=CNV` to `DEL` or `DUP` from the ALT, or on sex chromosomes from
     `CN` against the sample's expected ploidy. Stripping `SVTYPE` so that `SvReader` reads
     the ALT would also work, but the explicit rewrite is clearer.
2. **`DRAGEN:REF:` records** (ALT `.`, no `SVTYPE`) are non-events. Drop them; RBCeq2 would
   skip them anyway.
3. **The `cnvLength` filter** flags CNVs under 10 kb, the majority of records, and RBCeq2
   keeps only `FILTER=PASS` unless run with `--no_filter`.
   - **Fix (SPEC §5.4, Q1 resolved):** keep a `cnvLength` record that survives triage and has
     at least 3 bins, and rewrite its FILTER to PASS with the original in `INFO/SVFILTER`.
     `--no_filter` is never used, since it is global.

---

## How they combine

```mermaid
flowchart TD
    subgraph DRAGEN["DRAGEN 3.7.8 output (per sample)"]
        GVCF["1 · SNV gVCF<br/>recal_gvcf/&lt;sg&gt;.hard-filtered.recal.gvcf.gz<br/><i>&lt;NON_REF&gt; ref-blocks</i>"]
        SV["2 · SV VCF (DRAGEN SV caller)<br/>dragen_metrics/&lt;sg&gt;/&lt;sg&gt;.sv.vcf.gz<br/><i>SVTYPE=DEL/DUP/INS/BND · exact breakpoints</i>"]
        CNV["3 · CNV VCF (when cnv_metrics.csv shows the caller ran)<br/>dragen_metrics/&lt;sg&gt;/&lt;sg&gt;.cnv.vcf.gz<br/><i>every record SVTYPE=CNV · genomes: 1–2 kb bins, &lt;10 kb = cnvLength · exomes: capture targets vs panel of normals, events only</i>"]
        PLOIDY["ploidy_estimation_metrics.csv<br/><i>Ploidy estimation: XX / XY / other</i><br/>(read directly only if sg.meta['qc'] is absent)"]
    end
    SEX["sg.meta['qc'] from single_sample_qc_popgen<br/>ploidy_estimation · norm_x/y_coverage<br/>somalier f_stat_raw, x_het_rate, y_calls → karyotype_from_signals<br/>qc_checks_failed (ploidy vs reported sex)"]
    SEX -.->|"primary sex source<br/>+ somalier fallback"| CNVv
    SEX -.->|"recorded reported-sex<br/>failure → SEXCHECK"| QC

    GVCF -->|"convert → sites VCF<br/>split multiallelics, drop &lt;NON_REF&gt;/ref-blocks,<br/>restrict to bg_regions.GRCh38.bed<br/>(existing stage)"| SNVv["snv sites VCF"]
    SV -->|"drop BND; region-restrict;<br/>rename sample"| SVv["sv records"]
    CNV -->|"drop DRAGEN:REF records;<br/>rewrite SVTYPE=CNV → DEL/DUP from ALT;<br/>drop cnvLength with BC &lt; 3"| CNVv["cnv records"]
    PLOIDY -.->|"not XX/XY → drop or tag<br/>chrX/chrY cnv records<br/>(§5.6 karyotype gate)"| CNVv

    SVv --> TRIAGE
    CNVv --> TRIAGE
    TRIAGE["<b>Triage to database-relevant records (§5.5 rule 1)</b><br/>keep if (a) within SvMatcher tolerance of a db SV token,<br/>or (b) PASS deletion &lt;1 Mb spanning a defining SNV site;<br/>discard everything else<br/><i>150 genomes: 0–2 (a) + 0–3 (b) records per genome</i>"]

    TRIAGE --> ONE{"two records,<br/>same direction,<br/>recip overlap ≥ 0.5?"}
    ONE -->|"yes → keep the SV record<br/>(exact breakpoints), record<br/>partner in INFO/SVPARTNER"| TAG
    ONE -->|"no → keep as is"| TAG
    TAG["tag INFO/SVSRC = SV · CNV · CNV_LOWRES<br/>rewrite FILTER → PASS, original in INFO/SVFILTER<br/>assert no two records share CHROM/POS/END/SVTYPE<br/>assert GT is the first FORMAT key on every row"]

    SNVv --> MERGE["bcftools concat + sort<br/>→ one bgzipped, tabixed VCF"]
    TAG --> MERGE
    MERGE --> VCF["&lt;sg&gt;.rbceq2_input.vcf.gz"]

    VCF --> RBC["rbceq2 --vcf … --reference_genome GRCh38<br/>(never --no_filter, never --RH)"]
    RBC --> OUT["blood-group calls<br/>(geno / pheno TSVs)"]
    VCF -.->|"structural records +<br/>DRAGEN:REF tiling"| QC["FlagBloodGroupCallQc (SPEC §6)<br/>SVNOCOV · SVLOWRES · SVDEPTH · SVUNASSESSED<br/>SVDEL · KARYOTYPE · SEXCHECK"]
```

**Exomes** take the same path. With a production CNV VCF, the §5.4 step has no REF records to
drop, and the QC reads whether the region was assessed from the capture design rather than
from `DRAGEN:REF:` tiling (SPEC §6). Where a run produced no CNV VCF, as the test run did,
triage runs over the SV records alone and every 10 kb+ target is `SVUNASSESSED`.

**Why one record per event.** rbceq2 keeps the best db token per record, so a deletion that
arrives from both callers, offset by the CNV caller's bin snapping, is read as two different
alleles, a compound heterozygote. Its 2.4.4 tie error needs identical coordinates and never
fires on real pairs. The figure and evidence are in SPEC §5.5 and `sv_cnv_overlap_50_samples.md`.

### Combine step notes
- **Merge:** `bcftools concat` of the normalised VCFs, then sort, bgzip and tabix. All inputs
  must share the sample column name and contig naming (`chr`-prefixed hg38); RBCeq2 strips
  the `chr` prefix internally.
- **FORMAT fields (review asked, 2026-09-18):** nothing to harmonise. FORMAT is per record,
  so gVCF rows keep `GT:AD:DP:GQ:…`, SV rows `GT:FT:GQ:PL:PR:SR` and CNV rows
  `GT:SM:CN:BC:PE`. RBCeq2 requires only `GT` first on every row (`IO/vcf.py`), true of all
  three, and reads the SV geometry from INFO. The tags the SV and CNV headers share have
  identical Number and Type in DRAGEN 3.7.8, and a real concat produced no header warning
  (SPEC §5.5).
- **SV ∩ CNV overlap:** decided, SPEC §5.5 rule 3. Dedup at merge time, keeping the SV
  record and recording the CNV partner; `select_best_per_vcf` cannot do it because it keeps
  one db token per record, not per locus.
- **Adjacent records (review asked):** not merged. Across 150 genomes, adjacent
  same-direction CNV records were copy-number steps (3 of 254 pairs shared a `CN`), never a
  target split in two, and the SV caller produced one such pair in 17,175 records.
  `sv_cnv_overlap_50_samples.md` §13.
- **Sex chromosomes:** for the X-linked systems (XK/Kx, XG, CD99) the expected copy number
  depends on the sample's karyotype and on PAR versus non-PAR. The CNV direction logic takes
  expected ploidy from the sex sources above (SPEC §5.6), never a diploid assumption on
  chrX/chrY.

### Caveats carried from the research
- **RH and GYP hybrids are unreliable on short-read DRAGEN.** RHD/RHCE and GYPA/B/E are
  near-identical paralogs, so reads map with MAPQ0 and SV and CNV calls are noisy or absent.
  `--RH` is documented as long-read only. Large deletions and small indels are in reach;
  hybrid alleles are out of scope for this pipeline.
- **Phasing** is opt-in (`--phased`, default off) and never required for detection, so it
  stays off.
