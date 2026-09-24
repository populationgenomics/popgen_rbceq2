# Popgen_RBCeq2

## What Popgen_RBCeq2 is
A pipeline that produces an estimate of the blood type for each individual in the input cohort.
- __In__: Variant calls in gVCF format, for the sequencing groups of one or more Metamist cohorts.
- __Structure__: A [cpg-flow](https://github.com/populationgenomics/cpg-flow) workflow that converts each input into the correct format for RBCeq2, calls RBCeq2, flags the calls whose supporting sites were low quality, and concatenates the per-sample results per cohort.
- __Out__: Three TSV files per sequencing group, plus a QC TSV, and a combined copy of each per cohort:
  1. `<SGID>.geno.tsv` genotype calls,
  2. `<SGID>.pheno_alphanumeric.tsv` alphanumeric phenotype,
  3. `<SGID>.pheno_numeric.tsv` numeric phenotype,
  4. `<SGID>.qc.tsv` per-system quality flags, in the same one-column-per-system layout so it joins to the calls by column.

Each is registered in Metamist as an Analysis, with the calls or flags parsed into its `meta`.

## Why it exists

Popgen_RBCeq2 exists to provide a simple, automated way to estimate blood type from variant calls.
Blood donations are heavily relied on by healthcare institutions. Some blood groups are much rarer than
others among the general population. Lifeblood, the developer of RBCeq2, are interested in targeting
communities with high occurrences of rare blood groups to bolster blood donations. CPG has access
to blood samples from communities not commonly represented in genomic databases and are uniquely positioned
to provide insights into blood-group frequencies in these underrepresented populations to help direct
blood drive initiatives.

## Who is this for?
This pipeline is for PopGen team members at CPG to identify individual blood groups and
estimate frequency of blood groups across their datasets.

## Core Thesis
The value is a reproducible RBCeq2 wrapper that annotates CPG's underrepresented cohort gVCFs with estimated blood groups.

## Load-bearing design principles

**gVCF -> VCF conversion (`FilterAndConvertGvcfsForRbceq2`)**
- RBCeq2 cannot read gVCFs: its variant encoder rejects any ALT ending in `>`, so the `<NON_REF>` symbolic allele present on every gVCF record (both ref blocks and real variants like `G,<NON_REF>`) falls through to a fallback that mis-keys the variant and never matches the database. RBCeq2 does not error, it silently reverts affected systems to reference (observed: ABO -> Undetermined on HG00096). The conversion therefore splits multiallelics (`bcftools norm -m -any`) and drops `<NON_REF>` (`view -e 'ALT="<NON_REF>"'`) so only clean biallelic records reach RBCeq2, then trims ALT alleles no surviving genotype carries.
- The order matters: `norm -m -any` must precede the `<NON_REF>` exclusion. `ALT="<NON_REF>"` matches if *any* ALT matches, so filtering an unsplit gVCF would delete every variant in the file.

**Report low-quality sites; do not remove them**
- RBCeq2 has no depth or quality filtering of its own, and reads a blood-group site absent from its input as a *confident homozygous reference call*. Dropping a borderline genotype therefore does not produce a no-call — it manufactures a wild-type call at a site that defines an antigen.
- So genotypes are not filtered on DP or GQ. DRAGEN has already hard-filtered these gVCFs, keeping a failed record with its filter name, and rbceq2 excludes an allele whose defining variant is not `PASS`; post-hoc records are given `PASS` or DRAGEN's `LowDepth` in the merge so they are judged the same way, since HaplotypeCaller leaves FILTER `.` and rbceq2 would discard every recovered allele. `FlagBloodGroupCallQc` instead reports per-system flags (`LOWQ`, `DEL`, `NOCOV`, `NA`) naming the defining site and what the caller reported there, and the thresholds are recorded alongside the flags so a reader can tell what `LOWQ` meant on that run. See the [README](../README.md) for how to read a flag.
- This reverses an earlier decision to delete sub-threshold sites. Deletion was not neutral: it moved uncertain calls to false reference rather than removing them from analysis.

**Post-hoc recall of off-target sites, for exomes only (`PosthocGenotypeOffTargetSites`)**
- An exome gVCF is called against a capture-target BED with no padding, so DRAGEN stops emitting
  at the BED edges. ~165 of 1,599 non-HPA defining coordinates fall outside the capture per
  design, and they arrive as *no record at all* rather than a poor one — the exact input rbceq2
  reads as a confident homozygous reference call.
- The reads are usually there. ~105 of those sites sit within 100bp of a target edge with
  40-110x MAPQ>=20 depth in the CRAM, the FY GATA Duffy-null promoter sites among them. Only the
  caller stopped early, so we call those sites again from the CRAM with GATK HaplotypeCaller in
  `--dragen-mode`, streaming ~136kb of padded intervals rather than localising the CRAM. This
  takes assessable non-HPA coordinates from 91.7% to ~98% per design. Re-running DRAGEN over a
  cohort is not affordable, and sites >250bp off both designs are at ~0x and unrecoverable.
- **Silence is judged empirically; the fill is bounded by the capture design.** A post-hoc
  record is kept only where the DRAGEN gVCF has no record covering a defining site *and* that
  site lies outside the cohort's capture design. Silence is read per sample from the gVCF, never
  from metadata, so a design BED that misdescribes a sample's real footprint cannot overwrite a
  primary call or hide a hole — the design only narrows what may be filled. Bounding by the
  design keeps two findings apart that were being conflated: "the capture never targeted this
  site", which the recall exists to fix, and "the capture targeted it and DRAGEN still said
  nothing", which is a fact about that DRAGEN run and stays `NOCOV`. The design is configured
  per run, not defaulted, and the job fails if a DRAGEN reference block reaches a defining site
  outside it, which only a wrong file can produce; a DRAGEN deletion running off a capture edge
  is a carrier, and the site under it is covered and left unfilled. Which sites the
  design missed is subtracted once per design, in bedtools, and committed with the code,
  because that answer depends on no sample and no run, only on the design and the committed
  sites; judging each sample's silence stays per sample, in the conversion job, because it
  cannot be anything else.
- **The primary caller wins wherever both speak, and the design bound holds at every site a
  post-hoc record merely reaches.** A record kept for covering a hole is kept whole, so one
  anchored in a hole can extend over a neighbouring defining site it may not fill: one DRAGEN
  did call, or a hole inside the design. Where that record is a variant it is dropped, because
  keeping it would put two callers' alleles on one base in rbceq2's input with nothing to choose
  between them, or let the second caller decide a site the design targeted. The hole returns to
  `NOCOV`. A reference block is kept, since it asserts nothing rbceq2 sees and dropping it would
  lose the hole it was kept for; the QC then reads the same off-design BED the merge did and
  counts a post-hoc record only at a site the merge was allowed to fill, so an in-design hole
  the block spans stays `NOCOV`.
- **A call resting on a recovered site is always annotated as such.** The site's QC flag name
  carries `POSTHOC`, because the antigen then rests on a different caller, without the sample's
  DRAGstr model, over reads the capture design did not target. That is a fact a
  reviewer has to be able to see. The token does double duty: it equally says the system was
  **not typable from the primary caller alone**, since without the recall that site would have
  been `NOCOV` and there would be no call to report.
- **Provenance is not a severity level, and is not ranked against one.** Both are properties of
  the site, so a flag name carries both, joined with `+`: severity runs `NOCOV` > `DEL` > `LOWQ`
  and `POSTHOC` is appended to whichever applies, or stands alone when the recovered site is
  fine. Ranking them would let severity displace provenance, which on the validation cohorts
  would hide the recall's contribution to 234 of 681 reliant systems. `NOCOV` keeps one meaning:
  no record from either caller, and so no caller to name.

**Regions BED generation from RBCeq2's own `db.tsv`**
- RBCeq2 internally restricts the VCFs it handles to regions appearing in its `db.tsv` database, via a `build_intervals()` function. The biggest time-sink is the filter-and-convert step, because `bcftools` otherwise parses the entire genome.
- Clear time saving: we build a BED from `db.tsv` the same way RBCeq2 does and restrict every input gVCF to those regions. This takes run times from ~50min-1hr to a couple of minutes.
- The BED is generated by `scripts/gen_bg_resources.py` (almost an exact copy of RBCeq2's `build_intervals()`) and committed under `src/popgen_rbceq2/resources/`, so a run is reproducible against a known database version rather than whatever the tool ships today.
- Note: because this re-implements rather than calls RBCeq2, the BED can drift if upstream changes `build_intervals()`. It must stay a strict superset of every coordinate RBCeq2 reads, or calls go silently wrong. Regenerate it whenever the pinned RBCeq2 version moves.

**The DAG lives in one file**
- Stage classes carry no dependency information. `stages/pipeline.py` declares every edge and every Metamist registration, so the pipeline shape is readable in one place and a stage cannot quietly acquire an undeclared input. Tests fail if a stage reads an upstream it did not declare, or if a wired stage is unreachable from the requested set.

## Scope boundaries & ecosystem
### Scope
- Not a variant caller, with one deliberate exception: `PosthocGenotypeOffTargetSites` calls
  variants at blood-group defining sites an exome capture left uncalled. It exists to supply
  the primary caller's blind spots, not to replace it, and it never overrides a DRAGEN call.
- Not a fork or fix of RBCeq2
- Not the source/maintainer of blood-group data/knowledge
- No phasing (yet)
- No calling of structural-variant-defined alleles, which needs the DRAGEN SV and CNV VCFs rather than the gVCF alone. So ATP11C, CD99 and XG are blind today, and so are the null alleles of 16 further systems. Design and estimate in `docs/rbceq2_cnv_sv/`.

### Ecosystem
- **RBCeq2**: the genotyper
- **bcftools**: the gVCF conversion, the post-hoc merge, and the per-site DP/GQ extract
- **GATK HaplotypeCaller**: the post-hoc caller for exome off-target defining sites
- **cpg-flow**: CPG's workflow framework, over Hail Batch
- **analysis-runner**: how a run is launched, and where its config comes from
- **Metamist**: CPG's sample metadata system. Source of truth for inputs and where results are registered.

## Relevant external repositories and resources.
[RBCeq2 repo and source code](https://github.com/limcintyre/RBCeq2)

## Bets & open questions.
- Currently do not have phased data. This leads to multiple blood-group phenotypes assigned to individuals. All outputs are registered in the analysis object on metamist, regardless of this uncertainty.
- The QC thresholds (`min_depth = 10`, `min_gq = 20`) are deliberately quiet and gnomAD-aligned. They are there to catch the tail, not to reproduce a hard filter — raising them towards DP 20 / GQ 30 flags most of a typical sample's systems and makes the annotation useless. `min_gq` must sit on a gVCF GQ band edge to mean anything exact.
- The blood type is recorded in an Analysis `meta`, not on the SequencingGroup record. If it needs to be queryable directly on the SG, that is a SequencingGroup update mutation we have not written.
- Post-hoc recall runs without a DRAGstr model for the sample, so `--dragen-mode` is not fully
  DRAGEN-equivalent for STRs. Accepted to start with; revisit if post-hoc indel calls at
  STR-adjacent defining sites look discordant with DRAGEN's at sites both callers reach.
- Whether `POSTHOC` should count as flagged at all is a judgement, not a measurement. It is
  flagged today because a recovered site is a different kind of evidence and a reviewer should
  see it. If it turns out most exome systems carry one, the flag stops discriminating and the
  right answer is probably a separate provenance column rather than quietly demoting it to
  `PASS`.
- The recall is exome-only. A genome gVCF has no capture edge to stop at, so there is nothing
  for it to fix there — but that assumes any genome `NOCOV` site is genuinely unmappable rather
  than merely uncalled, which we have not tested.
- Post-hoc recall is validated on Twist VCGS and Agilent CREv2, 10 exomes each, and gains the
  same ~7.5 percentage points on both. It is untested on the other CREv2 variants in mackenzie
  (`SSXTLICREV2`, `SSQXTCREV2`, `SSQXTCRE`). The two designs measured behave alike and silence
  is still judged from the gVCF, so the expectation is that it generalises — but each new design
  now has to be named in config, and the run fails rather than guesses if the wrong one is.
- Post-hoc recall assumes the CRAM and the gVCF for a sequencing group are two outputs of the
  same DRAGEN run, and **the merge fails if their sample names disagree**. Merging mismatched
  inputs would splice another individual's genotypes into these calls at exactly the sites
  nothing else covers, where nothing downstream could catch it. An older block of mackenzie
  test CRAMs trips this benignly, carrying a retired sequencing-group ID for the same
  individual from an upstream reheadering bug; newer additions are clean. Those inputs are what
  should be fixed. Downgrading the check to a relabel-and-warn was tried for one commit and
  reverted, because from inside the job a real swap and a stale header are indistinguishable,
  and accommodating a test-bucket artefact in shipped behaviour weakens every cohort. Positive
  identity confirmation still belongs to somalier, whose output already sits beside these
  CRAMs.
- Which CRAM and gVCF `sequencing_group.cram`/`.gvcf` resolve to is Metamist's call, and this
  dataset registers several of each per sequencing group. cpg-flow keeps whichever the API
  returns last, with no ordering, so the pairing is not pinned by anything we control. It
  currently lands on the DRAGEN 3.7.8 CRAM and the matching recal gVCF, which is what we want,
  but registering another analysis would silently change it.

## The current slice.
An implementation of RBCeq2 as a cpg-flow workflow in CPG's infrastructure, ported from the
prototype stages in `ourdna_genomic_atlas`. Not yet running on production data from this repo.

## Domain terms
[See GLOSSARY.md](../GLOSSARY.md)
