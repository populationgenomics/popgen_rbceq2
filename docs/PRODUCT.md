# Popgen_RBCeq2

## What Popgen_RBCeq2 is
A pipeline that produces an estimate of the blood type for each individual in the input cohort.
- __In__: Variant calls in gVCF format.
- __Structure__: A Nextflow pipeline that converts the input into the correct format for RBCeq2, which it then calls.
- __Out__: Three TSV files:
  1. `<SGID>.geno.tsv` to genotype calls,
  2. `<SGID>.pheno_alphanumeric.tsv` to alphanumeric phenotype,
  3. `<SGID>.pheno_numeric.tsv` to numeric phenotype.

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
**gVCF -> VCF conversion + filtering (`filter_and_convert.nf`)**
- RBCeq2 cannot read gVCFs: Its variant encoder rejects any ALT ending in `>`, so the `<NON_REF>` symbolic allele present on every gVCF record (both ref blocks and real variants like G,`<NON_REF>`) falls through to a fallback that mis-keys the variant and never matches the database. RBCeq2 does not error, it silently reverts affected systems to reference (observed: ABO -> Undetermined on HG00096). `filter_and_convert.nf` therefore splits multiallelics (`bcftools norm -m -any`) and drops `<NON_REF>` (`view -e 'ALT="<NON_REF>"'`) so only clean biallelic records reach RBCeq2. RBCeq2 also has no depth/quality filtering, so depth/GQ filtering must happen upstream. Both are handled at conversion:
1. `bcftools norm -m -any`: splits multiallelic records to one ALT per row. RBCeq2 evaluates one allele at a time and misreads multiallelic sites, so this precedes filtering.
2. `bcftools view -e 'ALT="<NON_REF>" || FMT/DP="." || FMT/DP < min_depth || FMT/GQ="." || FMT/GQ < min_gq'`: drops `<NON_REF>` records and missing/sub-threshold depth or GQ. The explicit `="."` checks are required because a missing value does not reliably fail a `<` comparison in `bcftools`; without them, no-data sites pass the filter.
3. `--trim-alt-alleles`: drops ALT alleles not carried by any surviving genotype call.

**Regions BED generation from RBCeq2's own `db.tsv`**
- RBCeq2 internally restricts VCFs it handles to regions only appearing in its `db.tsv` database. It does this via a `build_intervals()` function. By far the biggest time-sink is the filter and conversion step using `bcftools` as it parses the entire genome. 
- Clear time saving: We build a BED file from the `db.tsv` the same way RBCeq2 does and restrict all input gVCF files to these regions. This DRAMATICALLY improves run times from ~50min-1hr to a couple of minutes.
- Pre-generated BED is located in root: `bg_regions.GRCh38.bed`. As well as the script to generate it `gen_bg_bed.py`. The script itself is almost an exact copy of `build_intervals()` used in RBCeq2.
- Note: no dependency on RBCeq2's implementation in this process means BED file can drift if upstream RBCeq2 changes `build_intervals()`.

## Scope boundaries & ecosystem
### Scope
- Not a variant caller
- Not a fork or fix of RBCeq2
- Not the source/maintainer of blood-group data/knowledge
- No phasing (yet)

### Ecosystem
- **RBCeq2**: the genotyper
- **bcftools**
- **Nextflow/Seqera**: workflow orchestrator and platform to manage workflows
- **Metamist**: CPG's sample metadata system. Source of truth for inputs and where results are registered.

## Relevant external repositories and resources.
[RBCeq2 repo and source code](https://github.com/limcintyre/RBCeq2)

## Bets & open questions.
- Currently do not have phased data. This leads to multiple blood-group phenotypes assigned to individuals. All outputs are registered in the analysis object on metamist, regardless of this uncertainty.
- Decision to remove sites that do not meet `min_depth`/`min_gq` instead of labelling them as `"./."` (no calls). Converting low quality sites to `"./."` would lead to false-positives at blood-group loci. Compared to the current implementation of deleting low quality sites leading to false-reference and therefore removed from analysis (except at Lane loci where they are added back in as hom-ref).

## The current slice.
An implementation of RBCeq2 using Nextflow in CPG's infrastructure. Not yet running on production data.

## Domain terms
[See GLOSSARY.md](../GLOSSARY.md)
