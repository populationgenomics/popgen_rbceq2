# popgen_rbceq2

Blood-group genotyping with [RBCeq2](https://www.rbceq.com/), as a
[cpg-flow](https://github.com/populationgenomics/cpg-flow) workflow.

RBCeq2 infers blood-group genotypes and phenotypes from short-read variant calls. It reads a
per-sample VCF, compares the calls at each system's allele-defining sites against its allele
database, and writes one column per system: the inferred genotype and the phenotype in numeric
and alphanumeric form. This workflow runs it over CPG's DRAGEN gVCFs and registers the results
in Metamist.

## Running it

```commandline
analysis-runner \
  --dataset <your-dataset> \
  --access-level test \
  --output-dir '' \
  --config src/popgen_rbceq2/config/popgen_rbceq2_default_config.toml \
  --config <your-cohort>.toml \
  --description "rbceq2 blood-group genotyping for <your-cohort>" \
  --image australia-southeast1-docker.pkg.dev/cpg-common/images/popgen_rbceq2:<image-tag> \
  popgen_rbceq2
```

Configs merge left to right, so your file overrides
[the defaults](src/popgen_rbceq2/config/popgen_rbceq2_default_config.toml). Set at least
`workflow.input_cohorts` and `workflow.sequencing_type`; see
[`config_template.toml`](src/popgen_rbceq2/config/config_template.toml). `analysis-runner`
requires `--output-dir`, but nothing reads it — stages derive their own paths.

To run part of the branch only, restrict the graph:

```toml
[workflow]
only_stages = ["FilterAndConvertGvcfsForRbceq2", "GenotypeBloodGroupsWithRbceq2"]
```

Sequencing groups without a gVCF are skipped, not failed.

## How the code is laid out

| path | what's in it |
|---|---|
| [`stages/pipeline.py`](src/popgen_rbceq2/stages/pipeline.py) | the DAG: what each stage depends on, what it registers in Metamist, and `REQUESTED_STAGES` |
| [`stages/blood_group_genotyping/`](src/popgen_rbceq2/stages/blood_group_genotyping) | one undecorated stage class per module |
| [`stage_support.py`](src/popgen_rbceq2/stage_support.py) | `wire()`, output prefixes, the config-section convention, job configuration |
| [`analysis_meta.py`](src/popgen_rbceq2/analysis_meta.py) | what each Analysis records in its `meta` |
| [`jobs/`](src/popgen_rbceq2/jobs) | scripts a stage runs in a batch job, invoked by path |
| [`scripts/`](src/popgen_rbceq2/scripts) | developer scripts, run by hand |
| [`resources/`](src/popgen_rbceq2/resources) | the committed blood-group site resources |

Dependencies are declared **only** in `stages/pipeline.py`. A stage class carries no `@stage`
decorator; `wire()` applies one, and adding a `requires=` there is what makes an
`inputs.as_path(...)` call legal. Read pipeline.py top to bottom to see the pipeline shape.

Each stage reads its own `[workflow.<section>]` config table, named after its class in
snake_case (`FlagBloodGroupCallQc` -> `[workflow.flag_blood_group_call_qc]`). Renaming a class
renames its section; `test_stage_support` fails if the two ever disagree.

## What the stages do

### `FilterAndConvertGvcfsForRbceq2` (per sequencing group)

Convert a gVCF into a VCF rbceq2 can read, using `bcftools`. This is the only stage that reads
the raw gVCF, 12-15Gb localised whole, so its resources are sized above the others. One pass
writes two outputs:

- `vcf`, the rbceq2 input. Multiallelics are split, the `<NON_REF>` symbolic allele is dropped
  because it breaks rbceq2, unused ALT alleles are trimmed, and a tabix index is written
  alongside because rbceq2 fetches blood-group regions by coordinate.
- `defining_sites`, holding FORMAT/GT, DP and GQ at every allele-defining coordinate for the
  QC stage. Do not derive this from the converted VCF: dropping `<NON_REF>` removes every
  reference block, and a reference block is what covers a defining site where no variant was
  called. GT tells the QC stage whether a deletion removed a defining base on one haplotype
  or both.

Two things to preserve when changing this stage:

- **Region restriction is unconditional**, using `resources/bg_regions.<genome>.bed`. This is
  what allows the single pass, and it beats reading the whole genome — minutes rather than the
  best part of an hour. It needs the gVCF `.tbi`, because `bcftools -R` jumps by index. Keep
  the BED a strict superset of every coordinate rbceq2 queries for `references.genome_build`,
  or calls go silently wrong.
- **Do not filter genotypes on DP or GQ.** rbceq2 treats a blood-group site missing from its
  input as a confident homozygous reference call, so dropping a borderline genotype invents a
  wild-type call instead of producing a no-call. DRAGEN has already hard-filtered these gVCFs.
  Report DP and GQ as a QC flag instead.

### `GenotypeBloodGroupsWithRbceq2` (per sequencing group)

Run rbceq2 in the pinned image and write its three TSVs (`geno`, `pheno_numeric`,
`pheno_alphanumeric`). Each sequencing group's calls are registered as a custom Metamist
Analysis, with the calls parsed into its `meta`. This does not write the blood type onto the
SequencingGroup record.

### `FlagBloodGroupCallQc` (per sequencing group)

Join the DP/GQ extract to the committed site-to-system map and write `<sg>.qc.tsv` in the same
one-column-per-system layout as rbceq2's own TSVs, so it joins to the calls by column. This
stage reads only the extract, never the gVCF.

A system is `PASS`, or carries the semicolon-joined flags of its defining sites:

| flag | meaning |
|---|---|
| `LOWQ` | DP or GQ below threshold, or missing |
| `DEL` | a deletion the sample carries removed the base the antigen is defined on |
| `NOCOV` | no gVCF record covers the site |
| `NA` | the system has no assessable defining site, so it was never checked |

Read a flag as two parts: the site the database defines, then what the caller reported there.
The first field in the parentheses is always the database's allele; every later field is
`key=value`, and which keys appear tells you what the numbers describe.

```
NOCOV:6:31992067(CGT>C)
LOWQ:9:133255766(T>C,DP=19,GQ=15)
DEL:1:3774964(A>G,del=CATGA>C,GT=0/1,DP=30,GQ=50)
LOWQ:9:133255766(T>C,block=4.2kb,DP=26,MIN_DP=19,GQ=15)
```

- The first has no covering record at all, so there are no metrics to report. That differs
  from a record that carries the field but leaves it empty, which renders as `DP=.`.
- The second is a low-quality call at the site itself.
- The third means rbceq2 saw no variant here and called the system reference, but the
  deletion in `del=` removed the base that call rests on. `DEL` outranks `LOWQ` however good
  the deletion's own DP and GQ, and `GT=` says whether it removed the base on one haplotype
  or both.
- The fourth has no call at that coordinate. Its numbers come from a 4.2kb reference block
  spanning it, so `DP` is that block's median depth and `MIN_DP` its shallowest base.

Set `min_depth` and `min_gq` in this stage's config section rather than the conversion stage's,
since this is the stage that reports them. Keep them quiet. They exist to catch the tail,
including apparent confident reference resting on very few reads, not to reproduce a hard
filter, which would flag most of a typical sample's systems.

### `CombineRbceq2OutputsPerCohort`

Concatenate the per-sequencing-group TSVs into cohort TSVs and register a cohort-level Analysis.

## Reference-block depth and GQ bands

A reference block covers many bases and reports one DP, one MIN_DP and one GQ for all of them.

- `DP` is the median depth across the block. Judge a site on this.
- `MIN_DP` is the block's minimum. Report it, do not filter on it, because one shallow base
  anywhere in the block sets it.

Ignore DRAGEN's v3.7 page, which calls `FORMAT/DP` the minimum across the band. `DP` exceeded
`MIN_DP` in 69% of 608,745 blocks, so the two cannot both be minima. The v3.10 page is right.

Read the GQ bands from the file, which declares them:

```
##GVCFBlock=minGQ=0(inclusive),maxGQ=10(exclusive)     ... 10-20, 20-30, 30-40 ...
##GVCFBlock=minGQ=40(inclusive),maxGQ=2147483647(exclusive)
```

We band at 10/20/30/40, coarser than DRAGEN's `1 10 20 30 40 60 80` default, to keep gVCF size
down. Some files carry no `##GVCFBlock` lines at all, so check rather than assume.

Blocks band on GQ, not on depth. A block's single GQ therefore answers a threshold sitting on a
band edge exactly: a block reporting GQ >= 20 holds no base below 20, and one reporting less
holds no base at or above it. **Keep `min_gq` on a band edge.** 10, 20, 30 and 40 all qualify;
25 does not, and would flag sites whose real GQ is 26 to 29.

Depth has no equivalent guarantee, but the GQ band still floors MIN_DP. Measured across the
blood-group regions of one sample, 608,745 blocks over 50.6 Mb:

| GQ band | blocks | median span | span p99 | MIN_DP min | MIN_DP p05 | MIN_DP median |
|---|---|---|---|---|---|---|
| `[0,10)` | 9.3% | 3 | 63 | **0** | 6 | 33 |
| `[10,20)` | 3.2% | 3 | 62 | **4** | 4 | 20 |
| `[20,30)` | 4.5% | 3 | 40 | **7** | 7 | 22 |
| `[30,40)` | 8.5% | 3 | 44 | **10** | 11 | 20 |
| `[40,inf)` | 74.5% | 37 | 1,323 | **14** | 17 | 30 |

A shallow base usually drops a block out of the top band and fragments it, so long blocks sit
in the top band and none reported MIN_DP below 14. The longest block in the region set, 10,202
bases at chr6:50,035,970, reports MIN_DP 25. Blocks run long where coverage is uniform.

Filtering on MIN_DP instead of DP would change one system on this sample: RHCE, from a 19bp
block reporting `DP=10,MIN_DP=9`. Across all blocks, MIN_DP fell below `min_depth` while DP did
not for 0.20% of them, never over 78 bases.

94% of defining sites resolve from a spanning block rather than a record of their own, which is
why the flag reports `block=<span>`.

## Committed resources

`resources/bg_*.<genome>.*` are generated from the rbceq2 allele database rather than fetched at
runtime, so a run is reproducible against a known database version:

| resource | contents |
|---|---|
| `bg_regions.<genome>.bed` | merged ±500kb intervals around the blood-group genes |
| `bg_defining_sites.<genome>.bed` | the allele-defining coordinates |
| `bg_site_systems.<genome>.tsv` | `chrom/pos/ref/alt/kind/system` rows |

When you bump the rbceq2 image, bump `RBCEQ2_VERSION` in
[`constants.py`](src/popgen_rbceq2/constants.py) and regenerate all three with
[`scripts/gen_bg_resources.py`](src/popgen_rbceq2/scripts/gen_bg_resources.py) against the
`db.tsv` from that same version. One parse writes all three, so they always cover the same
sites. A build with no committed resources for the configured reference fails at graph
construction rather than per sequencing group.

### Structural-variant entries

The database defines some alleles by a multi-kb event rather than a sequence, written
`<pos>_del_<N>kb`, `<pos>_DEL_<N>` or `<pos>_ins_<N>bp`. One base of DP and GQ cannot say
whether a sample carries a 21kb deletion, so `gen_bg_resources.py` excludes these from the
defining-sites BED and the site-to-system map: 60 of the database's 1828 site rows at v2.4.1.
Their positions stay in the regions BED, because rbceq2's own `parse_positions` keeps them and
the BED must stay a superset of what rbceq2 reads.

Structural variants are the only defining alleles for ABCC1, ATP11C and CD99, so those three
always report `NA` rather than `PASS`. `gen_bg_resources.py` names them when it runs, so the
list re-derives itself on a database bump.

ABCC1 is worse than unassessable. Its one row, `ABCC1*01N.01` on chr16, carries
`ABCC4*01N.01`'s coordinates verbatim in both builds (`GRCh37=95670705`, `GRCh38=95018451`),
positions that exist on chr16 in neither build. This is an upstream database error, and it
means rbceq2 cannot call the allele either. `test_committed_sites_fall_inside_their_contigs`
catches another one reaching the resources.

Calling these alleles needs the DRAGEN SV and CNV VCFs alongside the SNV calls, not the gVCF on
its own.

## Development

```commandline
uv sync --group dev      # install, including the package itself
uv run pytest            # the suite: no cloud, no Hail, no credentials
uv run ruff check .
uv run pyright
pre-commit install
```

CI runs those on every pull request, then builds the driver image and runs the suite again
inside it.

## Related

- [PRODUCT.md](docs/PRODUCT.md) — what this repo is for and the decisions that shape it
- [GLOSSARY.md](GLOSSARY.md) — domain terms
- [RBCeq2 source](https://github.com/limcintyre/RBCeq2)
