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
  --skip-repo-checkout \
  popgen_rbceq2
```

`--skip-repo-checkout` is required, not a convenience: analysis-runner otherwise clones this
repo into the job, and a dataset only permits repos on its allowed list in the infrastructure
config. The image already holds the code at `<image-tag>`, which is what runs either way.

On a first run, append `--dry_run` after `popgen_rbceq2`: it checks the config and builds the
stage graph without submitting any jobs, so a config mistake costs seconds rather than a batch.

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

Sequencing groups without a gVCF are skipped, not failed. Exome runs get one extra stage,
`PosthocGenotypeOffTargetSites` per sequencing group, and an exome `only_stages` list has to
name it alongside the conversion stage; genome runs are unaffected.

### An exome run must name its capture design

Exome runs recover defining sites the capture never targeted, so they have to say which capture
that was. Set `exome_design_bed` to the BED the cohort's gVCFs were called against:

```toml
[workflow]
input_cohorts = ['COH123']
sequencing_type = 'exome'
exome_design_bed = 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'
```

[`exome_example_config.toml`](src/popgen_rbceq2/config/exome_example_config.toml) is a complete
example for a CREv2 cohort, with a placeholder cohort ID.

The value is the design's key in the [references](https://github.com/populationgenomics/references)
repo, but the pipeline never opens the design BED. It selects a committed resource instead: the
defining sites outside that design, `resources/bg_off_design_sites.<key>.<genome>.bed`,
subtracted once per design with `bedtools` by
[`scripts/gen_off_design_sites.py`](src/popgen_rbceq2/scripts/gen_off_design_sites.py) (see
[committed resources](#committed-resources)). The key is deliberately not defaulted: every default
is the wrong design for some cohort, and a wrong design does not fail, it just recovers the wrong
set of sites. An exome run that omits the key, or names a design with no committed subtraction,
fails while the stage graph is built, before any job starts. Genome runs never read it.

The design is also a segment of every exome output path, directly under the release:
`.../rbceq2_<tool>_<release>/<design key>/<stage>/...`. cpg-flow reuses a stage whose outputs
already exist without asking how they were made, and every exome output from the conversion
onward is built from the holes the design leaves. So repointing `exome_design_bed` starts a
fresh tree for the whole run rather than reusing the previous design's genotypes and QC. Genome
paths have no design segment.

The keys for the designs seen so far, so you do not have to open the
[references](https://github.com/populationgenomics/references) repo:

| capture design | `exome_design_bed` value (prefix every one with `exome_probesets_hg38/`) | file |
|---|---|---|
| Twist Comprehensive Exome + VCGS custom content | `twist_vcgs_custom_exome_covered_targets_bed` | `Twist_VCGS_Exome_Covered_Targets_hg38.bed` |
| Agilent SureSelect Clinical Research Exome v2 | `agilent_sureselect_clinical_research_exome_v2_covered_by_probes_bed` | `S30409818_Covered.bed` |

These are the mackenzie designs, the only two the recall has been run against end to end, and
the only two with a committed subtraction. Other vendor panels are in the same
`exome_probesets_hg38` section of the references repo, under the same naming convention; for a
cohort on one of them, run the generator with its key and the BED's `gs://` path and commit the
two files it writes:

```commandline
uv run python -m popgen_rbceq2.scripts.gen_off_design_sites \
    exome_probesets_hg38/<design>_bed gs://cpg-common-main/references/exome-probesets/hg38/<file>.bed \
    GRCh38 src/popgen_rbceq2/resources
```

**Where a vendor ships both, take `Covered`, not `Regions`.** They are different files: `Regions`
is the intervals the design aims at, `Covered` is the footprint the probes actually reach, and
DRAGEN emits over the latter. On the CREv2 validation cohort, naming `Regions` would have
recovered 691 sites instead of 1,650. Some older designs were specified at exon resolution and
ship a `Regions` entry only; there the choice does not arise.

Naming the wrong file is caught rather than tolerated. DRAGEN emits reference blocks over exactly
the design with no padding, so a DRAGEN reference block reaching a defining site *outside* the
configured design means the named file is not the one the gVCF was called against, and the job
fails saying so. That is what the `Regions`-for-`Covered` mistake trips over: 309 such sites on
the CREv2 cohort. A DRAGEN *variant* proves nothing, because its REF can run past the target
edge; a deletion at a capture edge that covers an off-design site is a carrier, and that site is
simply not filled.

If you do not know a cohort's design, read `sequencing_library` from the sequencing-group meta
in Metamist, and confirm the cohort is one design rather than a mix.

## How the code is laid out

| path | what's in it |
|---|---|
| [`stages/pipeline.py`](src/popgen_rbceq2/stages/pipeline.py) | the DAG: what each stage depends on, what it registers in Metamist, and `REQUESTED_STAGES` |
| [`stages/blood_group_genotyping/`](src/popgen_rbceq2/stages/blood_group_genotyping) | one undecorated stage class per module |
| [`stages/blood_group_qc/`](src/popgen_rbceq2/stages/blood_group_qc) | one undecorated stage class per module |
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

### `PosthocGenotypeOffTargetSites` (per sequencing group, exomes only)

Call the blood-group defining sites again from the CRAM, to fill the blind spots the exome
capture leaves. Skipped for a genome, and for an exome with no CRAM or no gVCF.

DRAGEN calls an exome against the capture-target BED with no `--vc-target-bed-padding`, so gVCF
emission hard-stops on the BED edges. A defining coordinate outside the capture gets **no record
at all** — not a low-quality one, none — and rbceq2 reads a site missing from its input as a
confident homozygous reference call. `NOCOV` makes that visible, but the system stays
unassessable.

The reads are usually there. Measured over 20 DRAGEN 3.7.8 exomes (5 Twist VCGS, 15 Agilent
CREv2) in August 2026:

| | |
|---|---|
| non-HPA defining coordinates off-target per capture design | ~165 of 1,599 |
| of those, within 1-100bp of a target edge with 40-110x MAPQ>=20 depth | ~105 |
| assessable non-HPA coordinates, before and after | 91.7% -> ~98% |

Only the caller stopped early, so this stage runs GATK HaplotypeCaller in GVCF mode over the
padded defining-site intervals. Sites more than ~250bp off both capture designs sit at ~0x and
are not recoverable at any padding. Re-running DRAGEN across a cohort is not affordable.

The GATK pin lives in [`constants.py`](src/popgen_rbceq2/constants.py), in two constants that
have to move together: `GATK_IMAGE_TAG` selects the image, and `GATK_VERSION` is what reaches
the QC flags as `src=gatk-hc-<version>`. The image tag carries a CPG build suffix, so it is not
simply the version string.

Measured on 20 exomes, 10 Twist VCGS and 10 Agilent CREv2: the recall lifts defining
coordinates that have a covering record from 89.7% to 97.4% (Twist) and 89.8% to 97.2%
(CREv2), for ~130s and $0.004 per sample. A recovered site that passes is indistinguishable
from a DRAGEN-called one, with recovered depth tracking each cohort's own baseline (73 against
76, and 141 against 150) at an identical median GQ of 99.

The two designs gain the same amount but at **different sites** — hole sets agree at Jaccard
0.99 within a design and 0.33 across them — which is what confirms the recall tracks real
capture boundaries. Full numbers in
[`docs/rbceq2_posthoc_exome/RESULTS.md`](docs/rbceq2_posthoc_exome/RESULTS.md).

Three things to preserve when changing this stage:

- **The CRAM is streamed, not localised.** The `gs://` path goes straight to `-I` and GATK
  reads it over NIO, pulling only the blocks `-L` asks for. The interval list is 199 regions
  over 136kb, so this reads megabytes where localising would move gigabytes per sample.
- **The reference must be the one the CRAM was aligned against**, the DRAGEN masked assembly38
  at `references.broad.ref_fasta`. CRAM stores reads differentially against a reference, so
  decoding with a different assembly does not fail loudly, it yields wrong bases. The mackenzie
  DRAGEN 3.7.8 CRAMs were aligned to unmasked hg38, which shares all 3,366 contig names with
  the masked build; the 786 contigs whose checksums differ are all HLA, and every blood-group
  defining site is on a primary chromosome, so nothing this stage reads is affected.
- **CRAM 3.0 only.** GATK 4.6.2.0 rejects CRAM 3.1 outright (`CRAM version 3.1 is not
  supported`). CPG's CRAMs are 3.0 today, so this is a trap for the future rather than a
  current problem.
- **`--dragen-mode`, and no re-alignment.** The primary calls being supplemented come from
  DRAGEN 3.7.8, so the supplement is made as close to DRAGEN-equivalent as a re-call can be.
  Reads are used as aligned: no DRAGMAP. There is no DRAGstr model for the sample, so STR
  genotyping is not fully DRAGEN-equivalent — accepted for now, worth revisiting if post-hoc
  indel calls near STRs look discordant.

The stage calls every assessable defining site, not only the off-target ones, because which
sites are off-target depends on the sample's capture kit and nothing is saved by knowing:
the whole interval list is 136kb either way. Deciding what to **keep** is the conversion
stage's job.

The output gVCF goes to tmp and registers no Metamist Analysis. It is an intermediate the
conversion stage consumes through the cpg-flow graph, and what a reader needs — that a call
rests on a recovered site — reaches Metamist as a `POSTHOC` flag on the QC TSV instead.

### `FilterAndConvertGvcfsForRbceq2` (per sequencing group)

Convert a gVCF into a VCF rbceq2 can read, using `bcftools`. This is the only stage that reads
the raw gVCF, 12-15Gb localised whole, so its resources are sized above the others. One pass
writes two outputs:

- `vcf`, the rbceq2 input. Multiallelics are split, the `<NON_REF>` symbolic allele is dropped
  because it breaks rbceq2, unused ALT alleles are trimmed, and a tabix index is written
  alongside because rbceq2 fetches blood-group regions by coordinate.
- `defining_sites`, holding FORMAT/GT, DP and GQ at every allele-defining coordinate for the
  QC stage, plus `INFO/POSTHOC` naming the caller that supplied each record. Do not derive this
  from the converted VCF: dropping `<NON_REF>` removes every reference block, and a reference
  block is what covers a defining site where no variant was called. GT tells the QC stage
  whether a deletion removed a defining base on one haplotype or both.

The `POSTHOC` column is extracted on every run, genome and exome alike, where it reads `.` for
each record. One extract format everywhere means the parser needs no per-sequencing-type branch
to know how many columns to expect.

#### Merging the post-hoc calls (exomes only)

For an exome, this stage also merges in `PosthocGenotypeOffTargetSites`'s gVCF before the
extract and the conversion run. A post-hoc record survives at a defining site only if **both**
tests pass:

1. **The DRAGEN gVCF has no record covering the site.** Judged per sample from the gVCF itself,
   never from capture metadata, so a design BED that misdescribes a sample's real footprint
   cannot overwrite a DRAGEN call or hide a hole. DRAGEN wins wherever both speak.
2. **The site lies outside the cohort's capture design**, read from the committed subtraction
   for the design named by `exome_design_bed`. See
   [naming the capture design](#an-exome-run-must-name-its-capture-design) for the key, the
   values, and why a wrong one fails the job rather than quietly under-filling.

A hole *inside* the design is left alone and reaches the QC as `NOCOV`. Recalling it would use
a second caller to answer a question about that sample's DRAGEN run, which is a different
question from the one this feature exists to answer. On the validation cohorts this is 4 of
1,674 Twist recoveries and 7 of 1,657 CREv2 ones, in C4B, RHD, RHCE and A4GALT. Keeping that
promise takes two steps, because a post-hoc record kept for an off-design hole is kept whole and
can reach an in-design hole beside it: the merge drops a post-hoc *variant* that does, and the
QC, handed the same off-design BED, disregards a post-hoc *reference block* there.

The two tests are two subtractions, and they run in different places for a reason. Taking the
design's intervals out of the defining sites depends on nothing about any sample, or any run:
only on two fixed files. So it is done **once per design**, by `scripts/gen_off_design_sites.py`
in `bedtools intersect -v`, a standard tool built for exactly this operation, and the result is
committed under `resources/` beside the sites it was subtracted from. Taking the DRAGEN records'
spans out of what is left is per sample by nature, so it stays in awk inside the conversion job,
whose image has no bedtools.

The sites the first subtraction keeps and the second drops are the off-design sites DRAGEN
*did* cover. Where the covering record is a reference block the job fails, because blocks stop
at the target edge and one reaching an off-design site means the named design is not the one
the gVCF was called against. Where it is a variant, most often a deletion anchored inside the
design whose REF runs past the edge, the site is left covered and unfilled, and the QC reports
`DEL` from DRAGEN's record.

A record is selected for reaching a fillable hole and is selected *whole*, so one anchored in a
hole can extend over a neighbouring defining site it may not fill: one DRAGEN did call, or a hole
inside the design. Defining sites are dense enough that this is the ordinary case near a capture
edge, and what happens next depends on the record:

- **A post-hoc variant reaching a site it may not fill is dropped.** At a called base, keeping
  it would put two callers' alleles on one base in the file rbceq2 reads, with nothing to choose
  between them, and the QC would still report that base as `PASS` because `resolve_coverage`
  prefers the primary record. At an in-design hole, keeping it would let the second caller decide
  the genotype at a site the design targeted. Either way the hole the dropped record would have
  filled goes back to `NOCOV`.
- **A post-hoc reference block reaching such a site is kept whole.** It asserts nothing rbceq2
  sees, since the conversion drops every `<NON_REF>`-only record first, and dropping the block
  instead would throw away the hole it was kept for. At a called base two records then cover the
  site in the extract, and `resolve_coverage` prefers the one with no `INFO/POSTHOC`. At an
  in-design hole the block is the only record, so `FlagBloodGroupCallQc` reads the same
  off-design BED and counts a post-hoc record only at a site in it; the hole stays `NOCOV`.

Mechanically the drop is an `INFO/COVERED` mark from `bcftools annotate -m`, which matches on a
record's whole span rather than its POS, followed by removing every marked record whose alleles
are not `<NON_REF>` alone. Two details carry the weight. The mark's source is the defining sites
the merge may not fill, every site less the fillable holes, not the DRAGEN records' spans, so a
variant that merely clips the tail of a long reference block is left alone. And the mark and
drop run before `norm -m -any` splits a variant from its `<NON_REF>` allele, so a variant is
recognised by its alleles rather than by any INFO tag, and no split-off `<NON_REF>` twin is
around to pass as a reference block.

A genome sequencing group has no post-hoc input and never reads the design key, so nothing is
merged for it and the merge is a plain rename. Its command is not otherwise unchanged, though:
every run now gains the `INFO/POSTHOC` header line on the intermediate and a trailing `POSTHOC`
column in the extract, which is why the release version bumped.

Four more details that are easy to get wrong:

- **Post-hoc records are given a FILTER value, because HaplotypeCaller leaves it `.`.** rbceq2
  uses an allele only when its defining variant is literally `PASS`, and a DRAGEN gVCF keeps a
  failed record with its filter name rather than removing it, so rbceq2 already excludes
  DRAGEN's `LowDepth` and `DRAGENSnpHardQUAL` alleles today. Left as `.`, every recovered
  alternate allele would be discarded before genotyping and the site typed as reference by
  absence, under a `POSTHOC` flag saying the call rested on the recall. The merge runs
  `bcftools filter -s LowDepth -e 'FORMAT/DP<=1'`, which marks DRAGEN's depth rule and writes
  `PASS` on everything else. DRAGEN's QUAL rule is not copied: on these cohorts QUAL is
  ML-recalibrated by a second DRAGEN pass and its threshold is not comparable to
  HaplotypeCaller's scale, so a post-hoc variant is used at any QUAL and its DP and GQ reach the
  QC flags, exactly as a `PASS` DRAGEN variant's do.
- Hole-finding reads covered spans with `--targets-overlap 1`, not the `2` the extract uses.
  Mode 1 asks whether the **record** overlaps, which is what `%END` reports and what the QC
  counts as covering; mode 2 asks whether the **variant** does, and drops a deletion anchored on
  the site itself. Using mode 2 would make a covered site look like a hole and let a post-hoc
  record displace a DRAGEN call. Note that mode 2 does *not* have that property on an unsplit
  gVCF record: `<NON_REF>` makes bcftools match the whole span there, which is why selecting the
  supplement with mode 2 needs the trespass drop above.
- **Zero-depth post-hoc records are dropped.** Given `-L`, HaplotypeCaller emits a reference
  block across the whole interval, including stretches with no reads, as `DP=0,GQ=0`. Keeping
  those would put a record over every hole and retire `NOCOV` for exomes entirely: a site with
  no reads would read `LOWQ(DP=0)`, which says "poor data" where the truth is "no data".
- **`INFO/POSTHOC` is declared on the intermediate for every run**, genome included, even
  though only exome records ever carry a value. `bcftools query` fails outright on a tag the
  header does not define rather than rendering `.`, so without the unconditional declaration
  the extract aborts on every genome run.

The supplement is stripped to the fields the pipeline reads (GT, DP, GQ, MIN_DP, END), which
keeps `concat` from having to reconcile two callers' definitions of tags nothing reads.

**The job fails if the CRAM and the gVCF name different samples.** This is the first thing
the job does, read from the two input headers before anything is computed. The post-hoc caller takes
its sample name from the CRAM's read group and DRAGEN named the gVCF from the same run, so a
mismatch means the two files this sequencing group resolves to do not describe one individual.
Merging them would splice another person's genotypes into these calls at exactly the sites
nothing else covers, and the result would look like an ordinary recovery.

Relabelling the supplement to the gVCF's name would satisfy `concat`, which requires identical
sample sets, and would bury that. Some mackenzie DRAGEN 3.7.8 test CRAMs do trip this benignly,
carrying a retired sequencing-group ID for the same individual from an upstream test-set
reheadering bug since fixed for newer additions. That is a reason to fix those inputs, not to
weaken the check for every cohort: this is the only place the pipeline compares the two files
it was handed, and from here a real swap and a stale header look the same. Confirm identity
with somalier, then fix the input.

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
`pheno_alphanumeric`), plus rbceq2's own run log. Each sequencing group's calls are registered
as a custom Metamist Analysis, with the calls parsed into its `meta`. This does not write the
blood type onto the SequencingGroup record.

#### The rbceq2 run log

rbceq2 writes a log on every run whether or not `--debug` is passed; the flag only raises it
from INFO to DEBUG. At INFO it is a version banner and the arguments used. At DEBUG it adds the
per blood group account of which filters removed which candidate alleles, which is what makes it
possible to answer why a given antigen was called. Measured cost is tens of KB per sample
(2.7 KB at INFO against 36 to 55 KB at DEBUG), so `--debug` is unconditional rather than
configurable: a run that was not made verbose cannot be explained after the fact without
repeating it.

Two things about capturing it:

- rbceq2 names the log `<out>_<uuid>_log.txt`, using a `uuid4` generated at runtime, so the path
  cannot be declared in `expected_outputs` and the job renames the file. The uuid is not lost by
  renaming: rbceq2 also writes it into each TSV as the index header, which is what ties a log
  back to the calls it explains.
- The rename requires exactly one match. loguru is configured to rotate the log at 50 MB, and if
  it ever did, moving the live segment alone would silently discard everything written before
  the rotation. `nullglob` is set first so that no match expands to nothing rather than to the
  literal pattern, which would otherwise be counted as one file and fail later on a file
  named `*`.

The log is written to GCS but deliberately **not** registered in Metamist. Nothing downstream
reads it, and stages take each other's outputs by path through the cpg-flow graph rather than
through Metamist, so a record would only earn its keep if something outside this pipeline had to
find the log. Registering it would also need a stage of its own, because cpg-flow allows one
analysis type per stage and runs every `analysis_keys` entry through the same
`update_analysis_meta` callback, and this stage's callback parses the geno TSV. See
[`docs/rbceq2_debug_log/SPEC.md`](docs/rbceq2_debug_log/SPEC.md).

### `FlagBloodGroupCallQc` (per sequencing group)

Join the DP/GQ extract to the committed site-to-system map and write `<sg>.qc.tsv` in the same
one-column-per-system layout as rbceq2's own TSVs, so it joins to the calls by column. This
stage reads only the extract, never the gVCF.

A system is `PASS`, or carries the semicolon-joined flags of its defining sites:

| flag | meaning |
|---|---|
| `NOCOV` | no record from either caller covers the site |
| `DEL` | a deletion the sample carries removed the base the antigen is defined on |
| `LOWQ` | DP or GQ below threshold, or missing |
| `POSTHOC` | the post-hoc caller supplied this site; joined to any severity with `+` |
| `NA` | the system has no assessable defining site, so it was never checked |
| `NOT_REPORTED` | cohort TSVs only: rbceq2 emitted no column for this system for this sample, so there was no cell to copy |

Severity runs `NOCOV` > `DEL` > `LOWQ`, and that is the order the checks are applied in. A site
that clears both thresholds has no severity and is listed only if it was recovered.

#### `POSTHOC`: the call rests on a recovered site

A flag name states two independent findings about the site, joined by `+`:

- **severity** — `NOCOV`, `DEL` or `LOWQ`, or absent when the site clears both thresholds;
- **provenance** — `POSTHOC` when the post-hoc caller supplied the record, absent when the
  primary caller did.

```
POSTHOC:1:159204893(T>C,src=gatk-hc-4.6.2.0,DP=42,GQ=99)
LOWQ+POSTHOC:1:3774964(A>G,src=gatk-hc-4.6.2.0,block=91bp,DP=1,MIN_DP=1,GQ=3)
```

`POSTHOC` does double duty. It says the antigen rests on a different caller, without the
sample's DRAGstr model, over reads the capture design did not target. It equally says the
system **was not typable from the primary caller alone**: without the recall that site would
have been `NOCOV`, so the call would not exist.

The two halves are joined rather than ranked, because ranking makes one displace the other
and a recovered site that is also poor has to report both. Ranking provenance below `LOWQ`
would lose it wherever both apply, which on the validation cohorts is 234 of 681 reliant
systems. Grepping `POSTHOC` finds all of them.

A site that is neither poor nor recovered is not listed at all, which is what makes a bare
`PASS` cell mean "nothing to report". `NOCOV` never carries `POSTHOC`: no record from either
caller means there is no caller to name. A post-hoc record counts as covering a site only where
the merge was allowed to fill one, so a hole inside the capture design that a kept post-hoc
reference block happens to span is still `NOCOV`; the QC stage reads the same off-design BED
the merge did.

A cell can therefore say both things at once. Only `POSTHOC` names is a clean recovery; a
`LOWQ+POSTHOC` or a `POSTHOC` beside a `NOCOV` site is a recovery that is also compromised,
and counts with the quality problems rather than the recoveries.

Read a flag as two parts: the site the database defines, then what the caller reported there.
The first field in the parentheses is always the database's allele; every later field is
`key=value`, and which keys appear tells you what the numbers describe.

```
NOCOV:6:31992067(CGT>C)
LOWQ:9:133255766(T>C,DP=19,GQ=15)
DEL:1:3774964(A>G,del=CATGA>C,GT=0/1,DP=30,GQ=50)
LOWQ:9:133255766(T>C,block=4.2kb,DP=26,MIN_DP=19,GQ=15)
POSTHOC:1:159204893(T>C,src=gatk-hc-4.6.2.0,DP=42,GQ=99)
```

- The first has no covering record at all, so there are no metrics to report. That differs
  from a record that carries the field but leaves it empty, which renders as `DP=.`.
- The second is a low-quality call at the site itself.
- The third means rbceq2 saw no variant here and called the system reference, but the
  deletion in `del=` removed the base that call rests on. `DEL` outranks `LOWQ` however good
  the deletion's own DP and GQ, and `GT=` says whether it removed the base on one haplotype
  or both.
- The fourth has no call at that coordinate. Its numbers come from a 4.2kb reference block
  covering it, so `DP` is that block's median depth and `MIN_DP` its shallowest base. The
  block may start on the site or reach it from an earlier position; either way the numbers
  describe the band, not the site, which is why `block=` is there.
- The fifth clears both thresholds, but DRAGEN never reported this site: the exome capture
  stopped short of it and `src=` names the caller that filled it. It is listed despite passing
  precisely because it was recovered. `src=` appears on any flag whose numbers came from the
  post-hoc caller, so a failing recovered site reads
  `LOWQ+POSTHOC:...(T>C,src=gatk-hc-4.6.2.0,DP=6,GQ=12)`.

The job's log counts quality-flagged systems and post-hoc-only systems separately, so a
successful exome run — where recovering sites is the point — does not read as a cohort that got
worse.

Set `min_depth` and `min_gq` in this stage's config section rather than the conversion stage's,
since this is the stage that reports them. Keep them quiet. They exist to catch the tail,
including apparent confident reference resting on very few reads, not to reproduce a hard
filter, which would flag most of a typical sample's systems.

### `CombineRbceq2OutputsPerCohort`

Concatenate the per-sequencing-group TSVs into cohort TSVs and register a cohort-level Analysis.

The cohort TSVs differ from the per-sequencing-group files they are built from in two ways.
Column 1 carries the CPG sequencing-group ID rather than rbceq2's own row label (the name of
the intermediate VCF it read), so the rows join to Metamist. And the columns are the union of
every sample's systems, sorted: where a sample's own file had no column for a system, the cell
reads `NOT_REPORTED`. That differs from `NA`, which means the QC job checked the system but it
has no assessable defining site; `NOT_REPORTED` means rbceq2 emitted no column at all.

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

94% of defining sites resolve from a reference block rather than a call of their own, which is
why the flag reports `block=<span>`.

## Committed resources

`resources/bg_*.<genome>.*` are generated from the rbceq2 allele database rather than fetched at
runtime, so a run is reproducible against a known database version:

| resource | contents |
|---|---|
| `bg_regions.<genome>.bed` | merged ±500kb intervals around the blood-group genes |
| `bg_defining_sites.<genome>.bed` | the allele-defining coordinates |
| `bg_defining_sites_padded.<genome>.bed` | merged ±250bp around those, the post-hoc caller's intervals |
| `bg_site_systems.<genome>.tsv` | `chrom/pos/ref/alt/kind/system` rows |
| `bg_off_design_sites.<design>.<genome>.bed` | the defining sites outside one exome capture design, one file per design |
| `bg_off_design_sites.manifest.tsv` | what each off-design BED was built from: the design's path and MD5, the sites BED's MD5, the counts, the bedtools version |

When you bump the rbceq2 image, bump `RBCEQ2_VERSION` in
[`constants.py`](src/popgen_rbceq2/constants.py) and regenerate all four with
[`scripts/gen_bg_resources.py`](src/popgen_rbceq2/scripts/gen_bg_resources.py) against the
`db.tsv` from that same version. One parse writes all four, so they always cover the same
sites. A build with no committed resources for the configured reference fails at graph
construction rather than per sequencing group.

The off-design BEDs are derived from the defining sites in turn, by
[`scripts/gen_off_design_sites.py`](src/popgen_rbceq2/scripts/gen_off_design_sites.py), one run
per design (the command is under [naming the capture design](#an-exome-run-must-name-its-capture-design)).
Whenever the defining sites are regenerated, re-run it for every design in the manifest: the
manifest records the MD5 of the sites BED each subtraction came from, and both the test suite
and the graph build refuse a subtraction whose sites are not the shipped ones, so a stale file
is a red suite, never a wrong fill.

The padded BED is built from the same non-SV sites the QC assesses, so a site the QC checks is
always one the post-hoc caller could reach. Its ±250bp comes from the exome coverage analysis:
past that, off-target sites sit at ~0x and no amount of padding recovers them. Padding is baked
into the committed file rather than read from config — changing it means regenerating and
committing, like every other resource. Keep it small: at 199 intervals over 136kb the CRAM can
be streamed, and padding raised towards the ±500kb regions BED would quietly turn every exome
job into a whole-CRAM read.

### Structural-variant entries

The database defines some alleles by a multi-kb event rather than a sequence, written
`<pos>_del_<N>kb`, `<pos>_dup_<N>kb`, `<pos>_DEL_<N>` or `<pos>_ins_<N>bp`. One base of DP
and GQ cannot say whether a sample carries a 21kb deletion, so `gen_bg_resources.py` excludes
these from the defining-sites BED and the site-to-system map: 72 of the database's 1871 site
rows at db 2.5.1.
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
its own. The design for merging them, and an estimate of which systems it would open up, is in
`docs/rbceq2_cnv_sv/`. The ABCC1 coordinate error is still present in the 2.4.4 database, so
ABCC1 stays uncallable whatever the input.

## Development

```commandline
uv sync --group dev      # install, including the package itself
uv run pytest            # the suite: no cloud, no Hail, no credentials
uv run ruff check .
uv run pyright
pre-commit install
```

The tests that run the real merge shell and the real design subtraction need `bcftools`,
`bgzip`, `tabix` and `bedtools` on PATH and skip otherwise; `uv run pytest -rs` lists any skip.
Install them with `brew install bcftools htslib bedtools` on macOS (`bgzip` and `tabix` come
from `htslib`, not the `bcftools` formula) or `apt-get install bcftools tabix bedtools` on
Debian and Ubuntu, which is what CI does.
CI installs them and runs those checks on every pull request, then builds the driver image and
runs the suite again inside it.

## Related

- [PRODUCT.md](docs/PRODUCT.md) — what this repo is for and the decisions that shape it
- [GLOSSARY.md](GLOSSARY.md) — domain terms
- [RBCeq2 source](https://github.com/limcintyre/RBCeq2)
