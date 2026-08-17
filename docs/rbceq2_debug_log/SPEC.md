# Design Spec: capture the rbceq2 debug log

**Status:** AGREED. Ready to implement.
**Area:** rbceq2 blood-group genotyping pipeline (`GenotypeBloodGroupsWithRbceq2`), in `popgen_rbceq2`.
**Author:** Joshua Schmidt · **Reviewers:** (fill in)

Run rbceq2 with `--debug` and save the per sample log to GCS.

## Findings

- `--debug` creates no new file. rbceq2 always writes one log; the flag only raises it from INFO to DEBUG.
- We already discard that log. The stage captures only the three TSVs, so it dies with the job.
- It records the rbceq2 and database versions, the arguments, and per blood group which filters removed which candidate alleles.
- Measured sizes: 2.7 KB at INFO, 36 to 55 KB with `--debug`.
- The file name carries a random UUID, but the same UUID sits in each TSV header, so we can rename the log and still match it back.

## Proposal

- Add `--debug` to the command, always on. Collaborators want the verbose output on every run.
- Add a `log` key to `expected_outputs` on the genotype stage, kept out of `RBCEQ2_TSV_KEYS`. That constant drives the TSV resource group and the cohort combine job, so a log listed there would be treated as a fourth TSV to concatenate. The alternative is to add it to the constant and filter it out at each place the constant is used.
- Rename the log inside the job. Fail unless the glob matches exactly one file, which catches loguru rotating the log at 50 MB.
- Do not record the log in Metamist. Nothing downstream reads it, and stages inside this pipeline take each other's outputs by path through the cpg-flow graph rather than through Metamist, so a Metamist record would only earn its keep if something outside this pipeline had to find the log. Registering it is also awkward: cpg-flow allows one analysis type per stage and runs every registered output through the same meta callback, and ours parses the geno TSV, so the log would need a stage of its own whose only work was copying a file we had already written.
- Write the log to the main dataset prefix beside the TSVs, not tmp, so it is still there when someone asks why a call was made.

## Open questions

None outstanding. The Metamist question was settled by team consensus on 2026-08-17: the log is written but not registered.
