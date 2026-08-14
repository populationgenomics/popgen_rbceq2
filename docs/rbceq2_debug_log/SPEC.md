# Design Spec: capture the rbceq2 debug log

**Status:** DRAFT, for review. Do not implement until the design is approved.
**Area:** rbceq2 blood-group genotyping pipeline (`GenotypeBloodGroupsWithRbceq2`), in `popgen_rbceq2`.
**Author:** Joshua Schmidt · **Reviewers:** (fill in)

Run rbceq2 with `--debug` and save the per sample log to GCS and Metamist.

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
- Add a `RegisterRbceq2Log` stage typed `blood_group_log` that simply copies the log file, thus generating a registerable output. cpg_flow gives a stage one analysis type and runs every `analysis_keys` entry through the same meta callback, so registering the log from the genotype stage would type it as `blood_group_genotyping` and hand it to `blood_group_calls`, which parses the geno TSV and would fail. An alternative is to have a single cohort level stage that bundles every per sample log into a single registered file, which needs far fewer jobs but drops per sample lookup.
- Write the log to the main dataset prefix beside the TSVs, not tmp, because a registered Analysis cannot point at a file that gets cleaned up.

## Open questions

- One copy job per sequencing group is cheap, but is a job per sample acceptable at cohort scale?
