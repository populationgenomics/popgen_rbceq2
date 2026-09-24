# Release v5 against v4: every per-sample output change, accounted for

Run 2026-09-24. Three cohorts were rerun on the current image (`popgen_rbceq2:0.1.0-11`, rbceq2
2.4.4, db 2.5.1, release v5) and every cell of every per-sample table compared with the same
sample's output from the previous release (`popgen_rbceq2:0.1.0-8`, rbceq2 2.4.3, release v4 for
exomes and v2 for the genome cohort, which are identical for genomes). Tables compared: `geno`,
`pheno_alphanumeric`, `pheno_numeric` and `qc`, one cell per blood-group system. Phenotype cells
are compared as sets of alternatives, QC cells as sets of flags.

| cohort | data | samples compared | old tree | new tree |
|---|---|---|---|---|
| tob-wgs COH5103 | genome | 982 | `rbceq2_2_4_3_v2` | `rbceq2_2_4_4_v5` |
| mackenzie COH13065 | Twist exome | 400 | `rbceq2_2_4_3_v4` | `rbceq2_2_4_4_v5` |
| mackenzie COH13198 | CREv2 exome | 400 | `rbceq2_2_4_3_v4` | `rbceq2_2_4_4_v5` |

The genome cohort has 983 members; one has no v4 output because rbceq2 2.4.3 asserted on its
A4GALT\*0XN allele and killed the job. It completed under 2.4.4 with an ambiguous A4GALT
genotype (four candidate pairs) and a correct, if open, P1 status. All three v5 runs completed
every job with no failures, and the tob-wgs cohort tables were written for the first time.

Every changed cell falls into one of six causes. Nothing was left unexplained.

## 1. Single-copy chrX genotypes render as single copy (PR #20)

The conversion no longer rewrites DRAGEN's one-token chrX genotypes to two tokens, so a one-copy
sample's XK, GATA1 and ATP11C genotypes read `XK*01/-` where v4 read `XK*01/XK*01`. Every
change in these three systems has exactly that shape; no other shape occurred.

| cohort | samples whose XK, GATA1 and ATP11C all changed | reported male in Metamist | reported female |
|---|---|---|---|
| tob-wgs | 409 of 982 | 398 | 11 |
| Twist | 200 of 400 | 200 | 0 |
| CREv2 | 192 of 400 | 192 | 0 |

Conversely 20 tob-wgs and 8 CREv2 samples reported male in Metamist did not change, so DRAGEN
called them with two X copies. The 39 discordant samples read the same DRAGEN gVCF in both
runs, and DRAGEN's ploidy estimator is off (`--enable-ploidy-estimator false`), so the copy
number came from the sex the ICA pipeline was given. v4 hid this by rewriting every sample to
two tokens; v5 reports what DRAGEN called. **Follow-up outside this repo:** reconcile those 39
samples' reported sex with their DRAGEN karyotype.

## 2. Exome post-hoc fills gated out of single-copy chrX (PR #20)

For a one-copy exome sample the merge no longer fills off-design XK sites from the diploid
HaplotypeCaller call, so those sites return to `NOCOV`. Only samples in cause 1 are affected,
and only the XK QC cell.

| cohort | samples | QC flags removed | flags added |
|---|---|---|---|
| Twist | 200 | 1,783 `POSTHOC`, 52 `LOWQ+POSTHOC` | 1,835 `NOCOV` |
| CREv2 | 118 | 1 `POSTHOC`, 122 `LOWQ+POSTHOC` | 123 `NOCOV` |

Twist has ten off-design XK sites and CREv2 two, which is the size difference. On CREv2 the
fills were almost all sub-threshold anyway.

## 3. New rbceq2 filter removes an impossible GYPA alternative

rbceq2 2.4.4 adds `cant_pair_with_ref_cuz_shared_variant_has_too_few_copies`. A sample
heterozygous at all three GYPA\*01-defining sites used to be offered two pairs,
`GYPA*01/GYPA*02` and `GYPA*02/GYPA*08`; the second needs the shared c.59T>C on a copy that also
carries the other two variants, which a heterozygote does not have. The alternative, and its
`M+,N+,Mc+` phenotype, is gone. The remaining call is the one v4 also listed first.

| cohort | samples | transition |
|---|---|---|
| tob-wgs | 509 of 982 | `GYPA*01/GYPA*02,GYPA*02/GYPA*08` to `GYPA*01/GYPA*02` |
| Twist | 180 of 400 | same |
| CREv2 | 188 of 400 | same |

One further Twist sample went from `GYPA*02/GYPA*08` to `Undetermined`. Its c.72G site is
`NOCOV` and its c.71A>G and c.59T>C hets have GQ 3 and 5; 2.4.3 had dropped the c.71A>G het
from the variant pool and called a pair the data does not support. Both versions' QC cells
already flag every site involved, so the v5 result is the honest one.

## 4. Database edits in db 2.5.1

Six rows differ between the two databases. Four allele definitions changed (ABO\*A2.16 gained
chr9:133257521 T>TC and lost a duplicated variant; FY\*01N.09 and RHCE\*01.20.02.01 lost one
variant each; GCNT2\*01N.10 lost a self-contradicting `ref` at the site of its own variant),
the KEL\*02 phenotype row now reads VLAN−, VONG− and KYO−, and two XG rows changed only their
coding annotation. Output moved in three of those systems:

| system | cohort | samples | what changed |
|---|---|---|---|
| ABO | tob-wgs 31, Twist 6, CREv2 10 | 47 | the `A1,A2` phenotype alternative and its seven `A1.02/A2.16`-style pairs are gone; A2.16 now needs the insertion |
| ABO | tob-wgs | 1 | `A2 \| O` became `O` for the same reason |
| ABO | tob-wgs | 2 | a homozygous `O.01.83/O.01.83` pair dropped by the reworked `ensure_HET_SNP_used`; the het pairs remain (rbceq2 logic, not the database) |
| GCNT2 | Twist 6, CREv2 6 | 12 | `GCNT2*01/GCNT2*01.02` now offered beside `GCNT2*01/GCNT2*01` for c.G>C heterozygotes |
| GCNT2 | CREv2 | 3 | a `LOWQ` flag at the `ref` row of that site is gone, because the row left our site map with the database edit |
| KEL | tob-wgs 1, Twist 1, CREv2 2 | 4 | VLAN, VONG and KYO flipped to negative on the one KEL\*02 phenotype path that carried them |

No FY, RHCE or XG cell moved in any cohort.

## 5. Lewis phenotype now uses weak-secretor status

rbceq2 2.4.4 rewrote the FUT3 interpretation: with active FUT3, a weak secretor (`Se+w`) is now
reported `Le(a+b+)`. Eighteen samples gained that alternative, none lost one.

| cohort | samples | transitions |
|---|---|---|
| tob-wgs | 1 | `Le(a+b-)/Le(a-b+)` to `Le(a+b+)/Le(a+b-)/Le(a-b+)` |
| Twist | 8 | 6 `Le(a-b+)` to `Le(a+b+)`, 2 gained `Le(a+b+)` as an alternative |
| CREv2 | 9 | 5 `Le(a-b+)` to `Le(a+b+)`, 4 gained `Le(a+b+)` as an alternative |

## 6. Formatting only

KLF alternatives are listed in a different order in 2.4.4. Compared as sets they are identical,
so these do not appear in the counts above.

## What did not change

Every other system, in every table, in all 1,782 samples: no cell differs. In particular no QC
cell on a genome moved, no exome QC cell outside XK and GCNT2 moved, and no HPA, RHD or RHCE
column appeared or disappeared (rbceq2 writes none of them to these tables in either version).

## Known limitation carried forward

The HPA rows in db 2.5.1 still hold GRCh37 coordinates in the GRCh38 column, so our GRCh38 site
map places all 35 HPA systems at intergenic positions. HPA output is not in these tables and the
coverage statistics already exclude it; the fix is upstream.
