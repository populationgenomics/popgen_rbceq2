# Glossary

- RBCeq2: the genotyping tool this pipeline wraps; consumes VCF file and outputs blood group assignments per-sample.
- `db.tsv`: RBCeq2's bundled allele database; the source of truth for both calls and the regions BED.
- Blood group system/reported locus: RBCeq2 emits one genotype+phenotype call per gene/locus (88 systems in the v2.4.4 db, 35 of them HPA platelet systems; 86 without RHD and RHCE. Was 48 at v2.4.1, before HPA). These map to ISBT blood group systems but not 1:1, one system can involve several genes, and some reported loci are related transporters/regulators.
- Antigen vs phenotype vs genotype:
   - Genotype: The allele-pair call in ISBT nomenclature (`geno.tsv`), e.g. ABO*A1.01/ABO*O.01.05. RBCeq2 often lists many candidate pairs.
   - Phenotype (numeric): ISBT numeric antigen notation (`pheno_numeric.tsv`), e.g. ABCC1:1, CROM:1,2,-3 (system:antigen-number, sign = present/absent).
   - Phenotype (alphanumeric): conventional serological names (`pheno_alphanumeric.tsv`), e.g. P1+,Pk+, Fy(a+b−), Lan+. An antigen is the individual serological marker. The phenotype is the observed +/− pattern of antigens. The genotype is the  allele pair predicted to produce it.
- Antithetical antigens: a pair of antigens encoded by alternative alleles at the same locus, where having one usually means lacking the other.
Rare blood group: a phenotype for which antigen-compatible donor blood is hard to source, typically a rare antigen-negative combination; prevalence is population-specific.
- Lane variant: a position that's wildtype in the genomic reference but variant relative to the transcript (from Dr
Lane's paper); RBCeq2 adds the reference allele to complete the genotype.
- gVCF vs VCF: gVCFs contains information at every position in the genome, both reference and variant positions. VCFs contain information only at variant sites.
- `<NON_REF>`: sentinel symbolic ALT in gVCFs; breaks RBCeq2.
- `build_intervals` / regions / flank: how RBCeq2 (and the committed `bg_regions.<genome>.bed`) derive read regions from `db.tsv` (±500 kb).
- `Undetermined`: placeholder output value meaning RBCeq2 couldn't resolve a system, not read as "reference."
- `NOT_REPORTED`: cell value in the combined cohort TSVs for a system RBCeq2 emitted no column for in that sample's own file. Distinct from `NA` (system checked by QC but no assessable defining site) and `Undetermined` (RBCeq2 called the system but couldn't resolve it).
- Capture target BED / covered regions: the intervals an exome kit is designed to enrich, and
  which DRAGEN is given when calling an exome. With no `--vc-target-bed-padding` it stops
  emitting gVCF records at these edges, so anything outside them is silent rather than poorly
  called. Twist VCGS and Agilent CREv2 are the two designs in our cohorts; for CREv2 use the
  `Covered` BED, not `Regions`, because the gVCF boundaries match the covered footprint.
- Off-target site: a blood-group defining coordinate outside the capture target BED, so the
  exome gVCF has no record for it at all. Not the same as a low-quality site — there is nothing
  there to judge, which is why rbceq2 silently reads it as homozygous reference.
- Post-hoc calling / recall: calling defining sites again from the CRAM with GATK
  HaplotypeCaller, after the primary caller has run, to fill those off-target holes. Exome only.
  "Post-hoc" because it supplements calls already made rather than replacing the primary caller.
- `POSTHOC`: two related things. As an `INFO` field on a merged VCF record it names the caller
  that supplied that record where the primary gVCF was silent (`POSTHOC=gatk-hc-4.6.2.0`). As a
  QC flag name it marks a site the post-hoc caller supplied the record for. It is a provenance
  flag, not a quality one, and it is joined to the site's severity with `+` rather than ranked
  against it: a recovered site that passes the thresholds is `POSTHOC`, one that is also
  sub-threshold is `LOWQ+POSTHOC`. The same caller string appears inside the flag's metrics as
  `src=`.
- DRAGEN masked reference: `Homo_sapiens_assembly38_masked.fasta`, the assembly our CRAMs were
  aligned against (`references.broad.ref_fasta`). Required to decode a CRAM correctly — a
  different assembly yields wrong bases rather than an error.
- `--dragen-mode` / DRAGstr: GATK's DRAGEN-compatibility mode, and the per-sample STR model
  DRAGEN normally pairs with it. We run the mode without the model, so post-hoc STR genotyping
  is close to but not identical with DRAGEN's.
- geno / pheno_numeric / pheno_alphanumeric: the three output TSVs.
- SV VCF / CNV VCF: DRAGEN's two structural-variant outputs per sequencing group, beside the gVCF. The SV VCF is written by the DRAGEN SV caller, which integrates and extends Manta (breakpoint evidence, events of any size, proper `SVTYPE`; record IDs keep Manta's prefix); the CNV VCF is bin-based read-depth (dosage, `SVTYPE=CNV` on every record, events under 10 kb filtered as `cnvLength`). Genome runs produce both; the DRAGEN 3.7.8 exome runs produce only the SV VCF. Neither is an input to this pipeline yet; see `docs/rbceq2_cnv_sv/`.
- `DRAGEN:REF:` record: a CNV VCF record spanning an interval the CNV caller assessed and found at reference copy number, with its bin count (`BC`) and segment mean (`SM`). The CNV analogue of a gVCF reference block. The SV VCF has no equivalent: an absent SV call is silence.
- Structural-variant-defined allele: a db allele whose token names an event rather than a sequence change at a base (`<pos>_del_53kb`, or a spelled-out multi-kb deletion). rbceq2 matches these by fuzzy position and length through `SvReader`/`SvMatcher`, and only if the VCF record's `SVTYPE` equals the token's type. Excluded from the QC site map, so a system defined only by these reports `NA`.
- Hybrid allele: a gene-conversion product between paralogs (RHD/RHCE, GYPA/GYPB/GYPE), stored in the db as a paired deletion plus insertion or duplication. Not reliably detectable from short-read SV/CNV calls; long read only, and why `--RH` stays off.
- `Metamist`: CPG's sample-metadata system.
- `analysis-runner`: CPG's tool to launch workflows.
