# Writing your repo's GLOSSARY.md

- RBCeq2: the genotyping tool this pipeline wraps; consumes VCF file and outputs blood group assignments per-sample.
- `db.tsv`: RBCeq2's bundled allele database; the source of truth for both calls and the regions BED.
- Blood group system/reported locus: RBCeq2 emits one genotype+phenotype call per gene/locus (48 in v2.4.1: ABO, FY, VEL, GYPA, GYPB…). These map to ISBT blood group systems but not 1:1, one system can involve several genes, and some reported loci are related transporters/regulators.
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
- `build_intervals` / regions / flank: how RBCeq2 (and the pre-built BED) derive read regions from `db.tsv` (±500 kb).
- `Undetermined`: placeholder output value meaning RBCeq2 couldn't resolve a system, not read as "reference."
- geno / pheno_numeric / pheno_alphanumeric: the three output TSVs.
- `Metamist`: CPG's sample-metadata system.
- `analysis-runner`: CPG's tool to launch workflows.
