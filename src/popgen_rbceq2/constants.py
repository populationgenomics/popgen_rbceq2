"""Values shared by more than one stage, and the pins that tie a run to a tool version."""

# rbceq2 tool version. Recorded in the per-SG Analysis meta, and stage_support derives the
# `rbceq2_<version>_<release>` segment of every output prefix from it, so bumping it starts
# a fresh output tree. When it moves, regenerate resources/bg_*.<genome>.* against the
# db.tsv from the new version (scripts/gen_bg_resources.py).
RBCEQ2_VERSION = '2.4.3'
# The cpg-common image tag, a separate literal rather than derived from the version above:
# the build-number suffix is owned by the images CI, which increments it on any rebuild
# (2.4.3-1 -> 2.4.3-2) independently of the tool version. Check the registry for the
# current tag when bumping either.
RBCEQ2_IMAGE_TAG = '2.4.3-1'

# The `workflow.sequencing_type` value of an exome run. An exome gVCF is called against a
# capture-target BED, which is what puts defining sites outside it beyond reach; a genome gVCF
# has no such edge. Post-hoc calling exists for this type, and stage_support keys an exome
# run's output tree on the capture design that type makes load-bearing.
EXOME = 'exome'

# GATK, for the post-hoc exome caller (PosthocGenotypeOffTargetSites). Pinned as its own
# literal for the same reason as the rbceq2 image tag: the build-number suffix moves on any
# rebuild. Check the registry when bumping.
GATK_VERSION = '4.6.2.0'
GATK_IMAGE_TAG = '4.6.2.0-2'

# bedtools, for the once-per-run capture-design subtraction (SelectOffDesignDefiningSites).
# Only `intersect -v` is used, so the 2022 build is fine and nothing here needs a newer one.
# Same suffix caveat as the tags above. No matching version constant: unlike the post-hoc
# caller, this tool's version does not reach any output, only the set of sites it selects.
BEDTOOLS_IMAGE_TAG = '2.30.0-1'

# The value of INFO/POSTHOC on a record the post-hoc caller supplied, written by the merge in
# FilterAndConvertGvcfsForRbceq2 and surfaced by the QC flag. It names the caller and version
# rather than being a bare flag, so a QC TSV says which caller stood in for DRAGEN, and a
# version bump is visible in the output rather than only in the code that produced it.
POSTHOC_CALLER = f'gatk-hc-{GATK_VERSION}'

# rbceq2 emits one TSV per type per sample: <out>.geno.tsv etc. This drives rbceq2's own
# output resource group, so the QC TSV is deliberately NOT a member — it is produced by
# FlagBloodGroupCallQc, and adding it here would make rbceq2 expected to emit a file it
# never writes.
RBCEQ2_TSV_KEYS = ('geno', 'pheno_numeric', 'pheno_alphanumeric')
