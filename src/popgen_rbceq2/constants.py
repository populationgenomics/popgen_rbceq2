"""Values shared by more than one stage, and the pins that tie a run to a tool version."""

# rbceq2 tool version. Recorded in the per-SG Analysis meta, and stage_support derives the
# `rbceq2_<version>_<release>` segment of every output prefix from it, so bumping it starts
# a fresh output tree. When it moves, regenerate resources/bg_*.<genome>.* against the
# db.tsv the new version bundles (scripts/gen_bg_resources.py).
RBCEQ2_VERSION = '2.4.4'
# The allele database bundled with that version, versioned independently of the tool. The
# committed resources/bg_*.<genome>.* are built from this db, not from the tool version above,
# so it is what says whether they are current. Recorded for that; nothing reads it.
RBCEQ2_DB_VERSION = '2.5.1'
# The cpg-common image tag, a separate literal rather than derived from the version above:
# the build-number suffix is owned by the images CI, which increments it on any rebuild
# (2.4.4-1 -> 2.4.4-2) independently of the tool version. Check the registry for the
# current tag when bumping either.
RBCEQ2_IMAGE_TAG = '2.4.4-1'

# The single-copy stretch of chrX, 1-based inclusive: everything between PAR1 and PAR2. See
# README, "Sex-chromosome ploidy".
#
# Keyed by genome build because the boundaries move between builds, and absent rather than
# guessed for a build whose resources this repo does not ship. chrY needs no entry: no
# blood-group defining site or region is on it, so no chrY record ever reaches rbceq2.
NON_PAR_X: dict[str, tuple[int, int]] = {
    # PAR1 ends at 2,781,479 and PAR2 begins at 155,701,383.
    'GRCh38': (2_781_480, 155_701_382),
}


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
