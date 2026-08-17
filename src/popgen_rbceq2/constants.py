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

# rbceq2 emits one TSV per type per sample: <out>.geno.tsv etc. This drives rbceq2's own
# output resource group, so the QC TSV is deliberately NOT a member — it is produced by
# FlagBloodGroupCallQc, and adding it here would make rbceq2 expected to emit a file it
# never writes.
RBCEQ2_TSV_KEYS = ('geno', 'pheno_numeric', 'pheno_alphanumeric')
