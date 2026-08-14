"""Values shared by more than one stage, and the pins that tie a run to a tool version."""

# rbceq2 tool version, and the cpg-common image tag built from it (tool version plus the
# CPG image build number). The tool version is recorded in the per-SG Analysis meta and
# the tag resolves the image, so both come from one literal here rather than drifting.
# Bump when the image is bumped, and regenerate resources/bg_*.<genome>.* against the
# db.tsv from the new version (scripts/gen_bg_resources.py).
RBCEQ2_VERSION = '2.4.3'
RBCEQ2_IMAGE_TAG = f'{RBCEQ2_VERSION}-1'

# rbceq2 emits one TSV per type per sample: <out>.geno.tsv etc. This drives rbceq2's own
# output resource group, so the QC TSV is deliberately NOT a member — it is produced by
# FlagBloodGroupCallQc, and adding it here would make rbceq2 expected to emit a file it
# never writes.
RBCEQ2_TSV_KEYS = ('geno', 'pheno_numeric', 'pheno_alphanumeric')
