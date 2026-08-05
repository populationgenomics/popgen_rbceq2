# =============================================================================
# Driver image for the popgen_rbceq2 workflow.
#
# Base: cpg_hail (australia-southeast1-docker.pkg.dev/cpg-common/images), built by CPG.
# It provides Hail compiled from CPG's fork at the tag matching the image version — NOT the
# PyPI wheel — plus Python 3.11 (Debian bullseye), OpenJDK and Spark, with
# HAIL_QUERY_BACKEND=service preset.
#
# Tag layout: <hail-fork-tag>-<image-revision>, e.g. 0.2.138.cpg2-1 = hail fork tag
# 0.2.138.cpg2 plus image build revision -1.
#
# The image is the single source of truth for hail. Do NOT declare or pin hail in
# pyproject.toml — pip would risk reinstalling a PyPI wheel over this build. To bump hail,
# change the tag below.
#
# This image runs the cpg-flow driver and the jobs/ scripts. The compute images are pulled
# per stage from cpg-common (bcftools, rbceq2) and are not built here.
# =============================================================================

FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail:0.2.138.cpg2-1

# Read by .github/workflows/get_version.py to derive the next image tag, and kept in step with
# pyproject.toml's version by bumpversion.
ENV VERSION=0.1.0

COPY pyproject.toml ./
COPY README.md ./
COPY src src

RUN pip install .
