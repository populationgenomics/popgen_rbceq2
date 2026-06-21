#!/usr/bin/env bash
#
# Build and push the rbceq2 container image to the seqera-trial Artifact Registry.
#
# One-time setup (so Docker can auth to Artifact Registry via gcloud):
#   gcloud auth configure-docker australia-southeast1-docker.pkg.dev
#
# Usage:
#   ./build_image.sh           # builds & pushes tag 2.4.1 (matches pinned rbceq2 version)
#   ./build_image.sh 2.4.1-2   # override the tag, e.g. for an image rebuild
#
set -euo pipefail

# --- config -----------------------------------------------------------------
PROJECT_ID="seqera-trial"
REGION="australia-southeast1"
REPOSITORY="images"
IMAGE_NAME="rbceq2"
TAG="${1:-2.4.1}"                       # default tag = the rbceq2 version baked into the Dockerfile
CONTEXT="images/rbceq2"                 # build context (contains the Dockerfile)

IMAGE_PATH="${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPOSITORY}/${IMAGE_NAME}:${TAG}"

# --- build & push -----------------------------------------------------------
echo "Building ${IMAGE_PATH} (linux/amd64) from ${CONTEXT}/ ..."

docker build --platform=linux/amd64 -t "${IMAGE_PATH}" "${CONTEXT}/"

echo "Pushing ${IMAGE_PATH} ..."
docker push "${IMAGE_PATH}"

echo "Done. Reference this in nextflow.config as:"
echo "  rbceq2_container = '${IMAGE_PATH}'"
