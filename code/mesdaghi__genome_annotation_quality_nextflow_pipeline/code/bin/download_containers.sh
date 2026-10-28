#!/bin/bash
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Error: SIF directory must be specified." >&2
    echo "Usage: $0 /path/to/sif_dir" >&2
    exit 1
fi

if [ -z "${APPTAINER_TMPDIR:-}" ]; then
    echo "Error: APPTAINER_TMPDIR is not set. Please export it before running this script." >&2
    echo "  export APPTAINER_TMPDIR=/path/to/tmp" >&2
    exit 1
fi

SIF_DIR="$1"
REGISTRY="ghcr.io/mesdaghi/genome_annotation_quality_nextflow_pipeline"
TAG="sha-8bcab79"

mkdir -p "$SIF_DIR" || { echo "Error: Failed to create directory $SIF_DIR" >&2; exit 1; }

declare -A IMAGES=(
    ["GAQA_bash"]="bash"
    ["GAQA_python"]="python"
    ["GAQA_metapredict"]="metapredict"
    ["GAQA_protenix"]="protenix"
    ["GAQA_psauron"]="psauron"
    ["GAQA_interproscan"]="interproscan"
)

for SIF_NAME in "${!IMAGES[@]}"; do
    IMAGE_NAME="${IMAGES[$SIF_NAME]}"
    SOURCE="${REGISTRY}/${IMAGE_NAME}:${TAG}"
    DEST="${SIF_DIR}/${SIF_NAME}.sif"

    if [ -f "$DEST" ]; then
        echo "✓ Already exists, skipping: $DEST"
        continue
    fi

    echo "→ Building $DEST from docker://$SOURCE ..."
    apptainer build "$DEST" "docker://$SOURCE" || {
        echo "Error: Failed to build $DEST from docker://$SOURCE" >&2
        exit 1
    }
    echo "✓ Done: $DEST"
done

echo ""
echo "All images processed. SIFs are in: $SIF_DIR"
