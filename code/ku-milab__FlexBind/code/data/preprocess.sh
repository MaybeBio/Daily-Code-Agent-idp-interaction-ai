#!/bin/bash
set -euo pipefail

# Run from the directory containing this script
cd "$(dirname "$0")"

DRIVE_URL="https://drive.google.com/drive/folders/1yuIXYeA3iOvqufhN2uJCQiPea1Be8O3C"

mkdir -p raw processed

echo "0. Checking required tools..."

# psi-cd-hit.pl ships with the CD-HIT package and calls BLAST+ internally
for tool in gdown psi-cd-hit.pl blastp makeblastdb; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: '$tool' not found in PATH."
        exit 1
    fi
done

echo "1. Downloading raw datasets from Google Drive..."

# Download the entire Google Drive folder into raw/
gdown --folder "$DRIVE_URL" -O raw/ --continue

echo "2. Checking downloaded files..."

# Expected input filenames
for file in raw/dp81_raw.fasta raw/dp93_raw.fasta; do
    if [ ! -s "$file" ]; then
        echo "ERROR: Missing or empty file: $file"
        echo "Files downloaded:"
        find raw -type f
        echo "Check the actual filenames and update this script."
        exit 1
    fi
done

echo "3. Clustering DP81 (30% sequence identity)..."

# Outputs: processed/dp81_clustered.fasta (representatives)
#          processed/dp81_clustered.fasta.clstr (cluster membership)
psi-cd-hit.pl \
    -i raw/dp81_raw.fasta \
    -o processed/dp81_clustered.fasta \
    -c 0.3

echo "4. Clustering DP93 (25% sequence identity)..."

psi-cd-hit.pl \
    -i raw/dp93_raw.fasta \
    -o processed/dp93_clustered.fasta \
    -c 0.25

echo "5. Train/Validation/Test split generation..."

# Uncomment after implementing the splitting script:
# python scripts/generate_splits.py \
#     --dp81_clstr processed/dp81_clustered.fasta.clstr \
#     --dp93_clstr processed/dp93_clustered.fasta.clstr \
#     --out_dir processed/

echo "Download and clustering completed."
echo "Dataset splits have not been generated."
