#!/bin/bash
# Example: compare FlexBind against a baseline converted to the FlexBind prediction format.
# See baseline/baseline_instructions.md for the full protocol.
set -euo pipefail
cd "$(dirname "$0")/.."

FLEXBIND_PRED=${1:-checkpoints/DP81/seed42/test_predictions.pkl}
BASELINE_PRED=${2:-baseline/predictions/DP81/DisoFLAG.pkl}

python scripts/bootstrap_test.py \
    --pred_a "$FLEXBIND_PRED" \
    --pred_b "$BASELINE_PRED" \
    --n_iter 1000 \
    --out "results/$(basename "${BASELINE_PRED%.pkl}")_bootstrap.json"
