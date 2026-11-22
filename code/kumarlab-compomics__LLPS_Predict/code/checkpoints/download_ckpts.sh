#!/bin/bash

# absolute path to this script
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"

# change directory to script directory
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
cd "$SCRIPT_DIR"

# create model directories if not exist
FOUNDATION_MODELS_DIR="foundation_models"
TRAINED_MODELS_DIR="trained_models"

if [ ! -d "$FOUNDATION_MODELS_DIR" ]; then
    mkdir -p "$FOUNDATION_MODELS_DIR"
fi

if [ ! -d "$TRAINED_MODELS_DIR" ]; then
    mkdir -p "$TRAINED_MODELS_DIR"
fi

# change directory to trained models directory
cd "$TRAINED_MODELS_DIR"

# Use either wget or curl to download the checkpoints
if command -v wget &> /dev/null; then
    CMD="wget"
elif command -v curl &> /dev/null; then
    CMD="curl -L -O"
else
    echo "Please install wget or curl to download the checkpoints."
    exit 1
fi

# download model checkpoints
ENSEMBLE_URL="https://huggingface.co/mingjiezhao0508/LLPS_Predict/resolve/main/ensemble_baseline_h4_l1.pth"
ENSEMBLE_DRIVER_URL="https://huggingface.co/mingjiezhao0508/LLPS_Predict/resolve/main/ensemble_driver_baseline_h4_l1.pth"
ENSEMBLE_PARTNER_URL="https://huggingface.co/mingjiezhao0508/LLPS_Predict/resolve/main/ensemble_partner_baseline_h4_l1.pth"
ENSEMBLE_CRITICAL_REGION_URL="https://huggingface.co/mingjiezhao0508/LLPS_Predict/resolve/main/ensemble_critical_region_h4_l1.pth"

echo "Downloading first-step ensemble PS model checkpoint..."
$CMD $ENSEMBLE_URL
echo "Downloading second-step PS-Self model checkpoint..."
$CMD $ENSEMBLE_DRIVER_URL
echo "Downloading second-step PS-Part model checkpoint..."
$CMD $ENSEMBLE_PARTNER_URL
echo "Downloading critical region model checkpoint..."
$CMD $ENSEMBLE_CRITICAL_REGION_URL

echo "All checkpoints are downloaded successfully."
