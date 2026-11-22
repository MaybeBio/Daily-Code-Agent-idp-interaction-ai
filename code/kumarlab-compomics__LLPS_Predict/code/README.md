# LLPS_Predict

Our ensemble model integrates both protein sequence and structure information to predict the propensity of proteins to undergo liquid-liquid phase separation (LLPS). In the first step, the model distinguishes LLPS proteins from non-LLPS proteins. In the second step, proteins predicted to undergo phase separation are further classified as drivers or partners. Our method offers improved interpretability by pinpointing the sequence regions that drive phase separation. It also supports mutagenesis analysis and predicts the effects of point mutations on phase separation propensity.

## Installation
This tool needs to be installed before use. You can install this tool on a GPU machine. You may create a virtual environment. For example,
```
python -m venv llpslib
source llpslib/bin/activate
```

Clone the repository and install the requirements in `requreiments.txt`
```
git clone https://github.com/mingjie0508/LLPS_Predict.git
cd LLPS_Predict
pip install -r requirements.txt
```

## Model Checkpoints

Download the model weights by running
```
./checkpoints/download_ckpts.sh
```

This script creates two folders `foundation_models` and `trained_models`. It downloads project-specific model weights into `trained_models`:
- [Ensemble](https://huggingface.co/mingjiezhao0508/LLPS_Predict/blob/main/ensemble_baseline_h4_l1.pth): predicts the propensity of proteins to undergo LLPS
- [Ensemble_Driver](https://huggingface.co/mingjiezhao0508/LLPS_Predict/blob/main/ensemble_driver_baseline_h4_l1.pth): predicts the propensity of proteins to undergo homotypic LLPS
- [Ensemble_Partner](https://huggingface.co/mingjiezhao0508/LLPS_Predict/blob/main/ensemble_partner_baseline_h4_l1.pth): predicts the propensity of proteins to act as partners in LLPS
- [Ensemble Critical Region](https://huggingface.co/mingjiezhao0508/LLPS_Predict/blob/main/ensemble_critical_region_h4_l1.pth): predicts the LLPS profiles on the residue level using the ensemble model

Some foundation models may be cached in `foundation_models` while running subsequent scripts.

## Usage
To use as a library, refer to [tutorial.ipynb](https://github.com/kumarlab-compomics/LLPS_Predict/blob/HEAD/src/notebooks/tutorial.ipynb)
It illustrates how to:
- predict LLPS propensity scores on the protein level
- predict LLPS propensity scores on the residue level and identify critical regions
- perform mutagenesis analysis

## Project Structure
```
├── checkpoints           # Download pretrained weights here
│   ├── foundation_models # Some foundation models downloaded here automatically
│   └── trained_models    # Project-specific models
├── data                  # Processed data (more detail below)
└── src                   # Source code
    ├── experiments       # Configuration files, bash scripts, and output
    ├── lib               # High-level functions for users
    ├── loaders           # PyTorch dataset classes
    ├── models            # Model architecture
    ├── notebooks         # Jupyter notebook tutorials
    ├── scripts           # Train, test, and evaluate scripts
    └── utils             # Utils code (get structure, moving window, etc.)
```

## Dataset
The full datasets can be downloaded [here](https://drive.google.com/drive/folders/15QyVr6sMXyDO3L4jhZ0dn56vjVkGgS5S?usp=drive_link). Place all datasets in the ```data``` directory.
```
└── data                  # Processed data
    ├── clinvar
    ├── denovo_db
    ├── human_proteome
    ├── llps_driving_region
    ├── cellular_localization
    ├── test
    ├── training
    └── validation
```