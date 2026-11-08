# Installation and Setup

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [ProtTrans.tar.gz](ProtTrans.tar.gz)
- [install.sh](install.sh)
- [requirements.txt](requirements.txt)

</details>



## Purpose and Scope

This page provides comprehensive instructions for installing Disobind and setting up the required environment. It covers system requirements, dependency installation, environment configuration, and verification procedures. For guidance on running predictions after installation, see **1.2 Quick Start Guide**. For information about training custom models, see **4.5 Training Your Own Models**.

---

## System Requirements

### Operating System
Disobind is compatible with Linux and macOS systems. Windows users should use WSL (Windows Subsystem for Linux).

### Hardware Requirements
| Component | Minimum | Recommended |
|-----------|---------|-------------|
| RAM | 8 GB | 16 GB or more |
| Storage | 10 GB free | 50 GB or more (for datasets) |
| CPU Cores | 2 | 4 or more |
| GPU | None (CPU mode) | NVIDIA GPU with CUDA 11.8 support |

### Software Prerequisites
- **Python**: Version 3.9 (specified in the installation script [install.sh:15]())
- **Conda**: For environment management
- **CUDA Toolkit**: Version 11.8 (for GPU acceleration)
- **NVIDIA Drivers**: Latest compatible version (for GPU support)

Sources: [install.sh:15](), [requirements.txt:1-33]()

---

## Installation Workflow

```mermaid
flowchart TD
    Start["Start Installation"] --> CheckConda{"Conda<br/>Installed?"}
    CheckConda -->|No| InstallConda["Install Conda<br/>from conda.io"]
    CheckConda -->|Yes| CloneRepo
    InstallConda --> CloneRepo
    
    CloneRepo["Clone Repository<br/>git clone https://github.com/isblab/disobind.git"] --> EnterDir["Navigate to Directory<br/>cd disobind/"]
    
    EnterDir --> MakeExecutable["Make install.sh Executable<br/>chmod +x install.sh"]
    MakeExecutable --> RunInstall["Run Installation Script<br/>./install.sh"]
    
    RunInstall --> UnzipProt["Unzip ProtTrans dir<br/>tar -xzvf ./ProtTrans.tar.gz"]
    UnzipProt --> SetPath["Export PYTHONPATH<br/>Update ~/.bash_profile"]
    SetPath --> CreateEnv["Create Conda Environment 'diso'<br/>python=3.9"]
    CreateEnv --> InstallPkgs["Install Packages from<br/>requirements.txt"]
    
    InstallPkgs --> Verify["Verify Installation<br/>Test Import"] --> Complete["Installation Complete"]
    
    style Start fill:#f9f9f9
    style Complete fill:#f9f9f9
    style RunInstall fill:#e8e8e8
    style InstallPkgs fill:#e8e8e8
```

**Installation Workflow Diagram**: This flowchart illustrates the complete installation process defined in `install.sh`. The script automates the unzipping of language model assets, environment variable configuration, and dependency installation via `pip`.

Sources: [install.sh:1-27](), [ProtTrans.tar.gz:1-35]()

---

## Step-by-Step Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/isblab/disobind.git
cd disobind/
```

### Step 2: Run the Automated Installation Script

The repository includes a bash script `install.sh` that handles the entire setup process. 

```bash
chmod +x install.sh
./install.sh
```

The `install.sh` script performs the following operations:
1. **Asset Extraction**: Unzips `ProtTrans.tar.gz`, which contains essential components for the protein language model [install.sh:4]().
2. **Path Configuration**: Appends the current directory to `PYTHONPATH` in `~/.bash_profile` to ensure modules are importable [install.sh:7-12]().
3. **Environment Creation**: Creates a Conda environment named `diso` using Python 3.9 [install.sh:15]().
4. **Dependency Installation**: Activates the environment and installs all packages listed in `requirements.txt` [install.sh:18-21]().

Sources: [install.sh:1-27]()

---

## Core Dependencies

### Python Package Dependencies

The following diagram maps the major dependency categories to their specific packages:

```mermaid
graph TB
    subgraph "Deep Learning Framework"
        Torch["torch==2.0.1"]
        TorchVision["torchvision==0.15.2"]
        TorchAudio["torchaudio==2.0.2"]
        TorchMetrics["torchmetrics==1.3.1"]
    end
    
    subgraph "Protein Language Models"
        Transformers["transformers==4.33.1"]
        HFHub["huggingface-hub==0.15.1"]
        FairESM["fair-esm==2.0.0"]
        SentencePiece["sentencepiece==0.1.99"]
        Tokenizers["tokenizers==0.13.3"]
    end
    
    subgraph "Structural Biology"
        BioPython["biopython==1.81"]
    end
    
    subgraph "Data Processing"
        NumPy["numpy==1.24.3"]
        Pandas["pandas==1.5.0"]
        H5py["h5py==3.7.0"]
        Tables["tables==3.9.2"]
    end
    
    subgraph "Machine Learning Utilities"
        ScikitLearn["scikit-learn==1.3.0"]
        ImbalancedLearn["imbalanced-learn==0.11.0"]
        BetaCal["betacal==1.1.0"]
    end
    
    subgraph "Configuration & Utilities"
        OmegaConf["omegaconf==2.2.2"]
        TQDM["tqdm==4.64.0"]
        Requests["requests==2.29.0"]
        Wget["wget==3.2"]
        LXML["lxml==4.9.1"]
    end
```

**Dependency Categories**: Disobind relies on a specific stack of libraries. The `fair-esm` and `transformers` packages are critical for generating the sequence-based embeddings used as input to the `Epsilon_3` model.

Sources: [requirements.txt:1-33]()

---

## Post-Installation Verification

### Verify Python Environment

Activate the newly created Conda environment:

```bash
conda activate diso
```

Check that Python and key packages are correctly installed:

```bash
python --version  # Should show Python 3.9
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import transformers; print(f'Transformers: {transformers.__version__}')"
```

Sources: [install.sh:15](), [requirements.txt:15,23]()

### Verify Directory Structure

Ensure that `ProtTrans.tar.gz` was successfully extracted by the script:

```bash
ls -d ProtTrans/
```

The root directory should now also contain the `models/` and `params/` directories necessary for running `run_disobind.py`.

Sources: [install.sh:4]()

---

## Troubleshooting Common Installation Issues

### Issue: ProtTrans Extraction Failed

**Problem**: The `ProtTrans` directory is missing after running `install.sh`.

**Solution**: Manually extract the archive:
```bash
tar -xzvf ./ProtTrans.tar.gz
```
Ensure you have sufficient disk space, as language model assets can be large.

Sources: [install.sh:4]()

### Issue: PYTHONPATH Not Updated

**Problem**: Importing Disobind modules results in `ModuleNotFoundError`.

**Solution**: The `install.sh` script attempts to update `~/.bash_profile` [install.sh:7-12](). If you are using `zsh`, you may need to manually add the path to `~/.zshrc`:
```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)
```

### Issue: Dependency Conflicts

**Problem**: `pip install` fails due to version conflicts.

**Solution**: Disobind requires specific versions like `numpy==1.24.3` and `pandas==1.5.0`. Ensure you are installing into a clean `diso` environment as created by the script.

Sources: [requirements.txt:11,13](), [install.sh:15]()

---

## Next Steps

After successful installation:

1. **Quick Start**: Run a test prediction using the instructions in **1.2 Quick Start Guide**.
2. **GPU Setup**: If you have an NVIDIA GPU, verify that `torch.cuda.is_available()` returns `True` within the `diso` environment.
3. **Explore Data**: Review the `example/` folder to understand the input formats required by the system.

Sources: [install.sh:26-27]()

---