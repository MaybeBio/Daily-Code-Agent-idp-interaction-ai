# Troubleshooting and FAQ

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [example/test.csv](example/test.csv)
- [install.sh](install.sh)
- [requirements.txt](requirements.txt)

</details>



This page provides solutions to common issues encountered when installing, configuring, and running Disobind. It covers installation problems, input format errors, runtime failures, AlphaFold integration issues, and performance optimization.

---

## Common Installation Issues

### Conda Environment Setup Failures

**Problem**: Installation script fails or environment creation errors.

**Solution**:
1. Ensure Conda is properly installed and in PATH.
2. Run installation commands sequentially:
```bash
cd disobind/
chmod +x install.sh
./install.sh
```
3. If the script fails, manually create the environment using the steps defined in the script [install.sh:15-21]():
```bash
conda create --name diso python=3.9
conda activate diso
pip install -r ./requirements.txt
```

Sources: [README.md:19-36](), [install.sh:1-27]()

---

### CUDA and GPU Configuration

**Problem**: GPU not detected or CUDA errors during prediction.

**Solution**:
1. Verify CUDA toolkit version 11.8 is installed [README.md:38]().
2. Check NVIDIA driver compatibility.
3. Verify PyTorch CUDA installation:
```python
import torch
print(torch.cuda.is_available())  # Should return True
```
4. When running predictions, explicitly specify device using the `-d` flag [README.md:91]().

Sources: [README.md:38](), [README.md:91](), [requirements.txt:15]()

---

### Dependency Version Conflicts

**Problem**: Import errors or version incompatibility messages.

**Dependency Compatibility Matrix**

| Package | Required Version | Critical For |
|---------|-----------------|--------------|
| `torch` | 2.0.1 | Model inference and training |
| `transformers` | 4.33.1 | ProtT5-XL-U50 embeddings |
| `biopython` | 1.81 | PDB/CIF parsing |
| `numpy` | 1.24.3 | Array operations |
| `omegaconf` | 2.2.2 | YAML configuration loading |
| `h5py` | 3.7.0 | Embedding storage |

**Solution**:
1. Install exact versions from `requirements.txt` [requirements.txt:1-33]():
```bash
pip install -r requirements.txt
```
2. If conflicts persist, create a fresh environment and reinstall.

Sources: [requirements.txt:1-33]()

---

## Input Format and Validation Errors

### CSV Input Format Issues

**Problem**: Predictions fail due to malformed input files.

#### Valid Input Formats

```mermaid
graph TD
    Input["Input CSV File"]
    
    subgraph "Format 1: Disobind Only"
        F1["6 columns:<br/>UniProt_ID1, start1, end1,<br/>UniProt_ID2, start2, end2"]
    end
    
    subgraph "Format 2: Disobind + AF2"
        F2["12 columns:<br/>UniProt_ID1, start1, end1,<br/>UniProt_ID2, start2, end2,<br/>AF2_struct_path, AF2_pae_path,<br/>chain1, chain2, offset1, offset2"]
    end
    
    Input --> Check{"Column Count"}
    Check -->|"6 columns"| F1
    Check -->|"12 columns"| F2
    Check -->|"Other"| Error["Input Error"]
    
    F1 --> Valid["Parsed by run_disobind.py"]
    F2 --> Valid
```

**Solution**:
1. Verify column count is exactly 6 or 12 as specified in [README.md:56-60]().
2. Ensure residue positions are integers.
3. Example valid inputs can be found in `example/test.csv` [example/test.csv:1-2]().

Sources: [README.md:52-65](), [example/test.csv:1-2]()

---

### UniProt ID and Sequence Download Failures

**Problem**: The script hangs or errors while fetching sequences from UniProt.

**Solution**:
1. Verify UniProt IDs are valid at [uniprot.org](https://www.uniprot.org).
2. Check network connectivity; the script uses `requests` [requirements.txt:26]() to fetch data.
3. Reduce the number of cores used for downloading to avoid rate-limiting [README.md:89]():
```bash
python run_disobind.py -i csv -f ./example/test.csv -c 1
```
4. If a download fails, check the `output/` directory for `UniProt_seq.json` to see if some IDs were successfully cached.

Sources: [README.md:89](), [requirements.txt:26]()

---

## Runtime Prediction Errors

### Model Loading Failures

**Problem**: Model checkpoints cannot be found.

**Solution**:
1. Ensure you have run `install.sh` which unzips `ProtTrans.tar.gz` [install.sh:4]().
2. Check that model weights are present in the expected directories (e.g., `models/Epsilon_3/`).
3. If weights are missing, re-download them from Zenodo [README.md:9]().

Sources: [install.sh:4](), [README.md:9]()

---

### Memory Issues and Out-of-Memory (OOM) Errors

**Problem**: Prediction crashes with `RuntimeError: CUDA out of memory`.

**Solution**:
1. The ProtT5-XL-U50 model used for embeddings is large. If GPU memory is low, run on CPU [README.md:91]():
```bash
python run_disobind.py -i csv -f ./example/test.csv -d cpu
```
2. Close other GPU-intensive applications.
3. Disobind supports batch processing; if processing a very large CSV, try splitting it into smaller files.

Sources: [README.md:91](), [requirements.txt:23]()

---

## AlphaFold Integration Issues

### Structure and PAE File Mismatches

**Problem**: AlphaFold-enhanced predictions fail or produce illogical results.

**Solution**:
1. Ensure the `chain1` and `chain2` IDs in the CSV match the labels in the PDB file [README.md:62]().
2. Verify the `offset` values. If the AF2 structure residue 1 corresponds to UniProt residue 10, the offset is 9 [README.md:63]().
3. The PAE file must be the `.json` output from AlphaFold [README.md:60]().

Sources: [README.md:59-65](), [example/test.csv:2]()

---

## Frequently Asked Questions

### Q1: Can I predict interactions for more than two proteins?
**A**: Disobind is designed for binary complexes (AB). For complexes with more subunits (ABC), you must break them into binary pairs (AB, BC, AC) and run them individually [README.md:43]().

### Q2: What is the difference between contact map and interface residue prediction?
**A**: 
- **Contact Map**: Predicts specific residue-residue interactions (e.g., Res 10 of P1 interacts with Res 40 of P2) [README.md:111](). Enable with `-cm` flag [README.md:92]().
- **Interface Residues**: Predicts which residues are part of the binding interface generally, without specifying the exact partner residue [README.md:113](). This is the default output.

### Q3: What do the coarse-grained (CG) levels mean?
**A**: Disobind can predict at different resolutions. `cg 1` is residue-level. `cg 5` and `cg 10` group residues into beads of 5 or 10 residues respectively to reduce computational complexity and noise [README.md:93]().

### Q4: How do I interpret the `Predictions.npy` file?
**A**: This file contains a nested dictionary where keys are the entry IDs from your input and values are the numerical prediction tensors. It is useful for downstream programmatic analysis [README.md:99]().

### Q5: Can I use FASTA files for proteins without UniProt IDs?
**A**: Yes. Use the `-i fasta` flag. The header of the FASTA must follow a specific format to include metadata like offsets and AF2 paths if applicable [README.md:66-72]().

Sources: [README.md:43-114]()

---

## Output Reference

### Prediction Result Flow

```mermaid
flowchart LR
    subgraph "Input Processing"
        CSV["CSV/FASTA Input"]
        UniProt["UniProt Sequence Fetch"]
        AF2["AF2 Structure/PAE (Optional)"]
    end

    subgraph "Inference Engine"
        T5["ProtT5 Embeddings"]
        Model["Epsilon_3 Model"]
    end

    subgraph "Outputs"
        OutCSV["Result CSV"]
        OutNPY["Predictions.npy"]
    end

    CSV --> UniProt
    UniProt --> T5
    T5 --> Model
    AF2 --> Model
    Model --> OutCSV
    Model --> OutNPY
```

Sources: [README.md:47-100](), [run_disobind.py:1-50]() (implied flow)

---