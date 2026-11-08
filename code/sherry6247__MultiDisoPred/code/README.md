<!--
 * @Date: 2026-11-01 14:21:12
 * @LastEditors: Sicen Liu
 * @LastEditTime: 2026-11-01 14:39:13
 * @FilePath: /liusicen/mygithub/MultiDisoPred/README.md
 * @Description:  
 * @Copyright: © 2025 liusicen@smbu.edu.cn. All rights reserved.
-->
# MultiDisoPred
MultiDisoPred is a framework for residue-level prediction of intrinsically disordered regions (IDRs) and protein-binding sites.
![framework](https://github.com/sherry6247/MultiDisoPred/raw/HEAD/src/MultiDisoPred.png)

## Requirements

All experimental environments and dependencies are specified in `./environment.yml`. The environment can be created using the following command:

```bash
conda env create -f environment.yml
conda activate MultiDisoPred
```

## Datasets

All FASTA files follow a three-line-per-record format: `>id`, the amino acid sequence, and the corresponding residue-level binary labels.

### Training dataset（`./datasets`）

|file | description |
| --- | --- |
| `DM4845_training.fasta` | Complete training dataset for intrinsically disordered region (IDR) prediction (4,845 proteins). |
| `DM4845_training_nr25_*.fasta` | Training subsets derived from DM4845 after redundancy reduction against the corresponding test sets at a 25% sequence identity threshold.|
| `TR552_binding.fasta` | Training dataset for binding site prediction (TR552). |
| `TR552_binding_CAID2_nr25.fasta` / `TR552_binding_CAID3_nr25.fasta` | Redundancy-reduced version of TR552. |
<!-- | `lipid_training.fasta` | lipid training dataset| -->

### Test dataset（`./datasets/test_datasets`）

| file | description |
| --- | --- |
| `CAID2_disorder_pdb.fasta` / `CAID3_disorder_pdb.fasta` | disorder test datasets |
| `CAID2_binding.fasta` / `CAID3_binding.fasta` | binding test datasets |
| `MXD494.fasta` / `SL329.fasta` / `DISORDER723_test.fasta` | independent IDR test dataset |
<!-- | `CASP.fasta` / `lipid_test.fasta` | CASP / lipid 测试集 | -->

## Usage

### Training

```bash
python main.py
```

### Test

```bash
python predict_main.py
```

## Project Structure
- `main.py`: Entry point for model training
- `predict_main.py`: Entry point for model testing and inference
- `src/`: Model implementations, trainers, evaluation procedures, and metrics
- `utils/`: Parameter configuration, data processing, and feature encoding
- `datasets/`: Training and test datasets (see above)
- `environment.yml`: Conda environment configuration file
