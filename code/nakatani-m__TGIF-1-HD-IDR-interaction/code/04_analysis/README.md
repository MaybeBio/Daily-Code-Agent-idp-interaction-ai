# 04_analysis

This directory contains scripts used for analysis of mD-VcMD simulation trajectories and for preparation of subsequent mD-VcMD iterations.

---

## Analysis Workflow

The analysis workflow performed in this study was:

```text
mD-VcMD trajectories
        ↓
PCA
        ↓
DBSCAN clustering
        ↓
Representative structure extraction
```

---

## Directory Structure

```text
04_analysis/
├── PCA/
├── Qcano/
└── DBSCAN/
```

---

## PCA/

Principal component analysis (PCA) based on HD-IDR distance coordinates.

### Files

#### mda_pca_traj_selectatoms.py

Python script used to calculate HD-IDR distance coordinates from simulation trajectories.

#### pca.bash

Execution script used for PCA calculations.

#### run_ca.pdb

Reference structure used during PCA analysis.

### Description

The PCA calculations were performed using HD-IDR intermolecular distance coordinates as descriptors.

Part of the workflow was adapted from:

- [Omegagene Toolkit](https://github.com/kotakasahara/omegagene/tree/master/toolkit)

The generated PCA coordinates were subsequently used for DBSCAN clustering.

---

## Qcano/

Canonical distribution reconstruction and parameter updates used between mD-VcMD iterations.

### Files

#### merge_vcmd_qraw.py

Merges distributions obtained from all replicas.

#### power_params.py

Updates VcMD parameters for subsequent iterations.

#### prep_next_qraw.bash

Prepares input files required for Qcano calculations.

#### output/

Representative output files generated during Qcano calculations.

### Description

Qcano analysis was performed after each mD-VcMD iteration.

The reconstructed canonical distribution was used for:

- updating VcMD parameters
- assigning initial structures for the next iteration
- improving sampling efficiency in subsequent simulations

This workflow was performed using:

- [Omegagene Toolkit](https://github.com/kotakasahara/omegagene/tree/master/toolkit)
- [Omegagene Qcano (qcano_szfit)](https://github.com/kotakasahara/omegagene/tree/master/vcmd/qcano_szfit)

---

## DBSCAN/

Density-based clustering of PCA coordinates.

### pIDR-HD/

DBSCAN analysis for the phosphorylated model.

Files include:

- `dbscan.py`

Input:

- PCA coordinates generated from the PCA workflow (`pca_traj_P.txt`)

Representative parameters:

```python
eps = 0.013
min_samples = 650
```

---

### non-pIDR-HD/

DBSCAN analysis for the non-phosphorylated model.

Files include:

- `dbscan.py`

Input:

- PCA coordinates generated from the PCA workflow (`pca_traj_nonP.txt`)

Representative parameters:

```python
eps = 0.013
min_samples = 400
```

---

### Clustering Procedure

The DBSCAN workflow consisted of:

```text
PCA coordinates
        ↓
z-score normalization
        ↓
DBSCAN clustering
        ↓
Cluster relabeling by population
        ↓
Cluster statistics
        ↓
Visualization
```

Cluster labels were reassigned according to cluster population:

- C1: largest cluster
- C2: second largest cluster
- C3: third largest cluster
- ...

Noise points remain labeled as:

```text
-1
```

Representative structures of the major clusters were subsequently extracted and analyzed.

---

## Notes

The actual simulations consisted of 1000 independent replicas for each system and iteration.

For reproducibility, representative analysis scripts are provided in this repository.

The PCA coordinate files used as inputs for DBSCAN clustering analyses are available in the Zenodo dataset referenced in the repository root README.

Representative cluster structures extracted from these analyses are provided in:

```text
05_clustered_structure/
```
