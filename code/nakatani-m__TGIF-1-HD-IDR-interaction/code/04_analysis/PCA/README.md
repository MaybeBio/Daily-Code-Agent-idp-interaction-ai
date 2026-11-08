# PCA Analysis

Scripts for principal component analysis (PCA) of VcMD trajectories.

## Overview

PCA was performed using pairwise distance coordinates between HD residues and IDR residues.

Each frame was weighted according to the canonical distribution (`prob_mix.txt`).

The analysis was executed using:

- `pca.bash`
- `mda_pca_traj_selectatoms.py`

## Input Files

### `P_nonP_ca.trr`

Trajectory file used for PCA.

### `run_ca.pdb`

Reference structure.

### `prob_mix.txt`

Reweighted probabilities used as frame weights.

`P_nonP_ca.trr` and `prob_mix.txt` are not included in this repository because of their large file sizes.

## Scripts

### `pca.bash`

Main execution script.

### `mda_pca_traj_selectatoms.py`

PCA calculation script.

This script was originally developed by a former member of the laboratory.

This script:

- reads trajectory data
- applies probability-based weighting
- calculates principal components
- outputs PCA coordinates and eigenvectors

## External Dependencies

- kkstruct.py
- kkpdb.py
- kkmdconf.py

Available from:

https://github.com/kotakasahara/omegagene/tree/master/toolkit

## Workflow

```text
P_nonP_ca.trr
run_ca.pdb
prob_mix.txt
      ↓
pca.bash
      ↓
mda_pca_traj_selectatoms.py
      ↓
PCA coordinates
```
