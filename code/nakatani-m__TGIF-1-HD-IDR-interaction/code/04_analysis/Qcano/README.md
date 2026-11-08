# Qcano Reconstruction

Scripts for reconstructing the canonical distribution (Qcano) and generating the parameter set for the next VcMD iteration.

## Overview

This workflow combines qraw data obtained from multiple VcMD iterations, reconstructs the canonical distribution, and updates VcMD parameters for the next iteration.

Older iterations are assigned exponentially decaying weights:

*w_i* = 0.95^(*N*-*i*)

where:

- *N*: final iteration number
- *i*: iteration number
- *w_i*: weight assigned to iteration *i*

---

## Scripts

### `prep_next_qraw.bash`

Main workflow script.

Functions:

1. Collects `qraw_is.dat` and `qraw.dat` from multiple iterations.
2. Generates iteration-dependent weights (`qrawis_weight.txt`).
3. Reconstructs the canonical distribution.
4. Generates the parameter files for the next iteration.

---

### `merge_vcmd_qraw.py`

Merges qraw data and generates updated VcMD parameters.
Modified version of the original merge_vcmd_qraw.py script.

The original script is available from:
https://github.com/kotakasahara/omegagene/tree/master/toolkit

This version contains modifications used in the present study.

Outputs:

- `params_next.dat`
- `vcmd_next_qraw.dat`

Dependencies:

- kkstruct.py
- kkpdb.py
- kkmdconf.py

---

### `power_params.py`

Applies an additional weighting correction to the VcMD parameters.

Output:

- `params_next_power.dat`

The file `params_next_power.dat` is used as the input parameter file for the next VcMD iteration.

---

## Output Files

### `qrawis_list.txt`

List of all `qraw_is.dat` files used for Qcano reconstruction.

### `qraw_list.txt`

List of all `qraw.dat` files used for parameter updates.

### `qrawis_weight.txt`

Weights assigned to each iteration.

### `qraw_sum.dat`

Merged qraw distribution used for Qcano fitting.

### `qw_opt.dat`

Optimized Qcano distribution obtained by `qcano_szfit`.

### `vcmd_next_qraw.dat`

Merged qraw information for parameter updates.

### `params_next.dat`

VcMD parameter file generated from the merged qraw distributions.

### `params_next_power.dat`

Final VcMD parameter file after power-law correction.

**This file is used as the input for the next iteration.**

### `param_changes.log`

Log file summarizing parameter modifications.

---

## Notes

The files in the `output/` directory were generated from the phosphorylated model analyzed in this study.

The workflow itself can be applied to other systems by modifying the input files and iteration settings.

The iteration numbers and file paths in `prep_next_qraw.bash` correspond to the phosphorylated model used in this work.

---

## External Dependencies

The following external scripts are required but are not included in this repository.

### OmegaGene toolkit

- `vcmd_sum_qraw_is.py`
- `kkstruct.py`
- `kkpdb.py`
- `kkmdconf.py`

Available from:

https://github.com/kotakasahara/omegagene/tree/master/toolkit

### qcano_szfit

- `qcano_szfit`

Available from:

https://github.com/kotakasahara/omegagene/tree/master/vcmd/qcano_szfit

---

## Workflow

```
qraw_is.dat / qraw.dat
        ↓
prep_next_qraw.bash
        ↓
vcmd_sum_qraw_is.py
        ↓
qcano_szfit
        ↓
qw_opt.dat
        ↓
merge_vcmd_qraw.py
        ↓
params_next.dat
        ↓
power_params.py
        ↓
params_next_power.dat
        ↓
Next VcMD iteration
```
