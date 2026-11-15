# Data

This directory contains the training, validation, and test datasets used for the development and evaluation of DARUMA-binding and its associated models.

## Directory structure

- `training/`: Training datasets.
- `validation/`: Validation datasets.
- `test/`: Test datasets used for evaluation.

## Training data

- `train_binding.txt`: Training data for the DARUMA-binding binding-site prediction model.
- `train_plddt.txt`: Training data for the DAF-pLDDT model.
- `train_rsa.txt`: Training data for the DAF-RSA model.
- `train_ss.txt`: Training data for the DAF-SS model.

The large training datasets are archived separately on Zenodo:

https://doi.org/10.5281/zenodo.24127965

## Validation data

- `val_binding.txt`: Validation data for the DARUMA-binding binding-site prediction model.
- `val_plddt.txt`: Validation data for the DAF-pLDDT model.
- `val_rsa.txt`: Validation data for the DAF-RSA model.
- `val_ss.txt`: Validation data for the DAF-SS model.

## Test data

- `caid2_binding.fasta`: Test sequences for the CAID 2 Binding challenge.
- `caid3_binding.fasta`: Test sequences for the CAID 3 Binding challenge.
- `caid3_binding_idr.fasta`: Test sequences for the CAID 3 Binding-IDR challenge.

These test datasets were obtained from the CAID Web Portal.

Detailed procedures for dataset construction, preprocessing, redundancy reduction, and dataset splitting are described in the associated manuscript.