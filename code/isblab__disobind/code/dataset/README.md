# Dataset creation
Move to the `dataset` directory.  

1. Obtain all PDB IDs from databases of disordered proteins including DIBS, MFIB, FuzDB, PBDtot, PDBcdr, DisProt, IDEAL, MobiDB.
```
python 1_disobind_database.py -c 250
```

2. Download all PDB files, obtain relevant info and create binary complexes with the first protein being disordered.
```
python 2_create_database_dataset_files.py -c 250
```

3. Create merged binary complexes.
```
python 3_create_merged_binary_complexes.py -c 200
```

4. Create training and OOD test set.
```
python 4_create_non_redundant_dataset.py -c 100
```

5. Obtain embeddings for the dataset and split into Train:Dev:Test set.
```
python create_input_embeddings.py
```

# Comparison to competing methods 

Create input sequences to alternate methods (AIUPred, MORFchibi, DeepDISOBind).  
```
prep_other_methods_input.py
```
Use the respective web servers for [MORFchibi](https://gsponerlab.msl.ubc.ca/software/morf_chibi/mc2/) and 
[DeepDISOBind](https://www.csuligroup.com/DeepDISOBind/) to obtain predictions using the input FASTA files.  
For AIUPred, clone the Git repo in the disobind directory using:
```
git clone https://github.com/doszilab/AIUPred.git
```

# IDPPI dataset input for Disobind
We used the [IDPPI](https://doi.org/10.1038/s41598-018-28815-x) dataset for assessing Disobind on the PPI prediction task.
Use the following script to create input for running Disobind on the IDPPI dataset:
```
python prep_idppi_input2.py
```

# Case-specific comparisons 

Prepare contact maps from the PDB.
Use this script to create input files (entry_ids and contact maps) compatible with Disobind.  
This script was used to create the Misc dataset in Disobind.
```
prepare_entry_from_pdb.py
``` 

