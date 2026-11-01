![header](https://github.com/Floboysky/dvpILC/raw/HEAD/dev_ILC.png)

*Nfil3* during the overview of mouse ILC development (adapted from [Das et al. 2026](https://doi.org/10.1038/s41590-025-02414-1))

---

# NFIL3 triggers innate lymphoid cell specification through gene repression.

This repository contains all the scripts and resources used in a paper currently being published.

This source code enables the generation of graphs and data that help to characterize functional residues, analyze the structure-function relationship, and investigate potential protein partners of the NFIL3 transcription factor, to better understanding its mechanisms of action in the development of murine ILCs.

## Installation

### Requirement

For Windows:
- conda version 24.11.3
- Ubuntu version 24.04

### Steps

In an Anaconda terminal, run the following commands:

```bash
git clone https://github.com/Floboysky/dvpILC.git
conda env create -f environment.yml
conda activate NFIL3_ILC
```

## Description and Usage

The `example_data/` directory contains sample data that can be manipulated by the source code. In all the code examples provided below, the `path` and `name` variables can be changed by the user.

### Conservation Analysis

This chapter describes best practices for using the code files in the `Analysis/` directory.

The pLDDT score is extracted using the `select_plddt_CA.sh` script, which processes the CIF files obtained with [AlphaFold3](https://www.nature.com/articles/s41586-024-07487-w) (AF3) directly. It's used automatically in the files `CA_bfactors_plddt.ipynb` and `conservation_nfil3.ipynb`, and must be located in their parent folder.

- `CA_bfactors_plddt.ipynb` is used to compare and graphically analyze b-factor scores used in crystallography with the pLDDT score used by AF3. Example of input files for the pLDDT study (ouput of AF3) in `example_data/MeCP2_TBL1R/fold_mmecp2_mtbl1r/`, and for the b-factor study ([5NAF](https://www.rcsb.org/structure/5NAF)) in `example_data/MeCP2_TBL1R/Xp/`. The output graphs are in the folder `example_data/MeCP2_TBL1R/Plots/`.

- `conservation.ipynb` enables graphical analysis that aggregates data into a dataframe (in `example_data/NFIL3/Results/`), the results for disorder probability, reliability, alignment quality, and conservation pressure for each residue of a specific protein:
    - The input file for the intrinsic disorder analysis is a JSON file obtained from the [AIUPred](https://academic.oup.com/nar/article/52/W1/W176/7673484) server (`example_data/NFIL3/AIUPred/`).
    - For reliability (pLDDT), an AF3 output was used as an example (`example_data/NFIL3/fold_mnfil3/`).
    - The alignment statistics are obtained from the output files of the [Jalview](https://academic.oup.com/bioinformatics/article/25/9/1189/203460) tool following a global alignment of the protein sequences (`example_data/NFIL3/Jalview/`):

        - The "jalview_files.json" file is generated using the “Files/Output to Textbox” tab for sequence data.
        - The "jalview_annotation.csv" file is generated using the "Files/Export Annotations..." tab for the raw data.

    -The analysis of the holding pressure was obtained using the [EasyCodeML](https://onlinelibrary.wiley.com/doi/10.1002/ece3.5015) program (`example_data/NFIL3/CodeML/`). To run, EasyCodeML requires a protein alignment in FASTA format, the corresponding phylogenetic tree in NWK format, and the mapping of the aligned protein sequences to the genomic sequences (codons) in PAML format (obtained from the [PAL2NAL](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkl315) server). The `clean_fasta.py` script is used to clean the index of the various aligned genomic sequences before they are used in EasyCodeML. For each site extracted from EasyCodeML in the Model 8 (M8), obtained using the “Site Model” option, the script outputs its estimate the omega ratio, uncertainty, and the probability that the site belongs to the class with the highest conservation pressure (class 1).

### Analysis of protein-protein interactions

This chapter describes best practices for using the code contained in the `PPI/` folder.

As with the chapter [“Conservation Analysis”](#conservation-analysis), pLDDT extraction is performed automatically using the `select_plddt_CA.sh` script. It is used by `alphafold3_lis_contact.ipynb` and must be placed in its parent directory.

- `alphafold3_lis_contact.ipynb`: is a modified version of [Kim et al. 2024](http://biorxiv.org/lookup/doi/10.1101/2024.02.19.580970), for the analysis of PAE, iPTM, LIS, cLIS, and LIA scores of a two protein complex obtained by AF3. An example of input data has been placed in the folder `example_data/NFIL3_TBL1R/`. It can also generate a plot of the alpha-carbon pLDDTs for all residues, as well as plots of the LIS and cLIS values for the positive residues (maximum and average values) of each chain in the complex. The program's output is a CSV file containing the residues (filtered or unfiltered) with at least one positive LIS or cLIS value for each protein.

### BioGRID data

- `data_BioGRID.py` can merge CSV data files [BioGRID](https://thebiogrid.org/), to count the number of unique occurrences in a specific key. Example of the merging of TBL1 and TBL1R interactions in the file `data_BioGRID/TBL1X/` and `data_BioGRID/TBL1XR1/`.

### Useful scripts in PyMOL

Folder containing scripts to assist in the analysis of 3D protein models. For direct use in [PyMOL](https://www.pymol.org/).

- `helix_axis.py` is a program that draws a line passing as close as possible to the center of mass of each helix turn (protein or DNA).

- `PyMOL_color.py` is used to color the alpha carbons of proteins using the values contained in the output files of the `conservation.ipynb` (jalview_output.json) and `alphafold3_lis_contact.ipynb` (data_all.csv) scripts.

## Citation

```bibtex

```
