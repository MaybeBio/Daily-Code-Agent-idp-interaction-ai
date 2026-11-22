# DARUMA-binding

DARUMA-binding is a prediction tool for identifying binding sites in intrinsically disordered regions (IDRs) of proteins.

## Installation

Clone this repository:

```bash
git clone https://github.com/antepontem/DARUMA-binding.git
cd DARUMA-binding
```

PyTorch must be installed separately according to your operating system and CUDA environment.

After installing PyTorch, install the remaining dependencies:

```bash
pip install -r requirements.txt
```

## Usage

DARUMA-binding takes a FASTA file containing one or more amino acid sequences and predicts residue-level binding sites.

Example:

```bash
python src/daruma_binding/predict.py -i input.fasta -o output.txt
```

Arguments:

- `-i`: Input FASTA file containing one or more amino acid sequences.
- `-o`: Output file path.

See the `examples/` directory for example input files and execution examples.

## Input

Input sequences must be provided in FASTA format.

A FASTA file may contain one or more amino acid sequences.

## Output

For each input sequence, DARUMA-binding outputs the following information:

1. Accession / FASTA header
2. Amino acid sequence
3. Binding probability for each residue
4. Binary binding prediction for each residue (`1`: binding, `0`: non-binding; threshold = 0.26)
5. Predicted pLDDT
6. Predicted relative solvent accessibility (RSA)
7. Predicted three-state secondary structure (SS3; `H`: helix, `S`: strand, `C`: coil)
8. Predicted eight-state secondary structure (SS8)

In addition to binding predictions, DARUMA-binding outputs predicted structural features used by the model.

## Repository Structure

- `src/`: Source code for DARUMA-binding.
- `examples/`: Example inputs and execution scripts.
- `data/`: Validation and test datasets. Large training datasets are available separately from Zenodo.
- `reproduction/`: Results and scripts for reproducing the figures and tables in the paper.
- `caid/`: Code submitted to CAID 4.

## Reproduction

Scripts and processed results used to reproduce the figures and tables reported in the manuscript are available in `reproduction/`.

Each script corresponds to a figure or table in the manuscript and can be run using the processed results provided in the repository.

See `reproduction/README.md` for details.

## Data Availability

The large training datasets used in this study are available from Zenodo:

https://doi.org/10.5281/zenodo.24127965

Validation and test datasets used for evaluation are included in this repository where applicable.

## CAID 4 Submission

The code submitted to CAID 4 was previously located in the root directory of this repository and has been moved to [`caid/`](https://github.com/antepontem/DARUMA-binding/blob/HEAD/caid).

## Web Server

DARUMA-binding is available at:

https://antepontem.org/daruma_binding/

## Tested Environment

DARUMA-binding was tested in the following environment:

- Ubuntu 18.04.5 LTS
- Python 3.10.19
- PyTorch 1.12.1+cu102
- CUDA 10.2

PyTorch should be installed according to the user's operating system, GPU, and CUDA environment.

## License

See [`LICENSE`](https://github.com/antepontem/DARUMA-binding/blob/HEAD/LICENSE).

## Citation

The version of DARUMA-binding corresponding to the manuscript is archived on Zenodo:

https://doi.org/10.5281/zenodo.24272574
