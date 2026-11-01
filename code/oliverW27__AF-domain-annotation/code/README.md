# AF domain annotation
AlphaFold3 confidence plots annotated with InterProScan domains.

## Quick start
`plot_af3_domains.py` draws the predicted aligned error (PAE) and per-residue confidence (pLDDT) of AlphaFold3 models, with protein domains, coiled coils and disordered regions from InterProScan drawn on top. It works for single proteins and for complexes, where all chains are shown side by side. This makes it easy to see whether the annotated domains match what AlphaFold confidently built, and how confidently domains and chains are placed relative to each other.

You need:
- Python 3 with numpy and matplotlib
- AlphaFold Server downloads, unzipped (see [Input](#input))
- InterProScan TSV output for every protein in them
- This codebase

```bash
python plot_af3_domains.py --json-dir AF3_output --tsv interpro_input
```

Two graphs per model are written to `graphs/`, next to the script. See [Output](#output) for details.

## Input
### AlphaFold3 models (`--json-dir`)
A folder of unzipped [AlphaFold Server](https://alphafoldserver.com/) downloads, must contain the `full_data` JSON and its `model` CIF in the same folder. 

```
AF3_output/
├── fold_is26_tpase_dimer/
│   ├── fold_is26_tpase_dimer_full_data_0.json
│   ├── fold_is26_tpase_dimer_model_0.cif
│   └── ...
└── fold_notch1/
    └── ...
```

Only protein chains are plotted. DNA, RNA, ligands etc. are left out of both graphs.

### InterProScan output (`--tsv`)
One or more InterProScan TSV files, or folders of them (every `.tsv` inside is read). Proteins can be split over files however you like. Run InterProScan from the [InterPro website](https://www.ebi.ac.uk/interpro/search/sequence/) and download the TSV output.

## Options
There are two options, both required: `--json-dir` (folder of AlphaFold Server downloads) and `--tsv` (InterProScan TSV files or folders). Other settings are variables at the top of `plot_af3_domains.py`: `MAX_EVALUE` (1e-3), `TIERS` (databases used, in order of preference), `MAX_OVERLAP_AA` and `MAX_OVERLAP_FRAC` (10 aa and 25%), `PALETTE` (domain colours) and `OUT_DIR` (`graphs/`).

## How domains are chosen
InterProScan reports many overlapping hits for the same region, so one set of non-overlapping domains is chosen per protein from hits with an e-value ≤ `MAX_EVALUE`, taking Pfam first, then SMART and CDD, then SUPERFAMILY and CATH-Gene3D, then family-level databases. Each hit is kept only if it overlaps the domains already chosen by at most 10 aa and 25% of the shorter one.

## Output
Everything is written to `graphs/` next to the script. `NAME_pae.png` shows the PAE matrix with domains outlined on the diagonal and colour strips along both axes. In multimers, red lines separate the chains. `NAME_plddt.png` shows the pLDDT of each residue's C-alpha, coloured by AlphaFold's confidence bands, above a cartoon of the domains, coiled coils and disordered regions.

`NAME` is the InterProScan protein ID for a single-chain model, and the AlphaFold job name for a complex.

## Examples

### Notch1
<p align="center">
  <img src="https://github.com/oliverW27/AF-domain-annotation/raw/HEAD/graphs/Notch1_pae.png" width="520" alt="Notch1 PAE">
  <img src="https://github.com/oliverW27/AF-domain-annotation/raw/HEAD/graphs/Notch1_plddt.png" width="100%" alt="Notch1 pLDDT">
</p>

A very long (2,555 aa) multi-domain protein. The EGF repeats are only confident locally (along the diagonal), the ankyrin repeats and the negative regulatory region (LNR and NOTCH domains) form confident blocks. The C-terminal ~450 residues have very low pLDDT, and MobiDB-lite independently predicts much of this region as disordered.

### IS26 transposase dimer
<p align="center">
  <img src="https://github.com/oliverW27/AF-domain-annotation/raw/HEAD/graphs/is26_tpase_dimer_pae.png" width="520" alt="IS26 dimer PAE">
  <img src="https://github.com/oliverW27/AF-domain-annotation/raw/HEAD/graphs/is26_tpase_dimer_plddt.png" width="100%" alt="IS26 dimer pLDDT">
</p>

A homodimer. Both copies match the same InterProScan entry, and the red lines split the PAE into the two chains and the regions between them.

### DNA polymerase I
<p align="center">
  <img src="https://github.com/oliverW27/AF-domain-annotation/raw/HEAD/graphs/DNA_pol_I_pae.png" width="520" alt="DNA polymerase I PAE">
  <img src="https://github.com/oliverW27/AF-domain-annotation/raw/HEAD/graphs/DNA_pol_I_plddt.png" width="100%" alt="DNA polymerase I pLDDT">
</p>

The N-terminal 5'-3' exonuclease (residues ~1–300) is confident on its own but its position relative to the rest of the protein is not, and the linker between them has very low pLDDT. This is the flexible hinge at which the Klenow fragment is cut off.

The example graphs are derived from AlphaFold Server output and are subject to the [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms).

## Citations:
- Abramson J, et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500. https://doi.org/10.1038/s41586-024-07487-w
- Blum M, et al. (2025) InterPro: the protein sequence classification resource in 2025. *Nucleic Acids Research* 53:D444–D456. https://doi.org/10.1093/nar/gkae1082
