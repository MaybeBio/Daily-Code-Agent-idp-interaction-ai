# IDP MD-SAXS benchmarking

<img width="6112" height="2812" alt="H4SAXSMDgraphic" src="https://github.com/user-attachments/assets/328c5a72-924c-43e9-93a6-d43f74d66d93" />

**Overview** <br />
Repository containing tools to compare intrinsically disordered protein (IDP) structural ensembles generated from molecular dynamics (MD) with small-angle X-ray scattering (SAXS) data. <br />
An example case is provided of the H4 histone tail, as part of the paper below. If you happen to use the code in this repository, please cite this paper. <br />
**Tong Wang^, Egor S. Kolesnikov^, Lois Pollack, Alexey V. Onufriev. Combining atomistic simulation and solution scattering to reveal charge-driven histone H4 tail conformational changes. (2026, in submission)** <br />

**Contents** <br /> 
Regarding the ensemble optimization method (EOM): To calculate SAXS profiles from SAXS data, please consult the ATSAS manual, specifically GAJOE, which uses CRYSOL: https://www.embl-hamburg.de/biosaxs/manuals/atsas-2.4/eom.html <br /> 
The folder ExampleData contains an example GAJOE results output folder. <br />
Experimental SAXS profiles for wild type and K-to-Q (charge neutralized) H4 histone tail are provided in the folder SAXSData. <br />

Once GAJOE results are generated, you may run make2Dplot_fromEOMresults.m to generate a 2D histogram comparing a full MD ensemble with a SAXS-refined ensemble, as well as calculate the Jensen-Shannon divergence between the two distributions. <br />
You may also run plotEOMSAXScomparison.m to generate an averaged SAXS profile from an EOM/GAJOE result and plot and compare to experimental data. <br /> 
Both of these scripts are run by interacting with the "User-defined parameters" in the first code block. Mostly the paths just need to be matched. <br /> 

**Hardware Notes** <br />
The analysis should only require a standard computer that is able to run the commonly available software below, with enough memory to complete moderate intensity computing tasks. <br />

**Software Notes** <br />
EOM and GAJOE require the ATSAS SAXS modelling software. Using GAJOE v2.1 will be easiest. <br />
Molecular dynamics simulations are run using AMBER 24 <br />
MATLAB R2021 or later may be required <br />

**License** <br />
This project is covered under the Apache 2.0 License.

# EOM Tutorial (ATSAS)

This tutorial demonstrates the EOM workflow using 20 conformations of the H4 histone tail. For comprehensive analysis please consult the ATSAS manual.

After installing [ATSAS](https://www.embl-hamburg.de/biosaxs/software.html), open the directory `Example_tutorial` in a terminal.

The directory contains an experimental SAXS dataset (`H4WT_150Na.dat`), 20 example PDB structures, two Python plotting scripts (`graph.py`, `2d.py`), and two input files (`assignment.txt`, `protein.seq`) required by RANCH.

### 1. RANCH — Generate conformations

```bash
mkdir pool_20

ranch --repetitions 20 --prefix pool_20/model_ \
    --seed 12345 assignment.txt protein.seq
```

**Output:** `pool_20/` containing 20 randomly generated conformations in CIF format.

### 2. FFMAKER — Calculate SAXS profiles

```bash
ffmaker --output intensities_20.csv \
    --pool size_statistics_20.csv pool_20/*.cif
```

**Output:**
- `intensities_20.csv` — theoretical SAXS intensity profiles for all conformations.
- `size_statistics_20.csv` — structural properties of each conformation: \(R_g\), \(D_{\max}\), end-to-end distance (\(R_{ee}\)), and volume.

To analyze the provided PDB structures instead, replace `pool_20/*.cif` with the path to those PDB files.

### 3. GAJOE — Fit SAXS data

```bash
gajoe \
    --intensity intensities_20.csv \
    --size size_statistics_20.csv \
    --times 1 \
    --generation 100 \
    --ensemble 20 \
    --minimum 1 \
    --maximum 10 \
    --no-repeated \
    H4WT_150Na.dat
```

**GAJOE parameters:**

- `--intensity intensities_20.csv` — theoretical SAXS profiles calculated by FFMAKER.
- `--size size_statistics_20.csv` — structural properties of the conformations (\(R_g\), \(D_{\max}\), \(R_{ee}\), etc.).
- `--times 1` — number of independent genetic algorithm runs.
- `--generation 100` — number of generations per run.
- `--ensemble 20` — number of candidate ensembles evaluated in each generation.
- `--minimum 1` — minimum number of conformations in the selected ensemble.
- `--maximum 10` — maximum number of conformations in the selected ensemble.
- `--no-repeated` — prevents the same conformation from appearing multiple times within an ensemble.
- `H4WT_150Na.dat` — experimental SAXS data used for fitting.

**Note:** These settings are intended for a quick demonstration. For quantitative analysis, increase the number of generations and independent runs, and use a sufficiently large conformational pool.

**Output:** A directory named `GA001`, `GA002`, etc., containing:

- `profiles_*.fit` — experimental SAXS data and fitted ensemble intensity, including \(\chi^2\).
- `Rg_distr_*.txt` — \(R_g\) distributions of the full pool and selected ensemble.
- `Size_distr_*.txt` — corresponding \(D_{\max}\) distributions.
- `CaCa_distr_*.txt` — end-to-end distance distributions.
- `pdbs/` — conformations selected by GAJOE.
- `logFile_*.log` — optimization settings and results.

### 4. Plot the results

Install the required Python packages:

```bash
python -m pip install numpy matplotlib
```

Set the GAJOE output directory in both scripts, then run:

```bash
python graph.py
python 2d.py
```

- `graph.py` — experimental SAXS curve and EOM fit.
- `2d.py` — \(R_g\) versus \(R_{ee}\) scatter plot with selected conformations and marginal histograms.

**Note:** This small pool is intended only to demonstrate the workflow. Quantitative ensemble analysis requires substantially larger, adequately sampled conformational pools.
