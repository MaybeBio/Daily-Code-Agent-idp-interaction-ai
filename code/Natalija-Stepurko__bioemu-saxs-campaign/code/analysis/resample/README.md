# BioEmu v1.2 re-sample against the archive's ensembles

The design, written before any result, is [`docs/RESAMPLE_DESIGN.md`](https://github.com/Natalija-Stepurko/bioemu-saxs-campaign/blob/HEAD/docs/RESAMPLE_DESIGN.md): 100 samples
per protein from bioemu-v1.2 (notebook `notebooks/bioemu_resample.ipynb`, Colab T4) for the 51 proteins of
`notebooks/resample_list.csv`, compared with the archive's BioEmu ensembles on the same proteins by the Cα radius of
gyration against the measured Guinier radius, and by the raw χ² of the ensemble-average curve.

| Script | Environment | Does |
|---|---|---|
| `01_inventory.py` | main | per protein: npz batches and samples, the xtc files the notebook wrote and their frames, progress.json, class and measured Rg; extracts missing archive ensembles from `Predictions.tar.gz` in one pass (`--archive`) |
| `02_convert.py` | conversion (`requirements.txt`) | every sample recovered from the npz batches with BioEmu's own `save_pdb_and_xtc` (`filter_samples=False`); the notebook's conversion reproduced and checked frame by frame; BioEmu's physical filter per conformer on both sides |
| `03_rg.py` | main | Cα radius of gyration of every conformer, both sides |
| `04_saxs.py` | main | per conformer, both sides: PULCHRA back-mapping from the Cα trace, then Pepsi-SAXS with the measured profile; cached per conformer, resumable, three workers on three cores |
| `05_chi2.py` | main | raw χ² and NRMSD by the study's scoring for paths A (v1.2), B (archive, same steps) and C (archive's own curves), over all conformers and over those passing the filter; path C checked against `results/scores.csv` |
| `06_analyse.py` | main | per-protein table, paired statistics by class, E12 and E13, rejected shares, the figure |

`common.py` holds paths and settings; `rslib.py` the functions the scripts share (sample counting, RMSD, Rg, the raw
fit, paired statistics), tested in `tests/test_resample.py`. `run_all.sh` runs the steps in order.

Tracked outputs in `results/resample/`: `inventory.csv`, `conversion.csv`, `filter.csv` (per conformer and side),
`rg_conformers.csv`, `chi2.csv`, `curve_order.csv`, `per_protein.csv`, `summary.json` (every reported number) and
`versions.json` (software, inputs, wall times and code commit of each step). Recovered conformers and per-conformer
curves stay in `RESAMPLE_WORK`. The figure is `results/figures/fig_resample_rg.png`.

Pepsi-SAXS is distributed for academic use only and is not redistributed here; `analysis/simulations/02_fetch_sources.py`
downloads it and PULCHRA.
