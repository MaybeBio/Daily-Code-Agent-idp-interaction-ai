"""Paths, data sources and the fixed analysis settings (see docs/DESIGN.md)."""
import os
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("BSC_ROOT", REPO))
DATA = ROOT / "data"
RESULTS = ROOT / "results"

# Zenodo record of PeptoneBench (Invernizzi et al. 2025), CC-BY 4.0
ZENODO_RECORD = "17306061"
ZENODO_FILES = {
    "PeptoneDBs.tar.gz": "experimental SAXS profiles, sequences, pH and disorder scores",
    "Predictions.tar.gz": "generated ensembles and back-calculated curves for every model and entry",
}
# sha256 of each archive is recorded by `bsc fetch` into data/checksums.json on first download and
# verified on every later run; the Zenodo API's own md5 is checked at download time.
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"

SAXS_TABLE = DATA / "PeptoneDB-SAXS" / "PeptoneDB-SAXS.csv"
SASBDB_SUMMARY = DATA / "sasbdb_summary.json"      # entry summaries from the SASBDB REST API (bsc fetch)
SASBDB_API = "https://www.sasbdb.org/rest-api/entry/summary/{label}/?format=json"
OLIGOMER_MW_RATIO = 1.6          # measured mass over sequence mass above this: not a monomer in solution
SAXS_DIR = DATA / "PeptoneDB-SAXS" / "sasbdb-clean_data"
PRED_DIR = DATA / "Predictions" / "PeptoneDB-SAXS-expt"   # back-calculated curves per model

MODEL = "bioemu"                 # the model under study, as named in the archive
COMPARATORS = ["alphafold", "esmfold", "boltz2", "idpfold", "idpsam", "peptron"]   # as named in the archive
PREDICTOR = "Pepsi"              # forward model that produced the back-calculated curves

# disorder classes from the mean per-residue disorder score supplied with the table
DISORDER_CLASSES = [("folded", 0.0, 0.2), ("partly disordered", 0.2, 0.6), ("disordered", 0.6, 1.01)]
LENGTH_BINS = [0, 100, 200, 350, 10_000]

# reweighting: maximum-entropy weights w ∝ exp(-lambda·x) fit by minimising chi2/2 + theta·KL
THETA_GRID = 10.0 ** np.linspace(-2, 4, 25)
ESS_TARGET_BENCHMARK = 10        # PeptoneBench's fixed effective-sample-size target, for comparison
CHI2_TARGET = 1.0                # a fit at the experimental noise level

GUINIER_QRG_MAX = 1.3
AGGREGATION_UPTURN = 0.10        # relative rise of I(q) above the Guinier line at the lowest q
# a point below -3 sigma marks an over-subtracted buffer; noise around zero at high angle does not
NEGATIVE_SIGMA = 3.0

SEED = 42
N_JOBS = int(os.environ.get("BSC_N_JOBS", "4"))
