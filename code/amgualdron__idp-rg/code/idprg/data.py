"""The data tables (data/README.md), the evaluation splits and access to the MD frame store.

Splits, as (fit rows, evaluation rows) into the protein table:
    standard   fit = training split; evaluation = validation split (fold -1), or 3-fold cross-validation inside the
               training split grouped by sequence cluster (fold 0, 1, 2)
    test       fit = training split; evaluation = test split
    length     fit = training proteins with N <= 200; evaluation = training proteins with N >= 300
    sticky     fit = training proteins with contact energy <= its training 75th percentile; evaluation = >= 90th
For length and sticky, fit proteins that share a sequence cluster with an evaluation protein are dropped.
"""
import os
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SPLITS = ("standard", "test", "length", "sticky")
N_FOLDS = 3
FRAMES_ENV = "IDPRG_FRAMES"
FRAMES_SHAPE = (4743, 30000)
NU_F, LN_A_F = 0.57874, -1.56603               # Flory fit ln <Rg> = LN_A_F + NU_F ln N of the training proteins


def proteins():
    return pd.read_csv(DATA / "proteins.csv")


def md_fits():
    return pd.read_csv(DATA / "md_fits.csv")


def benchmark():
    return pd.read_csv(DATA / "benchmark.csv")


def split(table, name, fold=-1):
    """(fit index, evaluation index) into the protein table."""
    tr = np.where(table.split == "train")[0]
    if name == "standard":
        if fold < 0:
            return tr, np.where(table.split == "val")[0]
        cl = np.unique(table.cluster.iloc[tr])
        k = table.cluster.iloc[tr].map(dict(zip(np.random.default_rng(0).permutation(cl), np.arange(len(cl)) % N_FOLDS)))
        k = k.to_numpy()
        return tr[k != fold], tr[k == fold]
    if name == "test":
        return tr, np.where(table.split == "test")[0]
    t = table.iloc[tr]
    if name == "length":
        fit, ev = tr[(t.N <= 200).to_numpy()], tr[(t.N >= 300).to_numpy()]
    elif name == "sticky":
        q75, q90 = t.contact_energy.quantile(.75), t.contact_energy.quantile(.90)
        fit, ev = tr[(t.contact_energy <= q75).to_numpy()], tr[(t.contact_energy >= q90).to_numpy()]
    else:
        raise ValueError(f"unknown split {name!r} (expected one of {SPLITS})")
    fit = fit[~table.cluster.iloc[fit].isin(set(table.cluster.iloc[ev])).to_numpy()]
    return fit, ev


def frames(path=None):
    """The MD frame store, memory-mapped: float32 (4743, 30000), Rg in nm, row = the table's frame_row. Path from the
    argument or the IDPRG_FRAMES environment variable."""
    path = path or os.environ.get(FRAMES_ENV)
    if not path:
        raise FileNotFoundError(f"this step needs the MD frame store: pass --frames PATH or set {FRAMES_ENV}")
    F = np.load(path, mmap_mode="r")
    if F.shape != FRAMES_SHAPE or F.dtype != np.float32:
        raise ValueError(f"{path}: expected a float32 array of shape {FRAMES_SHAPE}, got {F.dtype} {F.shape}")
    return F


def compress(x):
    """Unique positive frame values and their counts (the stored Rg are rounded to 0.001 nm)."""
    x = np.asarray(x, float)
    u, c = np.unique(x[np.isfinite(x) & (x > 0)], return_counts=True)
    return u, c.astype(float)


def packed(F, rows):
    """All frames of the given frame-store rows, compressed per protein and flattened for the jax objectives:
    (ln Rg, count, protein index, number of proteins)."""
    comp = [compress(F[int(r)]) for r in rows]
    lr = np.concatenate([np.log(u) for u, _ in comp])
    c = np.concatenate([c for _, c in comp])
    seg = np.concatenate([np.full(len(u), i) for i, (u, _) in enumerate(comp)])
    return lr, c, seg, len(comp)
