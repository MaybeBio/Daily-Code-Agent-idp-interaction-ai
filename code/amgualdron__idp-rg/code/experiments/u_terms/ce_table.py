"""results/ce_threshold_cv.csv for the adopted form, from the fits in runs/adopt (idprg fit-closed-form --c1 ...)."""
import pathlib
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from idprg import data  # noqa: E402

CL = data.proteins().set_index("id").cluster
A = ROOT / "runs" / "adopt"
tag = lambda c: f"{c:g}".replace("-", "m").replace(".", "p")
cv = lambda c: pd.concat(pd.read_csv(A / f"c1_{tag(c)}" / f"fold{f}" / "per_protein.csv") for f in range(3)).set_index("id")
base = cv(0.0)
rows = []
for c in (-0.5, 0.0, 0.25, 0.5, 0.75, 1.0):
    d = cv(c)
    x = d.nll - base.nll.reindex(d.index)
    g = x.groupby(CL.reindex(x.index)).sum()
    r = dict(c1=c, cv_nll=d.nll.mean(), cv_nll_minus_constant=x.mean(), se=np.sqrt((g ** 2).sum() * len(g) / (len(g) - 1)) / len(x),
             cv_mean_err_pct=d.mean_err_pct.abs().mean())
    for s in ("length", "sticky"):
        p = A / f"c1_{tag(c)}" / s / "per_protein.csv"
        e = pd.read_csv(p) if p.exists() else None
        r[f"{s}_nll"] = e.nll.mean() if e is not None else np.nan
        r[f"{s}_mean_err_pct"] = e.mean_err_pct.abs().mean() if e is not None else np.nan
    rows.append(r)
pd.DataFrame(rows).to_csv(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results" / "ce_threshold_cv.csv", index=False)
