"""Stage `envelopes` (DESIGN §7, A3): an ab initio electron-density envelope of each worked example
from its experimental curve, with DENSS (Grant 2018), for the worked-examples figure; and the
helpers the report uses to fit each example and dock its highest-weight BioEmu-1 conformer into the
averaged map.

DENSS needs numpy < 2 and lives in its own environment (BSC_DENSS_ENV, default /scratch/.venv-denss);
when that environment is absent the reconstruction step prints a notice and is skipped, so `bsc
report` still runs elsewhere with the tracked figure. Per example: `denss-all` with N_MAPS
reconstructions in FAST mode and Dmax from the SASBDB record (pddf_dmax, nm, converted to Å), aligned
and averaged.

Output: results/denss/<label>/ (maps and logs, not tracked) and
results/denss_stats.json (fit statistics, tracked).
"""
from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import requests

from bsc import analyse, saxs, score
from bsc import config as C

DENSS_ENV = Path(os.environ.get("BSC_DENSS_ENV", "/scratch/.venv-denss"))
N_MAPS = 10
MODE = "FAST"
OUT = C.RESULTS / "denss"


def dmax_angstrom(label: str) -> float | None:
    """Dmax from the SASBDB record: the summary written by `bsc fetch`, else the REST API."""
    summ = json.load(open(C.SASBDB_SUMMARY)) if C.SASBDB_SUMMARY.exists() else {}
    v = summ.get(label, {}).get("pddf_dmax")
    if v is None:
        try:
            r = requests.get(C.SASBDB_API.format(label=label), timeout=60)
            r.raise_for_status()
            v = r.json().get("pddf_dmax")
        except (requests.RequestException, ValueError):
            v = None
        if v is not None and C.SASBDB_SUMMARY.exists():
            summ.setdefault(label, {})["pddf_dmax"] = v
            json.dump(summ, open(C.SASBDB_SUMMARY, "w"), indent=1)
    return 10.0 * float(v) if v is not None else None


def write_dat(label: str, folder: Path) -> Path:
    exp = score.read_experiment(label)
    exp = exp[exp["sigma"] > 0]
    p = folder / f"{label}.dat"
    exp.to_csv(p, sep=" ", header=False, index=False, float_format="%.6e")
    return p


def parse_final_log(p: Path) -> dict:
    txt = p.read_text()
    out = {}
    for key, pat in (("final_chi2", r"Final Chi2:\s*([0-9.eE+-]+)"),
                     ("final_rg", r"Final Rg:\s*([0-9.eE+-]+)"),
                     ("support_volume", r"Final Support Volume:\s*([0-9.eE+-]+)")):
        m = re.findall(pat, txt)
        if m:
            out[key] = float(m[-1])
    return out


def parse_average_log(p: Path) -> dict:
    txt = p.read_text()
    out = {}
    m = re.search(r"Resolution = ([0-9.]+) \+- ([0-9.]+)", txt)
    if m:
        out["resolution_A"], out["resolution_sd_A"] = float(m.group(1)), float(m.group(2))
    m = re.search(r"Mean of correlation scores: ([0-9.]+)", txt)
    if m:
        out["mean_correlation"] = float(m.group(1))
    m = re.search(r"Number of aligned maps accepted: (\d+)", txt)
    if m:
        out["maps_accepted"] = int(m.group(1))
    return out


def run_denss(label: str, dmax: float, cores: int) -> dict:
    folder = OUT / label
    folder.mkdir(parents=True, exist_ok=True)
    dat = write_dat(label, folder)
    cmd = [str(DENSS_ENV / "bin" / "denss-all"), "-f", dat.name, "-d", f"{dmax:.1f}", "-m", MODE,
           "-nm", str(N_MAPS), "-j", str(cores), "-o", label, "--plot_off"]
    with open(folder / "denss-all.log", "w") as log:
        subprocess.run(cmd, cwd=folder, stdout=log, stderr=subprocess.STDOUT, check=True)
    maps = folder / label
    logs = [p for p in sorted(maps.glob(f"{label}_[0-9]*.log")) if "_aligned" not in p.name]
    per_map = [parse_final_log(p) for p in logs]
    per_map = [m for m in per_map if "final_chi2" in m]
    rec = {"label": label, "dmax_A": dmax, "n_maps": N_MAPS, "mode": MODE,
           "avg_map": str((maps / f"{label}_avg.mrc").relative_to(C.RESULTS)),
           "chi2_per_map": [m["final_chi2"] for m in per_map],
           "chi2_mean": float(np.mean([m["final_chi2"] for m in per_map])),
           "chi2_median": float(np.median([m["final_chi2"] for m in per_map])),
           "rg_per_map_mean": float(np.mean([m["final_rg"] for m in per_map])),
           "support_volume_mean_A3": float(np.mean([m["support_volume"] for m in per_map]))}
    rec.update(parse_average_log(maps / f"{label}_final.log"))
    return rec


# ---- worked examples: fit at the operating point, docking ----------------------------------------
POROD_A3_PER_DA = 1.7          # protein volume from mass, for the isovalue that encloses the expected volume
RAW_FIT_CHI2 = 2.0             # the operating point of the kinds: the first point on the path with chi2 <= 2
# two density levels of the averaged DENSS map, each the isosurface enclosing a volume: the particle
# volume DENSS assigned (at least OUTER_MIN_FACTOR times the protein volume) and the protein's expected
# volume, which is the level used for docking and for the C-alpha counts
OUTER_MIN_FACTOR = 1.5
CA_BREAK_A = 4.2               # consecutive C-alpha atoms further apart than this: a break in the chain


def ensembles_dir() -> Path:
    return C.DATA / "Predictions" / "PeptoneDB-SAXS-ensembles"


def operating_point(path_chi2: list[float], raw_chi2: float, target: float = RAW_FIT_CHI2) -> int | None:
    """Where the kinds read the reweighting path: -1 when the raw ensemble already fits (chi2 <= target),
    else the index of the first point along the path (decreasing prior strength) with chi2 <= target,
    else None (the target is not reached)."""
    if raw_chi2 <= target:
        return -1
    for i, c in enumerate(path_chi2):
        if c <= target:
            return i
    return None


def example_fit(label: str) -> dict | None:
    """Data, the raw and reweighted BioEmu-1 curves and the AlphaFold2 curve for one entry, each scaled
    to the data. The reweighted curve is taken at the operating point (first point along the path with
    chi2 <= 2), or at the chi2 minimum when that target is never reached; the chi2 minimum is recorded
    either way, and its weights pick the conformer drawn in the envelope."""
    exp = score.read_experiment(label)
    keep = (exp["sigma"] > 0).to_numpy()
    q, I, s = (exp[c].to_numpy()[keep] for c in ("q", "I", "sigma"))
    curves = score.read_prediction(C.MODEL, label)
    af = score.read_prediction("alphafold", label)
    if curves is None:
        return None
    ok = np.isfinite(curves).all(axis=1)
    idx = np.flatnonzero(ok)
    curves = curves[ok][:, keep]
    raw = saxs.ensemble_curve(curves)
    path, ws = saxs.reweighting_curve(curves, I, s, C.THETA_GRID, return_weights=True)
    chi_path = [p["chi2"] for p in path]
    i_min = int(np.argmin(chi_path))
    i_op = operating_point(chi_path, saxs.chi2(raw, I, s))
    w_min = ws[i_min]
    # the reweighted curve shown: the operating point; the chi2 minimum when chi2 <= 2 is never reached
    w_op = (np.full(len(curves), 1.0 / len(curves)) if i_op == -1 else w_min if i_op is None else ws[i_op])
    op, mn = saxs.ensemble_curve(curves, w_op), saxs.ensemble_curve(curves, w_min)

    def rg(curve):
        return saxs.guinier(q, curve, 0.01 * curve + 1e-9, C.GUINIER_QRG_MAX)["rg"]

    out = {"label": label, "q": q, "I": I, "sigma": s, "curves": {},
           "operating_reached": i_op is not None, "operating_is_raw": i_op == -1,
           "phi_operating": saxs.kish_fraction(w_op), "phi_minimum": saxs.kish_fraction(w_min),
           "chi2_minimum": float(chi_path[i_min]),
           "top_conformer": int(idx[int(np.argmax(w_min))]), "top_weight": float(w_min.max()),
           "rg_raw": rg(raw), "rg_operating": rg(op), "rg_minimum": rg(mn),
           "rg_exp": saxs.guinier(q, I, s, C.GUINIER_QRG_MAX)["rg"]}
    for name, m in (("raw", raw), ("operating", op), ("alphafold", af[0][keep] if af is not None else None)):
        if m is None:
            continue
        c, b = saxs.scale_and_background(m, I, s)
        out["curves"][name] = {"fit": c * m + b, "chi2": saxs.chi2(m, I, s)}
    return out


def principal_frame(x: np.ndarray, weights: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Centroid and rotation (columns: principal axes, longest first, right-handed)."""
    wts = np.ones(len(x)) if weights is None else weights
    c = (wts[:, None] * x).sum(0) / wts.sum()
    cov = ((x - c) * wts[:, None]).T @ (x - c) / wts.sum()
    vals, vecs = np.linalg.eigh(cov)
    R = vecs[:, ::-1]
    if np.linalg.det(R) < 0:
        R[:, 2] *= -1
    return c, R


def load_map(mrc_path, upsample: int = 3) -> tuple[np.ndarray, float]:
    """The averaged DENSS map, cubic-upsampled and clipped at zero, with its voxel size in Å."""
    import mrcfile
    from scipy import ndimage
    with mrcfile.open(mrc_path) as m:
        rho = np.asarray(m.data, float)
        vox = float(m.voxel_size.x)
    return np.clip(ndimage.zoom(rho, upsample, order=3), 0, None), vox / upsample


def isovalue(rho: np.ndarray, vox: float, volume_A3: float) -> float:
    """The density level whose isosurface encloses `volume_A3`."""
    frac = min(0.95, volume_A3 / (rho.size * vox**3))
    return float(np.quantile(rho, 1 - frac))


class Docking:
    """A conformer superposed on an envelope by principal axes. The four proper axis-sign choices and
    the mirror image of the envelope (SAXS does not fix handedness) are tried; the one placing most
    C-alpha atoms inside the isosurface is kept. The displayed frame is the envelope's principal-axes
    frame (x longest, then y, then z); the conformer keeps its own handedness."""

    def __init__(self, ca: np.ndarray, rho: np.ndarray, vox: float, iso: float):
        self.rho, self.vox = rho, vox
        inside = np.argwhere(rho > iso) * vox
        self.c_map, self.R_map = principal_frame(inside, rho[rho > iso])
        self.c_ca, self.R_ca = principal_frame(ca)
        x_frame = (ca - self.c_ca) @ self.R_ca
        best = None
        for mirror in (1, -1):
            for sx, sy in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                S = np.diag([sx, sy, sx * sy])
                M = np.array([mirror, 1, 1])
                frac = self.inside_fraction(x_frame @ S, iso, M)
                if best is None or frac > best[0]:
                    best = (frac, S, M)
        self.inside, self.S, self.M = best

    def conformer(self, x: np.ndarray) -> np.ndarray:
        """Conformer coordinates (Å, any atoms) -> displayed frame."""
        return ((x - self.c_ca) @ self.R_ca) @ self.S

    def map_points(self, v: np.ndarray) -> np.ndarray:
        """Map-grid coordinates (Å) -> displayed frame."""
        return ((v - self.c_map) @ self.R_map) * self.M

    def grid(self, x: np.ndarray, M: np.ndarray | None = None) -> np.ndarray:
        M = self.M if M is None else M
        g = ((x * M) @ self.R_map.T + self.c_map) / self.vox
        return np.clip(np.round(g).astype(int), 0, np.array(self.rho.shape) - 1)

    def inside_fraction(self, x_display: np.ndarray, iso: float, M: np.ndarray | None = None) -> float:
        ijk = self.grid(x_display, M)
        return float(np.mean(self.rho[ijk[:, 0], ijk[:, 1], ijk[:, 2]] > iso))


def load_conformer(label: str, conformer: int):
    import mdtraj as md
    d = ensembles_dir() / C.MODEL
    return md.load(str(d / f"{label}.xtc"), top=str(d / f"{label}.pdb"))[conformer]


def density_volumes(volume: float, support: float) -> list[float]:
    """Enclosed volumes of the two density levels, outer first: the particle volume DENSS assigned (at
    least OUTER_MIN_FACTOR times the protein volume) and the protein volume."""
    return [max(support, OUTER_MIN_FACTOR * volume), volume]


def docked_example(label: str, conformer: int, avg_map: Path, volume: float, support: float) -> dict:
    """The conformer docked into the averaged map at the protein-volume level; for each of the two
    density levels the map points above it, and the C-alpha trace, all in the displayed frame; and the
    share of C-alpha atoms inside each level."""
    rho, vox = load_map(avg_map)
    vols = density_volumes(volume, support)
    isos = [isovalue(rho, vox, v) for v in vols]
    frame = load_conformer(label, conformer)
    ca = 10.0 * frame.xyz[0][frame.topology.select("name CA")]
    dock = Docking(ca, rho, vox, isos[1])
    ca_show = dock.conformer(ca)
    surfaces = []
    for name, iso, vol in zip(("particle", "protein"), isos, vols, strict=True):
        points = dock.map_points(np.argwhere(rho > iso) * vox)
        surfaces.append({"name": name, "points": points, "volume_A3": vol,
                         "iso_over_max": iso / float(rho.max()),
                         "ca_inside": dock.inside_fraction(ca_show, iso)})
    return {"ca": ca_show, "inside": dock.inside, "inside_outer": surfaces[0]["ca_inside"], "voxel_A": vox,
            "chain_breaks": int(np.sum(np.linalg.norm(np.diff(ca, axis=0), axis=1) > CA_BREAK_A)),
            "surfaces": surfaces}


def main(cores: int | None = None) -> None:
    stats_path = C.RESULTS / "denss_stats.json"
    stats = json.load(open(stats_path)) if stats_path.exists() else {}
    if not (DENSS_ENV / "bin" / "denss-all").exists():
        print(f"  DENSS environment {DENSS_ENV} not found; reconstructions skipped (set BSC_DENSS_ENV)",
              flush=True)
    else:
        _, d = analyse.load()
        cores = cores or min(C.N_JOBS, 5)
        for cls, label in analyse.choose_examples(d).items():
            if label in stats and (C.RESULTS / stats[label]["avg_map"]).exists():
                print(f"  {cls}: {label} already reconstructed", flush=True)
                continue
            dmax = dmax_angstrom(label)
            if dmax is None:
                print(f"  {cls}: {label} has no Dmax in its SASBDB record; skipped", flush=True)
                continue
            print(f"  {cls}: {label}, Dmax {dmax:.0f} Å, {N_MAPS} maps in {MODE} mode", flush=True)
            rec = run_denss(label, dmax, cores)
            rec["disorder_class"] = cls
            stats[label] = rec
            json.dump(stats, open(stats_path, "w"), indent=1)
            print(f"    mean chi2 of the maps {rec['chi2_mean']:.2f}; "
                  f"resolution {rec.get('resolution_A', float('nan')):.0f} Å", flush=True)


if __name__ == "__main__":
    main()
