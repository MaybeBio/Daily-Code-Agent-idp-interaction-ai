"""Functions shared by the simulation scripts (docs/SIMULATION_DESIGN.md): the input table (proteins,
sequences, conditions from the SASBDB record, run settings), the radius of gyration of simulated frames,
the convergence check and the ratio summaries. numpy and pandas only, so the tests run in the main
environment."""

from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

# ------------------------------------------------------------------------------------------- settings
DEFAULT_TEMP_K = 293.15
DEFAULT_PH = 7.0
DEFAULT_IONIC_M = 0.15
MIN_IONIC_M = 0.01
PBS_IONIC_M = 0.17  # Dulbecco's PBS: 137 mM NaCl, 2.7 mM KCl, 10 mM Na2HPO4, 1.8 mM KH2PO4
N_FRAMES = 1010
N_DISCARD = 10
SEED_BASE = 20261015


def run_settings(n: int) -> dict:
    """The single-chain protocol of the IDRome study (Tesei et al. 2024), as in the design §4: cubic box side
    (N - 1) * 0.38 + 4 nm; a frame every 7,000 steps for N < 150, else every ceil(3e-4 N^2) * 1,000 steps;
    1,010 frames of which the first 10 are discarded; 10 fs per step."""
    n_save = 7000 if n < 150 else int(math.ceil(3e-4 * n * n) * 1000)
    return {
        "box_nm": round((n - 1) * 0.38 + 4, 3),
        "n_save": n_save,
        "n_frames": N_FRAMES,
        "n_discard": N_DISCARD,
        "steps": N_FRAMES * n_save,
        "ns": N_FRAMES * n_save * 1e-5,
    }


# ------------------------------------------------------------------------------------------- buffer rule
_NUM = r"(?P<c>\d+(?:\.\d+)?)\s*(?P<u>mM|M|µM|μM|uM)?\s*"
# species that need a concentration in front of them; (?i:...) keeps the unit case-sensitive (M against mM)
_SPECIES = {
    "NaCl": r"(?i:nacl)\b",
    "KCl": r"(?i:kcl)\b",
    "NaNO3": r"(?i:nano3)\b",
    "LiCl": r"(?i:licl)\b",
    "CaCl2": r"(?i:cacl2)\b",
    "MgCl2": r"(?i:mgcl2)\b",
    "Na2SO4": r"(?i:na2so4)\b",
    "(NH4)2SO4": r"(?i:\(nh4\)2so4)",
    "urea": r"(?i:urea)\b",
    "GdmCl": r"(?i:gdmcl|gdnhcl|guhcl|gdn-?hcl|guanidin\w*[\s-]*(?:hydro)?chloride)",
    "Bis-Tris": r"(?i:bis[\s-]?tris)\b",
    "Tris": r"(?i:tris)\b(?!\s*\()",
    "HEPES": r"(?i:hepes)\b",
    "MES": r"(?i:mes)\b",
    "phosphate": r"(?i:sodium\s+phosphate|potassium\s+phosphate|na-?phosphate|napi|kpi"
    r"|nah2po4|na2hpo4|kh2po4|k2hpo4"
    r"|phosphate(?!\s+buffered))",
}
_PBS = r"(?i:\bpbs\b|phosphate[\s-]+buffered[\s-]+saline)"
ONE_ONE = {"NaCl", "KCl", "NaNO3", "LiCl", "GdmCl"}
TWO_ONE = {"CaCl2", "MgCl2", "Na2SO4", "(NH4)2SO4"}
PKA = {"Tris": 8.1, "Bis-Tris": 6.5, "HEPES": 7.5, "MES": 6.1, "phosphate": 7.2}
DENATURANTS = ("urea", "GdmCl")


def _norm(s) -> str:
    if not isinstance(s, str):
        return ""
    return s.replace("₂", "2").replace("₃", "3").replace("₄", "4").replace("/ ", "/")


def parse_buffer(text) -> dict:
    """Molar concentration of each recognised species in a buffer description, and whether PBS is named. A
    species is counted once, at its first mention; a number without a unit is read as mM."""
    t = _norm(text)
    out: dict = {"PBS": bool(re.search(_PBS, t))}
    for name, pat in _SPECIES.items():
        # the concentration must stand directly before the name, so "50 mM BisTris" is not read as Tris and
        # "0.5 mM tris(2-carboxyethyl)phosphine" (TCEP) is not read at all
        m = re.search(_NUM + pat, t)
        if m:
            u = m.group("u") or "mM"
            out[name] = float(m.group("c")) * {"M": 1.0, "mM": 1e-3}.get(u, 1e-6)
    return out


def ionic_strength(species: dict, ph: float) -> float:
    """Ionic strength (M) from parse_buffer's species by the design's rule (§4)."""
    i = PBS_IONIC_M if species.get("PBS") else 0.0
    for name, c in species.items():
        if name == "PBS":
            continue
        if name in ONE_ONE:
            i += c
        elif name in TWO_ONE:
            i += 3 * c
        elif name == "phosphate":
            a = 1.0 / (1.0 + 10 ** (PKA[name] - ph))  # fraction HPO4(2-)
            i += c * (1 + 2 * a)
        elif name in ("Tris", "Bis-Tris"):
            i += c / (1.0 + 10 ** (ph - PKA[name]))  # protonated (cationic) fraction
        elif name in ("HEPES", "MES"):
            i += c / (1.0 + 10 ** (PKA[name] - ph))  # deprotonated (anionic) fraction
    return i


def conditions(buffer, additive, temp_c, ph) -> dict:
    """Temperature (K), pH and ionic strength (M) of one entry with the source of each value; denaturants."""
    rec: dict = {}
    if temp_c is not None and np.isfinite(temp_c):
        rec["temp_k"], rec["temp_source"] = round(float(temp_c) + 273.15, 2), "recorded"
    else:
        rec["temp_k"], rec["temp_source"] = DEFAULT_TEMP_K, "default"
    if ph is not None and np.isfinite(ph):
        rec["ph"], rec["ph_source"] = float(ph), "recorded"
    else:
        rec["ph"], rec["ph_source"] = DEFAULT_PH, "default"
    sp = parse_buffer(buffer)
    for k, v in parse_buffer(additive).items():  # the additive field often repeats the buffer's salts
        if k == "PBS":
            sp["PBS"] = sp["PBS"] or v
        else:
            sp.setdefault(k, v)
    ionic = ionic_strength(sp, rec["ph"])
    if ionic > 0:
        rec["ionic_m"], rec["ionic_source"] = round(max(ionic, MIN_IONIC_M), 4), "parsed"
    else:
        rec["ionic_m"], rec["ionic_source"] = DEFAULT_IONIC_M, "default"
    rec["species"] = ";".join(
        [f"{k}={v:g}" for k, v in sp.items() if k != "PBS"] + (["PBS"] if sp["PBS"] else [])
    )
    rec["urea_m"] = sp.get("urea", 0.0)
    rec["gdmcl_m"] = sp.get("GdmCl", 0.0)
    rec["denaturant"] = bool(rec["urea_m"] > 0 or rec["gdmcl_m"] > 0)
    return rec


def prepare_inputs(resolv: pd.DataFrame, entries: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """The simulated proteins: clean disordered entries with a BioEmu-1 ensemble, their construct sequence,
    conditions and run settings, in label order; seed = SEED_BASE + position."""
    r = resolv[
        (resolv["model"] == "bioemu")
        & resolv["clean"].astype(bool)
        & (resolv["disorder_class"] == "disordered")
    ]
    e = entries.set_index("label")
    md = meta.set_index("label")
    rows = []
    for i, lab in enumerate(sorted(r["label"])):
        seq = e.at[lab, "sequence"]
        m = md.loc[lab]
        c = conditions(m.get("buffer"), m.get("buffer_additive"), m.get("temp"), m.get("ph"))
        rows.append(
            {
                "label": lab,
                "length": len(seq),
                "sequence": seq,
                "rg_exp_A": float(e.at[lab, "rg_exp"]),
                **c,
                **run_settings(len(seq)),
                "seed": SEED_BASE + i,
                "buffer_text": m.get("buffer"),
                "additive_text": m.get("buffer_additive"),
            }
        )
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------- analysis
def residue_masses(seq: str, mw: dict) -> np.ndarray:
    """Bead masses as CALVADOS assigns them: residue masses, +2 Da on the first, +16 Da on the last bead."""
    m = np.array([mw[a] for a in seq], dtype=float)
    m[0] += 2.0
    m[-1] += 16.0
    return m


def unwrap_chain(xyz: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Make each frame's chain whole: every bond vector by the minimum-image convention. xyz (F, N, 3)."""
    xyz = np.array(xyz, dtype=float, copy=True)
    box = np.asarray(box, dtype=float).reshape(-1, 1, 3) if np.ndim(box) > 1 else np.asarray(box, float)
    d = np.diff(xyz, axis=1)
    d -= box * np.round(d / box)
    xyz[:, 1:] = xyz[:, :1] + np.cumsum(d, axis=1)
    return xyz


def rg_frames(xyz: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Mass-weighted radius of gyration of each frame, in the units of xyz. xyz (F, N, 3)."""
    w = np.asarray(masses, float) / np.sum(masses)
    com = np.einsum("n,fnk->fk", w, xyz)
    r2 = ((xyz - com[:, None, :]) ** 2).sum(axis=2)
    return np.sqrt(r2 @ w)


def ensemble_rg(rg: np.ndarray) -> float:
    """Ensemble radius of gyration: the root mean square over frames (as for BioEmu-1's coordinate radii)."""
    return float(np.sqrt(np.mean(np.asarray(rg, float) ** 2)))


def convergence(rg: np.ndarray, n_blocks: int = 5, sem_max: float = 0.02, halves_max: float = 0.05) -> dict:
    """Block statistics of Rg^2 (design §5): block Rg values, the relative standard error of the ensemble Rg
    from the block means, the relative difference between the halves, and whether both are under their
    limits."""
    r2 = np.asarray(rg, float) ** 2
    blocks = np.array_split(r2, n_blocks)
    bm = np.array([b.mean() for b in blocks])
    mean = r2.mean()
    sem_r2 = bm.std(ddof=1) / np.sqrt(n_blocks)
    rel_sem = 0.5 * sem_r2 / mean  # d(sqrt x)/sqrt x = dx / 2x
    h1, h2 = np.array_split(r2, 2)
    halves = abs(np.sqrt(h1.mean()) - np.sqrt(h2.mean())) / np.sqrt(mean)
    return {
        "block_rg": [float(np.sqrt(x)) for x in bm],
        "rel_sem": float(rel_sem),
        "halves_rel_diff": float(halves),
        "converged": bool(rel_sem < sem_max and halves < halves_max),
    }


def guinier_bias_gaussian_chain(guinier, rg: float = 30.0, qrg_max: float = 1.3) -> float:
    """Guinier radius over true radius for the Debye function of a Gaussian chain, I(q) = 2(e^-x + x - 1)/x^2,
    x = (q Rg)^2, fitted by `guinier` (bsc.saxs.guinier) the way the archive's curves are fitted."""
    q = np.linspace(0.2, 6.0, 400) / rg
    x = (q * rg) ** 2
    i = 2 * (np.exp(-x) + x - 1) / x**2
    return float(guinier(q, i, 0.01 * i + 1e-12, qrg_max)["rg"] / rg)


def closer_to_one(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """True where ratio a is closer to 1 than ratio b on the log scale."""
    return np.abs(np.log(np.asarray(a, float))) < np.abs(np.log(np.asarray(b, float)))


# ------------------------------------------------------------------------------------------- compute plan
def fit_speed(n: np.ndarray, steps_per_second: np.ndarray) -> tuple[float, float]:
    """Power law steps/s = a * N^-b fitted by least squares on the logs; returns (a, b)."""
    slope, intercept = np.polyfit(
        np.log(np.asarray(n, float)), np.log(np.asarray(steps_per_second, float)), 1
    )
    return float(np.exp(intercept)), float(-slope)


def makespan(hours: list[float], workers: int) -> float:
    """Wall time of independent jobs on `workers` identical workers, longest first to the least loaded."""
    load = [0.0] * workers
    for h in sorted(hours, reverse=True):
        i = int(np.argmin(load))
        load[i] += h
    return max(load)


def compute_plan(
    t: pd.DataFrame,
    a: float,
    b: float,
    workers: int = 6,
    budget_h: float = 24.0,
    long_h: float = 12.0,
    reduced_frames: int = 505,
    reduced_discard: int = 5,
) -> pd.DataFrame:
    """The design's subset rule (§8). t has label, length, n_save, n_frames, n_discard and steps_done.
    Projected hours from the fitted speed; proteins over long_h run reduced_frames; if the makespan still
    exceeds the budget, proteins above the longest length that fits are left out."""
    t = t.copy()
    sps = a * t["length"].astype(float) ** -b
    t["projected_hours_full"] = (t["n_frames"] * t["n_save"] - t["steps_done"]).clip(lower=0) / sps / 3600
    reduce = t["projected_hours_full"] > long_h
    t["frames_planned"] = np.where(reduce, reduced_frames, t["n_frames"])
    t["discard_planned"] = np.where(reduce, reduced_discard, t["n_discard"])
    t["projected_hours"] = (t["frames_planned"] * t["n_save"] - t["steps_done"]).clip(lower=0) / sps / 3600
    t["included"] = True
    if makespan(t["projected_hours"].tolist(), workers) > budget_h:
        for cut in sorted(t["length"].unique(), reverse=True):
            keep = t["length"] <= cut
            if makespan(t.loc[keep, "projected_hours"].tolist(), workers) <= budget_h:
                t["included"] = keep
                break
    return t
