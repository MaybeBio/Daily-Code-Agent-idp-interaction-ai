"""Analysis of the CALVADOS 3 runs against the measured radii and BioEmu-1 (design §3, §5, §6).

    BSC_ROOT=<where bsc fetch put data/> SIM_WORK=... uv run python analysis/simulations/07_analyse.py

Runs in the main environment (bsc, mdtraj). Per protein: the mass-weighted bead Rg of every production
frame, the ensemble Rg (root mean square over frames), its ratio to the measured Guinier Rg, the block
convergence check, and, where the computed curves exist, the reduced chi2, NRMSD and Guinier Rg of the
uniform ensemble curve scored exactly as the archive's BioEmu-1 curves (bsc.saxs). Paired with BioEmu-1's
curve ratio (primary), its Calpha coordinate ratio and its raw chi2, over all simulated proteins and over
the subsets the design names. Writes results/simulations/per_protein.csv, convergence.csv, summary.json,
run_metadata.json and results/figures/fig_sim_rg.png."""

from __future__ import annotations

import json
import os
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd
from common import DL, OUT, RES, RUNS, SAXS
from scipy import stats
from simlib import (
    closer_to_one,
    convergence,
    ensemble_rg,
    guinier_bias_gaussian_chain,
    residue_masses,
    rg_frames,
    unwrap_chain,
)

from bsc import config as C
from bsc import report, robust, saxs

SEED = 20261015
N_BOOT = 2000
E10_BAND = (0.90, 1.10)
CALVADOS_COL = "#0E9F6E"  # checked against BioEmu-1's role colour with the dataviz palette validator
SAXS_DIR = Path(os.environ.get("BSC_ROOT", RES.parent)) / "data" / "PeptoneDB-SAXS" / "sasbdb-clean_data"


def wilson(k: int, n: int, z: float = 1.959964) -> list[float]:
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(c - h), float(c + h)]


def protein_rg(label: str, r: pd.Series, mw: dict) -> tuple[np.ndarray, dict]:
    run = RUNS / label
    chunks = sorted(run.glob("chunk_*.dcd"))
    t = md.load([str(c) for c in chunks], top=str(run / "top.pdb"))
    xyz = unwrap_chain(t.xyz, np.array([r["box_nm"]] * 3)) * 10.0  # nm -> Angstrom
    rg = rg_frames(xyz, residue_masses(r["sequence"], mw))[int(r["n_discard"]) :]
    return rg, {"n_frames_analysed": int(len(rg))}


def score_curves(label: str) -> dict | None:
    f = SAXS / label / "curves.csv"
    if not f.exists():
        return None
    curves = pd.read_csv(f, index_col=0).to_numpy(float)
    exp = pd.read_csv(
        SAXS_DIR / f"{label}-bift.dat", sep=r"\s+", comment="#", header=None, names=["q", "I", "sigma"]
    )
    q, I, s = exp["q"].to_numpy(), exp["I"].to_numpy(), exp["sigma"].to_numpy()
    keep = s > 0
    ok = np.isfinite(curves).all(axis=1)  # conformers Pepsi-SAXS failed on are dropped, as in score.py
    curves, q, I, s = curves[ok][:, keep], q[keep], I[keep], s[keep]
    raw = saxs.ensemble_curve(curves)
    g = saxs.guinier(q, raw, 0.01 * raw + 1e-9, C.GUINIER_QRG_MAX)  # as for the archive's curves (score.py)
    return {
        "n_curves": int(len(curves)),
        "n_curves_dropped": int((~ok).sum()),
        "chi2_sim": saxs.chi2(raw, I, s),
        "nrmsd_sim": saxs.nrmsd_log(raw, I, s),
        "rg_curve_sim_A": g["rg"],
    }


def wall_clock_span() -> float:
    """Hours from the start of the first chunk to the end of the last (file times, chunk durations)."""
    starts, ends = [], []
    for prog_f in RUNS.glob("*/progress.json"):
        prog = json.loads(prog_f.read_text())
        for c in prog["chunks"]:
            f = prog_f.parent / f"chunk_{c['chunk']:03d}.dcd"
            if f.exists():
                ends.append(f.stat().st_mtime)
                starts.append(f.stat().st_mtime - c["seconds"])
    return round((max(ends) - min(starts)) / 3600, 2) if ends else 0.0


def paired(a: np.ndarray, b: np.ndarray, seed: int) -> dict:
    """a: CALVADOS ratios, b: BioEmu-1 ratios on the same proteins."""
    n = len(a)
    closer = closer_to_one(a, b)
    k = int(closer.sum())
    out = {
        "n": n,
        "calvados_median": float(np.median(a)),
        "calvados_ci95": list(robust.bootstrap_ci(a, n_boot=N_BOOT, seed=seed)),
        "bioemu_median": float(np.median(b)),
        "bioemu_ci95": list(robust.bootstrap_ci(b, n_boot=N_BOOT, seed=seed + 1)),
        "calvados_share_within_0.9_1.1": float(np.mean((a >= E10_BAND[0]) & (a <= E10_BAND[1]))),
        "bioemu_share_within_0.9_1.1": float(np.mean((b >= E10_BAND[0]) & (b <= E10_BAND[1]))),
        "closer_to_1": {
            "k": k,
            "n": n,
            "share": k / n,
            "ci95": wilson(k, n),
            "sign_test_p_two_sided": float(stats.binomtest(k, n, 0.5).pvalue),
        },
    }
    if n >= 6:
        out["wilcoxon_log_ratio_difference_p"] = float(stats.wilcoxon(np.log(a) - np.log(b)).pvalue)
        out["wilcoxon_abs_log_ratio_p"] = float(stats.wilcoxon(np.abs(np.log(a)) - np.abs(np.log(b))).pvalue)
    out["E10_met"] = bool(E10_BAND[0] <= out["calvados_median"] <= E10_BAND[1])
    out["E10_interval_inside_band"] = bool(
        E10_BAND[0] <= out["calvados_ci95"][0] and out["calvados_ci95"][1] <= E10_BAND[1]
    )
    out["E11_met"] = bool(
        out["closer_to_1"]["share"] > 0.5 and out["closer_to_1"]["sign_test_p_two_sided"] < 0.05
    )
    return out


def figure(t: pd.DataFrame, S: dict, path: Path) -> None:
    """Paired ratios: BioEmu-1 (left) and CALVADOS 3 (right) for each protein, joined by a line; log axis."""
    import matplotlib.pyplot as plt

    report.style()
    fig, ax = report.figure(560, 400)
    ax.axhspan(*E10_BAND, color="#EEF1F2", lw=0, zorder=0)
    ax.axhline(1.0, color=report.MUTED, lw=0.8, ls="--", zorder=1)
    flag = t["flag_training"].to_numpy() | t["flag_validation"].to_numpy()
    jit = np.random.default_rng(SEED).uniform(-0.07, 0.07, len(t))  # the same offset on both sides of a pair
    for i, (_, r) in enumerate(t.iterrows()):
        ax.plot(
            [jit[i], 1 + jit[i]],
            [r["bioemu_rg_ratio"], r["calvados_rg_ratio"]],
            color="#B9C0C6",
            lw=0.8,
            zorder=2,
        )
    for x, col, key in (
        (0, report.ROLE_COL["bioemu"], "bioemu_rg_ratio"),
        (1, CALVADOS_COL, "calvados_rg_ratio"),
    ):
        v = t[key].to_numpy()
        ax.scatter(x + jit[~flag], v[~flag], s=34, color=col, alpha=0.75, lw=0, zorder=3)
        ax.scatter(x + jit[flag], v[flag], s=34, facecolor="white", edgecolor=col, lw=1.4, zorder=3)
    P = S["primary"]["all"]
    for x, med, ci in (
        (0, P["bioemu_median"], P["bioemu_ci95"]),
        (1, P["calvados_median"], P["calvados_ci95"]),
    ):
        ax.plot([x + 0.2] * 2, ci, color=report.INK, lw=1.6, zorder=4)
        ax.plot([x + 0.14, x + 0.26], [med] * 2, color=report.INK, lw=2.2, zorder=4)
        ax.text(x + 0.3, med, f"{med:.2f}", va="center", fontsize=9.5, color=report.INK)
    ax.set_yscale("log")
    ax.set_xlim(-0.35, 1.55)
    ax.set_xticks([0, 1], ["BioEmu-1\n(ensemble curve)", "CALVADOS 3\n(simulated beads)"])
    ax.set_ylabel("ensemble Rg / measured Rg")
    ticks = [0.5, 0.7, 0.9, 1.0, 1.1, 1.5, 2.0, 2.5]
    lo, hi = (
        min(t[["bioemu_rg_ratio", "calvados_rg_ratio"]].min()) * 0.92,
        max(t[["bioemu_rg_ratio", "calvados_rg_ratio"]].max()) * 1.08,
    )
    ax.set_yticks([v for v in ticks if lo <= v <= hi], [f"{v:g}" for v in ticks if lo <= v <= hi])
    ax.minorticks_off()
    ax.set_ylim(lo, hi)
    ax.scatter([], [], s=34, color=report.MUTED, lw=0, label="not in CALVADOS training or validation")
    ax.scatter(
        [],
        [],
        s=34,
        facecolor="white",
        edgecolor=report.MUTED,
        lw=1.4,
        label="≥ 90% identical to a CALVADOS training or validation protein",
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), frameon=False, ncol=1)
    fig.tight_layout()
    fig.savefig(path, dpi=report.DPI, facecolor="white")
    plt.close(fig)
    from PIL import Image

    Image.open(path).convert("RGB").save(path, optimize=True)


def main() -> None:
    inp = pd.read_csv(OUT / "inputs.csv")
    circ = pd.read_csv(OUT / "circularity.csv").set_index("label")
    res = pd.read_csv(RES / "resolvability.csv")
    bio = res[(res["model"] == "bioemu")].set_index("label")
    coord = pd.read_csv(RES / "coord_rg.csv").set_index("label")
    mw = pd.read_csv(DL / "residues_CALVADOS3.csv").set_index("one")["MW"].to_dict()
    plan = {p["label"]: p for p in json.loads((OUT / "plan.json").read_text())["proteins"]}
    rows, conv_rows, runs = [], [], []
    for _, r in inp.iterrows():
        lab = r["label"]
        r = r.copy()
        r["n_frames"], r["n_discard"] = plan[lab]["frames_planned"], plan[lab]["discard_planned"]
        prog_f = RUNS / lab / "progress.json"
        prog = json.loads(prog_f.read_text()) if prog_f.exists() else {"frames": 0}
        runs.append(
            {
                "label": lab,
                "length": int(r["length"]),
                "frames": prog["frames"],
                "target_frames": int(r["n_frames"]),
                "included": bool(plan[lab]["included"]),
                "wall_hours": round(prog.get("wall_seconds", 0) / 3600, 3),
                "steps_per_second": prog.get("steps_per_second"),
                "ns": round(prog["frames"] * int(r["n_save"]) * 1e-5, 2),
            }
        )
        if not plan[lab]["included"] or prog["frames"] < int(r["n_frames"]):
            continue
        rg, meta = protein_rg(lab, r, mw)
        cv = convergence(rg)
        conv_rows.append(
            {
                "label": lab,
                **{f"block_{i + 1}_rg_A": v for i, v in enumerate(cv["block_rg"])},
                "rel_sem": cv["rel_sem"],
                "halves_rel_diff": cv["halves_rel_diff"],
                "converged": cv["converged"],
            }
        )
        rg_sim = ensemble_rg(rg)
        rec = {
            "label": lab,
            "length": int(r["length"]),
            "rg_exp_A": r["rg_exp_A"],
            "rg_sim_A": rg_sim,
            "calvados_rg_ratio": rg_sim / r["rg_exp_A"],
            "bioemu_rg_ratio": float(bio.at[lab, "rg_ratio"]),
            "bioemu_coord_rg_ratio": float(coord.at[lab, "ratio_coord_exp"])
            if lab in coord.index
            else np.nan,
            "bioemu_chi2_raw": float(bio.at[lab, "chi2_raw"]),
            "bioemu_nrmsd_raw": float(bio.at[lab, "nrmsd_raw"]),
            "temp_k": r["temp_k"],
            "ionic_m": r["ionic_m"],
            "ph": r["ph"],
            "denaturant": bool(r["denaturant"]),
            "converged": cv["converged"],
            **meta,
            "flag_training": "training" in str(circ.at[lab, "flag_sets"]),
            "flag_validation": bool(circ.at[lab, "flagged"])
            and "training" not in str(circ.at[lab, "flag_sets"]),
            "best_identity_calvados_sets": float(circ.at[lab, "best_identity"]),
        }
        sc = score_curves(lab)
        if sc:
            rec.update(sc)
            rec["calvados_curve_rg_ratio"] = sc["rg_curve_sim_A"] / r["rg_exp_A"]
        rows.append(rec)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "per_protein.csv", index=False)
    pd.DataFrame(conv_rows).to_csv(OUT / "convergence.csv", index=False)
    R = pd.DataFrame(runs)
    meta = {
        "n_proteins": int(len(R)),
        "n_complete": int((R["included"] & (R["frames"] >= R["target_frames"])).sum()),
        # every simulation step run, the pilot of the protein the plan left out included
        "total_wall_hours_single_thread": round(float(R["wall_hours"].sum()), 2),
        "wall_hours_analysed_proteins": round(float(R.loc[R["included"], "wall_hours"].sum()), 2),
        "wall_hours_left_out_pilot": round(float(R.loc[~R["included"], "wall_hours"].sum()), 2),
        "median_wall_hours_per_analysed_protein": round(
            float(R.loc[R["included"], "wall_hours"].median()), 3
        ),
        "max_wall_hours_per_protein": round(float(R["wall_hours"].max()), 2),
        "total_ns": round(float(R["ns"].sum()), 1),
        "wall_clock_span_hours": wall_clock_span(),
        "per_protein": R.to_dict(orient="records"),
    }
    (OUT / "run_metadata.json").write_text(json.dumps(meta, indent=1) + "\n")
    if t.empty:
        print("no complete run yet")
        return
    subsets = {
        "all": np.ones(len(t), bool),
        "without_training_overlap": ~t["flag_training"].to_numpy(),
        "without_any_overlap": ~(t["flag_training"] | t["flag_validation"]).to_numpy(),
        "without_denaturant": ~t["denaturant"].to_numpy(),
        "converged_only": t["converged"].to_numpy(bool),
    }
    S = {
        "n": int(len(t)),
        "e10_band": list(E10_BAND),
        "n_boot": N_BOOT,
        "seed": SEED,
        # the archive's Guinier fit (q Rg <= 1.3) on the Debye function of a Gaussian chain: fitted / true Rg
        "guinier_over_true_rg_gaussian_chain": guinier_bias_gaussian_chain(
            saxs.guinier, qrg_max=C.GUINIER_QRG_MAX
        ),
        "n_flag_training": int(t["flag_training"].sum()),
        "n_flag_validation_only": int(t["flag_validation"].sum()),
        "n_denaturant": int(t["denaturant"].sum()),
        "n_converged": int(t["converged"].sum()),
        "primary": {
            k: paired(t.loc[m, "calvados_rg_ratio"].to_numpy(), t.loc[m, "bioemu_rg_ratio"].to_numpy(), SEED)
            for k, m in subsets.items()
            if m.sum() >= 3
        },
    }
    c = t.dropna(subset=["bioemu_coord_rg_ratio"])
    S["sensitivity_coordinates"] = paired(
        c["calvados_rg_ratio"].to_numpy(), c["bioemu_coord_rg_ratio"].to_numpy(), SEED + 10
    )
    if "calvados_curve_rg_ratio" in t:
        c = t.dropna(subset=["calvados_curve_rg_ratio"])
        S["sensitivity_curves"] = paired(
            c["calvados_curve_rg_ratio"].to_numpy(), c["bioemu_rg_ratio"].to_numpy(), SEED + 20
        )
        a, b = c["chi2_sim"].to_numpy(), c["bioemu_chi2_raw"].to_numpy()
        k = int((a < b).sum())
        S["chi2"] = {
            "n": int(len(c)),
            "calvados_median": float(np.median(a)),
            "bioemu_median": float(np.median(b)),
            "calvados_ci95": list(robust.bootstrap_ci(a, n_boot=N_BOOT, seed=SEED + 30)),
            "bioemu_ci95": list(robust.bootstrap_ci(b, n_boot=N_BOOT, seed=SEED + 31)),
            "calvados_share_le_2": float(np.mean(a <= 2)),
            "bioemu_share_le_2": float(np.mean(b <= 2)),
            "calvados_lower": {
                "k": k,
                "n": int(len(c)),
                "share": k / len(c),
                "ci95": wilson(k, len(c)),
                "sign_test_p_two_sided": float(stats.binomtest(k, len(c), 0.5).pvalue),
            },
            "wilcoxon_log_p": float(stats.wilcoxon(np.log(a) - np.log(b)).pvalue),
            "nrmsd_calvados_median": float(c["nrmsd_sim"].median()),
            "nrmsd_bioemu_median": float(c["bioemu_nrmsd_raw"].median()),
        }
    (OUT / "summary.json").write_text(json.dumps(S, indent=1) + "\n")
    figure(t, S, RES / "figures" / "fig_sim_rg.png")
    p = S["primary"]["all"]
    print(
        f"n = {p['n']}: CALVADOS median ratio {p['calvados_median']:.3f} {p['calvados_ci95']}, "
        f"BioEmu-1 {p['bioemu_median']:.3f}; "
        f"closer to 1 for {p['closer_to_1']['k']} of {p['n']}; E10 {p['E10_met']}, E11 {p['E11_met']}"
    )


if __name__ == "__main__":
    main()
