"""Read-only paired density audit; no dynamics-based selection of trajectories."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "density_audit_2026_10_09"
CASES = {"FUS": (280, 325), "A1": (230, 260), "TDP43": (260, 300)}
RADII = (3, 5, 7, 10)
WINDOWS = ((.1, 2), (.1, 5), (.1, 10), (.1, 20), (15, 20))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_branch(root, snap, temp, mass, edges):
    branch = root / "runs" / f"snapshot{snap:02d}_{temp}K"
    assert not list(branch.glob("*.incomplete")), branch
    parts = sorted(branch.glob("part[0-9][0-9][0-9][0-9]"))
    assert len(parts) == 40, branch
    identity = json.loads((branch / "identity.json").read_text())
    profiles, records = [], []
    for n, part in enumerate(parts, 1):
        meta = json.loads((part / "complete.json").read_text())
        assert meta["identity"] == identity
        assert meta["begin"]["step"] == (n - 1) * 50000
        assert meta["end"]["step"] == n * 50000
        assert (part / "trajectory.dcd").stat().st_size > 0
        assert sha(part / "density.npz") == meta["files"]["density.npz"]
        with np.load(part / "density.npz") as data:
            np.testing.assert_array_equal(data["edges_nm"], edges)
            assert data["profiles"].shape == (50, 560)
            profiles.append(data["profiles"])
            records.append(data["records"])
    profile, record = np.concatenate(profiles), np.concatenate(records)
    assert np.isfinite(profile).all() and np.isfinite(record).all()
    np.testing.assert_allclose(record[:, 0], np.arange(1, 2001) / 100, atol=1e-5)
    np.testing.assert_allclose(profile.sum(1) * 200 * 6.02214076e23 / 1e21,
                               mass, rtol=1e-12)
    z = (edges[:-1] + edges[1:]) / 2
    core = np.stack([profile[:, abs(z) < r].mean(1) for r in RADII], axis=1)
    np.testing.assert_allclose(core[:, :3], record[:, 1:4], atol=1e-12)
    return core


def main():
    OUT.mkdir(exist_ok=True)
    times = np.arange(1, 2001) / 100
    rows, blocks, overview, inventory = [], [], [], []
    plot_data = {}
    for case, temps in CASES.items():
        root = BASE / f"{case}_slab_temperature_jump"
        manifest = json.loads((root / "manifest.json").read_text())
        for snap in range(1, 9):
            initial_path = root / "inputs" / f"snapshot{snap:02d}.npz"
            assert sha(initial_path) == manifest["inputs"][initial_path.name]
            with np.load(initial_path) as data:
                edges = data["edges_nm"]
                z = (edges[:-1] + edges[1:]) / 2
                baseline = np.array([data["baseline_profile"][abs(z) < r].mean() for r in RADII])
                mass = data["masses_da"].sum()
            branches = [load_branch(root, snap, t, mass, edges) for t in temps]
            for start, end in WINDOWS:
                mask = (times > start + 1e-8) & (times <= end + 1e-8)
                low, high = [b[mask].mean(0) for b in branches]
                for k, radius in enumerate(RADII):
                    rows.append(dict(case=case, snapshot=snap, start_ns=start, end_ns=end,
                        radius_nm=radius, low_K=temps[0], high_K=temps[1],
                        baseline_g_ml=baseline[k], low_g_ml=low[k], high_g_ml=high[k],
                        paired_difference_pct=100 * (low[k] - high[k]) / baseline[k],
                        low_drift_pct=100 * (low[k] / baseline[k] - 1),
                        high_drift_pct=100 * (high[k] / baseline[k] - 1)))
            for start in np.arange(0, 20, .5):
                mask = (times > max(.1, start) + 1e-8) & (times <= start + .5 + 1e-8)
                low, high = [b[mask].mean(0) for b in branches]
                for k, radius in enumerate(RADII):
                    blocks.append(dict(case=case, snapshot=snap, time_ns=start + .25,
                        radius_nm=radius, low_g_ml=low[k], high_g_ml=high[k],
                        paired_difference_pct=100 * (low[k] - high[k]) / baseline[k]))
            inventory.append(dict(case=case, snapshot=snap, branches=2, frames_per_branch=2000,
                duration_ns=20, density_hashes_verified=80, initial_sha256=sha(initial_path)))
        plot_data[case] = [r for r in blocks if r["case"] == case and r["radius_nm"] == 5]
        for start, end in WINDOWS:
            for radius in RADII:
                selected = [r for r in rows if r["case"] == case and r["radius_nm"] == radius
                            and r["start_ns"] == start and r["end_ns"] == end]
                delta = np.array([r["paired_difference_pct"] for r in selected])
                overview.append(dict(case=case, start_ns=start, end_ns=end, radius_nm=radius,
                    mean_low_g_ml=np.mean([r["low_g_ml"] for r in selected]),
                    mean_high_g_ml=np.mean([r["high_g_ml"] for r in selected]),
                    mean_pair_difference_pct=delta.mean(), min_pair_difference_pct=delta.min(),
                    max_pair_difference_pct=delta.max(), pairs_mean_difference_within_5pct=int((abs(delta) <= 5).sum())))
        print(f"Verified {case}: 16 branches, 320 ns", flush=True)
    table("paired_windows.csv", rows)
    table("paired_half_ns_blocks.csv", blocks)
    table("system_summary.csv", overview)
    table("inventory.csv", inventory)
    plot(plot_data)
    notes = dict(density_definition="Mass density of all beads within centered |z| < radius, not whole-box density.",
        normalization="100 * (rho_low - rho_high) / pre-jump 0.5 ns baseline density, separately per snapshot.",
        windows="First 0.1 ns omitted consistently with existing protocol. Windows are descriptive, not selected by dynamics.",
        tolerance="5% counts apply only to window averages, not every block; not an equivalence test or original strict screen.",
        uncertainty="Eight starting snapshots from a source trajectory; no independence claim. Plot bands are snapshot ranges, not CIs.",
        validation="Initial and all 1920 density files hashed; step/time continuity, profile mass and stored core densities verified. DCD presence checked; DCD contents not rehashed here.",
        dynamics="Only FUS snapshot01 has existing full-10-ps contact/Rg observables. Full three-system paired dynamics not yet extracted.")
    (OUT / "audit_notes.json").write_text(json.dumps(notes, indent=2) + "\n")
    print(json.dumps([r for r in overview if r["radius_nm"] == 5], indent=2))


def plot(data):
    plt.rcParams.update({"font.size": 12, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True, layout="constrained")
    for col, (case, rows) in enumerate(data.items()):
        x = np.arange(.25, 20, .5)
        for key, color, label in (("low_g_ml", "#177b78", r"$T_{\mathrm{L}}$"),
                                  ("high_g_ml", "#b44562", r"$T_{\mathrm{H}}$")):
            values = np.array([[r[key] for r in rows if r["snapshot"] == s] for s in range(1, 9)])
            axes[0, col].fill_between(x, values.min(0), values.max(0), color=color, alpha=.12)
            axes[0, col].plot(x, values.mean(0), color=color, label=label, lw=2)
        axes[0, col].set_title(case)
        axes[0, col].legend(frameon=False)
        for s in range(1, 9):
            y = [r["paired_difference_pct"] for r in rows if r["snapshot"] == s]
            axes[1, col].plot(x, y, lw=1, alpha=.7, label=f"{s:02d}")
        axes[1, col].axhline(0, color=".4", lw=.8)
        for limit in (-5, 5):
            axes[1, col].axhline(limit, color=".6", ls=":", lw=1)
        axes[1, col].set(xlabel="Time after branching (ns)", xlim=(0, 20))
        axes[0, col].set_ylabel(r"Core density (g/mL), $|z|<5$ nm")
        axes[1, col].set_ylabel("Paired density difference (%)")
    axes[1, 2].legend(ncol=4, fontsize=9, title="Starting snapshot", frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"paired_density_check.{suffix}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
