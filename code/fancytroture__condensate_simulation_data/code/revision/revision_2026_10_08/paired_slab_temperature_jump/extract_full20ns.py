"""Extract all paired jump trajectories without modifying simulation outputs."""
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from MDAnalysis.coordinates.DCD import DCDReader
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
OUT = HERE / "analysis_full20ns_2026_10_09"
BOX = np.array([20., 20., 280.])
CASES = {"FUS": (163, 280, 325), "A1": (135, 230, 260), "TDP43": (148, 260, 300)}
CONTACT_SOURCE = ROOT / "revision_2026_09_30/interchain_recompute/compute.py"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pure_function(path, name, extra=None):
    node = next(n for n in ast.parse(path.read_text()).body
                if isinstance(n, ast.FunctionDef) and n.name == name)
    namespace = {"np": np, "cKDTree": cKDTree}
    namespace.update(extra or {})
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    return namespace[name]


original_contacts = pure_function(CONTACT_SOURCE, "partner_contact_counts")


def write_json(path, value):
    tmp = path.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def whole(x, nres):
    chains = x.reshape(200, nres, 3).copy()
    delta = np.diff(chains, axis=1)
    delta -= np.rint(delta / BOX) * BOX
    chains[:, 1:] = chains[:, :1] + np.cumsum(delta, axis=1)
    return chains


def contacts(x, nres, cohort):
    # A single periodic neighbor list yields the same directed partner counts
    # as querying every target bead; intrachain contacts remain unique pairs.
    pairs = cKDTree(x % BOX, boxsize=BOX).query_pairs(.66, output_type="ndarray")
    chain = pairs // nres
    inter = chain[:, 0] != chain[:, 1]
    counts = np.bincount(chain[inter, 0] * 200 + chain[inter, 1], minlength=40000).reshape(200, 200)
    counts += counts.T.copy()
    eligible = (~inter) & ((pairs[:, 1] - pairs[:, 0]) > 2)
    selected = pairs[eligible]
    delta = x[selected[:, 1]] - x[selected[:, 0]]
    delta -= np.rint(delta / BOX) * BOX
    strict = np.sum(delta**2, axis=1) < .66**2
    intra = np.bincount(selected[strict, 0] // nres, minlength=200)
    return counts[cohort].astype(np.int32), intra[cohort].astype(np.int32)


def check_contacts(x, nres, cohort, partner, intra):
    reference = original_contacts(x, BOX, nres, cohort, cutoff=.66)
    np.testing.assert_array_equal(partner, np.stack([reference[int(c)] for c in cohort]))
    i, j = np.triu_indices(nres, k=3)
    for index in (0, len(cohort) - 1):
        chain = x[cohort[index] * nres:(cohort[index] + 1) * nres]
        delta = chain[j] - chain[i]
        delta -= np.rint(delta / BOX) * BOX
        assert intra[index] == np.count_nonzero(np.sum(delta**2, axis=1) < .66**2)


def run_branch(job):
    case, snapshot, temperature = job
    start = time.monotonic()
    nres = CASES[case][0]
    root = BASE / f"{case}_slab_temperature_jump"
    branch = root / "runs" / f"snapshot{snapshot:02d}_{temperature}K"
    dest = OUT / "observables" / f"{case}_snapshot{snapshot:02d}_{temperature}K.npz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    initial_path = root / "inputs" / f"snapshot{snapshot:02d}.npz"
    manifest = json.loads((root / "manifest.json").read_text())
    assert sha(initial_path) == manifest["inputs"][initial_path.name]
    with np.load(initial_path) as initial:
        cohort = initial["initial_center_chain_ids"].astype(int)
        mass = initial["masses_da"].copy()
    parts = sorted(branch.glob("part[0-9][0-9][0-9][0-9]"))
    assert [p.name for p in parts] == [f"part{i:04d}" for i in range(1, 41)]
    assert not list(branch.glob("*.incomplete"))
    identity = dict(case=case, snapshot=snapshot, temperature_K=temperature,
        extractor_sha256=sha(__file__), contact_source_sha256=sha(CONTACT_SOURCE),
        initial_sha256=sha(initial_path), part_manifests=[sha(p / "complete.json") for p in parts],
        nres=nres, dt_ns=.01, cutoff_nm=.66, intra_sequence_exclusion=2,
        chain_ids=cohort.tolist(), rg="Equal bead weights; minimum-image bonded unwrapping",
        cohort="Initial center chain cohort held fixed; no outcome-dependent selection")
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    if dest.exists():
        with np.load(dest) as data:
            assert str(data["identity_sha256"]) == digest, f"Stale cache: {dest}"
        return dict(case=case, snapshot=snapshot, temperature_K=temperature, reused=True)
    rg, intra_all, partner_all, z_all, rho_all = [], [], [], [], []
    for index, part in enumerate(parts):
        meta = json.loads((part / "complete.json").read_text())
        assert meta["begin"]["step"] == index * 50000
        assert meta["end"]["step"] == (index + 1) * 50000
        assert sha(part / "trajectory.dcd") == meta["files"]["trajectory.dcd"]
        assert sha(part / "density.npz") == meta["files"]["density.npz"]
        with np.load(part / "density.npz") as density:
            rho_all.append(density["records"][:, 1:4])
        with DCDReader(str(part / "trajectory.dcd")) as reader:
            assert len(reader) == 50 and reader.n_atoms == 200 * nres
            np.testing.assert_allclose(reader.dt, 10., rtol=1e-6)
            for local, frame in enumerate(reader):
                tick = index * 50 + local + 1
                np.testing.assert_allclose(frame.time / 1000, tick * .01, atol=1e-5)
                np.testing.assert_allclose(frame.dimensions, [200, 200, 2800, 90, 90, 90], atol=.001)
                x = frame.positions.astype(float) / 10
                chain = whole(x, nres)[cohort]
                rg.append(np.sqrt(np.mean(np.sum((chain - chain.mean(1)[:, None])**2, axis=2), axis=1)))
                partner, intra = contacts(x, nres, cohort)
                if local == 0:
                    check_contacts(x, nres, cohort, partner, intra)
                intra_all.append(intra)
                partner_all.append(partner)
                center = np.angle(np.sum(mass * np.exp(2j * np.pi * x[:, 2] / BOX[2]))) * BOX[2] / (2 * np.pi)
                chain_mass = mass.reshape(200, nres)[cohort]
                chain_z = (chain[:, :, 2] * chain_mass).sum(1) / chain_mass.sum(1)
                z_all.append((chain_z - center + BOX[2] / 2) % BOX[2] - BOX[2] / 2)
        if (index + 1) % 10 == 0:
            print(f"{case} s{snapshot:02d} {temperature}K: {(index+1)*.5:g}/20 ns, {time.monotonic()-start:.0f} s", flush=True)
    data = dict(time_ns=np.arange(1, 2001) / 100, Rg=np.asarray(rg), N_intra=np.asarray(intra_all),
        partner_contacts=np.asarray(partner_all, dtype=np.int32), chain_ids=cohort,
        chain_center_z_nm=np.asarray(z_all), density_core_3_5_7=np.concatenate(rho_all),
        identity_sha256=np.asarray(digest))
    data["N_inter"] = data["partner_contacts"].sum(2)
    data["N_partner"] = (data["partner_contacts"] > 0).sum(2)
    assert np.isfinite(data["Rg"]).all() and np.all(data["Rg"] > 0)
    assert np.all(data["partner_contacts"][:, np.arange(len(cohort)), cohort] == 0)
    old_match = None
    if case == "FUS" and snapshot == 1:
        old_path = root / "acf_snapshot01_20ns" / f"observables_full10ps_{temperature}K.npz"
        with np.load(old_path) as old:
            for key in ("Rg", "N_intra", "N_inter", "partner_contacts", "chain_ids"):
                np.testing.assert_allclose(data[key], old[key], atol=1e-12, rtol=1e-12)
        old_match = True
    tmp = dest.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **data)
    tmp.replace(dest)
    record = dict(identity=identity, frames=2000, duration_ns=20, n_chains=len(cohort),
        direct_contact_checks=40, pilot_exact_match=old_match, elapsed_s=time.monotonic()-start,
        output_sha256=sha(dest))
    write_json(dest.with_suffix(".json"), record)
    print(f"DONE {case} s{snapshot:02d} {temperature}K ({record['elapsed_s']:.0f}s)", flush=True)
    return dict(case=case, snapshot=snapshot, temperature_K=temperature, reused=False,
                elapsed_s=record["elapsed_s"], chains=len(cohort))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    jobs = [(case, s, t) for s in range(1, 9) for case, (_, low, high) in CASES.items() for t in (low, high)]
    completed = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_branch, job): job for job in jobs}
        for future in as_completed(futures):
            completed.append(future.result())
            write_json(OUT / "extraction_progress.json", dict(completed=completed, total=48,
                       status="complete" if len(completed) == 48 else "running"))
    print("All 48 branches extracted and validated", flush=True)


if __name__ == "__main__":
    main()
