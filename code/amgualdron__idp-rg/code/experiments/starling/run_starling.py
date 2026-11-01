#!/usr/bin/env python
"""STARLING ensembles and R_g distributions for the experimental benchmark (data/benchmark.csv).

    mamba run -n starling python experiments/starling/run_starling.py --dry-run
    mamba run -n starling python experiments/starling/run_starling.py --out runs/starling/benchmark

STARLING generates residue-residue distance maps. Each map gives one R_g, so the distribution of R_g is a set of
independent samples, one per conformation; there is no fitted curve. Everything here derives from those samples.

Each protein is seeded from --seed and its name, so a result does not depend on which proteins ran before it. It is
written atomically as soon as it finishes. Rerunning with the same settings skips finished proteins and retries the
rest. Rerunning with different settings refuses to mix them in one directory.

Writes to --out:
    run_config.json     settings, versions, GPU and the weights' digests
    samples/<name>.npz  rg_nm (one R_g per conformation), bad_frame (bool, same length), meta (JSON string)
    rg_samples.npz      the same rg_nm arrays in one file, keyed by protein name
    summary.csv         one row for every input protein, whatever happened to it (see `status`)
    histograms.csv      name, bin_left_nm, bin_right_nm, count, density (Freedman-Diaconis bins for each protein)
    ensembles/          <name>.starling.xz, only with --save-ensembles
    failures.log        tracebacks, only if a protein failed

`status` in summary.csv is ok, too_long (over STARLING's length limit, never sampled), failed, or pending (not reached:
--names, --limit, or an interrupted run). bad_frame marks maps that STARLING's own check calls physically impossible,
a pair further apart than 4.81 A per bond. rg_mean_nm uses all frames, like STARLING's own mean; rg_mean_clean_nm
excludes the bad ones.
"""
import argparse
import gc
import hashlib
import json
import os
import pathlib
import platform
import re
import sys
import time
import traceback

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[2]
VALID_RESIDUES = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]+")
VALID_NAME = re.compile(r"[A-Za-z0-9_.\-]+")     # names become file names
IONIC_STRENGTHS = (20, 150, 300)                 # mM: the only values STARLING was trained on
MAX_CONSECUTIVE_FAILURES = 3                     # beyond this the GPU or environment is broken; stop rather than churn
# STARLING 2.0.2's DDIM schedule is range(0, 999, 1000 // steps) + 1. From 334 steps on it includes t = 999, where
# alpha_bar = 2.4e-10 and the x0 estimate is amplified 64,000 times: Tau_ht23's <R_g> went from 6.7 to 9.2 nm, silently.
DDIM_MAX_STEPS = 333
SAME_RUN_KEYS = ("conformations", "steps", "sampler", "ionic_strength_mM", "seed", "starling_version", "fast")
# SHA-256 of the release assets, from the GitHub API for idptools/starling tag v2.0.0. torch.hub does not check them.
WEIGHT_DIGESTS = {
    "STARLING_v2.0.0_ViT_VAE_2025_10_14.ckpt": "f373190466250030507c00bdb2bcf7757461c675097621a13aa4a92188949ec5",
    "STARLING_v2.0.0_ViT_DDPM_2025_10_14.ckpt": "d5642ad5ffce96c5c3f470e3314a8ab80d2d00387b148e2a23b769b9b942af41",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=pathlib.Path, default=ROOT / "data" / "benchmark.csv",
                   help="CSV with columns name and sequence")
    p.add_argument("--out", type=pathlib.Path, default=ROOT / "runs" / "starling" / "benchmark")
    p.add_argument("--conformations", type=int, default=2000)
    p.add_argument("--batch-size", type=int, default=25,
                   help="25 needs 5.5 GB of GPU memory and 100 does not fit in 8 GB; halved automatically if the GPU runs "
                        "out of memory")
    p.add_argument("--steps", type=int, default=45,
                   help=f"ddim denoising steps, at most {DDIM_MAX_STEPS}; STARLING's default is 30, which narrows the "
                        "distributions (see README). ddpm always uses its full schedule")
    p.add_argument("--sampler", choices=("ddim", "ddpm"), default="ddim")
    p.add_argument("--plain", action="store_true",
                   help="on the GPU, do not use the default float16 denoiser and torch.compile (see README)")
    p.add_argument("--ionic-strength", type=int, choices=IONIC_STRENGTHS, default=150, help="mM")
    p.add_argument("--device", default="cuda:0")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--names", nargs="+", help="run only these proteins")
    p.add_argument("--limit", type=int, help="run only the first K proteins that remain to be done")
    p.add_argument("--save-ensembles", action="store_true",
                   help="also keep the distance maps as compressed .starling files (0.1 A precision)")
    p.add_argument("--dry-run", action="store_true", help="validate the input, print the plan and exit; no GPU")
    args = p.parse_args()
    if args.conformations < 1 or args.batch_size < 1 or args.steps < 1:
        p.error("--conformations, --batch-size and --steps must be positive")
    if args.limit is not None and args.limit < 1:
        p.error("--limit must be positive")
    if args.sampler == "ddim" and args.steps > DDIM_MAX_STEPS:
        p.error(f"--steps above {DDIM_MAX_STEPS} reach t = 999, where STARLING 2.0.2's DDIM produces wrong ensembles")
    return args


def load_input(path):
    """The protein table, checked, with N and the reference R_g values (nm) that are available."""
    t = pd.read_csv(path)
    missing = {"name", "sequence"} - set(t.columns)
    if missing:
        sys.exit(f"{path}: missing columns {sorted(missing)}")
    bad = t[~t.sequence.map(lambda s: isinstance(s, str) and VALID_RESIDUES.fullmatch(s) is not None)]
    if len(bad):        # STARLING would silently delete non-standard residues
        sys.exit(f"{path}: sequences with non-standard residues: {bad.name.tolist()}")
    bad = t[~t.name.map(lambda s: VALID_NAME.fullmatch(str(s)) is not None)]
    if len(bad):
        sys.exit(f"{path}: names unusable as file names: {bad.name.tolist()}")
    if t.name.duplicated().any():
        sys.exit(f"{path}: duplicate names: {t.name[t.name.duplicated()].tolist()}")

    out = pd.DataFrame(dict(name=t.name, N=t.sequence.str.len(), sequence=t.sequence))
    # benchmark.csv has carried these in nm (committed) and in A (older compiled table).
    for new, nm, angstrom in (("saxs_rg_nm", "saxs_rg_nm", "SAXS Rg (A)"), ("hps1_rg_nm", "hps1_rg_nm", "HPS1 Rg (A) [sim]")):
        if nm in t:
            out[new] = t[nm]
        elif angstrom in t:
            out[new] = t[angstrom] / 10
    methods = path.parent / "benchmark_methods.csv"
    if methods.exists():    # STARLING was trained on Mpipi-GG simulations, so this is the column it should reproduce
        out = out.merge(pd.read_csv(methods)[["name", "mpipi_gg_rg_nm"]], on="name", how="left", validate="1:1")
    return out


def protein_seed(base, name):
    """Stable across runs and independent of processing order (unlike Python's hash)."""
    return int.from_bytes(hashlib.sha256(f"{base}:{name}".encode()).digest()[:8], "big") % 2**63


def atomic_write(path, write):
    """Call write(file object) on a temporary neighbour of path, then move it into place."""
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as f:
        write(f)
    os.replace(tmp, path)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_weights(torch, configs):
    """The two checkpoints, downloaded if needed and verified. Checkpoints are pickles: verify before loading them."""
    digests = {}
    for source in (configs.DEFAULT_ENCODER_WEIGHTS_PATH, configs.DEFAULT_DDPM_WEIGHTS_PATH):
        name = os.path.basename(source)
        expected = WEIGHT_DIGESTS.get(name)
        path = pathlib.Path(source)
        if source.startswith("http"):       # the location STARLING itself caches them
            path = pathlib.Path(torch.hub.get_dir()) / "checkpoints" / name
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                torch.hub.download_url_to_file(source, str(path), hash_prefix=expected)
        digest = sha256_file(path)
        if expected and digest != expected:
            sys.exit(f"{path}: SHA-256 {digest} does not match the release ({expected}). Delete it and rerun.")
        digests[name] = digest
    return digests


def check_or_write_config(out, config):
    path = out / "run_config.json"
    if path.exists():
        old = json.loads(path.read_text())
        diff = {k: (old.get(k), config[k]) for k in SAME_RUN_KEYS if old.get(k) != config[k]}
        if diff:
            sys.exit(f"{out} holds a run with other settings (old, new): {diff}. Use another --out.")
    atomic_write(path, lambda f: f.write(json.dumps(config, indent=2).encode()))


def load_sample(path):
    with np.load(path, allow_pickle=False) as z:
        return z["rg_nm"], z["bad_frame"], json.loads(str(z["meta"]))


def finished(out, name, sequence, conformations):
    """True if this protein has a complete result; a result for another sequence is an error, not a redo."""
    path = out / "samples" / f"{name}.npz"
    if not path.exists():
        return False
    rg, _, meta = load_sample(path)
    if meta["sequence"] != sequence:
        sys.exit(f"{path} was made for another sequence than the input's. Use another --out.")
    if len(rg) != conformations:
        sys.exit(f"{path} has {len(rg)} conformations, not {conformations}. Use another --out.")
    return True


def rg_from_ensemble(ensemble, sequence, conformations):
    """R_g per conformation (nm) and the bad-frame mask, after checking the ensemble is what was asked for."""
    if ensemble.sequence != sequence:
        raise ValueError("STARLING altered the sequence")
    L = len(sequence)
    maps = ensemble.distance_maps()
    if maps.shape != (conformations, L, L):
        raise ValueError(f"distance maps have shape {maps.shape}, expected {(conformations, L, L)}")
    if not np.isfinite(maps).all():
        raise ValueError("non-finite distances in the distance maps")

    rg_angstrom = np.asarray(ensemble.radius_of_gyration(), dtype=np.float64)     # the maps are in A
    # The same quantity from its definition, in double precision: guards against a change of meaning or unit upstream.
    check = np.array([np.sqrt(np.sum(np.square(m), dtype=np.float64) / (2 * L**2)) for m in maps])
    if not np.allclose(rg_angstrom, check, rtol=1e-4, atol=0):
        raise ValueError("STARLING's R_g disagrees with sqrt(sum d_ij^2 / 2N^2) from the distance maps")

    bad = np.zeros(conformations, dtype=bool)
    bad[ensemble.check_for_errors(remove_errors=False, verbose=False)] = True
    return rg_angstrom / 10, bad


def moments(rg, bad):
    good = rg[~bad]
    q = np.quantile(rg, [0.05, 0.25, 0.5, 0.75, 0.95])
    std = rg.std(ddof=1) if len(rg) > 1 else np.nan
    return dict(n_conformations=len(rg), rg_mean_nm=rg.mean(), rg_sem_nm=std / np.sqrt(len(rg)), rg_std_nm=std,
                rg_skew=stats.skew(rg), rg_exkurt=stats.kurtosis(rg), rg_min_nm=rg.min(), rg_q05_nm=q[0], rg_q25_nm=q[1],
                rg_median_nm=q[2], rg_q75_nm=q[3], rg_q95_nm=q[4], rg_max_nm=rg.max(), n_bad_frames=int(bad.sum()),
                rg_mean_clean_nm=good.mean() if len(good) else np.nan)


def histogram(name, rg):
    counts, edges = np.histogram(rg, bins=np.histogram_bin_edges(rg, bins="fd"))
    density = counts / (counts.sum() * np.diff(edges))
    return pd.DataFrame(dict(name=name, bin_left_nm=edges[:-1], bin_right_nm=edges[1:], count=counts, density=density))


def write_tables(out, table, max_length, failures):
    """summary.csv, histograms.csv and rg_samples.npz, rebuilt from samples/ so they always agree with it."""
    rows, hists, arrays = [], [], {}
    for r in table.itertuples():
        row = dict(name=r.name, N=r.N, status="pending", error="")
        path = out / "samples" / f"{r.name}.npz"
        if r.N > max_length:
            row["status"], row["error"] = "too_long", f"N > {max_length}"
        elif r.name in failures:
            row["status"], row["error"] = "failed", failures[r.name]
        elif path.exists():
            rg, bad, meta = load_sample(path)
            row.update(status="ok", **moments(rg, bad), batch_size=meta["batch_size"], seed=meta["seed"],
                       seconds=meta["seconds"], gpu_peak_mib=meta["gpu_peak_mib"])
            hists.append(histogram(r.name, rg))
            arrays[r.name] = rg
        rows.append(row)

    refs = table.drop(columns=["sequence", "N"])
    summary = pd.DataFrame(rows).merge(refs, on="name", how="left", validate="1:1")
    atomic_write(out / "summary.csv", lambda f: summary.to_csv(f, index=False))
    if hists:
        atomic_write(out / "histograms.csv", lambda f: pd.concat(hists).to_csv(f, index=False))
        atomic_write(out / "rg_samples.npz", lambda f: np.savez(f, **arrays))
    return summary


def main():
    args = parse_args()
    table = load_input(args.input)

    import starling     # slow, and pulls in torch: only after the arguments and the input have been checked
    import torch
    from starling import configs
    from starling.inference import generation

    too_long = table[table.N > configs.MAX_SEQUENCE_LENGTH]
    runnable = table[table.N <= configs.MAX_SEQUENCE_LENGTH]
    if args.names:
        unknown = set(args.names) - set(table.name)
        if unknown:
            sys.exit(f"--names not in {args.input}: {sorted(unknown)}")
        runnable = runnable[runnable.name.isin(args.names)]

    todo = [r for r in runnable.itertuples() if not finished(args.out, r.name, r.sequence, args.conformations)]
    n_done = len(runnable) - len(todo)
    if args.limit is not None:
        todo = todo[:args.limit]

    print(f"{len(table)} proteins in {args.input}; N {table.N.min()} to {table.N.max()}")
    print(f"STARLING {starling.__version__}: {args.conformations} conformations, {args.sampler}, {args.steps} steps, "
          f"{args.ionic_strength} mM, batch {args.batch_size}, device {args.device}")
    print(f"Output {args.out}: {n_done} finished already, {len(todo)} to run now")
    for r in too_long.itertuples():
        print(f"NOT SAMPLED: {r.name} has N = {r.N}, over STARLING's limit of {configs.MAX_SEQUENCE_LENGTH}")
    if args.dry_run:
        return 0

    if args.device.startswith("cuda") and not torch.cuda.is_available():
        sys.exit(f"--device {args.device} but torch sees no GPU. Pass --device cpu to run on the CPU on purpose.")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "samples").mkdir(exist_ok=True)
    if args.save_ensembles:
        (args.out / "ensembles").mkdir(exist_ok=True)

    use_gpu = args.device.startswith("cuda")
    fast = use_gpu and not args.plain
    config = dict(conformations=args.conformations, steps=args.steps, sampler=args.sampler,
                  ionic_strength_mM=args.ionic_strength, seed=args.seed, starling_version=starling.__version__, fast=fast,
                  requested_batch_size=args.batch_size, device=args.device, input=str(args.input),
                  torch_version=torch.__version__, torch_cuda=torch.version.cuda, python=platform.python_version(),
                  gpu=torch.cuda.get_device_name(args.device) if args.device.startswith("cuda") else None,
                  weights_sha256=ensure_weights(torch, configs))
    check_or_write_config(args.out, config)

    sampling = dict(steps=args.steps, sampler=args.sampler, ionic_strength=args.ionic_strength, device=args.device,
                    return_single_ensemble=True, show_progress_bar=False, show_per_step_progress_bar=False)

    def sample(sequences, conformations, batch):
        with torch.autocast("cuda", dtype=torch.float16, enabled=fast):
            return starling.generate(sequences, conformations=conformations, batch_size=batch, **sampling)

    t0 = time.perf_counter()
    if fast:    # not mode="reduce-overhead": after a GPU out-of-memory error its CUDA graphs fail for the whole process
        starling.set_compilation_options(enabled=True, mode="default", fullgraph=False)
    encoder, diffusion = generation.model_manager.get_models(device=torch.device(args.device))
    # STARLING 2.0.2 leaves its models in training mode. There, SequenceEncoder.forward zeroes the ionic strength for a
    # random 20% of its calls, and sampling calls it once per batch, so a fifth of the batches would ignore the salt.
    # This is the only effect of the mode here: the models have no dropout or batch norm.
    diffusion.eval()
    encoder.eval()
    if fast:    # float16 is for the denoiser only: squared distances in the decoded maps would overflow it
        decode = encoder.decode

        def decode_in_float32(x, *a, **k):
            with torch.autocast("cuda", enabled=False):
                return decode(x.float(), *a, **k)

        encoder.decode = decode_in_float32
    sample("GS" * 10, 1, 1)     # compiles; the first protein still pays for its own batch shape
    print(f"Model ready in {time.perf_counter() - t0:.1f} s", flush=True)

    batch_size = args.batch_size    # stays reduced after an out-of-memory error: the cost hardly depends on N
    failures, consecutive, times = {}, 0, []
    exit_code = 0
    try:
        for k, r in enumerate(todo, 1):
            seed = protein_seed(args.seed, r.name)
            try:
                while True:
                    out_of_memory = False
                    torch.manual_seed(seed)
                    if use_gpu:
                        torch.cuda.empty_cache()    # blocks cached since model loading add over 1 GB to the peak
                        torch.cuda.reset_peak_memory_stats(args.device)
                    t0 = time.perf_counter()
                    try:
                        ensemble = sample({r.name: r.sequence}, args.conformations, batch_size)
                    except torch.cuda.OutOfMemoryError:
                        if batch_size == 1:
                            raise
                        out_of_memory = True
                    if not out_of_memory:
                        break
                    # Outside the except block, so the traceback no longer pins the GPU tensors.
                    gc.collect()
                    torch.cuda.empty_cache()
                    batch_size = max(1, batch_size // 2)
                    print(f"GPU out of memory on {r.name}; retrying with batch size {batch_size}", flush=True)
                seconds = time.perf_counter() - t0
                peak = torch.cuda.max_memory_reserved(args.device) / 2**20 if use_gpu else float("nan")

                rg, bad = rg_from_ensemble(ensemble, r.sequence, args.conformations)
                if args.save_ensembles:     # before the sample file, which is what marks the protein as finished
                    ensemble.save(str(args.out / "ensembles" / r.name), compress=True, verbose=False)
                meta = dict(name=r.name, sequence=r.sequence, seed=seed, batch_size=batch_size, seconds=seconds,
                            gpu_peak_mib=peak)
                atomic_write(args.out / "samples" / f"{r.name}.npz", lambda f: np.savez(
                    f, rg_nm=rg, bad_frame=bad, meta=json.dumps(meta)))
                del ensemble
            except Exception as e:
                failures[r.name] = f"{type(e).__name__}: {e}".splitlines()[0]
                with open(args.out / "failures.log", "a") as f:
                    f.write(f"--- {r.name}\n{traceback.format_exc()}\n")
                print(f"[{k}/{len(todo)}] {r.name} FAILED: {failures[r.name]}", file=sys.stderr, flush=True)
                exit_code = 1
                consecutive += 1
                if consecutive >= MAX_CONSECUTIVE_FAILURES:
                    print(f"Stopping after {consecutive} consecutive failures; see failures.log", file=sys.stderr)
                    break
                continue

            consecutive = 0
            times.append(seconds)
            eta = np.mean(times) * (len(todo) - k) / 60
            ref = f"SAXS {r.saxs_rg_nm:.2f}" if hasattr(r, "saxs_rg_nm") else ""
            print(f"[{k}/{len(todo)}] {r.name} N={r.N}: <Rg> {rg.mean():.3f} +- {rg.std(ddof=1) / np.sqrt(len(rg)):.3f} "
                  f"nm, std {rg.std(ddof=1):.3f} {ref}; bad {int(bad.sum())}; {seconds:.1f} s, "
                  f"{peak:.0f} MiB; {eta:.1f} min left", flush=True)
    finally:    # also on Ctrl-C: what has finished is on disk, and the tables should say so
        summary = write_tables(args.out, table, configs.MAX_SEQUENCE_LENGTH, failures)
        print("Status: " + ", ".join(f"{s} {n}" for s, n in summary.status.value_counts().items()))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
