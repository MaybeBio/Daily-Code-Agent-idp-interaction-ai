"""Command line: idprg predict | evaluate | train | fit-closed-form | fit-md | ablations  (idprg <command> -h)."""
import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

from . import ablations, closed_form as cf, data, descriptors, distribution as dist, metrics, network, predict
from .data import ROOT

RESULTS = ROOT / "results"
SET_SPLIT = {"val": ("standard", -1), "test": ("test", -1), "length": ("length", -1), "sticky": ("sticky", -1),
             **{f"cv{k}": ("standard", k) for k in range(data.N_FOLDS)}}


def _sequences(args):
    names, seqs = [], []
    if args.fasta:
        for block in pathlib.Path(args.fasta).read_text().split(">")[1:]:
            head, *body = block.splitlines()
            names.append(head.split()[0] if head.split() else f"seq{len(names) + 1}")
            seqs.append("".join(body))
    for s in args.sequence:
        names.append(f"seq{len(names) + 1}")
        seqs.append(s)
    if not seqs:
        sys.exit("give sequences as arguments or with --fasta")
    return names, seqs


def cmd_predict(args):
    names, seqs = _sequences(args)
    out = predict(seqs, closed_form=args.closed_form)
    out.insert(0, "name", names)
    if args.density:
        rows = []
        for n, r in zip(names, out.itertuples()):
            lo, hi = dist.quantiles([1e-3, 1 - 1e-3], r.lnA, np.log(r.S_m), r.nu)
            rg = np.linspace(lo, hi, 200)
            rows.append(pd.DataFrame(dict(name=n, rg=rg, p=dist.pdf(rg, r.lnA, np.log(r.S_m), r.nu))))
        pd.concat(rows).to_csv(args.density, index=False)
    text = out.to_csv(index=False, float_format="%.5g")
    if args.out:
        pathlib.Path(args.out).write_text(text)
    else:
        print(text, end="")


def _model(kind, path):
    """The model's parameter function: descriptor table -> (ln A, ln S_m, u)."""
    if kind == "network":
        net, spec = network.load(*([path] if path else []))
        return lambda t: network.parameters(net, spec, t)
    coef = cf.load(*([path] if path else []))
    return lambda t: cf.parameters(coef, t)


def _update_metrics(out, kind, name, summary):
    p = out / "metrics.json"
    m = json.loads(p.read_text()) if p.exists() else {}
    m.setdefault(kind, {})[name] = summary
    p.write_text(json.dumps(m, indent=1))


def cmd_evaluate(args):
    kind = args.model.replace("-", "_")
    params = _model(kind, args.path)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.set == "benchmark":
        b = data.benchmark()
        lnA, lnSm, u = params(descriptors.table(b.sequence))
        nu = dist.nu_of(u)
        m, s, k, _ = dist.moments(lnA, lnSm, nu)
        d = b[["name", "family", "N", "saxs_rg_nm", "hps1_rg_nm"]].assign(lnA=lnA, lnSm=lnSm, nu=nu, mean=m, std=s, skew=k)
        e_s, e_h = 100 * (m / b.saxs_rg_nm - 1), 100 * (m / b.hps1_rg_nm - 1)
        summary = dict(n=len(d), mean_err_vs_hps1_pct=float(np.abs(e_h).mean()), bias_vs_hps1_pct=float(e_h.mean()),
                       mae_vs_hps1_nm=float(np.abs(m - b.hps1_rg_nm).mean()), mean_err_vs_saxs_pct=float(np.abs(e_s).mean()),
                       bias_vs_saxs_pct=float(e_s.mean()), mae_vs_saxs_nm=float(np.abs(m - b.saxs_rg_nm).mean()),
                       hps1_mae_vs_saxs_nm=float(np.abs(b.hps1_rg_nm - b.saxs_rg_nm).mean()))
    else:
        P = data.proteins()
        _, ev = data.split(P, *SET_SPLIT[args.set])
        d = metrics.score(*params(P.iloc[ev]), P.iloc[ev], data.frames(args.frames))
        summary = metrics.summary(d)
    d.to_csv(out / f"{kind}_{args.set}.csv", index=False)
    _update_metrics(out, kind, args.set, summary)
    print(json.dumps(summary, indent=1))


def cmd_train(args):
    P, M, F = data.proteins(), data.md_fits(), data.frames(args.frames)
    fit, ev = data.split(P, args.split, args.fold)
    inputs = args.inputs.split(",") if args.inputs else list(descriptors.INPUTS)
    t0 = time.time()
    net, spec = network.train(P, M, fit, F, args.split, inputs, args.width, args.depth, args.epochs, args.seed)
    spec["seconds"] = round(time.time() - t0, 1)
    network.save(net, spec, args.out)
    d = metrics.score(*network.parameters(net, spec, P.iloc[ev]), P.iloc[ev], F)
    d.to_csv(pathlib.Path(args.out) / "per_protein.csv", index=False)
    s = metrics.summary(d)
    (pathlib.Path(args.out) / "summary.json").write_text(json.dumps(s, indent=1))
    print(json.dumps(s, indent=1))


def cmd_fit_closed_form(args):
    P, M, F = data.proteins(), data.md_fits(), data.frames(args.frames)
    fit, ev = data.split(P, args.split, args.fold)
    coef, info = cf.fit(P, M, fit, F, c1=args.c1, drop=args.drop, interactions=not args.no_interactions, law=not args.free_lnA)
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "closed_form.json").write_text(json.dumps(dict(coef, **info), indent=1))
    d = metrics.score(*cf.parameters(coef, P.iloc[ev]), P.iloc[ev], F)
    d.to_csv(out / "per_protein.csv", index=False)
    s = metrics.summary(d)
    (out / "summary.json").write_text(json.dumps(s, indent=1))
    print(json.dumps(dict(info, eval=s), indent=1))


def _fit_md_one(task):
    from . import mdfit
    row, x, laws = task
    F = data.frames()
    fr = np.asarray(F[int(row)], float)
    lnA, lnSm, u, nll = mdfit.fit_free(fr)
    out = dict(lnA=lnA, lnSm=lnSm, u=u, nll=nll)
    for s, law in laws.items():
        out[f"lnSm_{s}"], out[f"u_{s}"] = mdfit.fit_law(fr, x["N"], x["mean_hydro"], x["scd_pos"], law, [lnSm, u])
    return out


def cmd_fit_md(args):
    import os
    from concurrent.futures import ProcessPoolExecutor
    if args.frames:
        os.environ[data.FRAMES_ENV] = args.frames
    data.frames()
    P, M = data.proteins(), data.md_fits()
    t = P.set_index("id").loc[M.id].reset_index()
    if args.sample:
        t = t.sample(args.sample, random_state=0).sort_index()
    # the law of each training split, from the free fits of its fit proteins (as in training)
    laws = {s: network.law_start(P, M, P.id.iloc[data.split(P, s)[0]]).tolist() for s in ("standard", "length", "sticky")}
    tasks = [(r.frame_row, dict(N=float(r.N), mean_hydro=r.mean_hydro, scd_pos=r.scd_pos), laws) for r in t.itertuples()]
    with ProcessPoolExecutor(args.jobs) as pool:
        res = list(pool.map(_fit_md_one, tasks, chunksize=8))
    out = pd.DataFrame(res, index=t.index)
    out.insert(0, "id", t.id.to_numpy())
    out.to_csv(args.out, index=False)
    ref = M.set_index("id").loc[out.id]
    for c in out.columns[1:]:
        print(f"{c:15s} max |refit - data/md_fits.csv| = {np.abs(out[c].to_numpy() - ref[c].to_numpy()).max():.1e}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="idprg", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    frames = dict(default=None, help=f"MD frame store (.npy); default: ${data.FRAMES_ENV}")

    p = sub.add_parser("predict", help="<Rg>, sigma, skewness and distribution parameters of sequences")
    p.add_argument("sequence", nargs="*")
    p.add_argument("--fasta")
    p.add_argument("--closed-form", action="store_true", help="the 21-parameter closed form instead of the network")
    p.add_argument("--density", help="also write p(Rg) on a 200-point grid per sequence to this CSV")
    p.add_argument("--out", help="CSV (default: stdout)")
    p.set_defaults(fn=cmd_predict)

    p = sub.add_parser("evaluate", help="score a model on a protein set against the MD frames (or the benchmark)")
    p.add_argument("--model", choices=["network", "closed-form"], default="network")
    p.add_argument("--path", help="model directory (network) or closed_form.json; default: the released model")
    p.add_argument("--set", choices=list(SET_SPLIT) + ["benchmark"], default="val",
                   help="cv0-cv2 are folds of the training split: score them with --path to a model fitted on the other folds")
    p.add_argument("--frames", **frames)
    p.add_argument("--out", default=str(RESULTS))
    p.set_defaults(fn=cmd_evaluate)

    p = sub.add_parser("train", help="train the network on a split and score its evaluation set")
    p.add_argument("--split", choices=["standard", "length", "sticky"], default="standard")
    p.add_argument("--fold", type=int, default=-1, help="standard split: cross-validation fold 0-2 (-1: validation)")
    p.add_argument("--inputs", help="comma-separated network inputs (default: the released model's)")
    p.add_argument("--width", type=int, default=512)
    p.add_argument("--depth", type=int, default=2)
    p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--frames", **frames)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_train)

    p = sub.add_parser("fit-closed-form", help="fit the closed form on a split and score its evaluation set")
    p.add_argument("--split", choices=["standard", "length", "sticky"], default="standard")
    p.add_argument("--fold", type=int, default=-1)
    p.add_argument("--c1", type=float, default=0.25, help="slope of the CE threshold in ln N")
    p.add_argument("--drop", nargs="*", default=[], choices=list(cf.DESC), help="leave these descriptors out")
    p.add_argument("--no-interactions", action="store_true", help="no descriptor x ln N terms in ln S_m")
    p.add_argument("--free-lnA", action="store_true", help="ln A linear in the descriptors instead of the law")
    p.add_argument("--frames", **frames)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_fit_closed_form)

    p = sub.add_parser("fit-md", help="refit the per-protein maximum-likelihood fits (data/md_fits.csv)")
    p.add_argument("--sample", type=int, help="only this many random proteins (a reproduction check)")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--frames", **frames)
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_fit_md)

    p = sub.add_parser("ablations", help="the ablation grid (results/ablations.csv); resumes where it stopped")
    p.add_argument("--only", choices=["network", "closed-form"])
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--frames", **frames)
    p.add_argument("--out", default=str(RESULTS / "ablations.csv"))
    p.set_defaults(fn=ablations.run)

    args = ap.parse_args(argv)
    try:
        args.fn(args)
    except (ValueError, FileNotFoundError) as e:
        sys.exit(f"idprg {args.cmd}: error: {e}")


if __name__ == "__main__":
    main()
