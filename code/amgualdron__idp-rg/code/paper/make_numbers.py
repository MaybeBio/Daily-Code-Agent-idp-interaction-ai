"""Every number the write-ups quote, computed from results/, model/ and data/ (no MD frames needed):

    paper/numbers.tex     one LaTeX macro per number: \\val{key} formatted, \\valraw{key} at full precision
    MODEL.md, README.md   the tables between <!-- numbers:NAME --> and <!-- /numbers --> markers

    python paper/make_numbers.py            # rewrite numbers.tex and the tables
    python paper/make_numbers.py --check    # exit 1 if any of them is out of date (run by tests/)

Keys (all listed in numbers.tex, grouped by source):
    data.*        protein counts, split sizes, CE quantiles, training ranges of the inputs
    net.* cf.*    the released network / closed form: constants, coefficients, and the scores in results/metrics.json
                  as <model>.<set>.<metric>, set = val, test, length, sticky, benchmark (closed form also cv)
    <model>.<set>.by.<group>.<metric>   per-protein breakdowns of val and test: N bins, CE bins, compact / not compact
    abl.*         results/ablations.csv; rows not yet run are "--" and abl.complete is "no"
    cethr.*       results/ce_threshold_cv.csv, the choice of the CE-threshold slope c1
    agree.*       the methods table: mean absolute difference in <R_g> on the benchmark between every method (the two
                  released models, data/benchmark_methods.csv, HPS1) and SAXS / each simulation, plain and per family
    signtest.*    is a predictor closer to SAXS than the simulation it was trained on, across families
                  (network vs HPS1, ALBATROSS vs Mpipi-GG)
    fig.*         paper/figure_metrics.json
Formatting: percentages to 2 decimals (1 decimal from 10%), signed for biases; NLL to 4 decimals; skewness errors to 3;
counts with thousands separators; everything else to 4 significant figures (coefficients 6).
"""
import argparse
import json
import math
import pathlib
import re
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS, MODEL, DATA = ROOT / "results", ROOT / "model", ROOT / "data"
TEX = HERE / "numbers.tex"
MODELS = {"net": "network", "cf": "closed_form"}
N_BINS = {"N_le50": (0, 50), "N_51_100": (50, 100), "N_101_200": (100, 200), "N_201_300": (200, 300), "N_gt300": (300, 1e9)}
INPUT_RANGES = ("N", "scd_pos", "scd_neg", "shd", "contact_energy", "fcr", "patch_neg", "patch_neu", "frac_sticky",
                "mean_hydro")
ABL_METRICS = ("n_parameters", "nll", "mean_err_pct", "bias_pct", "std_err_pct", "skew_mae", "w1_pct_median",
               "compact_mean_err_pct")
ABL_NETWORK = ["reference"] + [f"-{x}" for x in ("scd_pos", "scd_neg", "shd", "contact_energy", "fcr", "patch_neg",
                                                   "patch_neu", "frac_sticky")] + ["ln N only"]
ABL_ARCH = ("128x2", "1024x2", "512x1", "512x3")
ABL_CF = ("reference", "-scd_pos", "-scd_neg", "-mean_hydro", "-h", "no ln N interactions", "constant CE threshold",
          "free ln A")
ASYN = ("MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGA"
        "PQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA")                       # alpha-synuclein, the worked example
SET_NAMES = {"val": "validation", "test": "test", "cv": "cross-validation", "length": "length", "sticky": "sticky"}
METHODS = {"network": "network_this_work", "closed_form": "closed_form_formula_this_work", "albatross": "ALBATROSS",
           "idpgan": "idpGAN", "zheng": "Zheng", "afrc": "AFRC", "hps1": "HPS1", "mpipi_gg": "Mpipi_GG",
           "calvados2": "CALVADOS2"}                               # column -> key; the keys of the paper's table
REFERENCES = ("saxs", "hps1", "mpipi_gg", "calvados2")
STUDENTS = (("network", "hps1"), ("albatross", "mpipi_gg"))        # predictor, the simulation it was trained on


# ---------------------------------------------------------------- formatting
def slug(s):
    """An ablation name as a key: '-scd_pos' -> 'no_scd_pos', 'ln N only' -> 'ln_N_only', '128x2' -> 'w128x2'."""
    s = re.sub(r"^-", "no_", s).replace(" ", "_")
    return "w" + s if s[0].isdigit() else s


def fmt(key, v):
    """The display form of a number, chosen by its key; strings are escaped for TeX."""
    if isinstance(v, str):
        return v.replace("\\", "").replace("_", r"\_").replace("%", r"\%").replace("&", r"\&").replace("#", r"\#")
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "--"
    leaf = key.rsplit(".", 1)[-1]
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)) or ((leaf == "n" or leaf.startswith("n_")) and float(v).is_integer()):
        return f"{int(v):,}".replace(",", "{,}")
    v = float(v)
    if "bias" in leaf:
        s = f"{abs(v):.2f}" if abs(v) < 10 else f"{abs(v):.1f}"
        return ("\\ensuremath{+}" if v >= 0 else "\\ensuremath{-}") + s
    if leaf.endswith("pct") or "_pct_" in leaf or any(g.endswith("_pct") for g in key.split(".")[:-1]):
        s = f"{v:.2f}" if abs(v) < 10 else f"{v:.1f}"
    elif key.startswith("cethr") and leaf.endswith("nll"):
        s = f"{v:.5f}"
    elif leaf.startswith("nll") or leaf.endswith("_nll") or leaf == "objective_end":
        s = f"{v:.4f}"
    elif key.startswith("ex.") and leaf.endswith("_nm"):
        s = f"{v:.4f}"
    elif key.startswith(("agree.", "signtest.")) and leaf.endswith("_nm"):
        s = f"{v:.3f}"
    elif key.startswith("signtest.") and leaf.startswith("deviation_corr"):
        s = f"{v:.2f}"
    elif key.startswith("signtest.") and leaf.endswith("_p"):
        s = f"{v:.2g}"
    elif leaf.endswith("_nm"):
        s = f"{v:.3g}"
    elif "skew" in leaf and ("mae" in leaf or "err" in leaf):
        s = f"{v:.3f}"
    elif key.startswith("data.range."):                             # 3 significant figures: 3.80, -0.255, 490
        s = str(int(v)) if v.is_integer() else f"{v:#.3g}".rstrip(".")
    elif key.startswith("data.ce_q") or leaf.startswith("md_mean_"):
        s = f"{v:.3f}"
    elif key.startswith(("cf.coef", "cf.c0", "net.law", "net.off", "net.input", "data.nu_F", "data.ln_a_F", "ex.")):
        s = f"{v:.6g}"
    else:
        s = f"{v:.4g}"
    if "e" in s:                                                     # 1.572e-13 -> 1.572\times10^{-13}
        a, b = s.split("e")
        s = f"\\ensuremath{{{a}\\times10^{{{int(b)}}}}}"
    return s.replace("-", "\\ensuremath{-}", 1) if s.startswith("-") else s


def md(v, key=""):
    """The Markdown form: the TeX form without TeX."""
    return fmt(key, v).replace("{,}", ",").replace("\\ensuremath{-}", "−").replace("\\ensuremath{+}", "+").replace("\\_", "_")


# ---------------------------------------------------------------- the numbers
def per_protein_breakdown(N, kind, name, P, q75, q90):
    p = RESULTS / f"{kind}_{name}.csv"
    if not p.exists():
        return
    d = pd.read_csv(p).merge(P[["id", "contact_energy"]], on="id")
    groups = {k: (d.N > lo) & (d.N <= hi) for k, (lo, hi) in N_BINS.items()}
    groups.update(CE_le_q75=d.contact_energy <= q75, CE_q75_q90=(d.contact_energy > q75) & (d.contact_energy < q90),
                  CE_ge_q90=d.contact_energy >= q90, compact=d.compact, not_compact=~d.compact, all=d.N > 0)
    for g, m in groups.items():
        e = d[m]
        k = f"{name}.by.{g}"
        N[f"{k}.n"] = int(m.sum())
        N[f"{k}.nll"] = e.nll.mean()
        N[f"{k}.mean_err_pct"] = e.mean_err_pct.abs().mean()
        N[f"{k}.bias_pct"] = e.mean_err_pct.mean()
        N[f"{k}.std_err_pct"] = e.std_err_pct.abs().mean()
        N[f"{k}.std_bias_pct"] = e.std_err_pct.mean()
        N[f"{k}.skew_mae"] = e.skew_err.abs().mean()
        N[f"{k}.w1_pct_median"] = e.w1_pct.median()
    a = d.mean_err_pct.abs()
    N[f"{name}.mean_err_pct_median"] = a.median()
    N[f"{name}.mean_err_pct_p90"] = a.quantile(0.9)
    N[f"{name}.mean_err_pct_max"] = a.max()
    N[f"{name}.within5_pct"] = 100 * (a <= 5).mean()
    N[f"{name}.std_err_pct_median"] = d.std_err_pct.abs().median()
    N[f"{name}.std_err_pct_p90"] = d.std_err_pct.abs().quantile(0.9)
    for c in ("mean", "std", "skew"):
        N[f"{name}.r_{c}"] = float(np.corrcoef(d[c], d[f"md_{c}"])[0, 1])
    N[f"{name}.nu_min"], N[f"{name}.nu_median"], N[f"{name}.nu_max"] = d.nu.min(), d.nu.median(), d.nu.max()
    N[f"{name}.md_mean_min"], N[f"{name}.md_mean_max"] = d.md_mean.min(), d.md_mean.max()


def methods_comparison(N, b):
    """agree.<method>.*: mean absolute difference in <R_g> (nm) from SAXS and from each simulation, over the
    benchmark proteins and per family (each protein weighted by 1 / its family's size). signtest.<student>_<teacher>.*:
    absolute errors against SAXS averaged within each family, and a two-sided sign test across families (k = families
    where the predictor's error is the smaller); also the correlation of (predictor - simulation) with
    (simulation - SAXS), over proteins and over family means."""
    from scipy.stats import binomtest
    rg = pd.DataFrame(dict(saxs=b.saxs_rg_nm, hps1=b.hps1_rg_nm))
    for kind in ("network", "closed_form"):
        d = pd.read_csv(RESULTS / f"{kind}_benchmark.csv")
        assert d.name.tolist() == b.name.tolist(), f"{kind}_benchmark.csv rows do not match data/benchmark.csv"
        rg[kind] = d["mean"]
    o = pd.read_csv(DATA / "benchmark_methods.csv")
    assert o.name.tolist() == b.name.tolist(), "benchmark_methods.csv rows do not match data/benchmark.csv"
    for c in METHODS:
        if c not in rg:
            rg[c] = o[f"{c}_rg_nm"]
    for c, key in METHODS.items():
        k = f"agree.{key}"
        N[f"{k}.n"], N[f"{k}.n_families"] = int(rg[c].notna().sum()), b.family[rg[c].notna()].nunique()
        for ref in REFERENCES:
            if ref != c:
                e = (rg[c] - rg[ref]).abs()
                N[f"{k}.mae_vs_{ref}_nm"] = e.mean()
                N[f"{k}.family_weighted_mae_vs_{ref}_nm"] = e.groupby(b.family).mean().mean()
    for s, t in STUDENTS:
        k = f"signtest.{METHODS[s]}_{METHODS[t]}".replace("_this_work", "")
        g = pd.DataFrame(dict(family=b.family, ae_s=(rg[s] - rg.saxs).abs(), ae_t=(rg[t] - rg.saxs).abs(),
                              dev=rg[s] - rg[t], ff=rg[t] - rg.saxs)).dropna()
        fam = g.groupby("family").mean()
        n_fam, ties = len(fam), int((fam.ae_s == fam.ae_t).sum())
        kf, kp = int((fam.ae_s < fam.ae_t).sum()), int((g.ae_s < g.ae_t).sum())
        N[f"{k}.n_proteins"], N[f"{k}.n_families"] = len(g), n_fam
        N[f"{k}.mae_student_vs_saxs_nm"], N[f"{k}.mae_teacher_vs_saxs_nm"] = g.ae_s.mean(), g.ae_t.mean()
        N[f"{k}.family_weighted_mae_student_vs_saxs_nm"] = fam.ae_s.mean()
        N[f"{k}.family_weighted_mae_teacher_vs_saxs_nm"] = fam.ae_t.mean()
        N[f"{k}.families_student_closer"], N[f"{k}.sign_test_p"] = kf, binomtest(kf, n_fam - ties, 0.5).pvalue
        N[f"{k}.proteins_student_closer"] = kp
        N[f"{k}.protein_level_sign_test_p"] = binomtest(kp, len(g), 0.5).pvalue
        N[f"{k}.deviation_corr_proteins"] = float(np.corrcoef(g.ff, g.dev)[0, 1])
        N[f"{k}.deviation_corr_families"] = float(np.corrcoef(fam.ff, fam.dev)[0, 1])


def collect():
    N = {}
    P = pd.read_csv(DATA / "proteins.csv")
    tr = P[P.split == "train"]
    q75, q90 = tr.contact_energy.quantile(0.75), tr.contact_energy.quantile(0.90)
    sys.path.insert(0, str(ROOT))
    from idprg import data as D                                     # noqa: E402

    # data
    N["data.n_proteins"] = len(P)
    for s in ("train", "val", "test"):
        N[f"data.n_{s}"] = int((P.split == s).sum())
    N["data.n_clusters"] = P.cluster.nunique()
    b = pd.read_csv(DATA / "benchmark.csv")
    N["data.n_benchmark"], N["data.n_benchmark_families"] = len(b), b.family.nunique()
    N["data.benchmark.N_min"], N["data.benchmark.N_max"] = b.N.min(), b.N.max()
    e = 100 * (b.hps1_rg_nm / b.saxs_rg_nm - 1)
    N["data.benchmark.hps1_err_vs_saxs_pct"], N["data.benchmark.hps1_bias_vs_saxs_pct"] = e.abs().mean(), e.mean()
    N["data.benchmark.hps1_r_vs_saxs"] = float(np.corrcoef(b.hps1_rg_nm, b.saxs_rg_nm)[0, 1])
    N["data.ce_q75"], N["data.ce_q90"] = q75, q90
    N["data.nu_F"], N["data.ln_a_F"] = D.NU_F, D.LN_A_F
    for split in ("length", "sticky"):
        fit, ev = D.split(P, split)
        N[f"data.{split}.n_fit"], N[f"data.{split}.n_eval"] = len(fit), len(ev)
    fit, ev = D.split(P, "length")
    N["data.length.max_ratio"] = P.N.iloc[ev].max() / P.N.iloc[fit].max()
    for c in INPUT_RANGES:
        x = tr[c]
        N[f"data.range.{c}.min"], N[f"data.range.{c}.max"] = x.min(), x.max()
        N[f"data.range.{c}.p01"], N[f"data.range.{c}.p99"] = x.quantile(0.01), x.quantile(0.99)
    N["data.frac_scd_pos_pct"] = 100 * (tr.scd_pos > 0).mean()
    N["data.ce_slope_lnN"] = np.polyfit(tr.log_seq_length, tr.contact_energy, 1)[0]
    N["data.frac_scd_neg_pct"] = 100 * (tr.scd_neg < 0).mean()

    # the released models
    spec = json.loads((MODEL / "network.json").read_text())
    for k in ("width", "depth", "epochs", "seed", "n_fit", "n_parameters"):
        if k in spec:
            N[f"net.{k}"] = spec[k]
    N["net.n_inputs"] = len(spec["inputs"])
    N["net.off.b_S"], N["net.off.u_0"] = spec["b_S"], spec["u_0"]
    for k, mu, sd in zip(spec["inputs"], spec["input_mean"], spec["input_sd"]):
        N[f"net.input.{k}.mean"], N[f"net.input.{k}.sd"] = mu, sd
    for k, v in zip(spec["law_covariates"], spec["law"]):
        N[f"net.law.{k.replace('^', '').replace('/', '_over_')}"] = v
    N["net.sha256"] = spec["sha256"]["network.eqx"][:8]
    cf = json.loads((MODEL / "closed_form.json").read_text())
    N["cf.c0"], N["cf.c1"] = cf["ce_threshold"]["c0"], cf["ce_threshold"]["c1"]
    for r in ("ln_Sm", "u"):
        N[f"cf.coef.{r}.1"] = cf[r]["intercept"]
        for k, v in cf[r]["coef"].items():
            N[f"cf.coef.{r}.{k.replace('*', '.').replace('^2', '2').replace('/', '_over_')}"] = v
    for k, v in cf["ln_A"].items():
        N[f"cf.coef.ln_A.{k.replace('/', '_over_')}"] = v
    N["cf.n_parameters"], N["cf.n_coefficients"] = cf["n_parameters"], cf["n_coefficients"]
    N["cf.objective_end"] = cf["fit"]["objective_end"]

    # scores
    m = json.loads((RESULTS / "metrics.json").read_text())
    for short, kind in MODELS.items():
        for s, summary in m.get(kind, {}).items():
            for k, v in summary.items():
                N[f"{short}.{s}.{k}"] = v
        p = RESULTS / f"{kind}_benchmark.csv"
        if p.exists():
            d = pd.read_csv(p)
            e = (100 * (d["mean"] / d.hps1_rg_nm - 1)).abs()
            N[f"{short}.benchmark.mean_err_vs_hps1_pct_median"] = e.median()
            N[f"{short}.benchmark.mean_err_vs_hps1_pct_max"] = e.max()
            N[f"{short}.benchmark.r_vs_hps1"] = float(np.corrcoef(d["mean"], d.hps1_rg_nm)[0, 1])
            N[f"{short}.benchmark.r_vs_saxs"] = float(np.corrcoef(d["mean"], d.saxs_rg_nm)[0, 1])
        sub = {}
        for s in ("val", "test", "cv"):
            per_protein_breakdown(sub, kind, s, P, q75, q90)
        N.update({f"{short}.{k}": v for k, v in sub.items()})
    methods_comparison(N, b)

    # worked example: alpha-synuclein, through the closed form
    import idprg                                                    # noqa: E402
    from idprg import closed_form as CF, descriptors as DS, distribution as DI
    x = DS.table([ASYN])
    lnA, lnSm, u = (float(a[0]) for a in CF.parameters(cf, x))
    nu = float(DI.nu_of(u))
    g, dl = DI.exponents(nu)
    q, _ = DS.encode(ASYN)
    mean, std, skew, _ = (float(a) for a in DI.moments(lnA, lnSm, nu))
    ex = dict(N=int(x.N[0]), scd=DS.scd(q), scd_neg=x.scd_neg[0], mean_hydro=x.mean_hydro[0], ce=x.contact_energy[0],
              ce_threshold=N["cf.c0"] + N["cf.c1"] * x.log_seq_length[0], lnSm=lnSm, Sm_nm=np.exp(lnSm), u=u, nu=nu,
              delta=dl, g=g, lnA=lnA, mean_nm=mean, std_nm=std, skew=skew,
              hps1_mean_nm=float(b.hps1_rg_nm[b.sequence == ASYN].iloc[0]))
    for k, v in ex.items():
        N[f"ex.asyn.{k}"] = v
    for n_ in (20, 184, 490):
        N[f"cf.ce_threshold_at_{n_}"] = N["cf.c0"] + N["cf.c1"] * np.log(n_)

    # the CE-threshold choice
    c = pd.read_csv(RESULTS / "ce_threshold_cv.csv")
    for r in c.itertuples():
        k = "cethr.c1_" + f"{r.c1:g}".replace("-", "m").replace(".", "p")         # c1 = -0.5 -> cethr.c1_m0p5
        for col in c.columns[1:]:
            N[f"{k}.{col}"] = getattr(r, col)
        N[f"{k}.cv_nll_gain"] = -r.cv_nll_minus_constant                # how much lower than with c1 = 0

    # ablations
    a = pd.read_csv(RESULTS / "ablations.csv")
    complete = True
    rows = [("network", x, s) for s in ("val", "length", "sticky") for x in ABL_NETWORK]
    rows += [("network", x, "val") for x in ABL_ARCH] + [("closed_form", x, "val") for x in ABL_CF]
    for model, name, s in rows:
        k = f"abl.{'net' if model == 'network' else 'cf'}.{s}.{slug(name)}"
        r = a[(a.model == model) & (a.ablation == name) & (a.set == s)]
        seed = r[r.seed == (42 if model == "network" else 0)]
        complete &= len(seed) == 1
        for col in ABL_METRICS:
            N[f"{k}.{col}"] = float(seed[col].iloc[0]) if len(seed) else None
        if name == "reference" and model == "network":
            N[f"{k}.n_seeds"] = len(r)
            for col in ABL_METRICS[1:]:
                N[f"{k}.seed_mean.{col}"], N[f"{k}.seed_sd.{col}"] = r[col].mean(), r[col].std(ddof=1)
    N["abl.complete"] = bool(complete)
    N["abl.n_done"], N["abl.n_total"] = sum(1 for mdl, x, s in rows if ((a.model == mdl) & (a.ablation == x)
                                                                       & (a.set == s)).any()), len(rows)

    # figures
    def flat(prefix, v):
        if isinstance(v, dict):
            for k, w in v.items():
                flat(f"{prefix}.{k}", w)
        elif isinstance(v, list):
            for i, w in enumerate(v):
                flat(f"{prefix}.{i}", w)
        else:
            N[prefix] = v
    p = HERE / "figure_metrics.json"
    if p.exists():
        flat("fig", json.loads(p.read_text()))
    return N


# ---------------------------------------------------------------- numbers.tex
HEADER = r"""%% Generated by paper/make_numbers.py from results/, model/ and data/ -- do not edit; rerun the script instead.
%% \val{key}: the formatted number; \valraw{key}: full precision. An unknown key is a LaTeX error.
\makeatletter
\DeclareRobustCommand{\val}[1]{\ifcsname num@#1\endcsname\csname num@#1\endcsname\else
  \PackageError{numbers}{Undefined number key `\detokenize{#1}'}{Rerun paper/make_numbers.py or fix the key.}%
  \textbf{??\detokenize{#1}}\fi}
\DeclareRobustCommand{\valraw}[1]{\ifcsname numraw@#1\endcsname\csname numraw@#1\endcsname\else
  \PackageError{numbers}{Undefined number key `\detokenize{#1}'}{}\textbf{??\detokenize{#1}}\fi}
\newcommand{\num@set}[3]{\expandafter\def\csname num@#1\endcsname{#2}\expandafter\def\csname numraw@#1\endcsname{#3}}
"""


def tex(N):
    out, group = [HEADER], None
    for k in N:
        g = k.split(".")[0]
        if g != group:
            out.append(f"\n%% ---- {g}")
            group = g
        v = N[k]
        if isinstance(v, (int, np.integer)) and not isinstance(v, (bool, np.bool_)):
            raw = str(int(v))
        elif isinstance(v, (float, np.floating)) and math.isfinite(v):
            raw = repr(float(v))
        else:
            raw = fmt(k, v)
        out.append(f"\\num@set{{{k}}}{{{fmt(k, v)}}}{{{raw}}}")
    return "\n".join(out) + "\n\\makeatother\n"


# ---------------------------------------------------------------- Markdown tables
def table(head, rows):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(r) + " |" for r in rows])


def score_row(N, label, p):
    g = lambda k: md(N.get(f"{p}.{k}"), k)
    return [label, g("nll"), g("mean_err_pct") + "%", g("bias_pct") + "%", g("std_err_pct") + "%", g("skew_mae"),
            f"{g('w1_pct_median')}% [{g('w1_pct_p90')}]", g("compact_mean_err_pct") + "%"]


SCORE_HEAD = ["", "NLL", "⟨R_g⟩ err", "bias", "σ err", "skew err", "W₁ median [p90]", "compact ⟨R_g⟩ err"]


def by_rows(N, s, groups):
    rows = []
    for short, label in (("net", "network"), ("cf", "closed form")):
        rows.append([label] + [f"{md(N[f'{short}.{s}.by.{g}.mean_err_pct'], 'pct')}% / "
                               f"{float(N[f'{short}.{s}.by.{g}.std_err_pct']):.1f}%" for g in groups])
    return rows


def blocks(N):
    n = lambda k: md(N[k], k)
    B = {}
    B["accuracy"] = table(["", "proteins", "⟨R_g⟩ error", "σ error", "W₁", "NLL per frame"], [
        [f"{label}, {SET_NAMES[s]}", n(f"{short}.{s}.n"), n(f"{short}.{s}.mean_err_pct") + "%",
         n(f"{short}.{s}.std_err_pct") + "%", n(f"{short}.{s}.w1_pct_median") + "%", n(f"{short}.{s}.nll")]
        for s in ("val", "test") for short, label in (("net", "network"), ("cf", "closed form"))])
    for s in ("val", "test"):
        B[s] = table(SCORE_HEAD, [score_row(N, label, f"{short}.{s}") for short, label in
                                  (("net", "network"), ("cf", "closed form"))])
        B[f"{s}-length"] = table(
            ["N"] + [f"{lab} ({n(f'net.{s}.by.{g}.n')})" for g, lab in
                     zip(N_BINS, ("≤ 50", "51–100", "101–200", "201–300", "> 300"))], by_rows(N, s, N_BINS))
        B[f"{s}-ce"] = table(
            ["CE", f"≤ q75 = {N['data.ce_q75']:.3f} ({n(f'net.{s}.by.CE_le_q75.n')})",
             f"q75–q90 ({n(f'net.{s}.by.CE_q75_q90.n')})",
             f"≥ q90 = {N['data.ce_q90']:.3f} ({n(f'net.{s}.by.CE_ge_q90.n')})"],
            by_rows(N, s, ("CE_le_q75", "CE_q75_q90", "CE_ge_q90")))
    rows = []
    for s in ("length", "sticky"):
        rows.append([f"**{s}** ({n(f'net.{s}.n')} scored, {n(f'net.{s}.n_compact')} compact)"] + [""] * 7)
        rows += [score_row(N, label, f"{short}.{s}") for short, label in (("net", "network"), ("cf", "closed form"))]
    B["extrapolation"] = table(SCORE_HEAD, rows)
    b = lambda short, k: n(f"{short}.benchmark.{k}")
    B["benchmark"] = table(["⟨R_g⟩", "against HPS1: error / bias", "against SAXS: error / bias / MAE"], [
        [label, f"{b(short, 'mean_err_vs_hps1_pct')}% / {b(short, 'bias_vs_hps1_pct')}%",
         f"{b(short, 'mean_err_vs_saxs_pct')}% / {b(short, 'bias_vs_saxs_pct')}% / {b(short, 'mae_vs_saxs_nm')} nm"]
        for short, label in (("net", "network"), ("cf", "closed form"))] +
        [["HPS1 simulation itself", "—", f"{n('data.benchmark.hps1_err_vs_saxs_pct')}% / "
          f"{n('data.benchmark.hps1_bias_vs_saxs_pct')}% / {n('net.benchmark.hps1_mae_vs_saxs_nm')} nm"]])

    def abl(k, cols):
        out = []
        for c in cols:
            v = N.get(f"{k}.{c}")
            out.append("pending" if v is None else md(v, c) + ("%" if "pct" in c else ""))
        return out
    cols = ("nll", "mean_err_pct", "std_err_pct")
    rows = []
    ref = "abl.net.{}.reference.seed_{}.{}"
    rows.append(["reference (seeds 42–46, mean ± sd)", md(N["abl.net.val.reference.n_parameters"], "n_parameters")] +
                [f"{md(N[ref.format('val', 'mean', c)], c)} ± {float(N[ref.format('val', 'sd', c)]):.{4 if c == 'nll' else 2}f}"
                 + ("" if c == "nll" else "%") for c in cols] +
                [f"{md(N[ref.format(s, 'mean', 'mean_err_pct')], 'pct')} ± {float(N[ref.format(s, 'sd', 'mean_err_pct')]):.2f}%"
                 for s in ("length", "sticky")])
    for x in ABL_NETWORK[1:]:
        k = f"abl.net.{{}}.{slug(x)}"
        rows.append([f"`{x}`", abl(k.format("val"), ("n_parameters",))[0]] + abl(k.format("val"), cols) +
                    abl(k.format("length"), ("mean_err_pct",)) + abl(k.format("sticky"), ("mean_err_pct",)))
    for x in ABL_ARCH:
        k = f"abl.net.val.{slug(x)}"
        rows.append([f"width × depth {x}", abl(k, ("n_parameters",))[0]] + abl(k, cols) + ["", ""])
    B["ablations-network"] = table(["network, seed 42", "parameters", "val NLL", "val ⟨R_g⟩ err", "val σ err",
                                    "length ⟨R_g⟩ err", "sticky ⟨R_g⟩ err"], rows)
    B["ablations-closed-form"] = table(
        ["closed form", "parameters", "val NLL", "val ⟨R_g⟩ err", "val σ err", "val W₁ median"],
        [[f"`{x}`" if x.startswith("-") else x] + abl(f"abl.cf.val.{slug(x)}", ("n_parameters",) + cols + ("w1_pct_median",))
         for x in ABL_CF])
    names = {"N": "N", "scd_pos": "s₊", "scd_neg": "s₋", "shd": "SHD", "contact_energy": "CE", "fcr": "FCR",
             "mean_hydro": "⟨λ⟩"}
    rg = lambda c, a, b: f"{n(f'data.range.{c}.{a}')} to {n(f'data.range.{c}.{b}')}"
    B["domain"] = table(["", "range", "1–99%"], [[names.get(c, c), rg(c, "min", "max"), rg(c, "p01", "p99")]
                                                  for c in INPUT_RANGES])
    return B


def fill(text, B, path):
    def sub(mt):
        name = mt.group(1)
        if name not in B:
            raise KeyError(f"{path}: unknown numbers block {name!r}")
        return f"<!-- numbers:{name} -->\n{B[name]}\n<!-- /numbers -->"
    return re.sub(r"<!-- numbers:([\w-]+) -->.*?<!-- /numbers -->", sub, text, flags=re.S)


def outputs():
    """{path: the content make_numbers would write there}, and the numbers."""
    N = collect()
    B = blocks(N)
    new = {TEX: tex(N)}
    for p in (ROOT / "MODEL.md", ROOT / "README.md"):
        new[p] = fill(p.read_text(), B, p.name)
    return new, N


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="only check that the files are up to date")
    args = ap.parse_args()
    new, N = outputs()
    stale = [p for p, t in new.items() if not p.exists() or p.read_text() != t]
    if args.check:
        if stale:
            sys.exit("out of date (run python paper/make_numbers.py): " + ", ".join(p.name for p in stale))
        print("numbers up to date")
        return
    for p in stale:
        p.write_text(new[p])
        print(f"wrote {p.relative_to(ROOT)}")
    print(f"{len(N)} numbers; ablations {N['abl.n_done']}/{N['abl.n_total']} done")


if __name__ == "__main__":
    main()
