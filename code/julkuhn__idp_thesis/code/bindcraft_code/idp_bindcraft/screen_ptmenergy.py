#!/usr/bin/env python
"""Score / rank candidate binders against a target with pTMEnergy (BECraft, arXiv:2505.21241).

This is the *screening* use of pTMEnergy, as opposed to the design-time loss: given a target
and a pool of candidate binder sequences, predict each complex and rank by interface energy.
Lower energy = more likely a real binder. No labels, no co-crystal structure needed.

The energy is registered as a ColabDesign loss with weight 0, so AlphaFold computes and logs
it on every predict() without it influencing anything. That also means you can add it to an
existing ipTM run purely as a recorded metric.

Sequences come from a FASTA, or from any CSV with a 'Sequence' column -- including BindCraft's
own final_design_stats.csv / mpnn_design_stats.csv.

Example:
  python screen_ptmenergy.py \
    --target_pdb .../crops/APP_P05067_gapcrop0632_0700.pdb --chains A \
    --sequences .../output/APP_P05067_gapcrop0632_0700/final_design_stats.csv \
    --out /iopsstor/scratch/cscs/jkuhn/app_screen.csv
"""
import argparse, csv, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from colabdesign import mk_afdesign_model, clear_mem
from functions.colabdesign_utils import add_i_ptm_energy_loss

AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_sequences(path, seq_col="Sequence", name_col="Design", limit=None):
    """Return [(name, sequence)] from a FASTA or a CSV holding a sequence column."""
    out, seen = [], set()
    if path.lower().endswith((".fa", ".fasta", ".faa")):
        name, buf = None, []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if name and buf:
                        out.append((name, "".join(buf)))
                    name, buf = line[1:].split()[0] or f"seq{len(out)}", []
                elif line:
                    buf.append(line)
        if name and buf:
            out.append((name, "".join(buf)))
    else:
        with open(path) as fh:
            rdr = csv.DictReader(fh)
            if seq_col not in (rdr.fieldnames or []):
                sys.exit(f"error: no '{seq_col}' column in {path}\n"
                         f"       columns: {(rdr.fieldnames or [])[:12]} ...")
            for i, row in enumerate(rdr):
                s = (row.get(seq_col) or "").strip().upper()
                if not s:
                    continue
                out.append((row.get(name_col) or f"row{i}", s))

    clean = []
    for name, s in out:
        s = s.strip().upper()
        bad = set(s) - AA
        if bad:
            print(f"  skipping {name}: non-standard residues {sorted(bad)}")
            continue
        if s in seen:                      # identical sequences score identically
            continue
        seen.add(s)
        clean.append((name, s))
    return clean[:limit] if limit else clean


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--target_pdb", required=True)
    p.add_argument("--chains", default="A", help="target chain(s) in the pdb")
    p.add_argument("--sequences", required=True, help="FASTA, or CSV with a Sequence column")
    p.add_argument("--out", required=True, help="output CSV, ranked by energy ascending")
    p.add_argument("--af_params_dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "params"))
    p.add_argument("--tm_lambda", type=float, default=1.0,
                   help="1 = full TM-weighting (default), 0 = plain -logsumexp")
    p.add_argument("--models", type=int, nargs="+", default=[0],
                   help="AF2 model indices to average over (default: 0)")
    p.add_argument("--num_recycles", type=int, default=3)
    p.add_argument("--multimer", action="store_true", default=True)
    p.add_argument("--seq_col", default="Sequence")
    p.add_argument("--name_col", default="Design")
    p.add_argument("--limit", type=int, default=None, help="score only the first N candidates")
    a = p.parse_args()

    seqs = read_sequences(a.sequences, a.seq_col, a.name_col, a.limit)
    if not seqs:
        sys.exit("error: no usable sequences found")
    lengths = sorted({len(s) for _, s in seqs})
    print(f"loaded {len(seqs)} unique candidates; binder length(s): {lengths}")
    if len(lengths) > 1:
        print("  note: lengths differ, so a separate model is compiled per length")

    clear_mem()
    rows, models = [], {}
    t0 = time.time()
    for i, (name, seq) in enumerate(seqs, 1):
        L = len(seq)
        if L not in models:
            m = mk_afdesign_model(protocol="binder", data_dir=a.af_params_dir,
                                  use_multimer=a.multimer, num_recycles=a.num_recycles)
            m.prep_inputs(pdb_filename=a.target_pdb, chain=a.chains, binder_len=L)
            # weight 0 -> computed and logged, but contributes nothing
            add_i_ptm_energy_loss(m, weight=0.0, tm_lambda=a.tm_lambda)
            models[L] = m
            print(f"  compiled model for binder_len={L}")
        m = models[L]

        per_model = []
        for mn in a.models:
            m.predict(seq=seq, models=[mn], num_recycles=a.num_recycles, verbose=False)
            per_model.append(dict(m.aux["log"]))

        def avg(k):
            v = [d[k] for d in per_model if k in d and d[k] is not None]
            return float(np.mean(v)) if v else None

        if "i_ptm_energy" not in per_model[0]:
            sys.exit("error: i_ptm_energy absent from aux['log'] -- the weight-0 loss did not "
                     "register. Check add_i_ptm_energy_loss in functions/colabdesign_utils.py.")

        rows.append({"name": name, "length": L,
                     "i_ptm_energy": avg("i_ptm_energy"),
                     "i_ptm": avg("i_ptm"), "ptm": avg("ptm"),
                     "plddt": avg("plddt"), "i_pae": avg("i_pae"), "pae": avg("pae"),
                     "sequence": seq})
        if i % 10 == 0 or i == len(seqs):
            print(f"  {i}/{len(seqs)} scored ({(time.time()-t0)/i:.1f} s/candidate)")

    # lower energy = better binder
    rows.sort(key=lambda r: (r["i_ptm_energy"] is None, r["i_ptm_energy"]))
    for rank, r in enumerate(rows, 1):
        r["rank"] = rank

    cols = ["rank", "name", "length", "i_ptm_energy", "i_ptm", "ptm", "plddt",
            "i_pae", "pae", "sequence"]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})

    print(f"\nwrote {a.out}  ({len(rows)} candidates, ranked by pTMEnergy ascending)")
    print(f"{'rank':>4} {'name':<42} {'energy':>9} {'i_ptm':>7} {'plddt':>7}")
    for r in rows[:10]:
        print(f"{r['rank']:>4} {str(r['name'])[:42]:<42} {r['i_ptm_energy']:>9.3f} "
              f"{(r['i_ptm'] or 0):>7.3f} {(r['plddt'] or 0):>7.3f}")

    e = np.array([r["i_ptm_energy"] for r in rows if r["i_ptm_energy"] is not None])
    ip = np.array([r["i_ptm"] for r in rows if r["i_ptm"] is not None])
    if len(e) > 2 and len(e) == len(ip) and e.std() > 0 and ip.std() > 0:
        print(f"\nSpearman-free check: Pearson r(energy, i_ptm) = {np.corrcoef(e, ip)[0,1]:+.3f}")
        print("  (strongly negative is expected: low energy ~ high i_ptm)")


if __name__ == "__main__":
    main()
