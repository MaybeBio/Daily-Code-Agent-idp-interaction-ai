#!/usr/bin/env bash
#
# Fig. 5 panel (d): Rouse-chain simulations of the June toy model, whose
# end-to-end distance histogram is compared with exp(-U_tot)/Z.
#
# Run from the repo root:
#     bash figure_production/fig5/run_fig5_rouse.sh [ncpus]
#
# Eight runs -- two concentrations x four seeds -- of 5e7 steps each (spec:
# fig5/rouse_toy.jl). Each run holds its full distance trace in RAM (8 bytes/step
# = 0.4 GB), so ncpus is also the memory knob. Must divide 8. Default: picked from
# RAM and cores.
#
# No input files: the toy model's potential is analytic.
# Output: data/YYYYMMDD/fig5_rouse_toy/output_cpu*.json (gitignored), today's date.

set -euo pipefail

if [[ ! -f Project.toml || ! -f runJL ]]; then
    echo "ERROR: run this from the repo root (the directory holding Project.toml)." >&2
    exit 1
fi
command -v julia >/dev/null || { echo "ERROR: julia not found -- see figure_production/README.md, Fig. 4 step 1." >&2; exit 1; }

OUT=data/$(date +%Y%m%d)/fig5_rouse_toy

# ---- ncpus: from the argument, else the largest divisor of 8 that fits ------
if [[ $# -ge 1 ]]; then
    NCPUS=$1
else
    if [[ "$(uname)" == Darwin ]]; then
        CORES=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
        MEM_GB=$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))
    else
        CORES=$(nproc); MEM_GB=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 / 1024 ))
    fi
    NCPUS=1
    for n in 2 4 8; do
        # 0.4 GB per run, keep half the RAM free
        if (( n <= CORES && n * 4 <= MEM_GB * 10 / 2 )); then NCPUS=$n; fi
    done
    echo "==> $CORES performance cores, ${MEM_GB} GB RAM -> ncpus = $NCPUS"
fi
case $NCPUS in 1|2|4|8) ;; *) echo "ERROR: ncpus must divide 8 (got $NCPUS)" >&2; exit 1;; esac

echo "==> julia $(julia --version | awk '{print $3}'), instantiating the project"
julia --project=@. -e 'using Pkg; Pkg.instantiate()'

echo "==> 8 Rouse runs on $NCPUS CPUs -> $OUT   (started $(date +%H:%M))"
# bash, not ./runJL: file sync can drop runJL's executable bit
bash runJL sim/sweep.jl figure_production/fig5/rouse_toy.jl "$NCPUS" "$OUT"

# ---- verify ----------------------------------------------------------------
python3 - "$OUT" <<'PY'
import json, glob, sys, collections
out = sys.argv[1]
runs = []
for f in sorted(glob.glob(out + "/output_cpu*.json")):
    for d in json.load(open(f)).values():
        runs.extend(d.values())
by = collections.defaultdict(list)
for r in runs:
    by[r["alpha"]].append(r)
print()
for a, rs in sorted(by.items()):
    n = sum(len(r["distance_traj"]) for r in rs)
    print("  phi_inf = %.4g   %d runs   %8d distance samples   %.1f min/run"
          % (a, len(rs), n, sum(r["runtime_min"] for r in rs) / len(rs)))
bad = len(runs) != 8 or len(by) != 2 or any(len(r["distance_traj"]) == 0 for r in runs)
print("\n  " + ("PROBLEM: expected 8 runs over 2 concentrations, each with a distance trace."
                if bad else "All 8 runs present, both concentrations, distance traces kept."))
sys.exit(1 if bad else 0)
PY
echo
echo "Done ($(date +%H:%M)). Next: python figure_production/fig5/prepare_fig5_data.py"
