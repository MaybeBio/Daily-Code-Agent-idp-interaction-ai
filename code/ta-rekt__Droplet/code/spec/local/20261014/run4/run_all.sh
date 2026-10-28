#!/bin/bash
# Rouse escape/return FPTs vs s at fixed phi_inf (indices 13,14,15 of the 32-point
# sweep, around 0.45 c_sat).  One chain per s in PARAMS.rouse.s_values, ALL launched
# at once (each `( ... ) &` is a background subshell; `wait` blocks until every chain
# is done).  Each chain runs
#   bvp_potential_interp.jl -> fisher_rate_from_fpt_solver.jl   (cheap, s-dependent)
#   runJL ... rouse_mfpt_attractive.jl 3                         (3 concurrent Rouse runs)
# so length(s_values) x 3 Rouse processes run side by side.  Every value -- s list,
# output root, run length -- comes from param.jl / the Rouse script; this file only
# launches.
#
#   ./spec/local/20261009/run_all.sh smoke                 # local check: 3 s, ~5 min
#   nohup ./spec/local/20261009/run_all.sh production > run_all.log 2>&1 &
#
# Output, with OUT_ROOT = data/20261009 (production) or data/20261009/smoke:
#   OUT_ROOT/s<s>/bvp_potential_interp.json, fisher_rate_from_fpt_solver.json
#   OUT_ROOT/s<s>/rouse/output_cpu{1,2,3}.json
#   OUT_ROOT/s<s>/run.log
# Retrieve:  scp -r <cluster>:<repo>/data/20261009 data/

set -u
[ $# -eq 1 ] || { echo "usage: $0 smoke|production"; exit 1; }
export RUN_MODE=$1
export OPENBLAS_NUM_THREADS=1      # concurrent stages must not each grab every core

# This script's directory, from the launching command, and the repo root above it
# (the first ancestor holding runJL).  Work from the root -- runJL and sim/sweep.jl
# take repo-relative paths -- with `here` the script's directory relative to it.
here_abs=$(cd "$(dirname "$0")" && pwd)
repo=$here_abs
until [ -x "$repo/runJL" ]; do
    [ "$repo" = / ] && { echo "no runJL in any parent of $here_abs"; exit 1; }
    repo=$(dirname "$repo")
done
cd "$repo"
here=${here_abs#"$repo"/}

# The s list and output root, read from param.jl (ROUSE_S unset so it cannot reject).
read -r root svals <<< "$(env -u ROUSE_S julia --project=@. -e "include(\"$here/param.jl\");
    print(OUT_ROOT, \" \", join(Int.(PARAMS.rouse.s_values), \" \"))")"
[ -n "${svals:-}" ] || { echo "could not read s_values from $here/param.jl"; exit 1; }
echo "$(date '+%F %T')  RUN_MODE=$RUN_MODE  root=$root  s=[$svals]"

for s in $svals; do
    mkdir -p "$root/s$s"
    ( export ROUSE_S=$s
      julia --project=@. $here/bvp_potential_interp.jl &&
      julia --project=@. $here/fisher_rate_from_fpt_solver.jl &&
      ./runJL sim/sweep.jl $here/rouse_mfpt_attractive.jl 30 "$root/s$s/rouse"
    ) > "$root/s$s/run.log" 2>&1 \
        && echo "$(date '+%F %T')  s=$s done" \
        || echo "$(date '+%F %T')  s=$s FAILED -- see $root/s$s/run.log" &
done
wait
