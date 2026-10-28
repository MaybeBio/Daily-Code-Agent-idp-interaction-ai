#!/usr/bin/env bash
#
# Regenerate the BVP solver output behind Fig. 4 (panels B-G).
#
# Run from the repo root:
#     bash figure_production/fig4/run_fig4_solver.sh
#
# Writes eight files, one per run, each into the data/ directory its own
# param.jl names:
#
#     data/20261004/{small_R,medium_R,large_R}/demo_bvp_profiles.json   -> B, C, D, F
#     data/20261007/run{1,2,3,4,5}/demo_bvp_profiles.json               -> E, G
#
# Each run sweeps n_phi = 512 far-field volume fractions, solving the isolated-
# sphere Cahn-Hilliard BVP with branch tracking (thin/thick film, lower-F
# branch returned, the other reported). README step (g) times these at
# 10-20 s each, so the whole script is a couple of minutes.
#
# data/ is gitignored, so nothing here enters the repo.

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. Where are we
# ---------------------------------------------------------------------------
if [[ ! -f Project.toml || ! -d spec/local ]]; then
    echo "ERROR: run this from the repo root (the directory holding Project.toml)." >&2
    echo "       cd ~/NextCloudCNRS/Code/droplet/Droplet && bash figure_production/fig4/$(basename "$0")" >&2
    exit 1
fi

JULIA_VERSION=1.11.1        # Manifest.toml was resolved with exactly this

# ---------------------------------------------------------------------------
# 1. Julia
# ---------------------------------------------------------------------------
# Manifest.toml is pinned to 1.11.1, so that exact version is required. juliaup
# is the installer the README specifies. If `julia` is missing, install it:
#
#     curl -fsSL https://install.julialang.org | sh     # then open a NEW shell
#     juliaup add 1.11.1
#     juliaup default 1.11.1
#
# Pin 1.11.1 rather than the `release` channel: `juliaup update` moves `release`
# to the newest version and the Manifest would no longer resolve.

if ! command -v julia >/dev/null 2>&1; then
    cat >&2 <<'EOF'
ERROR: julia not found.

Install it (macOS/Linux), then open a NEW shell and re-run this script:

    curl -fsSL https://install.julialang.org | sh
    juliaup add 1.11.1
    juliaup default 1.11.1
    julia --version          # -> julia version 1.11.1
EOF
    exit 1
fi

have=$(julia --version | awk '{print $3}')
if [[ "$have" != "$JULIA_VERSION" ]]; then
    echo "WARNING: julia $have, but Manifest.toml was resolved with $JULIA_VERSION." >&2
    echo "         If instantiate or the runs fail, use:  juliaup add $JULIA_VERSION && juliaup default $JULIA_VERSION" >&2
    echo "         (or invoke this script's julia calls as 'julia +$JULIA_VERSION')" >&2
fi

echo "==> julia $have"
echo "==> instantiating the project (no-op if already resolved)"
julia --project=@. -e 'using Pkg; Pkg.instantiate()'

# ---------------------------------------------------------------------------
# 2. The eight runs
# ---------------------------------------------------------------------------
# Each script does include(joinpath(@__DIR__, "param.jl")) and writes to
# PARAMS.sim.output_dir. Nothing is passed on the command line; to change a run,
# edit its param.jl, never the script.

RUNS_R=(small_R medium_R large_R)                 # R = 66.66 / 166.66 / 10000  -> 20 / 50 / 3000 nm
RUNS_H=(run1 run2 run3 run4 run5)                 # omega_b = -0.05 .. -0.0345, run4 = h*

echo
echo "==> 20261004: varying sphere radius R  (panels B, C, D, F)"
for r in "${RUNS_R[@]}"; do
    echo "    -- $r"
    julia --project=@. "spec/local/20261004/$r/demo_bvp.jl"
done

echo
echo "==> 20261007: varying surface attraction h0  (panels E, G)"
for r in "${RUNS_H[@]}"; do
    echo "    -- $r"
    julia --project=@. "spec/local/20261007/$r/demo_bvp.jl"
done

# ---------------------------------------------------------------------------
# 3. Verify
# ---------------------------------------------------------------------------
# The 2026-08-25 run set has a known folder-ordering pitfall (README, step g):
# the spec folders were renamed after the run, so data/run1..run5 had to be
# rotated afterwards. The 2026-10-04 and 2026-10-07 sets do NOT need that --
# each param.jl's output_dir already matches its folder name. This check makes
# that visible rather than assumed: the printed values must be monotonic.

echo
echo "==> written files and their physical parameter"
python3 - <<'PY'
import json, os, sys

expect_R  = [("data/20261004/small_R",   66.66),
             ("data/20261004/medium_R", 166.66),
             ("data/20261004/large_R", 10000.0)]
expect_h  = [("data/20261007/run%d" % i, v) for i, v in
             zip(range(1, 6), [-0.05, -0.04, -0.0365, -0.0358, -0.0345])]

bad = 0
for d, want in expect_R + expect_h:
    p = os.path.join(d, "demo_bvp_profiles.json")
    if not os.path.exists(p):
        print("  MISSING  %s" % p); bad += 1; continue
    s = json.load(open(p))["params"]["surface"]
    got = s["R"] if "small_R" in d or "medium_R" in d or "large_R" in d else s["omega_b"]
    key = "R" if "_R" in d else "omega_b"
    ok = abs(got - want) < 1e-9
    bad += (not ok)
    print("  %-26s %-8s = %-10s %5.1f MB   %s"
          % (d, key, got, os.path.getsize(p) / 1e6, "ok" if ok else "MISMATCH (expected %s)" % want))

print()
if bad:
    print("  %d problem(s) above -- do not use this data until resolved." % bad); sys.exit(1)
print("  All eight runs present and in the expected order (no folder rotation needed).")
PY

echo
echo "Done. Now run figure_production/fig4/fig4.ipynb."
