#!/usr/bin/env bash
# Relaunch the closed-form u-term tests; see README.md here. From any directory, with IDPRG_FRAMES set.
#
#   run.sh smoke                     plumbing check on a few proteins, under a minute
#   run.sh check                     must pass before anything else counts
#   run.sh u                         step 1: dS, dH, hN, dSH, 3 folds each
#   run.sh drop MODEL PRODUCTS       step 2: one removal, e.g.  run.sh drop dS scd_pos   then  run.sh drop dS scd_pos,h
#   run.sh final MODEL [PRODUCTS]    step 3: length and contact-energy sets
#   run.sh validation MODEL [PRODUCTS]   step 4: the one validation score
#
# JOBS (default 2) = fits at once. A fold fit needs about 4 GB of memory and a full-training fit about 6 GB.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${IDPRG_FRAMES:?set IDPRG_FRAMES to the MD frame store}"
JOBS="${JOBS:-2}"
fit() { printf '%s\n' "$@" | xargs -P "$JOBS" -I{} python -u "$HERE/fit.py" {}; }
folds() { for f in 0 1 2; do echo "$1:standard:$f$2"; done; }
suffix() { [ -n "${1:-}" ] && echo ":$1" || true; }

check() {
    python -u "$HERE/fit.py" check:standard:0
    python - "$HERE" <<'EOF'
import json, pathlib, sys
here = pathlib.Path(sys.argv[1])
got = json.loads((here.parents[1] / "runs/u_terms/check.standard.0.json").read_text())["objective_end"]
want = json.loads((here / "reference/s_over_N.standard.0.json").read_text())["objective_end"]
print(f"check: objective {got:.6f}, reference {want:.6f}")
if abs(got - want) > 5e-5:
    sys.exit("check FAILED: fit.py does not reproduce the reference s+/N fit, so its other results cannot be trusted")
print("check passed")
EOF
}

case "${1:-}" in
    smoke)
        python -u "$HERE/fit.py" --smoke dS:standard:0 dS:standard:1 dS:standard:2 dH:length:-1 hN:sticky:-1 dSH:standard:0:scd_pos,h
        python "$HERE/compare.py" --dir "$HERE/../../runs/u_terms_smoke"
        ;;
    check) check ;;
    u)
        check
        fit $(for m in dS dH hN dSH; do folds "$m" ""; done)
        ;;
    drop) fit $(folds "$2" ":$3") ;;
    final)
        s="$(suffix "${3:-}")"
        fit "$2:length:-1$s" "$2:sticky:-1$s"
        ;;
    validation)
        s="$(suffix "${3:-}")"
        python -u "$HERE/fit.py" --validation "$2:standard:-1$s"
        ;;
    *) sed -n 2,11p "$0"; exit 1 ;;
esac
