#!/bin/bash
# test_faster_pipeline.sh - end-to-end test of `run.sh --faster` on a built system.
#
# Usage (GPU node, openmm_faster env):
#   bash utils/test_faster_pipeline.sh [system_dir] [gpu] [--keep]
#
#   system_dir  a folder where build.sh has run (step3_input.parm7/.rst7, restraints/,
#               step4.*_equilibration.inp, step5_production.inp); default: current folder
#   gpu         GPU index for run.sh (default 0)
#   --keep      keep the test folder even if everything passed
#
# The system folder is NOT touched: everything runs in a temporary copy
# (system_dir/test_faster_pipeline_XXXX, removed on success). The copy gets the
# utils/ of THIS repository and the system's own .inp files with only the
# lengths shortened (minimisation 200 steps, 2000 steps per equilibration, two
# production chunks of 2500 steps); every other setting is the system's own.
#
# It runs the real pipeline, run.sh --faster (minimisation, restrained
# equilibration, production chunk 1 from the .rst, chunk 2 from the .chk), then
# one more chunk WITHOUT --faster to test the engine switch, and checks:
#   - every step wrote its .rst (equilibration) / .dcd + .rst + .chk (production)
#   - every .out says which engine ran, and the faster runs say what was swapped
#   - temperatures finite and within 25 K of the target, no NaN
#   - chunk 2 continued from chunk 1's checkpoint (step counter carried on)
#   - the engine switch restarted from the portable .rst, not the .chk
# Exit code 0 = all PASS. A few minutes for a ~50k-atom system.
set -uo pipefail

SYS=.; GPU=0; KEEP=no
for arg in "$@"; do
    case "$arg" in
        --keep) KEEP=yes ;;
        [0-9]*) GPU="$arg" ;;
        *) SYS="$arg" ;;
    esac
done
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SYS="$(cd "$SYS" && pwd)" || { echo "ERROR: no folder $SYS"; exit 1; }
for f in step3_input.parm7 step3_input.rst7 step5_production.inp; do
    [ -f "$SYS/$f" ] || { echo "ERROR: $SYS/$f missing - run build.sh first"; exit 1; }
done
ls "$SYS"/step4.*_equilibration.inp >/dev/null 2>&1 || { echo "ERROR: no step4.*_equilibration.inp in $SYS"; exit 1; }
python -c "import nv_openmm_plugin" 2>/dev/null || { echo "ERROR: cuPotential plugin not found: conda activate openmm_faster"; exit 1; }

T=$(mktemp -d "$SYS/test_faster_pipeline_XXXX")
echo "test folder: $T"
cp "$SYS"/step3_input.* "$T"/
[ -d "$SYS/restraints" ] && cp -r "$SYS/restraints" "$T"/
cp "$SYS"/step4.*_equilibration.inp "$SYS"/step5_production.inp "$T"/
cp -r "$REPO/utils" "$T"/
cd "$T"

# shorten lengths only
sed -i -E 's/^(mini_nstep\s*=\s*)[0-9]+/\1200/; s/^(nstep\s*=\s*)[0-9]+/\12000/; s/^(nstout\s*=\s*)[0-9]+/\1500/; s/^(nstdcd\s*=\s*)[0-9]+/\11000/' step4.*_equilibration.inp
DT=$(awk '/^dt[[:space:]]*=/{print $3}' step5_production.inp)
TEMP=$(awk '/^temp[[:space:]]*=/{print $3}' step5_production.inp)
CHUNK_NS=$(python -c "print(2500 * $DT / 1000)")
sed -i -E "s/^(nstep\s*=\s*)[0-9]+/\12500/; s/^(nstout\s*=\s*)[0-9]+/\1500/; s/^(nstdcd\s*=\s*)[0-9]+/\11250/; s/^(Target_ns\s*=\s*)[0-9.]+/\1$(python -c "print(2 * $CHUNK_NS)")/" step5_production.inp

FAILS=0
check () {  # name, condition (0 = ok), detail
    if [ "$2" -eq 0 ]; then printf "PASS  %-52s %s\n" "$1" "$3"; else printf "FAIL  %-52s %s\n" "$1" "$3"; FAILS=$((FAILS + 1)); fi
}
last_T () { grep -E '^[0-9]' "$1" 2>/dev/null | tail -1 | cut -f5; }
first_step () { grep -E '^[0-9]' "$1" 2>/dev/null | head -1 | cut -f2; }
last_step () { grep -E '^[0-9]' "$1" 2>/dev/null | tail -1 | cut -f2; }
T_ok () { python -c "import sys, math; t = float(sys.argv[1]); sys.exit(0 if math.isfinite(t) and abs(t - $TEMP) < 25 else 1)" "$1" 2>/dev/null; }

echo "== run.sh $GPU --faster  (dt $DT ps, 2 production chunks of 2500 steps)"
bash utils/run.sh "$GPU" --faster > run_faster.log 2>&1
check "run.sh --faster finished" $? "log: run_faster.log"

for inp in step4.*_equilibration.inp; do
    s=${inp%.inp}
    [ -s "$s.rst" ]; check "$s: .rst written" $? ""
    grep -q "Engine: faster" "$s.out" 2>/dev/null; check "$s: ran with the faster engine" $? "$(grep -m1 -oE 'integrator.*' "$s.out" 2>/dev/null | head -1)"
    t=$(last_T "$s.out"); T_ok "$t"; check "$s: temperature sane" $? "T = $t K (target $TEMP)"
done
for c in 1 2; do
    s=step5_$c
    [ -s "$s.dcd" ] && [ -s "$s.rst" ] && [ -s "$s.chk" ]; check "$s: .dcd .rst .chk written" $? ""
    grep -q "cuPotential LangevinMiddleIntegrator\|stock LangevinIntegrator" "$s.out"; check "$s: components swapped as expected" $? \
        "$(grep -oE 'integrator (=|->) [^(]*' "$s.out" | head -1)"
    t=$(last_T "$s.out"); T_ok "$t"; check "$s: temperature sane" $? "T = $t K (target $TEMP)"
done
grep -qi "nan" step4.*.out step5_*.out; [ $? -ne 0 ]; check "no NaN in any .out" $? ""
f2=$(first_step step5_2.out); l1=$(last_step step5_1.out)
[ -n "$f2" ] && [ -n "$l1" ] && [ "$f2" -gt "$l1" ]; check "chunk 2 continued from chunk 1 checkpoint" $? "chunk 1 ended at step $l1, chunk 2 starts at $f2"
[ "$(cat run_engine.txt 2>/dev/null)" = faster ]; check "run_engine.txt records the engine" $? "$(cat run_engine.txt 2>/dev/null)"

echo "== engine switch: one more chunk WITHOUT --faster"
sed -i -E "s/^(Target_ns\s*=\s*)[0-9.]+/\1$(python -c "print(3 * $CHUNK_NS)")/" step5_production.inp
bash utils/run.sh "$GPU" > run_switch.log 2>&1
check "run.sh (classic) finished" $? "log: run_switch.log"
grep -q "Engine changed (faster -> classic): restarting from step5_2.rst" run_switch.log; check "switch restarted from .rst, not the .chk" $? ""
grep -q "Engine: classic" step5_3.out 2>/dev/null; check "step5_3 ran with the classic engine" $? ""
t=$(last_T step5_3.out); T_ok "$t"; check "step5_3: temperature sane" $? "T = $t K"

echo
if [ "$FAILS" -eq 0 ]; then
    echo "RESULT: run.sh --faster pipeline OK"
    [ "$KEEP" = yes ] && echo "kept: $T" || rm -rf "$T"
    exit 0
fi
echo "RESULT: $FAILS check(s) FAILED - test folder kept for inspection: $T"
exit 1
