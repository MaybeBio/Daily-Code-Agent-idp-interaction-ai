#!/usr/bin/env bash
# fetch_ped_capped.sh — like fetch_ped.sh, but with a hard cap on the total bytes downloaded.
# Each archive is streamed through `head -c <remaining budget + 1>`, so no transfer can go past the cap. A transfer that
# hits it is deleted and the run stops.
#   bash fetch_ped_capped.sh OUT_DIR PLAN_FILE CAP_MB
# PLAN_FILE lines: "<entry> <ensemble ids to try, comma-separated>", e.g. "PED00501 e002,e001"
set -u
OUT="${1:?output dir}"; PLAN="${2:?plan file}"; CAP=$(( ${3:?cap in MB} * 1000000 ))
mkdir -p "$OUT"
used=0
while read -r e list; do
  [ -z "$e" ] && continue
  f="$OUT/${e}_ensembles.tar.gz"
  if [ -e "$f" ]; then echo "$e: already present"; continue; fi
  ok=0
  for ens in ${list//,/ }; do
    remain=$(( CAP - used ))
    curl -s -f --max-time 900 "https://deposition.proteinensemble.org/v1/entries/$e/ensembles/$ens/ensemble-pdb" \
      | head -c $(( remain + 1 )) > "$f"
    st=("${PIPESTATUS[@]}")
    size=$(stat -f %z "$f" 2>/dev/null || stat -c %s "$f")
    if (( size > remain )); then
      rm -f "$f"; echo "$e $ens: STOPPED — would exceed the ${3} MB cap (used $((used / 1000000)) MB so far)"; exit 2
    fi
    if (( st[0] == 0 && size > 0 )); then
      used=$(( used + size )); ok=1
      mkdir -p "$OUT/$e" && tar -xzf "$f" -C "$OUT/$e" && echo "$e $ens: ok ($((size / 1000)) kB; total $((used / 1000000)) MB)"
      break
    fi
    rm -f "$f"; echo "$e $ens: not available (curl exit ${st[0]})"
  done
  (( ok )) || echo "$e: no ensemble downloaded"
done < "$PLAN"
echo "done: $((used / 1000000)) MB downloaded (cap ${3} MB)"
