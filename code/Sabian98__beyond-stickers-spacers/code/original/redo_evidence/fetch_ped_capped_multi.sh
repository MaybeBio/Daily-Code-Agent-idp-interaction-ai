#!/usr/bin/env bash
# fetch_ped_capped_multi.sh — like fetch_ped_capped.sh, but downloads EVERY listed ensemble of an entry (one archive per
# entry/ensemble pair) instead of the first available one. Hard cap on the total bytes, enforced with `head -c`.
#   bash fetch_ped_capped_multi.sh OUT_DIR PLAN_FILE CAP_MB
# PLAN_FILE lines: "<entry> <ensemble id>", e.g. "PED00142 e003". Archives go to OUT_DIR/<entry>/<entry><ens>.tar.gz and
# are extracted into OUT_DIR/<entry>/<ens>/.
set -u
OUT="${1:?output dir}"; PLAN="${2:?plan file}"; CAP=$(( ${3:?cap in MB} * 1000000 ))
used=0
while read -r e ens; do
  [ -z "$e" ] && continue
  d="$OUT/$e"; f="$d/${e}${ens}.tar.gz"
  mkdir -p "$d"
  if [ -s "$f" ]; then echo "$e $ens: already present"; continue; fi
  remain=$(( CAP - used ))
  curl -s -f --max-time 1800 "https://deposition.proteinensemble.org/v1/entries/$e/ensembles/$ens/ensemble-pdb" \
    | head -c $(( remain + 1 )) > "$f"
  st=("${PIPESTATUS[@]}")
  size=$(stat -f %z "$f" 2>/dev/null || stat -c %s "$f")
  if (( size > remain )); then
    rm -f "$f"; echo "$e $ens: STOPPED, would exceed the ${3} MB cap (used $((used / 1000000)) MB)"; exit 2
  fi
  if (( st[0] == 0 && size > 0 )); then
    used=$(( used + size ))
    mkdir -p "$d/$ens" && tar -xzf "$f" -C "$d/$ens" && echo "$e $ens: ok ($((size / 1000)) kB; total $((used / 1000000)) MB)"
  else
    rm -f "$f"; echo "$e $ens: not available (curl exit ${st[0]})"
  fi
done < "$PLAN"
echo "done: $((used / 1000000)) MB downloaded (cap ${3} MB)"
