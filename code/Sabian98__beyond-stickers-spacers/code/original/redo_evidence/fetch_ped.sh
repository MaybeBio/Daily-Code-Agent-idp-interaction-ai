#!/usr/bin/env bash
# fetch_ped.sh — download PED ensemble PDBs (one ensemble per entry) into OUT_DIR, using the endpoint documented at
# https://proteinensemble.org/api : https://deposition.proteinensemble.org/v1/entries/{id}/ensembles/{ens}/ensemble-pdb
# (returns a .tar.gz with the multi-model PDB; the older /api/v1/ path now returns 403).
#   bash fetch_ped.sh OUT_DIR PLAN_FILE        # PLAN_FILE lines: "<entry> <ensemble_id> <models>" (see ped_fetch_plan.txt)
set -u
OUT="${1:?output dir}"; PLAN="${2:?plan file}"
mkdir -p "$OUT"
while read -r e ens models; do
  [ -z "$e" ] && continue
  f="$OUT/${e}_ensembles.tar.gz"
  if [ -e "$f" ]; then echo "$e: already present"; continue; fi
  code=$(curl -s --max-time 900 -o "$f" -w "%{http_code}" "https://deposition.proteinensemble.org/v1/entries/$e/ensembles/$ens/ensemble-pdb")
  if [ "$code" = "200" ]; then mkdir -p "$OUT/$e" && tar -xzf "$f" -C "$OUT/$e" && echo "$e $ens: ok ($(du -h "$f" | cut -f1))"
  else rm -f "$f"; echo "$e $ens: HTTP $code"; fi
done < "$PLAN"
