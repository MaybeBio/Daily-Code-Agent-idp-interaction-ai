#!/usr/bin/env bash
# Download the CATH latest-release classification files (v4.4.0 when this analysis was run) and, for a
# sensitivity check, the S35 list of v4.3.0 (the release behind CATH-Gene3D 4.3.0 in InterPro). Files
# already present are kept; checksums and the access date are written when anything is downloaded.
set -euo pipefail
WORK=${COVERAGE_WORK:-$(cd "$(dirname "$0")" && pwd)/work}
OUT=${COVERAGE_DOWNLOADS:-$WORK/downloads}/cath
BASE=https://download.cathdb.info/cath/releases/latest-release/cath-classification-data
B43=https://download.cathdb.info/cath/releases/all-releases/v4_3_0/cath-classification-data
mkdir -p "$OUT"
fetched=0
for f in cath-domain-list.txt cath-domain-list-S35.txt cath-names.txt cath-superfamily-list.txt README-cath-list-file-format.txt README-cath-names-file-format.txt; do
  if [ ! -s "$OUT/$f" ]; then curl -sSfL --retry 3 -o "$OUT/$f" "$BASE/$f"; fetched=1; fi
done
for f in cath-domain-list-S35-v4_3_0.txt cath-superfamily-list-v4_3_0.txt; do
  if [ ! -s "$OUT/$f" ]; then curl -sSfL --retry 3 -o "$OUT/$f" "$B43/$f"; fetched=1; fi
done
if [ "$fetched" = 1 ]; then
  curl -sSfL "$BASE/" | sed 's/<[^>]*>/ /g' > "$OUT/listing.txt"
  ( cd "$OUT"; sha256sum cath-*.txt README-*.txt > SHA256SUMS; date -u +%FT%TZ > ACCESS_DATE )
fi
( cd "$OUT"; sha256sum -c --quiet SHA256SUMS )
head -3 "$OUT/cath-domain-list.txt"
