#!/usr/bin/env bash
# Swiss-Prot flat file (its release in reldate.txt) and the RCSB 30%-identity polymer-entity clusters.
# Files already present are kept, so the release recorded in reldate.txt stays the one downloaded.
set -euo pipefail
WORK=${COVERAGE_WORK:-$(cd "$(dirname "$0")" && pwd)/work}
D=${COVERAGE_DOWNLOADS:-$WORK/downloads}
mkdir -p "$D/uniprot" "$D/rcsb"
U=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete
if [ ! -s "$D/uniprot/uniprot_sprot.dat.gz" ]; then
  curl -sSfL "$U/reldate.txt" -o "$D/uniprot/reldate.txt"
  curl --http1.1 -sSfL --retry 3 -o "$D/uniprot/uniprot_sprot.dat.gz.part" "$U/uniprot_sprot.dat.gz"
  mv "$D/uniprot/uniprot_sprot.dat.gz.part" "$D/uniprot/uniprot_sprot.dat.gz"
  date -u +%FT%TZ > "$D/reference_sets.ACCESS_DATE"
fi
if [ ! -s "$D/rcsb/clusters-by-entity-30.txt" ]; then
  curl -sSfL --retry 3 -D "$D/rcsb/clusters30.headers" -o "$D/rcsb/clusters-by-entity-30.txt" \
   https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt
fi
gzip -t "$D/uniprot/uniprot_sprot.dat.gz"
echo "swissprot entries: $(zcat "$D/uniprot/uniprot_sprot.dat.gz" | /bin/grep -c '^ID ')"
echo "rcsb clusters: $(wc -l < "$D/rcsb/clusters-by-entity-30.txt")"
