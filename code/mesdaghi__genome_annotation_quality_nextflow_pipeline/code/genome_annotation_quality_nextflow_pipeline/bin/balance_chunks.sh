#!/usr/bin/env bash
# Deterministically split a protein FASTA into length-balanced chunks.
#
# Naive sequential/count-based chunking (N sequences per chunk, in file order)
# can produce chunks whose total residue count -- and therefore InterProScan
# runtime -- varies wildly. Large multi-domain genes (e.g. PKS/NRPS
# megasynthases) tend to sit next to each other in gene-sorted FASTA files
# (secondary metabolite clusters are physically adjacent in the genome), so a
# handful of chunks can end up with several giant outliers while others get
# only short sequences. That's what made one chunk take 6 hours while its
# neighbour finished almost instantly.
#
# This instead uses greedy longest-processing-time-first (LPT) bin-packing:
# sequences are sorted by length descending, and each is assigned to whichever
# chunk currently has the smallest total residue count. Ties are broken by
# header text and chunk index, so the same input and --chunks value always
# produce the exact same assignment (no randomness, no wall-clock or PID
# seeding) -- runs are reproducible and diffable.
set -euo pipefail

input=""
chunks=""
outdir=""
prefix="chunk_"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --chunks) chunks="$2"; shift 2 ;;
        --outdir) outdir="$2"; shift 2 ;;
        --prefix) prefix="$2"; shift 2 ;;
        *) echo "error: unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$input" ]] || { echo "error: --input is required" >&2; exit 1; }
[[ -n "$chunks" ]] || { echo "error: --chunks is required" >&2; exit 1; }
[[ -n "$outdir" ]] || { echo "error: --outdir is required" >&2; exit 1; }

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT
mkdir -p "$workdir/records"

# One file per record (rec dir/<idx>.fa, original header/sequence lines
# verbatim) plus a length<TAB>header<TAB>index table for sorting/packing.
awk -v recdir="$workdir/records" '
/^>/ {
    idx++
    header[idx] = $0
    len[idx] = 0
    file = recdir "/" idx ".fa"
    print $0 > file
    next
}
idx > 0 {
    file = recdir "/" idx ".fa"
    print $0 > file
    line = $0
    sub(/^[ \t\r]+/, "", line)
    sub(/[ \t\r]+$/, "", line)
    len[idx] += length(line)
}
END {
    for (i = 1; i <= idx; i++) print len[i] "\t" header[i] "\t" i
}
' "$input" >"$workdir/records.tsv"

n_records=$(wc -l <"$workdir/records.tsv")
if [[ "$n_records" -eq 0 ]]; then
    echo "error: no sequences found in $input" >&2
    exit 1
fi

n_chunks=$chunks
((n_chunks > n_records)) && n_chunks=$n_records
((n_chunks < 1)) && n_chunks=1

# Longest-first, deterministic tiebreak on header text.
sort -t $'\t' -k1,1nr -k2,2 "$workdir/records.tsv" >"$workdir/sorted.tsv"

totals=()
for ((b = 0; b < n_chunks; b++)); do totals[b]=0; done

: >"$workdir/assign.tsv"
while IFS=$'\t' read -r len header idx; do
    best=0
    for ((b = 1; b < n_chunks; b++)); do
        # smallest running total wins; ties -> lowest chunk index (deterministic)
        if ((totals[b] < totals[best])); then
            best=$b
        fi
    done
    totals[best]=$((totals[best] + len))
    printf '%s\t%s\n' "$idx" "$best" >>"$workdir/assign.tsv"
done <"$workdir/sorted.tsv"

mkdir -p "$outdir"
sort -t $'\t' -k2,2n -k1,1n "$workdir/assign.tsv" >"$workdir/assign_by_bin.tsv"

n_seqs=()
for ((b = 0; b < n_chunks; b++)); do n_seqs[b]=0; done

cur_bin=-1
out_file=""
while IFS=$'\t' read -r idx bin; do
    if [[ "$bin" != "$cur_bin" ]]; then
        cur_bin=$bin
        out_file=$(printf '%s/%s%05d.fa' "$outdir" "$prefix" "$((cur_bin + 1))")
        : >"$out_file"
    fi
    cat "$workdir/records/${idx}.fa" >>"$out_file"
    n_seqs[bin]=$((n_seqs[bin] + 1))
done <"$workdir/assign_by_bin.tsv"

# Bins that got zero records (more --chunks than sequences never happens here
# since n_chunks is capped at n_records, but keep this for safety) still need
# an (empty) output file.
for ((b = 0; b < n_chunks; b++)); do
    out_file=$(printf '%s/%s%05d.fa' "$outdir" "$prefix" "$((b + 1))")
    [[ -e "$out_file" ]] || : >"$out_file"
done

echo "Wrote $n_chunks chunk(s) to $outdir"
for ((b = 0; b < n_chunks; b++)); do
    echo "  $(printf '%s%05d.fa' "$prefix" "$((b + 1))"): ${n_seqs[b]:-0} sequences, ${totals[b]} residues"
done
