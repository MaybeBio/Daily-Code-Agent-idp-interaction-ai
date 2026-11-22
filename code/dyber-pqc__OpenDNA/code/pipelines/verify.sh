#!/bin/bash
# Does splitting a job across tasks and merging it back give the same answers?
#
# This runs what the Nextflow and Snakemake pipelines run -- check, split,
# predict each part, merge -- and compares the result against predicting the
# whole file in one go. It is the claim both pipelines rest on, and it is
# checkable without installing a workflow engine.
#
#   pipelines/verify.sh <models-dir> <fasta>
#
# The models directory should be one produced by `opendna stage`, so the run
# can use --offline the way a compute node would.
#
# Exit 0 if the merged result is byte-identical to the whole-file run.

set -u
MODELS="${1:?usage: verify.sh <models-dir> <fasta>}"
FASTA="${2:?usage: verify.sh <models-dir> <fasta>}"
OPENDNA="${OPENDNA:-opendna}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Batching groups sequences by length within whatever file they arrive in, so
# a different division of the input is a different set of batch-mates and a
# different amount of padding. 0 pads nothing, which is what makes the result
# independent of the part count -- and that independence is what this checks.
BATCH_TOKENS="${BATCH_TOKENS:-0}"
PARTS="${PARTS:-4}"

fail() { echo "FAILED: $*" >&2; exit 1; }

echo "== the installation verifies =="
"$OPENDNA" --models "$MODELS" models --check "$MODELS/models.lock" >/dev/null \
    || fail "models --check (exit $?)"
"$OPENDNA" --models "$MODELS" --offline selftest >/dev/null \
    || fail "selftest (exit $?)"
echo "   ok"

echo "== the whole file in one run =="
"$OPENDNA" --models "$MODELS" --offline --batch-tokens "$BATCH_TOKENS" \
    predict --task disorder -i "$FASTA" -o "$WORK/whole.csv" 2>/dev/null \
    || fail "whole-file predict (exit $?)"
echo "   $(grep -cv '^#' "$WORK/whole.csv") lines"

echo "== split into $PARTS parts, predict each, merge =="
"$OPENDNA" split -i "$FASTA" --out-dir "$WORK/parts" --parts "$PARTS" >/dev/null \
    || fail "split (exit $?)"
for part in "$WORK"/parts/part-*.fasta; do
    "$OPENDNA" --models "$MODELS" --offline --batch-tokens "$BATCH_TOKENS" \
        predict --task disorder -i "$part" -o "${part%.fasta}.csv" 2>/dev/null \
        || fail "predict on $part (exit $?)"
done
# Sorted, because the merge preserves the order it is given and the parts are
# named so that lexical order is numerical order.
# shellcheck disable=SC2046
"$OPENDNA" merge --output "$WORK/merged.csv" $(ls "$WORK"/parts/part-*.csv | sort) 2>/dev/null \
    || fail "merge (exit $?)"
echo "   $(grep -cv '^#' "$WORK/merged.csv") lines"

echo "== the two agree =="
# The provenance comments name the input file, which differs by construction,
# so the comparison is over the data rows. Those are what a result is.
grep -v '^#' "$WORK/whole.csv" | sort > "$WORK/a.txt"
grep -v '^#' "$WORK/merged.csv" | sort > "$WORK/b.txt"
if ! diff -q "$WORK/a.txt" "$WORK/b.txt" >/dev/null; then
    echo "   the merged result differs from the whole-file run:" >&2
    diff "$WORK/a.txt" "$WORK/b.txt" | head -10 >&2
    if [ "$BATCH_TOKENS" != "0" ]; then
        echo >&2
        echo "   BATCH_TOKENS is $BATCH_TOKENS. Padding depends on which sequences share" >&2
        echo "   a batch, so dividing the input differently changes the last digit of" >&2
        echo "   some values. This is expected above 0 and is why the pipelines set 0." >&2
    fi
    exit 1
fi
echo "   identical over $(wc -l < "$WORK/a.txt") rows"

echo "== a mixture is refused =="
# Half the parts relabelled as though they had run on another device: the
# merge must stop rather than produce a file nobody can describe.
mixed="$WORK/mixed"
mkdir -p "$mixed"
i=0
for csv in "$WORK"/parts/part-*.csv; do
    base="$(basename "${csv%.csv}")"
    cp "$csv" "$mixed/$base.csv"
    manifest="${csv%.csv}.manifest.json"
    if [ "$i" -eq 0 ]; then
        sed 's/"device": "[^"]*"/"device": "SOMETHING ELSE"/' "$manifest" > "$mixed/$base.manifest.json"
    else
        cp "$manifest" "$mixed/$base.manifest.json"
    fi
    i=$((i + 1))
done
# shellcheck disable=SC2046
if "$OPENDNA" merge --output "$WORK/mixed.csv" $(ls "$mixed"/part-*.csv | sort) >/dev/null 2>&1; then
    fail "merge accepted parts from two different devices"
fi
echo "   ok"

echo
echo "Chunking does not change the answers at --batch-tokens $BATCH_TOKENS."
