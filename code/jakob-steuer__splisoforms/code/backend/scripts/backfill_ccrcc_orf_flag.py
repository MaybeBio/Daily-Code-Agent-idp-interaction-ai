#!/usr/bin/env python3
"""
ccRCC folded-ORF QC flag (integrity + transparency).

The AF3 structure of a ccRCC isoform depicts the author-curated folded ORF
(af3_targets.fa = forAF). For a few % of transcripts that ORF's protein identity
diverges from the transcript's genomic-locus gene — e.g. a GAPDH-locus transcript
whose folded ORF best-matches UniProt A4UCT1, not GAPDH. We keep the genomic gene as the objective anchor but
FLAG the divergence and record the folded protein's best UniProt match, so a
structure is never silently presented as a gene it does not depict.

Signal (computed against the FOLDED protein, so it is correct whether or not the
DB sequence has been re-ingested to the folded ORF):
    fold_id = containment 5-mer similarity of af3_targets[structure_md5] to the
              gene's canonical sequence = shared / min(|folded|,|canonical|).
    mismatch = fold_id < THRESHOLD  (default 0.50; matches the classifier's
               'trusted' boundary). Isoforms whose gene has no canonical are left
               unflagged (cannot assess) but still get the folded_best_ref record.

Writes transcript_metadata.folded_orf_mismatch / folded_best_ref /
folded_best_ref_identity, and dumps the flagged list to logs/ for author review.

Usage:  python backend/scripts/backfill_ccrcc_orf_flag.py [--threshold 0.5]
"""
import argparse
import csv
import os
import re
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from _backfill_io import setup_logging, print_summary
from _sources_registry import SOURCES_DIR
from sqlalchemy import text
from app.database import SessionLocal

_HDR = re.compile(r"Gene=(?P<gene>\S+).*?BestRef=(?P<ref>\S+).*?Identity=(?P<id>[\d.]+|no_hit)")


def load_af3_targets(path):
    """md5 -> (seq, gene, best_ref, identity_pct_or_None)."""
    out = {}
    h = meta = None
    s = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h is not None:
                out[h] = ("".join(s),) + meta
            body = line[1:]
            h = body.split()[0]
            m = _HDR.search(body)
            if m:
                ident = None if m.group("id") == "no_hit" else float(m.group("id"))
                meta = (m.group("gene"), m.group("ref"), ident)
            else:
                meta = ("?", "?", None)
            s = []
        else:
            s.append(line)
    if h is not None:
        out[h] = ("".join(s),) + meta
    return out


def kset(seq, k=5):
    return set(seq[i:i + k] for i in range(len(seq) - k + 1)) if len(seq) >= k else set()


def containment(a_k, b_k):
    if not a_k or not b_k:
        return 0.0
    return len(a_k & b_k) / min(len(a_k), len(b_k))


def main():
    ap = argparse.ArgumentParser(description="ccRCC folded-ORF QC flag")
    ap.add_argument("--threshold", type=float, default=0.50,
                    help="fold_id below this flags a mismatch (default 0.50)")
    args = ap.parse_args()

    logfile = setup_logging("backfill_ccrcc_orf_flag")
    af3 = load_af3_targets(str(SOURCES_DIR / "ccrcc" / "af3_targets.fa"))
    print(f"[ccrcc_orf_flag] loaded {len(af3):,} folded proteins from af3_targets.fa")

    db = SessionLocal()
    rows = db.execute(text(
        "SELECT i.id, substring(st.af3_cif_url from '([0-9a-f]{32})') AS md5, "
        "       g.symbol, g.canonical_sequence "
        "FROM isoforms i "
        "JOIN structures st ON st.isoform_id = i.id "
        "JOIN genes g ON g.id = i.gene_id "
        "WHERE i.source_id = 2"
    )).fetchall()

    canon_kcache = {}
    flagged = []
    stats = {"total": 0, "mismatch": 0, "trusted": 0, "no_ref": 0, "no_fold": 0, "updated": 0}

    for iso_id, md5, symbol, canon in rows:
        stats["total"] += 1
        entry = af3.get(md5)
        if not entry:
            stats["no_fold"] += 1
            continue
        folded, _agene, best_ref, best_id = entry
        if best_ref in ("?", "none", "", None):
            best_ref = None

        if not canon:
            stats["no_ref"] += 1
            mismatch = False
            fold_id = None
        else:
            ck = canon_kcache.get(symbol)
            if ck is None:
                ck = kset(canon); canon_kcache[symbol] = ck
            fold_id = containment(kset(folded), ck)
            mismatch = fold_id < args.threshold
            stats["mismatch" if mismatch else "trusted"] += 1

        db.execute(text(
            "UPDATE transcript_metadata SET folded_orf_mismatch=:m, "
            "folded_best_ref=:r, folded_best_ref_identity=:i WHERE isoform_id=:iid"
        ), {"m": mismatch, "r": best_ref, "i": best_id, "iid": iso_id})
        stats["updated"] += 1

        if mismatch:
            flagged.append((iso_id, symbol, best_ref or "", best_id or "",
                            f"{fold_id:.3f}", len(folded)))

    db.commit()
    db.close()

    # Dump the flagged list for author review.
    out_tsv = os.path.join(os.path.dirname(logfile), "ccrcc_folded_orf_flagged.tsv")
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["isoform_id", "gene_symbol", "folded_best_ref", "folded_best_ref_identity",
                    "fold_id_to_gene", "folded_len"])
        w.writerows(sorted(flagged, key=lambda r: r[4]))
    print(f"[ccrcc_orf_flag] flagged list → {out_tsv}")

    print_summary(
        "ccRCC FOLDED-ORF QC FLAG COMPLETE",
        counts=[
            ("ccRCC structured isoforms", stats["total"]),
            ("rows updated",              stats["updated"]),
            ("✓ trusted (fold matches gene)", stats["trusted"]),
            ("⚑ mismatch (flagged)",      stats["mismatch"]),
            ("– no canonical ref",        stats["no_ref"]),
            ("– folded protein missing",  stats["no_fold"]),
        ],
        total=stats["total"], logfile=logfile,
    )


if __name__ == "__main__":
    main()
