#!/usr/bin/env python3
"""
backfill_canonical_designation.py — enforce exactly ONE canonical isoform per gene.

The source metadata's `is_canonical` column is not one-per-gene (GENCODE marks
several transcripts principal; some novel long-read rows arrive flagged too), so
the DB ended up with 5,852 genes carrying >1 `is_canonical` isoform and 11,871
novel PB.* rows flagged canonical. This pass re-designates the single canonical
per gene and clears the flag on every other isoform.

Selection rule (agreed 2026-07):
  1. MANE Select isoform, if the gene has one.
  2. else, among the gene's REFERENCE (ENST) isoforms: the GENCODE-canonical one
     (is_canonical flag); tie-break / fallback by longest protein.
  3. else (novel-only locus, no ENST): the longest isoform.
Ties broken by longest protein, then lowest id (deterministic).

Also sets `gene.canonical_transcript_id` to the chosen isoform for consistency.
Does NOT touch `is_mane_select` (source truth) or `gene.canonical_sequence`
(UniProt-derived reference, a separate concept).

Usage:
    python backend/scripts/backfill_canonical_designation.py            # dry-run
    python backend/scripts/backfill_canonical_designation.py --apply
    python backend/scripts/backfill_canonical_designation.py --gene UGT2A1
"""
import os
import sys
import argparse
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import func
from app.database import SessionLocal
from app.models.bio import Gene, Isoform


class Iso:
    __slots__ = ("id", "gene_id", "mane", "canon", "is_ref", "seqlen")

    def __init__(self, id, gene_id, mane, canon, is_ref, seqlen):
        self.id, self.gene_id = id, gene_id
        self.mane, self.canon, self.is_ref = bool(mane), bool(canon), is_ref
        self.seqlen = seqlen or 0


def _best(pool):
    return max(pool, key=lambda i: (i.seqlen, -i.id))


def choose(isos):
    mane = [i for i in isos if i.mane]
    if mane:
        return _best(mane)
    refs = [i for i in isos if i.is_ref]
    if refs:
        flagged = [i for i in refs if i.canon]
        return _best(flagged or refs)
    return _best(isos) if isos else None


def main():
    ap = argparse.ArgumentParser(description="Enforce one canonical isoform per gene")
    ap.add_argument("--apply", action="store_true", help="write changes (default: dry-run)")
    ap.add_argument("--gene", type=str, help="limit to one gene symbol (debug)")
    args = ap.parse_args()

    db = SessionLocal()

    q = db.query(
        Isoform.id, Isoform.gene_id, Isoform.is_mane_select, Isoform.is_canonical,
        Isoform.ensembl_transcript_id.isnot(None), func.length(Isoform.sequence),
    )
    if args.gene:
        gid = db.query(Gene.id).filter(Gene.symbol == args.gene).scalar()
        if not gid:
            sys.exit(f"gene {args.gene} not found")
        q = q.filter(Isoform.gene_id == gid)

    by_gene: dict[int, list[Iso]] = defaultdict(list)
    for row in q.all():
        by_gene[row[1]].append(Iso(*row))

    iso_updates, gene_updates = [], []
    flipped_on = flipped_off = genes_changed = 0

    for gene_id, isos in by_gene.items():
        chosen = choose(isos)
        chosen_id = chosen.id if chosen else None
        gene_changed = False
        for i in isos:
            should = (i.id == chosen_id)
            if should != i.canon:
                iso_updates.append({"id": i.id, "is_canonical": should})
                if should:
                    flipped_on += 1
                else:
                    flipped_off += 1
                gene_changed = True
        gene_updates.append({"id": gene_id, "canonical_transcript_id": chosen_id})
        if gene_changed:
            genes_changed += 1

    print(f"genes examined            : {len(by_gene)}")
    print(f"genes with is_canonical fix: {genes_changed}")
    print(f"  flags set   (→canonical) : {flipped_on}")
    print(f"  flags clear (→ not)      : {flipped_off}")

    if not args.apply:
        # show a few multi-canonical fixes for sanity
        sample = [gid for gid, isos in by_gene.items()
                  if sum(1 for i in isos if i.canon) > 1][:5]
        if sample:
            print("\nsample genes currently carrying >1 canonical (will be fixed):")
            for gid in sample:
                g = db.query(Gene.symbol).filter(Gene.id == gid).scalar()
                ch = choose(by_gene[gid])
                print(f"  {g}: {sum(1 for i in by_gene[gid] if i.canon)} canonical "
                      f"→ keep isoform id {ch.id}")
        print("\n(dry-run — re-run with --apply to write)")
        db.close()
        return

    print("\napplying …", flush=True)
    for i in range(0, len(iso_updates), 5000):
        db.bulk_update_mappings(Isoform, iso_updates[i:i+5000]); db.commit()
    for i in range(0, len(gene_updates), 5000):
        db.bulk_update_mappings(Gene, gene_updates[i:i+5000]); db.commit()

    # verify
    dupes = db.query(func.count()).select_from(
        db.query(Isoform.gene_id).filter(Isoform.is_canonical.is_(True))
        .group_by(Isoform.gene_id).having(func.count() > 1).subquery()).scalar()
    print(f"done. genes with >1 canonical remaining: {dupes}")
    db.close()


if __name__ == "__main__":
    main()
