"""
Recompute ONLY `ptm_accessibility.structural_context` from current TED domains.

Why this exists
---------------
`structural_context` (intra-domain / disrupted / inter-domain) is the sole
PtmAccessibility column derived from TED domains. After ingesting Chainsaw
domains for the GENCODE-novel + BC isoforms (generate_novel_ted.py), the
context for every isoform that already had PTM rows is stale — it was computed
when those isoforms had no domains, so every site fell to "inter-domain".

The PTM sites themselves (positions, types, SASA, pLDDT) come from UniProt /
structures and have NOT changed, so a full `backfill_ptm.py --force` re-crawl
of UniProt is wasted work (~20h). Everything needed to recompute the context is
already stored on each row:
    residue_position         = canonical position (can_pos)
    novel_residue_position   = mapped novel position (nov_pos; NULL = site lost)
    canonical_isoform_id     = the canonical used (NULL = no canonical / identity)

This reuses the AUTHORITATIVE domain-membership functions from
app.analysis.ptm_accessibility (is_in_ted_domain / get_ted_status) so the result
is identical to what analyze_isoform() would produce, with zero network I/O.

Isoforms that have NO ptm rows at all are NOT touched here (they need the real
backfill_ptm.py run to create rows).

Usage:
    python backend/scripts/recompute_ptm_structural_context.py [--apply]
                                     [--isoforms 1,2,3] [--report-only]
Default is a dry run (counts changes, writes nothing). Pass --apply to persist.
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import text
from app.database import SessionLocal
from app.analysis.ptm_accessibility import get_ted_status  # authoritative logic
from app.models.bio import TedDomain


def _context_for(can_pos, nov_pos, canonical_id, can_domains, nov_domains):
    """Mirror of ptm_accessibility.py lines 264-272, fed from stored values."""
    can_in_domain = get_ted_status(can_pos, can_domains) if canonical_id else False
    nov_in_domain = get_ted_status(nov_pos, nov_domains) if (nov_pos is not None) else False
    if nov_in_domain:
        return "intra-domain"
    if can_in_domain and not nov_in_domain:
        return "disrupted"
    return "inter-domain"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="Persist updates (default: dry run)")
    ap.add_argument("--isoforms", type=str, default=None,
                    help="Comma-separated isoform_ids to restrict to (validation)")
    ap.add_argument("--report-only", action="store_true",
                    help="Only print the resulting context distribution, no per-row work")
    args = ap.parse_args()

    db = SessionLocal()

    # 1. Load ALL ted domains into memory: isoform_id -> [TedDomain, ...]
    #    (get_ted_status only reads .chopping / .seq_start / .seq_end)
    dom_by_iso = defaultdict(list)
    for d in db.query(TedDomain.isoform_id, TedDomain.chopping,
                      TedDomain.seq_start, TedDomain.seq_end).all():
        # lightweight shim objects with the attributes get_ted_status needs
        class _D:  # noqa: E306
            pass
        o = _D()
        o.isoform_id, o.chopping, o.seq_start, o.seq_end = d
        dom_by_iso[d[0]].append(o)
    print(f"[*] Loaded domains for {len(dom_by_iso)} isoforms.")

    # 2. Stream ptm rows.
    where = ""
    params = {}
    if args.isoforms:
        ids = [int(x) for x in args.isoforms.split(",") if x.strip()]
        where = "WHERE isoform_id = ANY(:ids)"
        params["ids"] = ids

    rows = db.execute(text(
        f"SELECT id, isoform_id, residue_position, novel_residue_position, "
        f"canonical_isoform_id, structural_context FROM ptm_accessibility {where}"
    ), params).fetchall()
    print(f"[*] Scanning {len(rows):,} ptm rows...")

    updates = []          # (id, new_context)
    new_dist = defaultdict(int)
    changed_dist = defaultdict(int)  # old -> new transitions
    for rid, iso_id, can_pos, nov_pos, can_iso, old_ctx in rows:
        new_ctx = _context_for(
            can_pos, nov_pos, can_iso,
            dom_by_iso.get(can_iso, []),
            dom_by_iso.get(iso_id, []),
        )
        new_dist[new_ctx] += 1
        if new_ctx != old_ctx:
            changed_dist[f"{old_ctx} -> {new_ctx}"] += 1
            updates.append({"id": rid, "ctx": new_ctx})

    print("\n=== resulting structural_context distribution ===")
    for k in sorted(new_dist):
        print(f"  {k:14s} {new_dist[k]:>10,}")
    print(f"\n=== rows that would change: {len(updates):,} ===")
    for k in sorted(changed_dist, key=lambda x: -changed_dist[x]):
        print(f"  {k:28s} {changed_dist[k]:>10,}")

    if args.report_only:
        db.close(); return

    if not args.apply:
        print("\n[dry run] pass --apply to persist.")
        db.close(); return

    print(f"\n[*] Applying {len(updates):,} updates...")
    B = 5000
    for i in range(0, len(updates), B):
        batch = updates[i:i + B]
        db.execute(
            text("UPDATE ptm_accessibility SET structural_context = :ctx WHERE id = :id"),
            batch,
        )
        db.commit()
        print(f"  committed {min(i + B, len(updates)):,}/{len(updates):,}")
    print("[*] Done.")
    db.close()


if __name__ == "__main__":
    main()
