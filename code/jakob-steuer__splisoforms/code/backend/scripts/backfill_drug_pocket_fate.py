#!/usr/bin/env python3
"""
backfill_drug_pocket_fate.py — our own Axis-1: carry each drug complex's contact
pocket (defined on the co-folded/drugged transcript, backfill_drug_contacts.py)
onto every sibling isoform of the same gene, and score residue retention.

Anchoring is on the DRUGGED TRANSCRIPT, not a gene-level "canonical", so this is
unaffected by the is_canonical-vs-MANE ambiguity that complicates gene-level anchoring.

Per (drug_complex × sibling isoform):
  - align drugged-transcript sequence → isoform sequence (global, BLOSUM62,
    gap −10/−0.5); the alignment is cached per (transcript, isoform) pair and
    reused across all drugs on that transcript.
  - retained = contact residues that map to a non-gap (aligned) position, i.e.
    NOT deleted by splicing.
  - frac_retained → fate: lost <0.25, disrupted <0.75, preserved ≥0.75.
  - DECIDABILITY GATE: if the pocket is partially retained but the isoform's own
    per-residue pLDDT over the retained pocket positions averages < 70, emit
    `no_call` (a reported category) instead of trusting a low-confidence region.
    A fully deleted pocket (0 retained) is a confident `lost` — deletion is a
    sequence fact, not pLDDT-dependent.

Usage:
    python backend/scripts/backfill_drug_pocket_fate.py [--workers N] [--gene SYM]
                                                        [--limit N] [--force]
"""
import os
import sys
import json
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Bio import Align
from Bio.Align import substitution_matrices
from app.database import SessionLocal
from app.models.bio import Drug, DrugComplex, DrugPocketFate, Isoform, Gene, StructuralMetrics

PLDDT_GATE = 70.0
LOST, DISRUPTED, PRESERVED = 0.25, 0.75, 0.75

_BLOSUM = substitution_matrices.load("BLOSUM62")
_ALPHABET = set(_BLOSUM.alphabet)


def _clean(seq: str) -> str:
    seq = (seq or "").rstrip("*").upper()
    return "".join(c if c in _ALPHABET else "X" for c in seq)


def align_blocks(t_seq: str, u_seq: str):
    """Return alignment.aligned as ([t_blocks],[u_blocks]) or None if empty."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = _BLOSUM
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"
    try:
        aln = next(iter(aligner.align(t_seq, u_seq)))
    except (StopIteration, Exception):  # noqa: BLE001
        return None
    t_blocks, u_blocks = aln.aligned
    return ([[int(a), int(b)] for a, b in t_blocks],
            [[int(a), int(b)] for a, b in u_blocks])


def _worker(args):
    tid, uid, t_seq, u_seq = args
    return tid, uid, align_blocks(t_seq, u_seq)


def build_tidx_to_uidx(blocks):
    """Map 0-based target index -> 0-based query index over aligned blocks."""
    m = {}
    if not blocks:
        return m
    t_blocks, u_blocks = blocks
    for (ts, te), (us, _ue) in zip(t_blocks, u_blocks):
        for k in range(te - ts):
            m[ts + k] = us + k
    return m


def classify(frac):
    if frac < LOST:
        return "lost", "domain_deletion"
    if frac < PRESERVED:
        return "disrupted", "pocket_disruption"
    return "preserved", "preserved"


def main():
    ap = argparse.ArgumentParser(description="Carry drug pockets onto sibling isoforms")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--gene", type=str, help="limit to one gene symbol (debug)")
    ap.add_argument("--limit", type=int, help="cap drugged transcripts processed")
    ap.add_argument("--force", action="store_true", help="delete + recompute")
    args = ap.parse_args()

    db = SessionLocal()

    # complexes with a pocket + a drugged transcript
    cq = db.query(DrugComplex.id, DrugComplex.target_isoform_id, DrugComplex.gene_id,
                  DrugComplex.contact_residues_json).filter(
        DrugComplex.target_isoform_id.isnot(None),
        DrugComplex.n_contact_residues > 0,
    )
    if args.gene:
        gid = db.query(Gene.id).filter(Gene.symbol == args.gene).scalar()
        if not gid:
            sys.exit(f"gene {args.gene} not found")
        cq = cq.filter(DrugComplex.gene_id == gid)
    complexes = cq.all()

    # group complexes by drugged transcript
    by_tx = defaultdict(list)
    gene_of_tx = {}
    for cid, tid, gid, contacts in complexes:
        by_tx[tid].append((cid, json.loads(contacts)))
        gene_of_tx[tid] = gid
    tx_ids = list(by_tx)
    if args.limit:
        tx_ids = tx_ids[:args.limit]
    print(f"{len(complexes)} complexes across {len(tx_ids)} drugged transcripts", flush=True)

    # sibling isoforms per gene; sequences + plddt series for everything involved
    gene_ids = {gene_of_tx[t] for t in tx_ids}
    siblings = defaultdict(list)  # gene_id -> [isoform_id]
    seqs, plddt = {}, {}
    for iso_id, gid, seq in db.query(Isoform.id, Isoform.gene_id, Isoform.sequence).filter(
            Isoform.gene_id.in_(gene_ids)).all():
        siblings[gid].append(iso_id)
        seqs[iso_id] = _clean(seq)
    for iso_id, series in db.query(StructuralMetrics.isoform_id,
                                   StructuralMetrics.plddt_series).filter(
            StructuralMetrics.isoform_id.in_(list(seqs))).all():
        plddt[iso_id] = series

    if args.force:
        cids = [c for t in tx_ids for c, _ in by_tx[t]]
        for i in range(0, len(cids), 5000):
            db.query(DrugPocketFate).filter(
                DrugPocketFate.drug_complex_id.in_(cids[i:i+5000])).delete(
                synchronize_session=False)
        db.commit()

    # unique (tx, sibling) alignment tasks
    tasks = []
    for tid in tx_ids:
        gid = gene_of_tx[tid]
        if tid not in seqs:
            continue
        for uid in siblings[gid]:
            if uid == tid or uid not in seqs:
                continue
            tasks.append((tid, uid, seqs[tid], seqs[uid]))
    print(f"aligning {len(tasks)} (transcript, sibling) pairs on {args.workers} workers…",
          flush=True)

    amap = {}  # (tid,uid) -> tidx_to_uidx
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(_worker, t) for t in tasks]
        for fut in as_completed(futs):
            tid, uid, blocks = fut.result()
            amap[(tid, uid)] = build_tidx_to_uidx(blocks)
            done += 1
            if done % 3000 == 0:
                print(f"  aligned {done}/{len(tasks)}", flush=True)

    # compute fate rows
    print("scoring fate rows…", flush=True)
    rows = []
    n = 0
    for tid in tx_ids:
        gid = gene_of_tx[tid]
        for uid in siblings[gid]:
            if uid == tid or (tid, uid) not in amap:
                continue
            t2u = amap[(tid, uid)]
            u_plddt = plddt.get(uid)
            for cid, contacts in by_tx[tid]:
                n_pocket = len(contacts)
                if n_pocket == 0:
                    continue
                mapped = [t2u[p - 1] for p in contacts if (p - 1) in t2u]
                n_ret = len(mapped)
                frac = n_ret / n_pocket
                # decidability gate on retained-pocket pLDDT in the isoform
                decidable, fate, mech = True, None, None
                if n_ret == 0:
                    fate, mech = "lost", "domain_deletion"
                else:
                    pv = [u_plddt[i] for i in mapped if u_plddt and i < len(u_plddt)]
                    pocket_plddt = (sum(pv) / len(pv)) if pv else None
                    if pocket_plddt is not None and pocket_plddt < PLDDT_GATE:
                        decidable, fate, mech = False, "no_call", None
                    else:
                        fate, mech = classify(frac)
                pv = [u_plddt[i] for i in mapped if u_plddt and i < len(u_plddt)] if n_ret else []
                rows.append({
                    "drug_complex_id": cid, "isoform_id": uid,
                    "n_pocket_residues": n_pocket, "n_retained": n_ret,
                    "frac_retained": round(frac, 4), "fate": fate,
                    "mechanism_class": mech, "decidable": decidable,
                    "pocket_plddt": round(sum(pv) / len(pv), 2) if pv else None,
                    "retained_residues_json": json.dumps([m + 1 for m in mapped]),
                })
                n += 1
                if len(rows) >= 5000:
                    db.bulk_insert_mappings(DrugPocketFate, rows); db.commit(); rows = []
    if rows:
        db.bulk_insert_mappings(DrugPocketFate, rows); db.commit()

    total = db.query(DrugPocketFate).count()
    print(f"done: wrote {n} fate rows; {total} total in drug_pocket_fate")
    db.close()


if __name__ == "__main__":
    main()
