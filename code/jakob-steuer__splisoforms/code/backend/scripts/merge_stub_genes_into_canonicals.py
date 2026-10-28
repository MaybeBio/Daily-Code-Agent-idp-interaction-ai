#!/usr/bin/env python3
"""
merge_stub_genes_into_canonicals.py
====================================
Reconcile the ~7,600 stub Gene rows that the ccRCC/BC ingest created with
their canonical GENCODE counterpart by way of shared UniProt accession.

Why this exists
---------------
Phase 6a (`backfill_gene_uniprot.py`) resolved Ensembl gene IDs → UniProt
accessions for every GENCODE canonical gene. In ~7,679 cases the resolved
UniProt was already on a *stub* Gene row carved out by Phase 2 / Phase 3
ingest because the long-read pipeline saw a protein with no canonical
gene record at ingest time. The unique constraint on `genes.uniprot_id`
forced us to skip those assignments — the canonical gene was left
without its UniProt, and the stub kept a UniProt that biologically
belongs to the canonical.

This script closes the loop: for each canonical-without-UniProt that
resolves to an accession already held by a stub, we merge the stub
into the canonical:

    1. Reassign every Isoform.gene_id from stub → canonical.
       (Source memberships and structures FK on Isoform.id and follow
       automatically.)
    2. Clear the stub's uniprot_id (frees the unique constraint).
    3. Set the canonical's uniprot_id to the resolved accession.
    4. DELETE the now-orphan stub Gene row.

Usage
-----
    .venv/bin/python backend/scripts/merge_stub_genes_into_canonicals.py --dry-run
    .venv/bin/python backend/scripts/merge_stub_genes_into_canonicals.py --apply
    .venv/bin/python backend/scripts/merge_stub_genes_into_canonicals.py --apply --limit 50
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.database import SessionLocal
from app.models.bio import Gene, Isoform, CellTypeExpression
from sqlalchemy import func

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("merge_stub_genes")

API_BASE = "https://rest.uniprot.org/idmapping"
POLL_INTERVAL_SEC = 3
POLL_MAX_SEC = 600
BATCH_SIZE = 1000
RESULT_PAGE_SIZE = 500


def _strip_version(ensg: str) -> str:
    return ensg.split(".", 1)[0] if "." in ensg else ensg


def _resolve_ensg_to_swissprot(ids: list[str]) -> dict[str, str]:
    """Batch ENSG → reviewed SwissProt accession (one per gene): {ensg: acc}.

    Uses the current UniProt idmapping flow (to=UniProtKB-Swiss-Prot, plain
    /idmapping/results/), matching backfill_gene_uniprot.py. The older
    to=UniProtKB + /idmapping/uniprotkb/results/ path the recovered script used
    now returns HTTP 400.
    """
    out: dict[str, str] = {}
    for i in range(0, len(ids), BATCH_SIZE):
        chunk = ids[i:i + BATCH_SIZE]
        r = requests.post(f"{API_BASE}/run",
                          data={"from": "Ensembl", "to": "UniProtKB-Swiss-Prot",
                                "ids": ",".join(chunk)}, timeout=60)
        r.raise_for_status()
        job = r.json()["jobId"]
        t0 = time.time()
        while True:
            s = requests.get(f"{API_BASE}/status/{job}", timeout=60).json()
            if s.get("jobStatus") in ("RUNNING", "NEW"):
                if time.time() - t0 > POLL_MAX_SEC:
                    raise TimeoutError(f"job {job} stuck")
                time.sleep(POLL_INTERVAL_SEC)
                continue
            break
        url = (f"{API_BASE}/results/{job}"
               f"?format=tsv&fields=accession&size={RESULT_PAGE_SIZE}")
        while url:
            resp = requests.get(url, timeout=120)
            resp.raise_for_status()
            for line in resp.text.splitlines()[1:]:       # skip "From\tTo" header
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0] not in out:  # first hit wins
                    out[parts[0]] = parts[1]
            url = None
            m = re.search(r'<([^>]+)>\s*;\s*rel="next"', resp.headers.get("Link", ""))
            if m:
                url = m.group(1)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--dry-run", action="store_true",
                   help="Resolve mappings + identify merges; no DB writes")
    g.add_argument("--apply", action="store_true",
                   help="Perform the merges")
    parser.add_argument("--limit", type=int,
                        help="Cap canonical-genes-to-resolve (testing)")
    args = parser.parse_args()

    db = SessionLocal()

    # Eligible canonical genes: have ensembl_gene_id but NO uniprot_id
    # (i.e. the cases where Phase 6a couldn't assign because of a stub
    # collision OR no UniProt mapping was found).
    q = (db.query(Gene)
           .filter(Gene.ensembl_gene_id.isnot(None))
           .filter(Gene.uniprot_id.is_(None)))
    if args.limit:
        q = q.limit(args.limit)
    canonicals = q.all()
    log.info("canonical genes lacking uniprot_id: %d", len(canonicals))

    if not canonicals:
        log.info("nothing to do.")
        db.close()
        return 0

    # Resolve their UniProt accessions in batches
    ensg_list = [_strip_version(g.ensembl_gene_id) for g in canonicals]
    ensg_to_canonical = {_strip_version(g.ensembl_gene_id): g for g in canonicals}

    log.info("resolving via UniProt idmapping (Ensembl→Swiss-Prot) ...")
    ensg_to_acc = _resolve_ensg_to_swissprot(ensg_list)
    log.info("resolved %d ENSG → UniProt mappings", len(ensg_to_acc))

    # Index existing uniprot_id → Gene (these are the stub holders)
    log.info("indexing existing uniprot_id holders ...")
    acc_to_holder: dict[str, Gene] = {}
    for g in db.query(Gene).filter(Gene.uniprot_id.isnot(None)).all():
        acc_to_holder[g.uniprot_id] = g

    # Plan merges
    merges: list[tuple[Gene, Gene, str]] = []  # (canonical, stub, acc)
    no_resolution = 0
    no_holder = 0
    same_gene = 0
    shared_accession = 0
    # Each UniProt accession can be held by exactly one Gene (UNIQUE). When
    # several distinct canonical genes resolve to the SAME accession they are
    # paralogs / readthroughs sharing one SwissProt entry (e.g. AMY1A/AMY1C =
    # P0DTE7, FPGT-TNNI3K/TNNI3K) — NOT duplicates. Merge only the first; skip
    # the rest cleanly rather than letting them collide on the UNIQUE constraint.
    claimed_accs: set[str] = set()
    for ensg, canonical in ensg_to_canonical.items():
        acc = ensg_to_acc.get(ensg)
        if not acc:
            no_resolution += 1
            continue
        holder = acc_to_holder.get(acc)
        if not holder:
            # No stub holds it — Phase 6a's eligibility filter excluded
            # this canonical for some other reason. Skip; the trivial
            # case is to just assign the accession.
            no_holder += 1
            continue
        if holder.id == canonical.id:
            same_gene += 1
            continue
        if acc in claimed_accs:
            shared_accession += 1   # paralog/readthrough sharing one accession
            continue
        claimed_accs.add(acc)
        merges.append((canonical, holder, acc))

    log.info("planned merges:               %d", len(merges))
    log.info("  canonical missing resolution:  %d", no_resolution)
    log.info("  no holder (trivial assignment): %d  (could just be set)", no_holder)
    log.info("  same gene (already canonical):  %d", same_gene)
    log.info("  shared accession (paralog, skipped): %d", shared_accession)

    if args.dry_run:
        # Sample summary of what would happen
        for canonical, stub, acc in merges[:5]:
            n_iso = db.query(Isoform).filter(Isoform.gene_id == stub.id).count()
            log.info("  preview: stub id=%d symbol=%r (%d isoforms) → canonical id=%d symbol=%r [acc=%s]",
                     stub.id, stub.symbol, n_iso, canonical.id, canonical.symbol, acc)
        log.info("[dry-run] no changes written.")
        db.close()
        return 0

    # --- Apply merges ---
    log.info("applying %d merges ...", len(merges))
    n_merged = 0
    n_iso_moved = 0
    n_errors = 0
    for canonical, stub, acc in merges:
        try:
            # Reassign isoforms and update delta_length
            if canonical.canonical_length is not None:
                iso_count = (db.query(Isoform)
                             .filter(Isoform.gene_id == stub.id)
                             .update({"gene_id": canonical.id,
                                      "delta_length": func.length(Isoform.sequence) - canonical.canonical_length},
                                     synchronize_session=False))
            else:
                iso_count = (db.query(Isoform)
                             .filter(Isoform.gene_id == stub.id)
                             .update({"gene_id": canonical.id},
                                     synchronize_session=False))
            n_iso_moved += iso_count

            # If the stub was the target of any canonical_transcript_id
            # FK, clear it (the canonical gene already has its own).
            if stub.canonical_transcript_id is not None:
                stub.canonical_transcript_id = None

            # Drop the stub's gene-level HPA cell-type rows before deletion.
            # CellTypeExpression.gene_id is NOT NULL and the ORM would otherwise
            # try to NULL it on delete (relationship has no delete-orphan cascade),
            # raising NotNullViolation. These rows are redundant — the canonical
            # gene carries its own HPA data keyed by its ensembl_gene_id.
            db.query(CellTypeExpression).filter(
                CellTypeExpression.gene_id == stub.id
            ).delete(synchronize_session=False)

            # Clear stub's uniprot_id to free the unique constraint, then
            # assign to canonical.
            stub_acc = stub.uniprot_id
            stub.uniprot_id = None
            db.flush()
            canonical.uniprot_id = stub_acc

            # Delete the stub
            db.delete(stub)
            db.commit()
            n_merged += 1
            if n_merged % 100 == 0:
                log.info("  ... merged %d (isoforms moved: %d)",
                         n_merged, n_iso_moved)
        except Exception as e:
            db.rollback()
            n_errors += 1
            if n_errors <= 10:
                log.error("  merge failed for stub id=%d → canonical id=%d: %r",
                          stub.id, canonical.id, e)

    log.info("✅ summary:")
    log.info("  merges completed:  %d", n_merged)
    log.info("  isoforms moved:    %d", n_iso_moved)
    log.info("  errors:            %d", n_errors)
    db.close()
    return 0 if n_errors == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
