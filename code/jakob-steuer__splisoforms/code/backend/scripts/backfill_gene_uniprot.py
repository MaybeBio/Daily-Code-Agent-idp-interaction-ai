#!/usr/bin/env python3
"""
backfill_gene_uniprot.py — map GENCODE genes to a reviewed UniProt accession.

The unified ingest only sets Gene.uniprot_id from long-read BLAST best-hits, so
GENCODE-reference genes (which carry an ensembl_gene_id but no best-hit) are left
unmapped — limiting homodimer / STRING / cell-type / descriptions coverage.

This fills the gap via UniProt's batch ID-mapping API:
    from = "Ensembl"  (ENSG, version-stripped)
    to   = "UniProtKB-Swiss-Prot"   (reviewed canonical only — one accession/gene)

Only genes with `ensembl_gene_id` set and `uniprot_id` NULL are touched (gap-fill,
never overwrite an existing mapping). Gene.uniprot_id is UNIQUE, so an accession
already used by another gene is skipped.

Usage:
    python backend/scripts/backfill_gene_uniprot.py
    python backend/scripts/backfill_gene_uniprot.py --limit 200   # smoke test
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import requests

from app.database import SessionLocal
from app.models.bio import Gene
from _backfill_io import Progress, print_summary, setup_logging

log = logging.getLogger("backfill_gene_uniprot")

UNIPROT = "https://rest.uniprot.org"
CHUNK = 5000           # ENSGs per idmapping job
POLL_INTERVAL = 3
POLL_MAX = 200


def _strip_version(ensg: str) -> str:
    return ensg.split(".")[0]


def _run_mapping(ids: List[str]) -> Dict[str, str]:
    """Submit one idmapping job (Ensembl→Swiss-Prot), return {ENSG: accession}."""
    r = requests.post(f"{UNIPROT}/idmapping/run",
                      data={"from": "Ensembl", "to": "UniProtKB-Swiss-Prot",
                            "ids": ",".join(ids)}, timeout=60)
    r.raise_for_status()
    job = r.json()["jobId"]

    for _ in range(POLL_MAX):
        s = requests.get(f"{UNIPROT}/idmapping/status/{job}", timeout=60).json()
        status = s.get("jobStatus")
        if status in ("RUNNING", "NEW"):
            time.sleep(POLL_INTERVAL)
            continue
        break

    out: Dict[str, str] = {}
    url = (f"{UNIPROT}/idmapping/results/{job}"
           f"?format=tsv&fields=accession&size=500")
    while url:
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        for line in resp.text.splitlines()[1:]:        # skip header "From\tTo"
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0] not in out:  # first reviewed hit wins
                out[parts[0]] = parts[1]
        # follow pagination via Link: <url>; rel="next"
        url = None
        link = resp.headers.get("Link", "")
        if 'rel="next"' in link:
            url = link[link.index("<") + 1:link.index(">")]
    return out


def run(limit: int | None) -> None:
    setup_logging("backfill_gene_uniprot")
    db = SessionLocal()
    try:
        q = (db.query(Gene)
             .filter(Gene.uniprot_id.is_(None), Gene.ensembl_gene_id.isnot(None)))
        if limit:
            q = q.limit(limit)
        genes = q.all()
        log.info("genes to map (ensembl set, uniprot null): %d", len(genes))

        # ENSG (version-stripped) → Gene ; accessions already in use (unique col)
        ensg_to_gene: Dict[str, Gene] = {}
        for g in genes:
            ensg_to_gene.setdefault(_strip_version(g.ensembl_gene_id), g)
        used = set(x[0] for x in
                   db.query(Gene.uniprot_id).filter(Gene.uniprot_id.isnot(None)).all())

        ensgs = list(ensg_to_gene.keys())
        stats = {"mapped": 0, "no_swissprot": 0, "accession_taken": 0}

        with Progress(total=len(ensgs), desc="gene→uniprot", unit="gene") as p:
            for i in range(0, len(ensgs), CHUNK):
                chunk = ensgs[i:i + CHUNK]
                try:
                    result = _run_mapping(chunk)
                except Exception as e:  # noqa: BLE001
                    log.warning("idmapping chunk failed (%d ids): %s", len(chunk), e)
                    p.update(len(chunk))
                    continue
                for ensg in chunk:
                    acc = result.get(ensg)
                    g = ensg_to_gene[ensg]
                    if not acc:
                        stats["no_swissprot"] += 1
                    elif acc in used:
                        stats["accession_taken"] += 1
                    else:
                        g.uniprot_id = acc
                        used.add(acc)
                        stats["mapped"] += 1
                db.commit()
                p.update(len(chunk), status=f"mapped={stats['mapped']}")

        print_summary("backfill_gene_uniprot",
                      [("genes considered", len(genes)),
                       ("newly mapped", stats["mapped"]),
                       ("no reviewed entry", stats["no_swissprot"]),
                       ("accession already used", stats["accession_taken"])])
    finally:
        db.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true",
                    help="(reserved) overwrite existing uniprot_id — not yet enabled")
    args = ap.parse_args()
    run(args.limit)


if __name__ == "__main__":
    main()
