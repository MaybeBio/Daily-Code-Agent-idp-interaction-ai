#!/usr/bin/env python3
"""
QC: compare the designated canonical isoform against the UniProt reviewed
(SwissProt) canonical sequence.

RETIRED FROM THE BACKFILL CHAIN (2026-10, review ING-3). `Gene.canonical_sequence`
is now written by `finalize_reattribution.py` from `gene.canonical_transcript_id`
(the designated canonical isoform: MANE → APPRIS/GENCODE-canonical → longest).
That is the baseline Pfam/NMD/neoantigen compare against, and it must be a
sequence that exists as an isoform in the DB. A UniProt canonical frequently does
not (older sequence version, a TrEMBL entry, or — for 119 genes — it was being
overwritten by a novel long-read sequence), so this script no longer writes it.

Default mode is REPORT ONLY: for each gene it fetches the UniProt canonical and
logs whether it is identical to, longer/shorter than, or absent for the
designated canonical. `--write-missing` fills `canonical_sequence` only for genes
that have none (no canonical isoform at all); nothing is ever overwritten.

Usage:
    python backend/scripts/backfill_canonical_uniprot.py                    # QC report, all genes
    python backend/scripts/backfill_canonical_uniprot.py --genes ENO1,VHL   # specific genes
    python backend/scripts/backfill_canonical_uniprot.py --limit 100        # test run
    python backend/scripts/backfill_canonical_uniprot.py --write-missing    # fill NULLs only
"""
import os
import sys
import time
import argparse
import logging
import requests
from typing import Optional, Tuple

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal
from app.models.bio import Gene

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

# SwissProt accession pattern: single letter + digit + 3 alphanum + digit (6 chars)
# or the newer 10-char format: single letter + digit + 8 alphanum + digit
import re
SWISSPROT_PATTERN = re.compile(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$')


def is_swissprot(accession: str) -> bool:
    """Check if a UniProt accession looks like a reviewed SwissProt entry."""
    return bool(SWISSPROT_PATTERN.match(accession))


def fetch_uniprot_sequence(accession: str) -> Optional[Tuple[str, str, bool]]:
    """
    Fetch protein sequence from UniProt REST API.

    Returns:
        Tuple of (sequence, accession, is_reviewed) or None on failure.
    """
    try:
        resp = requests.get(
            f"{UNIPROT_BASE}/{accession}.json",
            headers={"Accept": "application/json"},
            timeout=10,
        )
        if resp.status_code == 404:
            return None
        resp.raise_for_status()

        data = resp.json()
        seq = data.get("sequence", {}).get("value")
        reviewed = data.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"
        acc = data.get("primaryAccession", accession)

        return (seq, acc, reviewed) if seq else None

    except Exception as e:
        logger.warning(f"Failed to fetch {accession}: {e}")
        return None


def search_reviewed_canonical(gene_symbol: str) -> Optional[Tuple[str, str]]:
    """
    Search UniProt for the reviewed human canonical entry for a gene symbol.

    Returns:
        Tuple of (sequence, accession) or None.
    """
    try:
        query = f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true"
        resp = requests.get(
            f"{UNIPROT_BASE}/search",
            params={
                "query": query,
                "format": "json",
                "fields": "accession,sequence,length",
                "size": 5,
            },
            headers={"Accept": "application/json"},
            timeout=10,
        )
        resp.raise_for_status()

        results = resp.json().get("results", [])
        if not results:
            return None

        # Pick the longest reviewed entry (canonical is usually the longest)
        best = None
        best_len = 0
        for entry in results:
            seq = entry.get("sequence", {}).get("value", "")
            acc = entry.get("primaryAccession", "")
            if len(seq) > best_len:
                best = (seq, acc)
                best_len = len(seq)

        return best

    except Exception as e:
        logger.warning(f"Failed to search for {gene_symbol}: {e}")
        return None


def backfill(gene_symbols=None, limit=None, write_missing=False):
    db = SessionLocal()

    query = db.query(Gene).filter(Gene.uniprot_id.isnot(None), Gene.uniprot_id != "")

    if gene_symbols:
        query = query.filter(Gene.symbol.in_(gene_symbols))

    if write_missing:
        query = query.filter(Gene.canonical_sequence.is_(None))

    if limit:
        query = query.limit(limit)

    genes = query.all()
    total = len(genes)
    logger.info(f"Processing {total} genes (write_missing={write_missing})")

    t_start = time.time()
    fetched = 0
    upgraded = 0  # TrEMBL → SwissProt upgrades
    failed = 0
    already_reviewed = 0
    identical = longer = shorter = written = 0

    def _apply(gene: Gene, seq: str, acc: str) -> None:
        """Compare the UniProt canonical with the designated canonical; never overwrite."""
        nonlocal identical, longer, shorter, written
        if gene.canonical_sequence:
            if gene.canonical_sequence == seq:
                identical += 1
            elif len(seq) > len(gene.canonical_sequence):
                longer += 1
                logger.info(f"  {gene.symbol}: UniProt {acc} ({len(seq)} aa) is LONGER than the "
                            f"designated canonical ({len(gene.canonical_sequence)} aa)")
            else:
                shorter += 1
                logger.info(f"  {gene.symbol}: UniProt {acc} ({len(seq)} aa) differs from / is "
                            f"shorter than the designated canonical ({len(gene.canonical_sequence)} aa)")
        elif write_missing:
            gene.canonical_sequence = seq
            gene.canonical_length = len(seq)
            written += 1

    for i, gene in enumerate(genes):
        try:
            # Step 1: Try fetching the stored uniprot_id directly
            result = fetch_uniprot_sequence(gene.uniprot_id)

            if result:
                seq, acc, is_reviewed = result

                if is_reviewed:
                    _apply(gene, seq, acc)
                    already_reviewed += 1
                    fetched += 1
                else:
                    # Stored ID is TrEMBL — try to find the reviewed entry
                    reviewed_result = search_reviewed_canonical(gene.symbol)
                    if reviewed_result:
                        rev_seq, rev_acc = reviewed_result
                        _apply(gene, rev_seq, rev_acc)
                        upgraded += 1
                        fetched += 1
                        logger.info(
                            f"  {gene.symbol}: reviewed entry for {gene.uniprot_id} (TrEMBL, {len(seq)}aa) "
                            f"is {rev_acc} (SwissProt, {len(rev_seq)}aa)"
                        )
                    else:
                        _apply(gene, seq, acc)
                        fetched += 1
            else:
                # Stored uniprot_id not found — try gene symbol search
                reviewed_result = search_reviewed_canonical(gene.symbol)
                if reviewed_result:
                    rev_seq, rev_acc = reviewed_result
                    _apply(gene, rev_seq, rev_acc)
                    upgraded += 1
                    fetched += 1
                else:
                    failed += 1
                    logger.debug(f"  {gene.symbol}: no sequence found for {gene.uniprot_id}")

            db.commit()

            if (i + 1) % 50 == 0:
                elapsed = time.time() - t_start
                rate = (i + 1) / elapsed
                logger.info(
                    f"[{i+1}/{total}] {fetched} fetched, {upgraded} upgraded, "
                    f"{failed} failed ({rate:.1f} genes/sec)"
                )

            # Rate limit: UniProt allows ~100 req/sec but be polite
            time.sleep(0.1)

        except Exception as e:
            db.rollback()
            failed += 1
            logger.error(f"Error on {gene.symbol} ({gene.uniprot_id}): {e}")

    elapsed = time.time() - t_start
    logger.info(
        f"\nDone in {elapsed:.0f}s — {fetched} fetched ({already_reviewed} already reviewed, "
        f"{upgraded} TrEMBL→SwissProt), {failed} failed | designated canonical vs UniProt: "
        f"{identical} identical, {longer} UniProt longer, {shorter} UniProt differs/shorter, "
        f"{written} written (missing only)"
    )
    db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="QC: designated canonical isoform vs UniProt reviewed canonical (report only)"
    )
    parser.add_argument("--genes", type=str, help="Comma-separated gene symbols")
    parser.add_argument("--limit", type=int, help="Max genes to process")
    parser.add_argument("--write-missing", action="store_true",
                        help="Fill canonical_sequence ONLY for genes that have none (never overwrite)")
    args = parser.parse_args()

    symbols = args.genes.split(",") if args.genes else None
    backfill(gene_symbols=symbols, limit=args.limit, write_missing=args.write_missing)
