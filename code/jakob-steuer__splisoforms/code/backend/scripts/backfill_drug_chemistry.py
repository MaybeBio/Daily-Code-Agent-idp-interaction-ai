#!/usr/bin/env python3
"""
backfill_drug_chemistry.py — fill drug chemistry from PubChem.

The drug layer is ingested from the drug-dataset manifest, which carries only pubchem_cid,
broad_id and name. `smiles`, `inchikey`, `molecular_formula`, `molecular_weight` and
`iupac_name` all arrive empty, so the drug page has no chemistry to display.

This fetches them from PubChem PUG-REST, keyed on pubchem_cid.

BATCHED: PubChem accepts many CIDs in one property request, so all ~3,725 drugs resolve in
~40 requests rather than 3,725. Keep CHUNK modest — a URL with too many CIDs is rejected.

NOTE the property is `SMILES`, not `CanonicalSMILES`: PubChem deprecated the latter, and it
silently returns an EMPTY string inside a multi-property batch rather than erroring, which
is an easy way to backfill blank SMILES without noticing.

Dependencies: network (PubChem), DB table `drugs`.
Idempotent — only fetches drugs missing chemistry unless --force.

Usage:
    python backend/scripts/backfill_drug_chemistry.py [--dry-run] [--force]
"""
import argparse
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sqlalchemy as sa  # noqa: E402

from _backfill_io import Progress, print_summary, setup_logging  # noqa: E402
from app.database import SessionLocal  # noqa: E402
from app.models.bio import Drug  # noqa: E402

log = logging.getLogger("backfill_drug_chemistry")

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
PROPS = "MolecularFormula,MolecularWeight,SMILES,InChIKey,IUPACName"
CHUNK = 100          # CIDs per request
THROTTLE_S = 0.25    # PubChem asks for <= 5 requests/second
TIMEOUT_S = 30


def fetch_chunk(cids: list[str]) -> dict[str, dict]:
    """CID -> property dict for one batch. Returns {} on a failed chunk (logged, not fatal)."""
    url = f"{PUBCHEM}/{','.join(cids)}/property/{PROPS}/JSON"
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT_S) as resp:
            data = json.load(resp)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
        log.warning(f"chunk of {len(cids)} failed: {exc}")
        return {}

    out: dict[str, dict] = {}
    for p in data.get("PropertyTable", {}).get("Properties", []):
        cid = str(p.get("CID"))
        if cid:
            out[cid] = p
    return out


def run(dry_run: bool, force: bool) -> None:
    db = SessionLocal()
    try:
        q = db.query(Drug).filter(Drug.pubchem_cid.isnot(None))
        if not force:
            q = q.filter(
                sa.or_(
                    Drug.molecular_formula.is_(None),
                    Drug.smiles.is_(None),
                )
            )
        drugs = q.all()
        by_cid = {d.pubchem_cid: d for d in drugs}
        cids = list(by_cid)

        log.info(f"{len(cids):,} drugs to enrich  ({(len(cids) + CHUNK - 1) // CHUNK} requests)")

        ok = missing = 0
        with Progress(len(cids), "PubChem", unit="drug") as pbar:
            for i in range(0, len(cids), CHUNK):
                batch = cids[i:i + CHUNK]
                props = fetch_chunk(batch)

                for cid in batch:
                    p = props.get(cid)
                    drug = by_cid[cid]
                    if not p:
                        missing += 1
                    else:
                        drug.molecular_formula = p.get("MolecularFormula")
                        mw = p.get("MolecularWeight")
                        drug.molecular_weight = float(mw) if mw not in (None, "") else None
                        drug.smiles = p.get("SMILES") or None
                        drug.inchikey = p.get("InChIKey") or None
                        drug.iupac_name = p.get("IUPACName") or None
                        ok += 1
                    pbar.update(1, status=f"ok={ok} miss={missing}")

                if not dry_run:
                    db.commit()
                time.sleep(THROTTLE_S)

        if dry_run:
            db.rollback()
            log.info("--dry-run: rolled back.")

        print_summary("DRUG CHEMISTRY BACKFILL", [
            ("✓ Enriched", ok),
            ("– Not found in PubChem", missing),
        ], total=len(cids))
    finally:
        db.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="Fetch but do not write.")
    ap.add_argument("--force", action="store_true", help="Re-fetch drugs that already have chemistry.")
    args = ap.parse_args()

    setup_logging("backfill_drug_chemistry")
    run(args.dry_run, args.force)


if __name__ == "__main__":
    main()
