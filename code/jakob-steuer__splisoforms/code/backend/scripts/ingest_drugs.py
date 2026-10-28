#!/usr/bin/env python3
"""
ingest_drugs.py — load the drug-complex metadata TSVs (from build_drug_metadata.py)
into the DB: populate `drugs` and `drug_complexes`. Idempotent (upsert on the
natural keys pubchem_cid/slug for drugs, job_name for complexes).

Structure-derived metrics (iptm, ligand_plddt, contacts, …) are NOT filled here —
backfill_drug_contacts.py owns everything read from the CIFs.

Usage:
    python backend/scripts/ingest_drugs.py [--drugs-dir data/drugs] [--limit N]
"""
import os
import sys
import csv
import argparse
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal
from app.models.bio import Drug, DrugComplex

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def _int(x):
    x = (x or "").strip()
    return int(x) if x.isdigit() else None


def main():
    ap = argparse.ArgumentParser(description="Ingest drugs + drug_complexes from TSVs")
    ap.add_argument("--drugs-dir", type=Path, default=REPO_ROOT / "data" / "drugs")
    ap.add_argument("--limit", type=int, help="Cap complexes ingested (debug)")
    args = ap.parse_args()

    drugs_tsv = args.drugs_dir / "drugs.tsv"
    complexes_tsv = args.drugs_dir / "drug_complexes.tsv"
    for p in (drugs_tsv, complexes_tsv):
        if not p.exists():
            sys.exit(f"missing {p} — run build_drug_metadata.py first")

    db = SessionLocal()

    # ── drugs ────────────────────────────────────────────────────────────────
    # get-or-create keyed on pubchem_cid when present, else slug.
    existing_by_cid = {d.pubchem_cid: d for d in db.query(Drug).filter(Drug.pubchem_cid.isnot(None))}
    existing_by_slug = {d.slug: d for d in db.query(Drug).filter(Drug.pubchem_cid.is_(None))}
    drug_key_to_id: dict[str, int] = {}
    n_new_drugs = 0

    with open(drugs_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row["pubchem_cid"] or "").strip() or None
            slug = row["slug"].strip()
            d = existing_by_cid.get(cid) if cid else existing_by_slug.get(slug)
            if d is None:
                d = Drug(pubchem_cid=cid, slug=slug, name=row["name"].strip(),
                         broad_id=(row["broad_id"] or "").strip() or None)
                db.add(d)
                db.flush()  # assign id
                n_new_drugs += 1
                (existing_by_cid if cid else existing_by_slug)[cid or slug] = d
            drug_key_to_id[row["drug_key"]] = d.id
    db.commit()
    print(f"drugs: {n_new_drugs} new, {len(drug_key_to_id)} total keys", flush=True)

    # ── drug_complexes ───────────────────────────────────────────────────────
    existing_jobs = {c.job_name: c for c in db.query(DrugComplex).all()}
    n_new_c = n_upd_c = 0

    with open(complexes_tsv, newline="") as fh:
        for i, row in enumerate(csv.DictReader(fh, delimiter="\t")):
            if args.limit and i >= args.limit:
                break
            job = row["job_name"]
            drug_id = drug_key_to_id.get(row["drug_key"])
            if drug_id is None:
                continue  # drug row absent (shouldn't happen)
            fields = dict(
                drug_id=drug_id,
                gene_id=_int(row["gene_id"]),
                target_isoform_id=_int(row["target_isoform_id"]),
                target_uniprot_id=(row["target_uniprot_id"] or "").strip() or None,
                target_enst=(row["target_enst"] or "").strip() or None,
                seq_md5=(row["seq_md5"] or "").strip() or None,
                structure_dir=(row["structure_dir"] or "").strip() or None,
                target_length=_int(row["length"]),
                has_structure=bool(_int(row["has_structure"])),
            )
            c = existing_jobs.get(job)
            if c is None:
                db.add(DrugComplex(job_name=job, **fields))
                n_new_c += 1
            else:
                for k, v in fields.items():
                    setattr(c, k, v)
                n_upd_c += 1
            if (n_new_c + n_upd_c) % 5000 == 0:
                db.commit()
    db.commit()

    total = db.query(DrugComplex).count()
    with_struct = db.query(DrugComplex).filter(DrugComplex.has_structure.is_(True)).count()
    print(f"drug_complexes: {n_new_c} new, {n_upd_c} updated; "
          f"{total} total ({with_struct} with structure)")
    db.close()


if __name__ == "__main__":
    main()
