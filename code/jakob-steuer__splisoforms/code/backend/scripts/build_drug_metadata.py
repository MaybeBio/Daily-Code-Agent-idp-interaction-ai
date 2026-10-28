#!/usr/bin/env python3
"""
build_drug_metadata.py — join the AF3 co-folded drug–protein complex set to
the live SPLISOFORMS DB and emit an ingest-ready metadata artifact.

This is the drug analogue of build_source_metadata.py, but drugs are NOT a new
isoform *source*: each co-fold is a drug docked against a protein that already
exists in the DB. So this builder resolves every complex to an existing gene
(via UniProt) and, where possible, to the exact drugged isoform (via protein-seq
MD5, the strongest anchor — ENST is a weaker fallback because of transcript
renames). Nothing is written to the DB; this only produces:

    data/drugs/drugs.tsv            one row per unique drug (keyed on PubChem CID)
    data/drugs/drug_complexes.tsv   one row per co-fold, with resolved gene_id /
                                    target_isoform_id (or empty when unresolved)

and prints a coverage report. ingest_drugs.py consumes these TSVs.

Usage:
    python backend/scripts/build_drug_metadata.py \
        [--drug-root /path/to/drug_dataset] \
        [--out data/drugs]
"""
import os
import sys
import csv
import argparse
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal
from app.models.bio import Gene, Isoform

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
# Root of the co-folded drug dataset (af3_jobs_v2/manifest.csv + structures/).
# Set SPLISOFORMS_DRUG_DATASET_ROOT or pass --drug-root.
DEFAULT_DRUG_ROOT = Path(
    os.environ.get("SPLISOFORMS_DRUG_DATASET_ROOT")
    or "/path/to/drug_dataset"
)


def slugify(name: str) -> str:
    return "".join(c if c.isalnum() else "-" for c in name.strip().lower()).strip("-")


def main():
    ap = argparse.ArgumentParser(description="Build drug-complex ingest metadata")
    ap.add_argument("--drug-root", dest="drug_root", type=Path, default=DEFAULT_DRUG_ROOT,
                    help="Root of the co-folded drug dataset (manifest + structures)")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "data" / "drugs",
                    help="Output directory for the metadata TSVs")
    args = ap.parse_args()

    manifest = args.drug_root / "af3_jobs_v2" / "manifest.csv"
    structures_root = args.drug_root / "structures"
    if not manifest.exists():
        sys.exit(f"manifest not found: {manifest}")
    args.out.mkdir(parents=True, exist_ok=True)

    # ── Load DB join maps into memory (one pass each) ────────────────────────
    db = SessionLocal()
    print("loading DB join maps …", flush=True)
    uniprot_to_gene = dict(
        db.query(Gene.uniprot_id, Gene.id).filter(Gene.uniprot_id.isnot(None)).all()
    )
    # md5 -> (isoform_id, gene_id); many isoforms may share an md5, keep the first
    # (a MANE/canonical row is preferred so the drugged transcript resolves to the
    # reference copy where possible).
    md5_to_iso: dict[str, tuple[int, int]] = {}
    for iso_id, gene_id, md5, mane, canon in db.query(
        Isoform.id, Isoform.gene_id, Isoform.protein_md5,
        Isoform.is_mane_select, Isoform.is_canonical,
    ).filter(Isoform.protein_md5.isnot(None)).all():
        prev = md5_to_iso.get(md5)
        if prev is None:
            md5_to_iso[md5] = (iso_id, gene_id, bool(mane), bool(canon))
        else:
            # upgrade to a MANE/canonical row if we find one
            _, _, pm, pc = prev
            if (mane and not pm) or (canon and not pc and not pm):
                md5_to_iso[md5] = (iso_id, gene_id, bool(mane), bool(canon))
    enst_to_iso = dict(
        db.query(Isoform.ensembl_transcript_id, Isoform.id)
        .filter(Isoform.ensembl_transcript_id.isnot(None)).all()
    )
    db.close()
    print(f"  genes(uniprot)={len(uniprot_to_gene)}  isoforms(md5)={len(md5_to_iso)}  "
          f"isoforms(enst)={len(enst_to_iso)}", flush=True)

    # ── Walk the manifest ────────────────────────────────────────────────────
    drugs: dict[str, dict] = {}          # cid|slug -> drug record
    complexes: list[dict] = []
    stats = {
        "rows": 0, "struct_missing": 0,
        "gene_hit": 0, "iso_md5_hit": 0, "iso_enst_hit": 0, "iso_any_hit": 0,
        "grade_from_summary": 0,
    }

    with open(manifest, newline="") as fh:
        for row in csv.DictReader(fh):
            stats["rows"] += 1
            job = row["job_name"]
            cid = (row.get("cid") or "").strip()
            name = (row.get("drug_name") or "").strip()
            uniprot = (row.get("uniprot_id") or "").strip()
            enst = (row.get("transcript_id") or "").strip()
            md5 = (row.get("seq_md5") or "").strip()

            # unique-drug key: CID when present, else slug of the name
            drug_key = f"cid:{cid}" if cid else f"name:{slugify(name)}"
            if drug_key not in drugs:
                drugs[drug_key] = {
                    "drug_key": drug_key,
                    "pubchem_cid": cid,
                    "broad_id": (row.get("broad_id") or "").strip(),
                    "name": name,
                    "slug": slugify(name),
                    "n_complexes": 0,
                }
            drugs[drug_key]["n_complexes"] += 1

            # structure presence
            sdir = structures_root / job
            cif = sdir / f"{job}_model.cif"
            has_struct = cif.exists()
            if not has_struct:
                stats["struct_missing"] += 1

            # DB resolution
            gene_id = uniprot_to_gene.get(uniprot, "")
            if gene_id:
                stats["gene_hit"] += 1
            iso_hit = md5_to_iso.get(md5)
            target_iso_id = ""
            resolved_via = ""
            if iso_hit:
                target_iso_id = iso_hit[0]
                resolved_via = "md5"
                stats["iso_md5_hit"] += 1
                if not gene_id:  # backfill gene from the isoform when uniprot missed
                    gene_id = iso_hit[1]
            elif enst in enst_to_iso:
                target_iso_id = enst_to_iso[enst]
                resolved_via = "enst"
                stats["iso_enst_hit"] += 1
            if target_iso_id:
                stats["iso_any_hit"] += 1

            complexes.append({
                "job_name": job,
                "drug_key": drug_key,
                "pubchem_cid": cid,
                "drug_name": name,
                "target_uniprot_id": uniprot,
                "target_enst": enst,
                "seq_md5": md5,
                "length": (row.get("length") or "").strip(),
                "largemem": (row.get("largemem") or "").strip(),
                "gene_id": gene_id,
                "target_isoform_id": target_iso_id,
                "resolved_via": resolved_via,
                "structure_dir": str(sdir),
                "has_structure": int(has_struct),
            })

    # ── Write artifacts ──────────────────────────────────────────────────────
    drugs_tsv = args.out / "drugs.tsv"
    with open(drugs_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(next(iter(drugs.values())).keys()),
                           delimiter="\t")
        w.writeheader()
        for d in sorted(drugs.values(), key=lambda r: -r["n_complexes"]):
            w.writerow(d)

    complexes_tsv = args.out / "drug_complexes.tsv"
    with open(complexes_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(complexes[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(complexes)

    # ── Report ───────────────────────────────────────────────────────────────
    # gene/isoform resolution is counted from the final complex records (which
    # include the md5 gene-backfill), not the raw uniprot-hit tally.
    gene_resolved = sum(1 for c in complexes if c["gene_id"] != "")
    iso_resolved = sum(1 for c in complexes if c["target_isoform_id"] != "")
    n = stats["rows"]
    pct = lambda x: f"{100*x/n:.1f}%"
    print("\n─── drug-complex ingest metadata ─────────────────────────────")
    print(f"manifest rows                : {n}")
    print(f"unique drugs                 : {len(drugs)}")
    print(f"structures present on disk   : {n - stats['struct_missing']} ({pct(n-stats['struct_missing'])})")
    print(f"resolved to a gene           : {gene_resolved} ({pct(gene_resolved)})")
    print(f"resolved to an isoform (any) : {iso_resolved} ({pct(iso_resolved)})")
    print(f"    via protein_md5          : {stats['iso_md5_hit']} ({pct(stats['iso_md5_hit'])})")
    print(f"    via ENST fallback        : {stats['iso_enst_hit']} ({pct(stats['iso_enst_hit'])})")
    print(f"UNRESOLVED (no gene)         : {n - gene_resolved} ({pct(n - gene_resolved)})")
    print(f"\nwrote:\n  {drugs_tsv}\n  {complexes_tsv}")


if __name__ == "__main__":
    main()
