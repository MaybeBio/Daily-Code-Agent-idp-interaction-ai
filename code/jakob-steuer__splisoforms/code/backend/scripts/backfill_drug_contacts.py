#!/usr/bin/env python3
"""
backfill_drug_contacts.py — the one script that owns everything read from the
The co-folded drug–protein CIFs. For each drug_complex with a structure it fills:

  - iptm, ptm, interface_pae_min   (from <job>_summary_confidences.json)
  - holo_grade                     (A ≥0.75 / B ≥0.50 / C, from iptm)
  - contact_residues_json, n_contact_residues
        protein residues with any atom within `--cutoff` Å (default 4.5) of any
        ligand atom — the ground-truth drug pocket, defined by real contacts
  - ligand_plddt                   (mean per-atom pLDDT over ligand atoms)
  - pocket_plddt                   (mean per-atom pLDDT over contact-residue atoms)

Structure layout (AF3): chain A = protein (standard amino acids, B-factor column =
per-atom pLDDT); ligand = the non-amino-acid residue(s). Residue numbering is the
1-based protein sequence position, so contacts carry directly onto isoform
sequence alignments downstream (backfill_drug_pocket_fate.py).

Usage:
    python backend/scripts/backfill_drug_contacts.py [--workers N] [--limit N]
                                                     [--cutoff 4.5] [--force]
"""
import os
import sys
import json
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import gemmi
from app.database import SessionLocal
from app.models.bio import DrugComplex

CUTOFF = 4.5


def _is_amino(resname: str) -> bool:
    tr = gemmi.find_tabulated_residue(resname)
    return bool(tr and tr.is_amino_acid())


def compute_metrics(job_name: str, structure_dir: str, cutoff: float) -> dict | None:
    """Pure worker: read CIF + summary JSON, return metrics dict. No DB access."""
    cif = os.path.join(structure_dir, f"{job_name}_model.cif")
    summ = os.path.join(structure_dir, f"{job_name}_summary_confidences.json")
    if not os.path.exists(cif):
        return None
    try:
        st = gemmi.read_structure(cif)
    except Exception as e:  # noqa: BLE001
        return {"job_name": job_name, "error": str(e)}
    model = st[0]

    lig_atoms, prot_res_atoms = [], {}  # (chain,seqid) -> list[atom]
    for ch in model:
        for res in ch:
            if _is_amino(res.name):
                prot_res_atoms[(ch.name, res.seqid.num)] = list(res)
            elif res.name != "HOH":
                lig_atoms.extend(res)
    if not lig_atoms or not prot_res_atoms:
        return {"job_name": job_name, "error": "no ligand or no protein"}

    # neighbour search over the whole model; query each ligand atom
    ns = gemmi.NeighborSearch(model, st.cell, cutoff + 1.0).populate()
    contact_res: set[tuple[str, int]] = set()
    for la in lig_atoms:
        for mark in ns.find_atoms(la.pos, "\0", radius=cutoff):
            cra = mark.to_cra(model)
            if _is_amino(cra.residue.name):
                contact_res.add((cra.chain.name, cra.residue.seqid.num))

    # residue numbers (protein sequence position); if multiple protein chains,
    # keep the numbers — single-chain is the norm for these co-folds.
    contact_nums = sorted({num for (_c, num) in contact_res})
    lig_plddt = sum(a.b_iso for a in lig_atoms) / len(lig_atoms)
    pocket_atoms = [a for key in contact_res for a in prot_res_atoms.get(key, [])]
    pocket_plddt = (sum(a.b_iso for a in pocket_atoms) / len(pocket_atoms)) if pocket_atoms else None

    iptm = ptm = pae_min = None
    if os.path.exists(summ):
        try:
            d = json.load(open(summ))
            iptm = d.get("iptm")
            ptm = d.get("ptm")
            cpm = d.get("chain_pair_pae_min")
            if cpm and len(cpm) >= 2:
                pae_min = cpm[0][1]  # protein(0) ↔ ligand(1)
        except Exception:  # noqa: BLE001
            pass
    grade = None
    if iptm is not None:
        grade = "A" if iptm >= 0.75 else "B" if iptm >= 0.50 else "C"

    return {
        "job_name": job_name,
        "iptm": iptm, "ptm": ptm, "interface_pae_min": pae_min,
        "holo_grade": grade,
        "ligand_plddt": round(lig_plddt, 2),
        "pocket_plddt": round(pocket_plddt, 2) if pocket_plddt is not None else None,
        "n_contact_residues": len(contact_nums),
        "contact_residues_json": json.dumps(contact_nums),
    }


def main():
    ap = argparse.ArgumentParser(description="Extract drug-contact pockets + confidence")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--cutoff", type=float, default=CUTOFF)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--force", action="store_true",
                    help="recompute even if n_contact_residues already set")
    args = ap.parse_args()

    db = SessionLocal()
    q = db.query(DrugComplex.job_name, DrugComplex.structure_dir).filter(
        DrugComplex.has_structure.is_(True))
    if not args.force:
        q = q.filter(DrugComplex.n_contact_residues.is_(None))
    jobs = q.all()
    if args.limit:
        jobs = jobs[:args.limit]
    db.close()
    total = len(jobs)
    print(f"computing contacts for {total} complexes on {args.workers} workers "
          f"(cutoff {args.cutoff} Å)…", flush=True)

    done = err = 0
    batch: list[dict] = []
    db = SessionLocal()

    def flush(rows):
        if not rows:
            return
        db.bulk_update_mappings(
            DrugComplex,
            [{"job_name": r["job_name"], **{k: v for k, v in r.items()
                                            if k not in ("job_name", "error")}}
             for r in rows if "error" not in r],
        )
        db.commit()

    # bulk_update_mappings needs the PK; map job_name -> id up front
    jn_to_id = dict(db.query(DrugComplex.job_name, DrugComplex.id).all())

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(compute_metrics, jn, sd, args.cutoff) for jn, sd in jobs]
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            if r is None or "error" in r:
                err += 1
            else:
                r["id"] = jn_to_id[r["job_name"]]
                batch.append({k: v for k, v in r.items() if k != "job_name"})
            if len(batch) >= 2000:
                db.bulk_update_mappings(DrugComplex, batch); db.commit(); batch = []
            if done % 2000 == 0:
                print(f"  {done}/{total}  ({err} errors)", flush=True)
    if batch:
        db.bulk_update_mappings(DrugComplex, batch); db.commit()

    filled = db.query(DrugComplex).filter(DrugComplex.n_contact_residues.isnot(None)).count()
    print(f"done: {done} processed, {err} errors; {filled} complexes now have contacts")
    db.close()


if __name__ == "__main__":
    main()
