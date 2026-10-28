#!/usr/bin/env python3
"""
ingest_from_metadata.py — source-agnostic ingester for SPLISOFORMS.

Reads the unified `data/sources/<src>/metadata.tsv` files (produced by
`scripts/build_source_metadata.py`, schema in docs/unified_metadata_schema.md) and
builds the core DB entities. Replaces the per-dataset phase scripts
(ingest_gencode_phase1 / ingest_ccrcc_phase2 / ingest_bc_phase3 / rename_to_gencode
/ backfill_source_memberships).

Identity & dedup
----------------
`protein_md5` is the universal identity key. Sources are ingested in registry
`order` (GENCODE first = the reference / dedup anchor). Within a source the rows
are ingested in PRIORITY order, not file order, so that when several transcripts
encode the same protein the surviving Isoform record is the most authoritative one
(review ING-1/ING-2: first-row-wins used to drop 2,276 MANE flags and file the
CDC42 MANE protein under a readthrough gene):

    MANE Select  >  canonical/APPRIS  >  not a readthrough transcript
    >  gene has a real symbol (not a bare ENSG)  >  has a structure  >  file order

For each row (in that order):

  * md5 already an Isoform → add/merge an IsoformSourceMembership (OR-merge tissue
    flags onto both the membership and the Isoform fast-lookup flags). No new row.
  * else, INCLUSION POLICY: create a new Isoform only if the source is_reference OR
    the row has_structure. Otherwise skip (structureless non-reference rows — afdb /
    already_in_db / rejected — carry no analyzable signal).

`transcript_id` is globally UNIQUE; on collision it is namespaced `<src_key>:<tid>`.

Optional metadata columns (schema 30–31): `sample_source` and `gencode_tags` are
stored on `transcript_metadata` when present; older 29-column files still ingest.

Second pass (after all requested sources)
-----------------------------------------
  * Gene canonical: for each gene, pick its is_canonical (MANE > APPRIS) isoform and
    set Gene.canonical_transcript_id / canonical_sequence / canonical_length.
  * delta_length: per isoform, protein_length − gene.canonical_length (when known).

`canonical_match_id` / `gencode_match_identity` for novel isoforms are NOT set here;
they are a separate post-ingest step (`backfill_canonical_match.py`).

Usage
-----
    python backend/scripts/ingest_from_metadata.py --source gencode   # one source
    python backend/scripts/ingest_from_metadata.py --all              # all, in order
    python backend/scripts/ingest_from_metadata.py --all --dry-run
    python backend/scripts/ingest_from_metadata.py --finalize-only     # 2nd pass only

Idempotent: re-running upserts memberships and skips existing isoforms.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

# backend/ root + scripts/ dir on path (mirrors reset_db.py / backfill convention)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sqlalchemy.orm import Session

from app.database import SessionLocal
from app.models.bio import (
    Gene, Isoform, IsoformSourceMembership, Source, Structure,
    TranscriptMetadata, SourceTypeEnum,
)
from _sources_registry import SourceDef, ordered_sources, SOURCES_BY_KEY

csv.field_size_limit(10_000_000)

REQUIRED_COLUMNS = {
    "protein_md5", "source", "transcript_id", "orig_transcript_id", "gene_symbol",
    "gene_id", "seq_len", "sequence", "chrom", "strand", "genomic_start",
    "genomic_end", "exon_count", "biotype", "structural_category", "subcategory",
    "is_canonical", "is_mane_select", "appris_tag", "fate", "structure_provenance",
    "has_structure", "exons", "in_tumor", "in_normal", "in_metastasis",
    "in_fibroblast", "best_hit_acc", "best_hit_identity",
}

OPTIONAL_COLUMNS = {"sample_source", "gencode_tags"}

BATCH = 2000


def _b(v: str) -> bool:
    return str(v).strip().lower() == "true"


def row_priority(row: Dict[str, str], idx: int) -> tuple:
    """Sort key: lower sorts first = preferred survivor of an md5 collision."""
    tags = row.get("gencode_tags") or ""
    sym = (row.get("gene_symbol") or "").strip()
    return (
        0 if _b(row.get("is_mane_select", "")) else 1,
        0 if _b(row.get("is_canonical", "")) else 1,
        1 if "readthrough_transcript" in tags else 0,
        1 if (not sym or sym.startswith("ENSG")) else 0,
        0 if _b(row.get("has_structure", "")) else 1,
        idx,
    )


def load_rows_by_priority(path: Path) -> list:
    """Read a metadata.tsv fully and return its rows in ingest-priority order."""
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    order = sorted(range(len(rows)), key=lambda i: row_priority(rows[i], i))
    return [rows[i] for i in order]


def _int(v: str) -> Optional[int]:
    v = (v or "").strip()
    try:
        return int(v)
    except (ValueError, TypeError):
        return None


# ── Gene resolution (adapted from ingest_ccrcc_phase2.resolve_or_create_gene) ────
def resolve_or_create_gene(db: Session, *, gene_symbol: str, gene_id: str,
                           uniprot_fallback: Optional[str], md5: str,
                           chrom: str, strand: str, ensembl_gene_id: str) -> Gene:
    """Find or create a Gene. Priority: ensembl_gene_id → symbol → uniprot →
    source gene_id → stub `UNK_<md5[:10]>`. Each candidate is re-checked before
    INSERT to avoid UNIQUE violations from concurrent stub collapse."""
    if ensembl_gene_id:
        g = db.query(Gene).filter_by(ensembl_gene_id=ensembl_gene_id).first()
        if g:
            return g
    if gene_symbol:
        g = db.query(Gene).filter_by(symbol=gene_symbol).first()
        if g:
            return g
    if uniprot_fallback:
        g = db.query(Gene).filter_by(uniprot_id=uniprot_fallback).first()
        # Only accept the UniProt best-hit gene if it is on the same chromosome as
        # this transcript. Otherwise a distant paralog (sharing the same UniProt
        # best hit, e.g. a UGT2-family member) captures transcripts from an
        # unrelated locus — the cross-chromosome mis-attribution bug. Fall through
        # to source gene_id / stub when the best-hit gene sits on another chromosome.
        if g and (not chrom or not g.chromosome or g.chromosome == chrom):
            return g
    if gene_id:
        g = db.query(Gene).filter_by(symbol=gene_id).first()
        if g:
            return g

    symbol = gene_symbol or gene_id or uniprot_fallback or f"UNK_{md5[:10]}"
    g = db.query(Gene).filter_by(symbol=symbol).first()
    if g:
        return g
    # Don't claim a UniProt accession another gene already holds. We only reach
    # here for the best-hit accession when the cross-chromosome guard above
    # rejected the gene that owns it (a different locus) — reusing it would
    # violate ix_genes_uniprot_id. Leave it NULL; gene-attribution backfills
    # (merge_stub_genes_into_canonicals) reconcile these gene-duplication cases.
    uid = uniprot_fallback or None
    if uid and db.query(Gene).filter_by(uniprot_id=uid).first():
        uid = None
    g = Gene(
        symbol=symbol,
        uniprot_id=uid,
        ensembl_gene_id=ensembl_gene_id or None,
        chromosome=chrom or None,
        strand=strand or None,
    )
    db.add(g)
    db.flush()
    return g


def _tissue_flags(row: Dict[str, str]) -> Dict[str, bool]:
    # metadata `in_tumor` → model `in_primary_tumor`; others 1:1
    return {
        "in_normal": _b(row["in_normal"]),
        "in_primary_tumor": _b(row["in_tumor"]),
        "in_metastasis": _b(row["in_metastasis"]),
        "in_fibroblast": _b(row["in_fibroblast"]),
    }


def _validate_columns(path: Path) -> None:
    with path.open() as fh:
        header = set(next(csv.reader(fh, delimiter="\t")))
    missing = REQUIRED_COLUMNS - header
    if missing:
        raise SystemExit(f"{path}: metadata is missing required columns: "
                         f"{sorted(missing)}")


def ensure_source(db: Session, sd: SourceDef) -> Source:
    src = db.query(Source).filter_by(name=sd.name).first()
    if src is None:
        src = Source(
            name=sd.name,
            source_type=SourceTypeEnum(sd.source_type),
            cancer_type=sd.cancer_type or None,
            tissue_type=sd.tissue_type or None,
            sequencing_platform=sd.sequencing_platform or None,
            notes=sd.notes or None,
        )
        db.add(src)
        db.flush()
    return src


# ── Per-source ingest ────────────────────────────────────────────────────────────
def ingest_source(db: Session, sd: SourceDef, taken_tids: Set[str],
                  seen_memberships: Set[Tuple[int, int]], dry_run: bool) -> Dict:
    stats = {"rows": 0, "new_isoform": 0, "membership": 0, "membership_same_source": 0,
             "skip_no_structure": 0, "tid_namespaced": 0, "structures": 0}
    _validate_columns(sd.metadata_tsv)
    src = None if dry_run else ensure_source(db, sd)

    rows = load_rows_by_priority(sd.metadata_tsv)
    if True:
        for row in rows:
            stats["rows"] += 1
            md5 = row["protein_md5"].strip()
            if not md5:
                continue
            flags = _tissue_flags(row)

            existing = db.query(Isoform).filter_by(protein_md5=md5).first()
            if existing is not None:
                stats["membership"] += 1
                if not dry_run and existing.source_id == src.id:
                    stats["membership_same_source"] += 1
                if dry_run:
                    continue
                key = (existing.id, src.id)
                if key not in seen_memberships:
                    db.add(IsoformSourceMembership(isoform_id=existing.id,
                                                   source_id=src.id, **flags))
                    seen_memberships.add(key)
                else:
                    m = (db.query(IsoformSourceMembership)
                         .filter_by(isoform_id=existing.id, source_id=src.id).first())
                    if m:
                        for k, v in flags.items():
                            if v and not getattr(m, k):
                                setattr(m, k, True)
                # OR-merge onto the Isoform fast-lookup flags
                for k, v in flags.items():
                    if v and not getattr(existing, k):
                        setattr(existing, k, True)
                _commit_batch(db, stats["rows"])
                continue

            # INCLUSION POLICY
            if not (sd.is_reference or _b(row["has_structure"])):
                stats["skip_no_structure"] += 1
                continue
            stats["new_isoform"] += 1
            if dry_run:
                continue

            # transcript_id collision guard
            tid = row["transcript_id"].strip()
            if tid in taken_tids:
                tid = f"{sd.key}:{tid}"
                stats["tid_namespaced"] += 1
            taken_tids.add(tid)

            is_ref = sd.is_reference
            ensembl_tid = row["orig_transcript_id"] if is_ref else None
            gene = resolve_or_create_gene(
                db, gene_symbol=row["gene_symbol"].strip(),
                gene_id=row["gene_id"].strip(),
                uniprot_fallback=(row["best_hit_acc"].strip() or None),
                md5=md5, chrom=row["chrom"], strand=row["strand"],
                ensembl_gene_id=(row["gene_id"].strip() if is_ref else ""),
            )

            iso = Isoform(
                transcript_id=tid,
                sequence=row["sequence"],
                gene_id=gene.id,
                source_id=src.id,
                protein_md5=md5,
                legacy_transcript_id=row["transcript_id"].strip(),
                ensembl_transcript_id=ensembl_tid,
                gencode_version=("v47" if is_ref else None),
                is_canonical=_b(row["is_canonical"]),
                is_mane_select=_b(row["is_mane_select"]),
                structure_provenance=(row["structure_provenance"] or None),
                exons_json=(row["exons"] or None),
                **flags,
            )
            db.add(iso)
            db.flush()

            db.add(IsoformSourceMembership(isoform_id=iso.id, source_id=src.id, **flags))
            seen_memberships.add((iso.id, src.id))

            if _b(row["has_structure"]):
                db.add(Structure(
                    isoform_id=iso.id,
                    af3_cif_url=f"{sd.struct_url_prefix}/{md5}/{md5}_model.cif",
                    af3_pae_url=f"{sd.struct_url_prefix}/{md5}/{md5}_confidences.json",
                ))
                stats["structures"] += 1

            db.add(TranscriptMetadata(
                isoform_id=iso.id,
                transcript_id=tid,
                metadata_kind=("gencode" if is_ref else "sqanti3"),
                chrom=(row["chrom"] or None),
                strand=(row["strand"] or None),
                length=_int(row["seq_len"]),
                protein_length_aa=(float(_int(row["seq_len"])) if _int(row["seq_len"]) else None),
                biotype=(row["biotype"] or None),
                structural_category=(row["structural_category"] or None),
                subcategory=(row["subcategory"] or None),
                appris_tag=(row["appris_tag"] or None),
                is_mane_select=_b(row["is_mane_select"]),
                enst_id=(ensembl_tid if is_ref else None),
                uniprot_id=(row["best_hit_acc"].strip() or None),
                sample_source=((row.get("sample_source") or "").strip() or None),
                # ARRAY(String) column: metadata carries a comma-joined list
                gencode_tags=([t for t in (row.get("gencode_tags") or "").split(",") if t]
                              or None),
            ))
            _commit_batch(db, stats["rows"])

    if not dry_run:
        db.commit()
    return stats


def _commit_batch(db: Session, n: int) -> None:
    if n % BATCH == 0:
        db.commit()


# ── Second pass: gene canonical + delta_length ──────────────────────────────────
def finalize(db: Session) -> Dict:
    stats = {"genes_canonical": 0, "delta_set": 0}
    # 1. canonical per gene: prefer is_mane_select, then is_canonical, then longest
    genes = db.query(Gene).all()
    for g in genes:
        isos = [i for i in g.isoforms]
        if not isos:
            continue
        cands = [i for i in isos if i.is_canonical] or isos
        winner = max(cands, key=lambda i: (
            1 if i.is_mane_select else 0,
            len(i.sequence or ""),
        ))
        g.canonical_transcript_id = winner.id
        g.canonical_sequence = winner.sequence
        g.canonical_length = len(winner.sequence) if winner.sequence else None
        stats["genes_canonical"] += 1
    db.commit()

    # 2. delta_length per isoform
    for g in db.query(Gene).all():
        if not g.canonical_length:
            continue
        for i in g.isoforms:
            ln = len(i.sequence) if i.sequence else None
            if ln is not None:
                i.delta_length = ln - g.canonical_length
                stats["delta_set"] += 1
    db.commit()
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=[s.key for s in ordered_sources()])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--finalize-only", action="store_true",
                    help="Run only the gene-canonical + delta_length pass.")
    ap.add_argument("--no-finalize", action="store_true",
                    help="Skip the second pass (use when ingesting sources one at "
                         "a time; run --finalize-only after the last).")
    args = ap.parse_args()

    db = SessionLocal()
    try:
        if args.finalize_only:
            print("Finalize pass (gene canonical + delta_length) ...")
            print("  ", finalize(db))
            return

        if args.all:
            targets = ordered_sources()
        elif args.source:
            targets = [SOURCES_BY_KEY[args.source]]
        else:
            ap.error("specify --source <key>, --all, or --finalize-only")

        # pre-seed idempotency caches from existing DB state
        taken_tids: Set[str] = set(
            x[0] for x in db.query(Isoform.transcript_id).all())
        seen_memberships: Set[Tuple[int, int]] = set(
            (m.isoform_id, m.source_id)
            for m in db.query(IsoformSourceMembership).all())

        for sd in targets:
            print(f"\n=== Ingesting {sd.key} ({sd.name}) ===")
            st = ingest_source(db, sd, taken_tids, seen_memberships, args.dry_run)
            print(f"  rows={st['rows']} new_isoform={st['new_isoform']} "
                  f"membership={st['membership']} "
                  f"(same-source collapse={st['membership_same_source']}) "
                  f"skip_no_structure={st['skip_no_structure']} "
                  f"structures={st['structures']} tid_namespaced={st['tid_namespaced']}")

        if args.all and not args.dry_run and not args.no_finalize:
            print("\nFinalize pass (gene canonical + delta_length) ...")
            print("  ", finalize(db))
    finally:
        db.close()


if __name__ == "__main__":
    main()
