"""
Backfill per-sample expression metrics (ccRCC) from the SQANTI3 CustomFiltered
classification files and (optionally) the per-sample sparse cell matrices.

Sources
-------
1. FL counts (bulk-equivalent)
   ``data/sources/ccrcc/expression/<sample>/
       scisoseq_classification.filtered_lite_custom_classification.txt``
   Column ``FL``: raw full-length read count per isoform per sample.

2. Cell prevalence (single-cell)
   ``data/sources/ccrcc/expression/<sample>/isoforms_seurat/
       matrix.mtx, genes.tsv, barcodes.tsv``
   10x-style sparse matrix: rows=isoforms, cols=cells.

Identity (review EXP-1)
-----------------------
PacBio ``PB.*`` ids are numbered independently per sample run, so the same id
can name different transcripts in different samples. An isoform's reads must
therefore be looked up ONLY in the sample(s) its transcript(s) were assembled in.
That provenance lives in ``data/sources/ccrcc/metadata.tsv`` (``sample_source``,
one row per (pbid, sample)); several metadata rows can share a ``protein_md5``
(identical protein assembled in several samples) and collapse onto one DB
isoform, so the DB isoform's expression is the union over all of them:

    isoform  ──protein_md5──▶  {(pbid, sample), …}  ──▶  FL in that sample's file

This also covers GENCODE isoforms that carry a ccRCC source membership (same
protein seen in the cohort), which the old transcript_id join never reached.

Derived metrics
---------------
  cpm          = fl_count / total_fl_in_sample * 1_000_000
  isoform_frac = fl_count / sum(fl_count for all DB isoforms of same gene in sample)
  n_cells      = number of cells with ≥1 UMI (from sparse matrix, if available)
  cell_frac    = n_cells / total_cells_in_sample

Usage
-----
    python backend/scripts/backfill_expression.py [--force] [--gene SYMBOL] [--no-cells]
"""
import argparse
import csv
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from sqlalchemy import text

from _backfill_io import exit_status
from _sources_registry import SOURCES_BY_KEY
from app.database import SessionLocal
from app.models.bio import IsoformExpression

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

csv.field_size_limit(10_000_000)

# ── Paths ─────────────────────────────────────────────────────────────────────
_BACKEND_DIR = Path(__file__).resolve().parent.parent
_DATA_ROOT   = _BACKEND_DIR.parent / "data" / "sources" / "ccrcc" / "expression"
_METADATA    = SOURCES_BY_KEY["ccrcc"].metadata_tsv

SAMPLES = {
    "Normal": "normal",
    "ccRCC2": "tumor",
    "ccRCC3": "tumor",
    "ccRCC4": "tumor",
    "ccRCC5": "tumor",
}


# ── FL count loader ───────────────────────────────────────────────────────────
def load_fl_counts(sample: str) -> dict[str, int]:
    """{pbid: fl_count} for *sample* (rows with FL NA or 0 skipped)."""
    path = _DATA_ROOT / sample / "scisoseq_classification.filtered_lite_custom_classification.txt"
    if not path.exists():
        log.warning(f"Classification file not found: {path}")
        return {}
    counts: dict[str, int] = {}
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        fl_col = header.index("FL")
        for line in fh:
            parts = line.rstrip().split("\t")
            fl_str = parts[fl_col]
            if fl_str == "NA":
                continue
            try:
                fl = int(float(fl_str))
            except ValueError:
                continue
            if fl > 0:
                counts[parts[0]] = fl
    log.info(f"  {sample}: {len(counts):,} transcripts with FL > 0")
    return counts


# ── Sparse matrix loader ──────────────────────────────────────────────────────
def load_cell_counts(sample: str) -> tuple[dict[str, int], int]:
    """({pbid: n_cells_with_umi}, total_cells) for *sample*."""
    seurat_dir    = _DATA_ROOT / sample / "isoforms_seurat"
    genes_path    = seurat_dir / "genes.tsv"
    barcodes_path = seurat_dir / "barcodes.tsv"
    matrix_path   = seurat_dir / "matrix.mtx"
    for p in (genes_path, barcodes_path, matrix_path):
        if not p.exists():
            log.warning(f"  {sample}: sparse matrix file missing: {p}")
            return {}, 0

    genes: list[str] = []
    with open(genes_path) as fh:
        for line in fh:
            genes.append(line.rstrip().split("\t")[0].split(":")[0])
    total_cells = sum(1 for _ in open(barcodes_path))

    cell_counts: dict[int, int] = defaultdict(int)
    with open(matrix_path) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) == 3 and float(parts[2]) > 0:
                cell_counts[int(parts[0])] += 1

    result = {genes[i - 1]: n for i, n in cell_counts.items() if 1 <= i <= len(genes)}
    log.info(f"  {sample}: {len(result):,} transcripts with cell data ({total_cells} cells)")
    return result, total_cells


# ── Provenance: protein_md5 → [(pbid, sample)] from the ccRCC metadata ────────
def load_md5_provenance() -> dict[str, list[tuple[str, str]]]:
    prov: dict[str, list[tuple[str, str]]] = defaultdict(list)
    if not _METADATA.exists():
        log.error(f"ccRCC metadata not found: {_METADATA}")
        return prov
    with _METADATA.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "sample_source" not in (reader.fieldnames or []):
            log.error("ccRCC metadata.tsv has no `sample_source` column — rebuild it with "
                      "scripts/build_source_metadata.py --source ccrcc (schema ≥31 cols)")
            return prov
        for r in reader:
            sample = (r.get("sample_source") or "").strip()
            pbid = (r.get("orig_transcript_id") or r.get("transcript_id") or "").split("@")[0]
            if sample in SAMPLES and pbid and r["protein_md5"]:
                prov[r["protein_md5"]].append((pbid, sample))
    log.info(f"provenance: {len(prov):,} proteins → {sum(len(v) for v in prov.values()):,} (pbid, sample) pairs")
    return prov


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Backfill per-sample isoform expression (ccRCC)")
    parser.add_argument("--force",    action="store_true", help="Delete existing records and recompute")
    parser.add_argument("--gene",     type=str,            help="Limit to a specific gene symbol")
    parser.add_argument("--no-cells", action="store_true", help="Skip sparse cell matrix loading")
    args = parser.parse_args()

    db = SessionLocal()

    log.info("Loading FL counts from CustomFiltered SQANTI3 files...")
    fl_by_sample:   dict[str, dict[str, int]] = {s: load_fl_counts(s) for s in SAMPLES}
    cell_by_sample: dict[str, tuple[dict[str, int], int]] = {}
    if not args.no_cells:
        cell_by_sample = {s: load_cell_counts(s) for s in SAMPLES}
    total_fl = {s: sum(fl_by_sample[s].values()) for s in SAMPLES}
    log.info(f"Total FL reads per sample: { {s: f'{v:,}' for s, v in total_fl.items()} }")

    provenance = load_md5_provenance()
    if not provenance:
        db.close()
        sys.exit(1)

    # Every DB isoform with a ccRCC membership (source_id of the ccRCC source), by md5.
    ccrcc_src_id = db.execute(text("SELECT id FROM sources WHERE name = :n"),
                              {"n": SOURCES_BY_KEY["ccrcc"].name}).scalar()
    sql = """
        SELECT DISTINCT i.id, i.gene_id, i.protein_md5
          FROM isoforms i
          JOIN isoform_source_memberships m ON m.isoform_id = i.id AND m.source_id = :src
    """
    params = {"src": ccrcc_src_id}
    if args.gene:
        sql += " JOIN genes g ON g.id = i.gene_id WHERE g.symbol = :sym"
        params["sym"] = args.gene
    isoforms = db.execute(text(sql), params).fetchall()
    log.info(f"Processing {len(isoforms):,} isoforms with a ccRCC membership...")

    # Per (isoform, sample) FL = sum over that isoform's transcripts assembled in that sample.
    iso_fl: dict[int, dict[str, int]] = {}
    iso_cells: dict[int, dict[str, int]] = {}
    gene_fl: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    n_unmatched = 0
    for iso_id, gene_id, md5 in isoforms:
        pairs = provenance.get(md5, [])
        if not pairs:
            n_unmatched += 1
            continue
        per_sample: dict[str, int] = defaultdict(int)
        per_sample_cells: dict[str, int] = defaultdict(int)
        for pbid, sample in pairs:
            fl = fl_by_sample[sample].get(pbid)
            if fl:
                per_sample[sample] += fl
            if sample in cell_by_sample:
                nc = cell_by_sample[sample][0].get(pbid)
                if nc:
                    per_sample_cells[sample] += nc
        if per_sample:
            iso_fl[iso_id] = dict(per_sample)
            for sample, fl in per_sample.items():
                gene_fl[gene_id][sample] += fl
        if per_sample_cells:
            iso_cells[iso_id] = dict(per_sample_cells)

    inserted = 0
    for iso_id, gene_id, _md5 in isoforms:
        if args.force:
            db.query(IsoformExpression).filter(IsoformExpression.isoform_id == iso_id).delete()
        for sample, fl in iso_fl.get(iso_id, {}).items():
            if not args.force:
                if db.query(IsoformExpression.id).filter(
                        IsoformExpression.isoform_id == iso_id,
                        IsoformExpression.sample == sample).first():
                    continue
            cpm = (fl / total_fl[sample] * 1_000_000) if total_fl[sample] > 0 else None
            gene_total = gene_fl[gene_id].get(sample, 0)
            isoform_frac = (fl / gene_total) if gene_total > 0 else None
            n_cells = iso_cells.get(iso_id, {}).get(sample)
            total_cells = cell_by_sample[sample][1] if sample in cell_by_sample else 0
            cell_frac = (n_cells / total_cells) if (n_cells is not None and total_cells > 0) else None
            db.add(IsoformExpression(
                isoform_id=iso_id, sample=sample, condition=SAMPLES[sample],
                fl_count=fl,
                cpm=round(cpm, 4) if cpm is not None else None,
                isoform_frac=round(isoform_frac, 5) if isoform_frac is not None else None,
                n_cells=n_cells,
                cell_frac=round(cell_frac, 5) if cell_frac is not None else None,
            ))
            inserted += 1
            if inserted % 10_000 == 0:
                db.commit()
                log.info(f"  Committed {inserted:,} records so far...")

    db.commit()
    db.close()

    log.info(f"\nDone: {inserted:,} expression records inserted over "
             f"{len(iso_fl):,} isoforms; {n_unmatched:,} isoforms had no ccRCC "
             f"transcript in the metadata; "
             f"{len(isoforms) - len(iso_fl) - n_unmatched:,} had transcripts but FL=0/NA")
    sys.exit(exit_status(len(isoforms), 0 if inserted or not isoforms else len(isoforms)))


if __name__ == "__main__":
    main()
