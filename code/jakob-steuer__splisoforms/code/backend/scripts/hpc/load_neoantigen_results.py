#!/usr/bin/env python3
"""
Load HPC neoantigen prediction results into the local database.

Reads the Parquet output from hpc_predict_neoantigens.py (raw binding scores),
joins with the local DB to compute annotation columns (pLDDT, PTM impact,
consensus rank, agreement tier), and bulk-inserts NeoantigenPrediction rows.

Usage (from project root):
    python backend/scripts/hpc/load_neoantigen_results.py \\
        --input predictions.parquet

    # Dry run (no DB writes):
    python backend/scripts/hpc/load_neoantigen_results.py \\
        --input predictions.parquet --dry-run
"""
import argparse
import json
import logging
import os
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from sqlalchemy import func
from app.database import SessionLocal
from app.models.bio import Gene, Isoform, IsoformExpression, NeoantigenPrediction, TranscriptMetadata
from app.analysis.ptm_data import get_ptm_sites, PTMSite
from app.analysis import neoantigen_tiers as tiers

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("load_neoantigen")

# ── Constants (must match neoantigen.py) ─────────────────────────────────────

STRONG_BINDER_THRESHOLD = 50.0
WEAK_BINDER_THRESHOLD = 500.0
PRESENTED_THRESHOLD = 0.5

_BA_RESCUE_MIN_CYS = 2
_BA_RESCUE_RANK_THRESH = 1.0
_EL_LOOSER_MULTIPLIER = 2.0
_MHCFLURRY_STRONG_RANK = 0.5
_NETMHCPAN_STRONG_RANK = 0.5
_BIGMHC_STRONG_SCORE = 0.5


# ── Annotation functions (mirrored from neoantigen.py) ───────────────────────

def _compute_consensus_rank(
    mhcflurry_percentile: Optional[float],
    netmhcpan_rank_el: Optional[float],
    bigmhc_el_score: Optional[float],
) -> Optional[float]:
    ranks: List[float] = []
    if mhcflurry_percentile is not None:
        ranks.append(mhcflurry_percentile)
    if netmhcpan_rank_el is not None:
        ranks.append(netmhcpan_rank_el)
    if bigmhc_el_score is not None:
        ranks.append(max(0.0, min(100.0, (1.0 - bigmhc_el_score) * 100.0)))
    if not ranks:
        return None
    return round(statistics.median(ranks), 3)


def _flag_ba_rescue(
    peptide: str,
    netmhcpan_rank_el: Optional[float],
    netmhcpan_rank_ba: Optional[float],
) -> bool:
    if netmhcpan_rank_ba is None or netmhcpan_rank_ba >= _BA_RESCUE_RANK_THRESH:
        return False
    if peptide.count("C") < _BA_RESCUE_MIN_CYS:
        return False
    if netmhcpan_rank_el is None:
        return False
    return netmhcpan_rank_el > _EL_LOOSER_MULTIPLIER * netmhcpan_rank_ba


def _compute_agreement(
    mhcflurry_percentile: Optional[float],
    netmhcpan_rank_el: Optional[float],
    bigmhc_el_score: Optional[float],
) -> Tuple[int, str]:
    n = 0
    if mhcflurry_percentile is not None and mhcflurry_percentile < _MHCFLURRY_STRONG_RANK:
        n += 1
    if netmhcpan_rank_el is not None and netmhcpan_rank_el < _NETMHCPAN_STRONG_RANK:
        n += 1
    if bigmhc_el_score is not None and bigmhc_el_score > _BIGMHC_STRONG_SCORE:
        n += 1
    tier = "high" if n >= 3 else "medium" if n == 2 else "low"
    return n, tier


# ── DB prefetch helpers ──────────────────────────────────────────────────────

def _build_canonical_map(db) -> Dict[int, int]:
    """gene_id → canonical isoform_id."""
    genes = db.query(Gene.id, Gene.canonical_transcript_id).filter(Gene.canonical_transcript_id.isnot(None)).all()
    return {g.id: g.canonical_transcript_id for g in genes}


def _build_uniprot_map(db) -> Dict[int, str]:
    """gene_id → uniprot_id, using Gene.uniprot_id with TranscriptMetadata fallback."""
    result: Dict[int, str] = {}

    for gene_id, uid in db.query(Gene.id, Gene.uniprot_id).filter(Gene.uniprot_id.isnot(None)):
        result[gene_id] = uid

    fallback = (
        db.query(Isoform.gene_id, TranscriptMetadata.uniprot_id)
        .join(TranscriptMetadata, Isoform.id == TranscriptMetadata.isoform_id)
        .filter(
            TranscriptMetadata.uniprot_id.isnot(None),
            ~Isoform.gene_id.in_(list(result.keys())) if result else True,
        )
        .distinct()
        .all()
    )
    for gene_id, uid in fallback:
        if gene_id not in result:
            result[gene_id] = uid

    return result


def _build_structural_metrics_map(db, isoform_ids: set) -> Dict[int, tuple]:
    """isoform_id → (plddt_series, pae_blocks_parsed).

    Returns only isoforms that have structural_metrics populated.
    """
    from app.models.bio import StructuralMetrics
    result: Dict[int, tuple] = {}

    batch_size = 5000
    id_list = sorted(isoform_ids)
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        rows = (
            db.query(StructuralMetrics)
            .filter(StructuralMetrics.isoform_id.in_(batch))
            .all()
        )
        for sm in rows:
            plddt = sm.plddt_series or []
            pae_blocks = []
            if sm.pae_blocks:
                try:
                    pae_blocks = json.loads(sm.pae_blocks)
                except Exception:
                    pass
            result[sm.isoform_id] = (plddt, pae_blocks)

    return result


def _build_gene_id_map(db, isoform_ids: set) -> Dict[int, int]:
    """isoform_id → gene_id."""
    result: Dict[int, int] = {}
    batch_size = 5000
    id_list = sorted(isoform_ids)
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        rows = db.query(Isoform.id, Isoform.gene_id).filter(Isoform.id.in_(batch)).all()
        for iso_id, gene_id in rows:
            result[iso_id] = gene_id
    return result


def _build_expressed_set(db, isoform_ids: set) -> set:
    """isoform_ids detected in expression data (any sample with fl_count > 0).
    Decorating flag only (ranking v2) — expression is currently ccRCC long-read."""
    result = set()
    id_list = sorted(isoform_ids)
    for i in range(0, len(id_list), 5000):
        batch = id_list[i:i + 5000]
        rows = (db.query(IsoformExpression.isoform_id)
                .filter(IsoformExpression.isoform_id.in_(batch),
                        IsoformExpression.fl_count > 0)
                .distinct().all())
        result.update(r[0] for r in rows)
    return result


def _build_sequence_map(db, isoform_ids: set) -> Dict[int, str]:
    """isoform_id → sequence."""
    result: Dict[int, str] = {}
    batch_size = 5000
    id_list = sorted(isoform_ids)
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        rows = db.query(Isoform.id, Isoform.sequence).filter(Isoform.id.in_(batch)).all()
        for iso_id, seq in rows:
            result[iso_id] = seq or ""
    return result


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Load HPC neoantigen results into local DB")
    parser.add_argument("--input", required=True, nargs="+",
                        help="Parquet file(s) from hpc_predict_neoantigens.py (accepts glob)")
    parser.add_argument("--dry-run", action="store_true", help="Compute annotations but skip DB writes")
    parser.add_argument("--batch-size", type=int, default=10000, help="DB insert batch size")
    parser.add_argument("--skip-ptm", action="store_true",
                        help="Skip the per-accession live UniProt/EBI PTM lookup that fills the "
                             "decorating ptm_impact_status / is_ptm_sensitive fields (~2 network "
                             "calls per unique accession). priority_tier and all presentation/"
                             "immunogenicity fields are unaffected (they come from the parquet). "
                             "Use this for a fast load, then decorate PTM separately with "
                             "backfill_neoantigen_ptm.py.")
    parser.add_argument("--replace", action="store_true",
                        help="Incremental batch: delete existing predictions for the input's "
                             "isoform_ids and reload them (instead of aborting on a non-empty table). "
                             "Safe re-run for a second/Nth HPC batch.")
    parser.add_argument("--id-map", type=str, default=None,
                        help="TSV with columns old_isoform_id, new_isoform_id (header optional) "
                             "produced by build_isoform_id_map.py. Parquet rows are keyed by the "
                             "numeric isoform_id of the DB the HPC input was staged from; after a "
                             "rebuild those ids are renumbered. Rows whose old id is absent from "
                             "the map are dropped (the protein no longer exists, is now the gene "
                             "canonical, or its canonical changed → must be re-predicted).")
    args = parser.parse_args()

    import glob
    import pyarrow.parquet as pq

    # Expand globs and read all input files
    input_files = []
    for pattern in args.input:
        input_files.extend(sorted(glob.glob(pattern)))
    if not input_files:
        log.error(f"No files matched: {args.input}")
        sys.exit(1)

    rows: list = []
    for f in input_files:
        table = pq.read_table(f)
        rows.extend(table.to_pylist())
        log.info(f"  {f}: {table.num_rows:,} rows")
    log.info(f"Read {len(rows):,} total prediction rows from {len(input_files)} file(s)")

    # Cross-file dedup: a resumed/re-chunked shard can recompute isoforms that another
    # shard already holds (2026-10-27: a sparse `--array=1,4` resubmission). The
    # delete-then-insert below works on the union of isoform_ids, so duplicate rows
    # would be inserted twice. Keep the first occurrence per natural key.
    seen: set = set()
    deduped: list = []
    for r in rows:
        key = (r["isoform_id"], r["peptide"], r["peptide_start"], r["hla_allele"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(r)
    if len(deduped) != len(rows):
        log.warning(f"Dropped {len(rows) - len(deduped):,} duplicate rows "
                    f"(same isoform/peptide/start/allele across input files)")
    rows = deduped
    del seen

    if not rows:
        log.info("No rows to load.")
        return

    if args.id_map:
        id_map: Dict[int, int] = {}
        with open(args.id_map) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2 or not parts[0].strip().lstrip("-").isdigit():
                    continue
                id_map[int(parts[0])] = int(parts[1])
        before = len(rows)
        old_ids = {r["isoform_id"] for r in rows}
        kept = []
        for r in rows:
            new_id = id_map.get(r["isoform_id"])
            if new_id is None:
                continue
            r["isoform_id"] = new_id
            kept.append(r)
        rows = kept
        log.info(f"--id-map {args.id_map}: {len(id_map):,} mappings; "
                 f"{before - len(rows):,} rows dropped (old isoform_id not mapped: "
                 f"{len(old_ids - set(id_map)):,} of {len(old_ids):,} ids); {len(rows):,} rows kept")
        if not rows:
            log.info("No rows left after id mapping.")
            return

    # Collect unique isoform_ids
    isoform_ids = {r["isoform_id"] for r in rows}
    log.info(f"Unique isoform_ids: {len(isoform_ids):,}")

    db = SessionLocal()

    # Check for existing predictions (a dry run never aborts or deletes).
    existing_count = db.query(func.count(NeoantigenPrediction.id)).scalar()
    if existing_count > 0 and not args.dry_run:
        if not args.replace:
            log.warning(f"neoantigen_predictions already has {existing_count:,} rows!")
            log.warning("Pass --replace to reload only this batch's isoform_ids, "
                        "or TRUNCATE the table for a full fresh load. Aborting.")
            db.close()
            sys.exit(1)
        # Incremental batch: delete existing rows for exactly the isoform_ids in this
        # input, then re-insert them. Idempotent; leaves all other isoforms untouched.
        deleted = (
            db.query(NeoantigenPrediction)
            .filter(NeoantigenPrediction.isoform_id.in_(list(isoform_ids)))
            .delete(synchronize_session=False)
        )
        db.commit()
        log.info(f"--replace: deleted {deleted:,} existing rows for {len(isoform_ids):,} "
                 f"isoform_ids before reload")

    # Prefetch all annotation data
    log.info("Prefetching annotation data from DB …")
    gene_id_map = _build_gene_id_map(db, isoform_ids)
    canonical_map = _build_canonical_map(db)
    uniprot_map = _build_uniprot_map(db)
    structural_map = _build_structural_metrics_map(db, isoform_ids)
    sequence_map = _build_sequence_map(db, isoform_ids)
    expressed_set = _build_expressed_set(db, isoform_ids)  # decorating flag (ranking v2)

    log.info(f"  gene_id_map:      {len(gene_id_map):,} isoforms")
    log.info(f"  canonical_map:    {len(canonical_map):,} genes with canonical")
    log.info(f"  uniprot_map:      {len(uniprot_map):,} genes with uniprot_id")
    log.info(f"  structural_map:   {len(structural_map):,} isoforms with structural_metrics")

    # PTM cache: uniprot_id → List[PTMSite]
    ptm_cache: Dict[str, List[PTMSite]] = {}

    # Annotate and build NeoantigenPrediction objects
    if args.skip_ptm:
        log.info("--skip-ptm: NOT fetching canonical PTM sites from UniProt; "
                 "ptm_impact_status defaults to 'Standard' (decorate later with backfill_neoantigen_ptm.py)")
    log.info("Annotating predictions …")
    records: List[NeoantigenPrediction] = []
    stats = defaultdict(int)

    for row in rows:
        iso_id = row["isoform_id"]
        gene_id = gene_id_map.get(iso_id)
        if gene_id is None:
            stats["missing_gene_id"] += 1
            continue

        # Canonical isoform ID
        canonical_id = canonical_map.get(gene_id)

        # Consensus rank
        consensus = _compute_consensus_rank(
            row.get("mhcflurry_percentile"),
            row.get("netmhcpan_rank_el"),
            row.get("bigmhc_el_score"),
        )

        # Agreement
        tools_agreeing, agreement_tier = _compute_agreement(
            row.get("mhcflurry_percentile"),
            row.get("netmhcpan_rank_el"),
            row.get("bigmhc_el_score"),
        )

        # BA rescue
        ba_rescue = _flag_ba_rescue(
            row["peptide"],
            row.get("netmhcpan_rank_el"),
            row.get("netmhcpan_rank_ba"),
        )

        # ── ranking v2: presentation tier × immunogenicity → priority_tier ──
        v2 = tiers.classify(
            mhcflurry_percentile=row.get("mhcflurry_percentile"),
            netmhcpan_rank_el=row.get("netmhcpan_rank_el"),
            bigmhc_el_score=row.get("bigmhc_el_score"),
            bigmhc_im_score=row.get("bigmhc_im_score"),
        )

        # Structural confidence (pLDDT window)
        start = row["peptide_start"]
        end = row["peptide_end"]
        plddt_series, pae_blocks = structural_map.get(iso_id, ([], []))

        window_plddt = plddt_series[start - 1: end] if plddt_series and len(plddt_series) >= end else []
        avg_plddt = sum(window_plddt) / len(window_plddt) if window_plddt else None
        min_plddt = min(window_plddt) if window_plddt else None

        is_rigid = any(
            block.get("start", 0) <= start and end <= block.get("end", 0)
            for block in pae_blocks
        )

        # PTM impact
        impact_status = "Standard"
        is_sensitive = False
        ptm_canonical_residue_lost = False

        uniprot_id = uniprot_map.get(gene_id)
        if uniprot_id and not args.skip_ptm:
            if uniprot_id not in ptm_cache:
                ptm_cache[uniprot_id] = get_ptm_sites(uniprot_id)
            canonical_ptms = ptm_cache[uniprot_id]

            isoform_seq = sequence_map.get(iso_id, "")

            if canonical_ptms:
                direct_ptms = [p for p in canonical_ptms if start <= p.residue_position <= end]
                flanking_ptms = [
                    p for p in canonical_ptms
                    if (start - 3 <= p.residue_position < start)
                    or (end < p.residue_position <= end + 3)
                ]

                if direct_ptms:
                    any_lost = False
                    for p in direct_ptms:
                        if 1 <= p.residue_position <= len(isoform_seq):
                            novel_aa = isoform_seq[p.residue_position - 1]
                            if novel_aa != p.residue_type:
                                any_lost = True
                        else:
                            any_lost = True
                    impact_status = (
                        "Direct PTM-Neoantigen (residue lost)"
                        if any_lost else
                        "Direct PTM-Neoantigen (residue retained)"
                    )
                    is_sensitive = True
                    ptm_canonical_residue_lost = any_lost
                elif flanking_ptms:
                    found_flanking = False
                    lost_flanking = False
                    for p in flanking_ptms:
                        if 1 <= p.residue_position <= len(isoform_seq):
                            if isoform_seq[p.residue_position - 1] == p.residue_type:
                                found_flanking = True
                            else:
                                lost_flanking = True
                        else:
                            lost_flanking = True

                    if found_flanking:
                        impact_status = "Processing-Blocked"
                        is_sensitive = True
                    elif lost_flanking:
                        # ranking v2: keep the mechanistic label only. The old code
                        # baked pLDDT/rigidity into the status ("Gold Medal"); that is
                        # now the orthogonal priority_tier + high-fidelity-fold flag.
                        impact_status = "Splice-Enabled"
                        is_sensitive = True

        record = NeoantigenPrediction(
            isoform_id=iso_id,
            peptide_sequence=row["peptide"],
            peptide_start=start,
            peptide_end=end,
            hla_allele=row["hla_allele"],
            binding_affinity_nm=row["mhcflurry_affinity_nm"],
            binding_category=row["binding_category"],
            percentile_rank=row.get("mhcflurry_percentile"),
            junction_type=row["junction_type"],
            canonical_isoform_id=canonical_id,
            presentation_score=row.get("presentation_score"),
            processing_score=row.get("processing_score"),
            ptm_impact_status=impact_status,
            avg_plddt=round(avg_plddt, 2) if avg_plddt is not None else None,
            min_plddt=round(min_plddt, 2) if min_plddt is not None else None,
            is_ptm_sensitive=is_sensitive,
            netmhcpan_affinity_nm=row.get("netmhcpan_affinity_nm"),
            netmhcpan_rank_el=row.get("netmhcpan_rank_el"),
            netmhcpan_rank_ba=row.get("netmhcpan_rank_ba"),
            bigmhc_el_score=row.get("bigmhc_el_score"),
            bigmhc_im_score=row.get("bigmhc_im_score"),
            consensus_rank=consensus,
            is_ba_rescue=ba_rescue,
            tools_agreeing=tools_agreeing,
            agreement_tier=agreement_tier,
            ptm_canonical_residue_lost=ptm_canonical_residue_lost,
            presentation_tier=v2["presentation_tier"],
            priority_tier=v2["priority_tier"],
            is_immunogenic=v2["is_immunogenic"],
            is_expressed=iso_id in expressed_set,
        )
        records.append(record)
        stats["annotated"] += 1

    log.info(f"Annotated {stats['annotated']:,} rows ({stats.get('missing_gene_id', 0)} skipped: missing gene_id)")

    if args.dry_run:
        log.info("Dry run — skipping DB writes.")
        _print_summary(records)
        db.close()
        return

    # Bulk insert
    log.info(f"Inserting {len(records):,} rows into neoantigen_predictions …")
    for i in range(0, len(records), args.batch_size):
        batch = records[i:i + args.batch_size]
        db.bulk_save_objects(batch)
        db.flush()
        log.info(f"  inserted {min(i + args.batch_size, len(records)):,}/{len(records):,}")

    db.commit()
    db.close()

    _print_summary(records)


def _print_summary(records: List[NeoantigenPrediction]):
    total = len(records)
    strong = sum(1 for r in records if r.binding_category == "strong_binder")
    presented = sum(1 for r in records if r.presentation_score is not None and r.presentation_score >= PRESENTED_THRESHOLD)
    with_plddt = sum(1 for r in records if r.avg_plddt is not None)
    ptm_sensitive = sum(1 for r in records if r.is_ptm_sensitive)
    gold = sum(1 for r in records if r.ptm_impact_status == "Splice-Enabled (Gold Medal)")
    high_agree = sum(1 for r in records if r.agreement_tier == "high")

    print(f"\n{'═' * 60}")
    print(f"  NEOANTIGEN LOAD COMPLETE")
    print(f"{'─' * 60}")
    print(f"  Total rows:              {total:,}")
    print(f"  Strong binders:          {strong:,}")
    print(f"  Presented (score≥0.5):   {presented:,}")
    print(f"  With pLDDT annotation:   {with_plddt:,}")
    print(f"  PTM-sensitive:           {ptm_sensitive:,}")
    print(f"  Gold Medal:              {gold:,}")
    print(f"  High agreement (3/3):    {high_agree:,}")
    print(f"{'═' * 60}")


if __name__ == "__main__":
    main()
