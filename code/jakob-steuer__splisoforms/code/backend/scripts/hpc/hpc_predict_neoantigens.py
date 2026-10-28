#!/usr/bin/env python3
"""
HPC neoantigen binding prediction — self-contained, no database required.

Reads the staged input TSV (from stage_hpc_input.py), runs the three binding
predictors (MHCflurry + NetMHCpan-4.2 + BigMHC) on novel junction peptides,
and writes raw prediction scores to Parquet.

The heavy compute (GPU for BigMHC, licensed binary for NetMHCpan) happens here.
Annotation (pLDDT, PTM, consensus ranking) happens locally in load_neoantigen_results.py.

Required files co-located in the same directory (or on PYTHONPATH):
  - splice_aware_align.py   (pure Python, from backend/app/analysis/)
  - netmhcpan_runner.py     (subprocess wrapper, from backend/app/analysis/)
  - bigmhc_runner.py        (subprocess wrapper, from backend/app/analysis/)

Required env vars:
  - NETMHCPAN_BIN: path to netMHCpan binary
  - BIGMHC_DIR:    path to BigMHC repo root (containing src/ and models/)

Usage:
    python hpc_predict_neoantigens.py \\
        --input hpc_input.tsv \\
        --output predictions.parquet \\
        --workers 4 --batch-size 20

    # Validation batch:
    python hpc_predict_neoantigens.py \\
        --input hpc_input.tsv --output test.parquet --limit 100 --workers 2
"""
import argparse
import contextlib
import csv
import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

# ── Ensure co-located modules are importable ─────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

logging.basicConfig(
    level=logging.WARNING,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("hpc_predict")
log.setLevel(logging.INFO)
for _lib in ("netmhcpan_runner", "bigmhc_runner"):
    logging.getLogger(_lib).setLevel(logging.INFO)

# ── Import pure-Python analysis helpers ──────────────────────────────────────
from splice_aware_align import splice_aware_align

# ── Constants (mirrored from neoantigen.py to avoid app dependency) ──────────
PEPTIDE_LENGTHS = [8, 9, 10, 11]
FLANK_LEN = 10
STRONG_BINDER_THRESHOLD = 50.0
WEAK_BINDER_THRESHOLD = 500.0
PRESENTED_THRESHOLD = 0.5

DEFAULT_HLA_ALLELES = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01",
    "HLA-A*23:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-A*29:02",
    "HLA-A*30:01", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01",
    "HLA-A*68:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*13:02", "HLA-B*15:01",
    "HLA-B*15:02", "HLA-B*15:03", "HLA-B*18:01", "HLA-B*27:05",
    "HLA-B*35:01", "HLA-B*38:01", "HLA-B*40:01", "HLA-B*44:02",
    "HLA-B*44:03", "HLA-B*51:01", "HLA-B*53:01", "HLA-B*58:01",
    "HLA-C*03:04", "HLA-C*04:01", "HLA-C*06:02", "HLA-C*07:01",
    "HLA-C*07:02", "HLA-C*12:03",
]

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ── Pure-function copies from neoantigen.py (no DB, no ORM) ─────────────────

@dataclass
class NovelRegion:
    start: int
    end: int
    region_type: str


def find_novel_regions(novel_sequence: str, canonical_sequence: str) -> List[NovelRegion]:
    if not novel_sequence or not canonical_sequence:
        return []

    can_aligned, nov_aligned = splice_aware_align(canonical_sequence, novel_sequence)
    if not can_aligned or not nov_aligned:
        return []

    regions: List[NovelRegion] = []
    try:
        novel_pos = 0
        current_region_start = None
        current_region_type = None

        for i in range(min(len(can_aligned), len(nov_aligned))):
            can_char = can_aligned[i]
            nov_char = nov_aligned[i]

            if nov_char == "-":
                if current_region_start is not None:
                    regions.append(NovelRegion(
                        start=current_region_start + 1,
                        end=novel_pos,
                        region_type=current_region_type,
                    ))
                    current_region_start = None
                continue

            novel_pos += 1

            if can_char == "-":
                if current_region_start is None:
                    current_region_start = novel_pos - 1
                    current_region_type = "novel_insertion"
            elif can_char != nov_char:
                if current_region_start is None:
                    current_region_start = novel_pos - 1
                    current_region_type = "novel_substitution"
            else:
                if current_region_start is not None:
                    regions.append(NovelRegion(
                        start=current_region_start + 1,
                        end=novel_pos - 1,
                        region_type=current_region_type,
                    ))
                    current_region_start = None

        if current_region_start is not None:
            regions.append(NovelRegion(
                start=current_region_start + 1,
                end=novel_pos,
                region_type=current_region_type,
            ))
    except Exception as e:
        log.error(f"Alignment parsing failed: {e}")

    if len(novel_sequence) < len(canonical_sequence) * 0.7:
        for region in regions:
            if region.end >= len(novel_sequence) * 0.9:
                region.region_type = "frameshift"

    return regions


def extract_junction_peptides(
    sequence: str,
    novel_regions: List[NovelRegion],
) -> List[Tuple[str, int, int, str, str, str]]:
    seq_len = len(sequence)
    peptides: List[Tuple[str, int, int, str, str, str]] = []
    seen: Set[Tuple[str, str]] = set()

    for region in novel_regions:
        for k in PEPTIDE_LENGTHS:
            window_start = max(0, region.start - k)
            window_end = min(seq_len - k, region.end - 1)

            for pos in range(window_start, window_end + 1):
                peptide = sequence[pos: pos + k]
                if len(peptide) < k:
                    continue
                if not all(c in STANDARD_AA for c in peptide):
                    continue

                key = (peptide, region.region_type)
                if key in seen:
                    continue
                seen.add(key)

                n_flank = sequence[max(0, pos - FLANK_LEN): pos]
                c_flank = sequence[pos + k: min(seq_len, pos + k + FLANK_LEN)]

                peptides.append((peptide, pos + 1, pos + k, region.region_type, n_flank, c_flank))

    return peptides


# ── MHCflurry wrapper (same as neoantigen.py predict_presentation) ──────────

@contextlib.contextmanager
def _silence_output():
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    saved_out_fd = os.dup(1)
    saved_err_fd = os.dup(2)
    saved_py_out = sys.stdout
    saved_py_err = sys.stderr
    os.dup2(devnull_fd, 1)
    os.dup2(devnull_fd, 2)
    sys.stdout = open(os.devnull, "w")
    sys.stderr = open(os.devnull, "w")
    try:
        yield
    finally:
        sys.stdout.close()
        sys.stderr.close()
        os.dup2(saved_out_fd, 1)
        os.dup2(saved_err_fd, 2)
        os.close(saved_out_fd)
        os.close(saved_err_fd)
        os.close(devnull_fd)
        sys.stdout = saved_py_out
        sys.stderr = saved_py_err


PredResult = Tuple[float, Optional[float], Optional[float], Optional[float]]

_cached_presentation_predictor = None
_cached_affinity_predictor = None
_mhcflurry_loaded = False


def _get_mhcflurry_predictors():
    global _cached_presentation_predictor, _cached_affinity_predictor, _mhcflurry_loaded
    if _mhcflurry_loaded:
        return _cached_presentation_predictor, _cached_affinity_predictor
    _mhcflurry_loaded = True
    try:
        from mhcflurry import Class1PresentationPredictor
        with _silence_output():
            _cached_presentation_predictor = Class1PresentationPredictor.load()
        log.info("MHCflurry: loaded Class1PresentationPredictor")
    except (ImportError, Exception) as e:
        log.warning(f"Class1PresentationPredictor unavailable ({e}), falling back to affinity-only")
    try:
        from mhcflurry import Class1AffinityPredictor
        with _silence_output():
            _cached_affinity_predictor = Class1AffinityPredictor.load()
    except (ImportError, Exception) as e:
        log.error(f"Class1AffinityPredictor also unavailable: {e}")
    return _cached_presentation_predictor, _cached_affinity_predictor


def predict_presentation(
    peptides: List[str],
    alleles: List[str],
    n_flanks: Dict[str, str],
    c_flanks: Dict[str, str],
) -> Dict[Tuple[str, str], PredResult]:
    if not peptides:
        return {}

    results: Dict[Tuple[str, str], PredResult] = {}
    pep_n_flanks = [n_flanks.get(p, "") for p in peptides]
    pep_c_flanks = [c_flanks.get(p, "") for p in peptides]

    presentation_predictor, affinity_predictor = _get_mhcflurry_predictors()
    if presentation_predictor is None and affinity_predictor is None:
        log.error("MHCflurry not installed. Run: pip install mhcflurry && mhcflurry-downloads fetch")
        return {}

    for allele in alleles:
        if presentation_predictor is not None:
            try:
                _predict_fn = getattr(
                    presentation_predictor, "predict",
                    getattr(presentation_predictor, "predict_to_dataframe", None),
                )
                if _predict_fn is None:
                    raise AttributeError("No predict method found")

                with _silence_output():
                    df = _predict_fn(
                        peptides=peptides,
                        alleles=[allele],
                        n_flanks=pep_n_flanks,
                        c_flanks=pep_c_flanks,
                    )

                for _, row in df.iterrows():
                    pep = row["peptide"]
                    affinity = float(row.get("affinity",
                                  row.get("mhcflurry_affinity",
                                  row.get("prediction", 999999))))
                    percentile = row.get("presentation_percentile",
                                  row.get("affinity_percentile",
                                  row.get("mhcflurry_affinity_percentile",
                                  row.get("prediction_percentile", None))))
                    pres_score = row.get("presentation_score", None)
                    proc_score = row.get("processing_score", None)
                    results[(pep, allele)] = (
                        affinity,
                        float(percentile) if percentile is not None else None,
                        float(pres_score) if pres_score is not None else None,
                        float(proc_score) if proc_score is not None else None,
                    )
                continue
            except Exception as e:
                log.debug(f"Presentation predictor failed for {allele}: {e}")

        if affinity_predictor is not None:
            try:
                with _silence_output():
                    df = affinity_predictor.predict_to_dataframe(
                        peptides=peptides,
                        alleles=[allele] * len(peptides),
                    )
                for _, row in df.iterrows():
                    pep = row["peptide"]
                    affinity = float(row.get("mhcflurry_affinity",
                                              row.get("prediction", 999999)))
                    percentile = row.get("mhcflurry_affinity_percentile",
                                         row.get("prediction_percentile", None))
                    results[(pep, allele)] = (
                        affinity,
                        float(percentile) if percentile is not None else None,
                        None, None,
                    )
            except Exception as e:
                log.warning(f"Affinity prediction also failed for {allele}: {e}")

    return results


def classify_binding(affinity_nm: float) -> str:
    if affinity_nm <= STRONG_BINDER_THRESHOLD:
        return "strong_binder"
    elif affinity_nm <= WEAK_BINDER_THRESHOLD:
        return "weak_binder"
    else:
        return "non_binder"


# ── Per-isoform CPU phase (MHCflurry + NetMHCpan, no BigMHC) ─────────────────

@dataclass
class IsoformPartial:
    isoform_id: int
    transcript_id: str
    unique_peptides: List[str]
    binding_results: Dict[Tuple[str, str], PredResult]
    netmhcpan_results: dict
    peptide_info_map: Dict[str, List[Tuple[int, int, str]]]


def _prepare_isoform(
    isoform_id: int,
    transcript_id: str,
    sequence: str,
    canonical_sequence: str,
    alleles: List[str],
) -> Optional[IsoformPartial]:
    """CPU-only phase: alignment + MHCflurry + NetMHCpan. BigMHC runs separately."""

    novel_regions = find_novel_regions(sequence, canonical_sequence)
    if not novel_regions:
        return None

    peptides_info = extract_junction_peptides(sequence, novel_regions)
    peptides_info = [e for e in peptides_info if e[0] not in canonical_sequence]
    if not peptides_info:
        return None

    n_flanks: Dict[str, str] = {}
    c_flanks: Dict[str, str] = {}
    for pep, _, _, _, n_flank, c_flank in peptides_info:
        if pep not in n_flanks:
            n_flanks[pep] = n_flank
            c_flanks[pep] = c_flank

    unique_peptides = list(n_flanks.keys())

    binding_results = predict_presentation(unique_peptides, alleles, n_flanks, c_flanks)
    if not binding_results:
        return None

    netmhcpan_results: dict = {}
    try:
        from netmhcpan_runner import is_available as _nmhc_ok, run_netmhcpan_multi_allele
        if _nmhc_ok():
            netmhcpan_results = run_netmhcpan_multi_allele(unique_peptides, alleles)
    except Exception as _e:
        log.warning(f"  NetMHCpan FAILED for {transcript_id}: {_e}")

    peptide_info_map: Dict[str, List[Tuple[int, int, str]]] = {}
    for pep, start, end, jtype, _, _ in peptides_info:
        peptide_info_map.setdefault(pep, []).append((start, end, jtype))

    return IsoformPartial(
        isoform_id=isoform_id,
        transcript_id=transcript_id,
        unique_peptides=unique_peptides,
        binding_results=binding_results,
        netmhcpan_results=netmhcpan_results,
        peptide_info_map=peptide_info_map,
    )


def _assemble_rows(partial: IsoformPartial, bigmhc_results: dict) -> List[dict]:
    """Combine CPU-phase results with batched BigMHC results into final rows."""
    rows: List[dict] = []
    seen_keys: Set[Tuple[str, str, int]] = set()

    for (peptide, allele), (affinity, percentile, pres_score, proc_score) in partial.binding_results.items():
        category = classify_binding(affinity)
        is_presented = (pres_score is not None and pres_score >= PRESENTED_THRESHOLD)
        if category == "non_binder" and not is_presented:
            continue

        _nmhc = partial.netmhcpan_results.get((peptide, allele))
        _bmhc = bigmhc_results.get((peptide, allele))

        for start, end, junction_type in partial.peptide_info_map.get(peptide, []):
            dedup_key = (peptide, allele, start)
            if dedup_key in seen_keys:
                continue
            seen_keys.add(dedup_key)

            rows.append({
                "isoform_id": partial.isoform_id,
                "peptide": peptide,
                "peptide_start": start,
                "peptide_end": end,
                "junction_type": junction_type,
                "hla_allele": allele,
                "binding_category": category,
                "mhcflurry_affinity_nm": round(affinity, 2),
                "mhcflurry_percentile": round(percentile, 4) if percentile is not None else None,
                "presentation_score": round(pres_score, 4) if pres_score is not None else None,
                "processing_score": round(proc_score, 4) if proc_score is not None else None,
                "netmhcpan_affinity_nm": round(_nmhc.affinity_nm, 2) if _nmhc else None,
                "netmhcpan_rank_el": round(_nmhc.rank_el, 4) if _nmhc else None,
                "netmhcpan_rank_ba": round(_nmhc.rank_ba, 4) if _nmhc else None,
                "bigmhc_el_score": round(_bmhc.el_score, 4) if (_bmhc and _bmhc.el_score is not None) else None,
                "bigmhc_im_score": round(_bmhc.im_score, 4) if (_bmhc and _bmhc.im_score is not None) else None,
            })

    return rows


# ── Batch processing with batched BigMHC ─────────────────────────────────────

def process_batch(batch: List[tuple], alleles: List[str]) -> Tuple[int, int, int, List[dict], List[str]]:
    """Process a batch: CPU predictors per-isoform, then one BigMHC call for the batch.

    Returns (success, skipped, failed, result_rows, errors).
    """
    success = 0
    skipped = 0
    failed = 0
    all_rows: List[dict] = []
    errors: List[str] = []
    partials: List[IsoformPartial] = []

    for isoform_id, transcript_id, sequence, canonical_sequence in batch:
        try:
            partial = _prepare_isoform(isoform_id, transcript_id, sequence, canonical_sequence, alleles)
            if partial:
                partials.append(partial)
            else:
                skipped += 1
        except Exception as e:
            failed += 1
            errors.append(f"{transcript_id} (id={isoform_id}): {type(e).__name__}: {e}")

    if not partials:
        return success, skipped, failed, all_rows, errors

    all_unique_peptides: Set[str] = set()
    for p in partials:
        all_unique_peptides.update(p.unique_peptides)

    bigmhc_results: dict = {}
    try:
        from bigmhc_runner import is_available as _bmhc_ok, run_bigmhc
        if _bmhc_ok():
            pep_list = list(all_unique_peptides)
            log.info(
                f"  BigMHC batch: {len(pep_list)} unique peptides × {len(alleles)} alleles "
                f"= {len(pep_list) * len(alleles):,} pairs from {len(partials)} isoforms"
            )
            bigmhc_results = run_bigmhc(pep_list, alleles)
    except Exception as _e:
        log.warning(f"  BigMHC FAILED for batch: {_e}")

    for partial in partials:
        try:
            rows = _assemble_rows(partial, bigmhc_results)
            if rows:
                all_rows.extend(rows)
                success += 1
            else:
                skipped += 1
        except Exception as e:
            failed += 1
            errors.append(f"{partial.transcript_id} (id={partial.isoform_id}): assembly: {e}")

    return success, skipped, failed, all_rows, errors


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="HPC neoantigen binding prediction (no DB required)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input", required=True, help="Input TSV from stage_hpc_input.py")
    parser.add_argument("--output", required=True, help="Output Parquet path")
    parser.add_argument("--workers", type=int, default=1,
                        help="Parallel workers within this task (default: 1, recommended for GPU jobs)")
    parser.add_argument("--batch-size", type=int, default=100,
                        help="Isoforms per BigMHC batch (default: 100)")
    parser.add_argument("--checkpoint-every", type=int, default=500,
                        help="Write parquet checkpoint every N isoforms (default: 500)")
    parser.add_argument("--limit", type=int, help="Limit isoforms to process")
    parser.add_argument("--allow-missing-tools", action="store_true",
                        help="Run even if NetMHCpan/BigMHC are unavailable (reduced ensemble; "
                             "the loader will see NULL scores for the missing tool)")
    parser.add_argument("--chunk", type=int, default=None,
                        help="0-based chunk index for SLURM array (auto-read from SLURM_ARRAY_TASK_ID if unset)")
    parser.add_argument("--num-chunks", type=int, default=None,
                        help="Total number of chunks (auto-read from SLURM_ARRAY_TASK_COUNT if unset)")
    args = parser.parse_args()

    # Auto-detect SLURM array environment
    chunk = args.chunk
    num_chunks = args.num_chunks
    if chunk is None and "SLURM_ARRAY_TASK_ID" in os.environ:
        chunk = int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1  # SLURM is 1-based
    if num_chunks is None and "SLURM_ARRAY_TASK_COUNT" in os.environ:
        num_chunks = int(os.environ["SLURM_ARRAY_TASK_COUNT"])
        # Chunks are interleaved by row index (i % num_chunks), so the chunk count
        # MUST equal the array's full width. A sparse resubmission such as
        # `sbatch --array=1,4` reports SLURM_ARRAY_TASK_COUNT=2 and would silently
        # re-chunk the input (2026-10-27: task 1 spent 18 h recomputing other shards'
        # isoforms and never reached its own remainder). Refuse unless --num-chunks
        # pins the original width.
        task_max = int(os.environ.get("SLURM_ARRAY_TASK_MAX", num_chunks))
        task_min = int(os.environ.get("SLURM_ARRAY_TASK_MIN", 1))
        if task_min != 1 or task_max != num_chunks:
            log.error(
                f"Sparse SLURM array (tasks {task_min}..{task_max}, count {num_chunks}): "
                f"cannot infer the chunk width. Pass --num-chunks <original array width> "
                f"(the sbatch for this input pins it)."
            )
            sys.exit(2)

    # Read input
    log.info(f"Reading input from {args.input}")
    records: List[tuple] = []
    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            records.append((
                int(row["isoform_id"]),
                row["transcript_id"],
                row["sequence"],
                row["canonical_sequence"],
            ))

    # Slice to this chunk if running as part of an array
    if chunk is not None and num_chunks is not None:
        records = [r for i, r in enumerate(records) if i % num_chunks == chunk]
        log.info(f"SLURM array: chunk {chunk + 1}/{num_chunks} → {len(records):,} isoforms")

    # Handle resume: skip isoform_ids already processed (with or without results).
    # The parquet only contains isoforms that had binder rows. A companion
    # .done file tracks ALL processed ids so skipped isoforms aren't re-run.
    done_ids: Set[int] = set()
    done_file = args.output.replace(".parquet", ".done")
    if os.path.exists(done_file):
        try:
            with open(done_file) as f:
                done_ids = {int(line.strip()) for line in f if line.strip()}
        except Exception:
            pass
    if os.path.exists(args.output):
        try:
            import pyarrow.parquet as pq
            existing = pq.read_table(args.output, columns=["isoform_id"])
            done_ids.update(existing.column("isoform_id").to_pylist())
        except Exception:
            pass
    if done_ids:
        log.info(f"Resuming: {len(done_ids)} isoform_ids already processed, skipping")

    records = [r for r in records if r[0] not in done_ids]

    if args.limit:
        records = records[:args.limit]

    total = len(records)
    log.info(f"Processing {total:,} isoforms | {args.workers} workers | batch {args.batch_size}")

    if total == 0:
        log.info("Nothing to process.")
        return

    # ── Tool availability check (runs once, up-front) ──────────────────────
    tools_ok = []
    tools_missing = []

    try:
        from mhcflurry import Class1PresentationPredictor
        tools_ok.append("MHCflurry")
    except ImportError:
        tools_missing.append("MHCflurry (pip install mhcflurry && mhcflurry-downloads fetch)")

    try:
        from netmhcpan_runner import is_available as _nmhc_check
        if _nmhc_check():
            tools_ok.append("NetMHCpan")
        else:
            tools_missing.append("NetMHCpan (set NETMHCPAN_BIN or add netMHCpan to PATH)")
    except ImportError:
        tools_missing.append("NetMHCpan (netmhcpan_runner.py not found)")

    try:
        from bigmhc_runner import is_available as _bmhc_check
        if _bmhc_check():
            tools_ok.append("BigMHC")
        else:
            tools_missing.append("BigMHC (set BIGMHC_DIR or place bigmhc at tools/bigmhc/)")
    except ImportError:
        tools_missing.append("BigMHC (bigmhc_runner.py not found)")

    print(f"  Tools available:  {', '.join(tools_ok) or 'NONE'}")
    if tools_missing:
        print(f"  Tools MISSING:    {', '.join(tools_missing)}")
        log.warning(f"Missing tools: {tools_missing}")
    if not tools_ok:
        log.error("No prediction tools available — aborting.")
        sys.exit(1)
    if "MHCflurry" not in [t for t in tools_ok]:
        log.error("MHCflurry is required (primary predictor) — aborting.")
        sys.exit(1)
    # The ranking is a 3-tool consensus. A missing NetMHCpan/BigMHC silently turns
    # it into a 2-of-2 vote (2026-07 run: shards 23–40 ran without NetMHCpan because
    # the wrapper's NMHOME pointed at a moved directory — review NEO-1). Abort unless
    # the caller explicitly accepts a reduced ensemble.
    if tools_missing and not getattr(args, "allow_missing_tools", False):
        log.error("Declared ensemble member(s) missing: %s — aborting. "
                  "Fix the tool (NETMHCPAN_BIN / BIGMHC_DIR) or pass --allow-missing-tools.",
                  tools_missing)
        sys.exit(2)
    print()

    alleles = DEFAULT_HLA_ALLELES
    batches = [records[i:i + args.batch_size] for i in range(0, total, args.batch_size)]

    # ── Checkpointing helpers ─────────────────────────────────────────────
    import pyarrow as pa
    import pyarrow.parquet as pq

    pending_rows: List[dict] = []
    total_rows_written = 0
    if done_ids and os.path.exists(args.output):
        try:
            total_rows_written = pq.read_metadata(args.output).num_rows
        except Exception:
            pass

    pending_done_ids: Set[int] = set()

    def _checkpoint(rows: List[dict], batch_done_ids: Set[int]) -> int:
        """Append rows to parquet and all processed ids to .done file."""
        nonlocal total_rows_written
        if rows:
            if os.path.exists(args.output):
                try:
                    existing = pq.read_table(args.output)
                    combined = pa.concat_tables([existing, pa.Table.from_pylist(rows)])
                except Exception:
                    combined = pa.Table.from_pylist(rows)
            else:
                combined = pa.Table.from_pylist(rows)
            pq.write_table(combined, args.output, compression="zstd")
            total_rows_written = len(combined)
        if batch_done_ids:
            with open(done_file, "a") as f:
                f.writelines(f"{iid}\n" for iid in batch_done_ids)
        return len(rows)

    success_total = skipped_total = failed_total = 0
    all_errors: List[str] = []

    try:
        from tqdm import tqdm
        pbar = tqdm(total=total, desc="Neoantigen", unit="iso")
    except ImportError:
        pbar = None

    t0 = time.time()
    processed = 0
    since_last_checkpoint = 0

    if args.workers <= 1:
        for batch in batches:
            s, sk, f, rows, errs = process_batch(batch, alleles)
            success_total += s
            skipped_total += sk
            failed_total += f
            pending_rows.extend(rows)
            pending_done_ids.update(iid for iid, _, _, _ in batch)
            all_errors.extend(errs)
            processed += len(batch)
            since_last_checkpoint += len(batch)
            if pbar:
                pbar.update(len(batch))
                pbar.set_postfix(ok=success_total, skip=skipped_total, err=failed_total,
                                 rows=total_rows_written + len(pending_rows))
            elif processed % 100 < len(batch):
                elapsed_so_far = time.time() - t0
                rate = processed / elapsed_so_far if elapsed_so_far > 0 else 0
                eta_min = (total - processed) / rate / 60 if rate > 0 else 0
                log.info(
                    f"Progress: {processed:,}/{total:,} ({processed*100//total}%) — "
                    f"{total_rows_written + len(pending_rows):,} rows — "
                    f"{rate:.1f} iso/s — ETA {eta_min:.0f} min"
                )
            if since_last_checkpoint >= args.checkpoint_every:
                n = _checkpoint(pending_rows, pending_done_ids)
                log.info(f"Checkpoint: +{n:,} rows → {args.output} (total {total_rows_written:,})")
                pending_rows.clear()
                pending_done_ids.clear()
                since_last_checkpoint = 0
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(process_batch, batch, alleles): batch for batch in batches}
            for future in as_completed(futures):
                batch = futures[future]
                try:
                    s, sk, f, rows, errs = future.result()
                except Exception as exc:
                    s, sk, f, rows, errs = 0, 0, len(batch), [], [str(exc)]
                success_total += s
                skipped_total += sk
                failed_total += f
                pending_rows.extend(rows)
                pending_done_ids.update(iid for iid, _, _, _ in batch)
                all_errors.extend(errs)
                processed += len(batch)
                since_last_checkpoint += len(batch)
                if pbar:
                    pbar.update(len(batch))
                    pbar.set_postfix(ok=success_total, skip=skipped_total, err=failed_total)
                if since_last_checkpoint >= args.checkpoint_every:
                    n = _checkpoint(pending_rows, pending_done_ids)
                    log.info(f"Checkpoint: +{n:,} rows → {args.output} (total {total_rows_written:,})")
                    pending_rows.clear()
                    pending_done_ids.clear()
                    since_last_checkpoint = 0

    # Final flush
    if pending_rows or pending_done_ids:
        _checkpoint(pending_rows, pending_done_ids)
        pending_rows.clear()
        pending_done_ids.clear()

    if pbar:
        pbar.close()

    elapsed = time.time() - t0
    log.info(f"Wrote {total_rows_written:,} total prediction rows → {args.output}")

    # Summary — read back from parquet for tool coverage stats
    final_rows_count = 0
    has_netmhcpan = 0
    has_bigmhc = 0
    if os.path.exists(args.output):
        try:
            final_table = pq.read_table(args.output)
            final_rows_count = len(final_table)
            if "netmhcpan_rank_el" in final_table.column_names:
                has_netmhcpan = final_table.column("netmhcpan_rank_el").drop_null().length()
            if "bigmhc_el_score" in final_table.column_names:
                has_bigmhc = final_table.column("bigmhc_el_score").drop_null().length()
        except Exception:
            final_rows_count = total_rows_written

    print(f"\n{'═' * 60}")
    print(f"  NEOANTIGEN HPC PREDICTION COMPLETE")
    print(f"{'─' * 60}")
    print(f"  Isoforms processed:     {success_total + skipped_total + failed_total:,}")
    print(f"  With binder rows:       {success_total:,}")
    print(f"  Skipped (no binders):   {skipped_total:,}")
    print(f"  Errors:                 {failed_total:,}")
    print(f"  Prediction rows:        {final_rows_count:,}")
    print(f"{'─' * 60}")
    print(f"  Tool coverage in output rows:")
    print(f"    MHCflurry:   {final_rows_count:,}/{final_rows_count:,} (always)")
    print(f"    NetMHCpan:   {has_netmhcpan:,}/{final_rows_count:,}")
    print(f"    BigMHC:      {has_bigmhc:,}/{final_rows_count:,}")
    if final_rows_count > 0 and (has_netmhcpan == 0 or has_bigmhc == 0):
        print(f"  WARNING: A tool produced 0 rows — check availability above")
    print(f"{'─' * 60}")
    print(f"  Wall time:              {elapsed/60:.1f} min")
    print(f"{'═' * 60}")

    if all_errors:
        print(f"\n  First {min(20, len(all_errors))} errors:")
        for e in all_errors[:20]:
            print(f"    {e}")


if __name__ == "__main__":
    main()
