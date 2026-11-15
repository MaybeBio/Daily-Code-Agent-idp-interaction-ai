#!/usr/bin/env python3
"""Select IDP target crops from the regions outside the accepted structured domains."""

import argparse
import json
import os
import tempfile
import concurrent.futures
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import seq1
from tqdm.auto import tqdm

# ── defaults ──────────────────────────────────────────────────────────────────
DEFAULT_DOMAIN_STATS_PATH = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Statistics" / "domain_stats.csv"
DEFAULT_FILTERED_DOMAIN_STATS_PATH = Path(__file__).resolve().parent.parent / "Data" / "Filtered" / "filtered_domain_stats.csv"
DEFAULT_CROPS_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "PDB_IDP_Crops"
DEFAULT_SETTINGS_DIR = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Settings_IDP"
DEFAULT_IDP_CROPS_STATS_PATH = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "Statistics" / "idp_crops_stats.csv"
DEFAULT_BINDCRAFT_JOBS_FILE = Path(__file__).resolve().parent.parent / "Data" / "AFDB_v6" / "idp_bindcraft_jobs.txt"
DEFAULT_DSSP_BINARY_PATH = Path(__file__).resolve().parent.parent / "BindCraft" / "functions" / "dssp"
DEFAULT_LOCAL_DATA_BASE_DIR = Path(__file__).resolve().parent.parent / "Data"
DEFAULT_REMOTE_DATA_BASE_DIR = Path("/leonardo_work/AIFAC_S03_099/config/scripts/Bindome/Data")

# Crop generation parameters
DEFAULT_MIN_CROP_LENGTH = 30
DEFAULT_MAX_CROP_LENGTH = 300
DEFAULT_CROP_LENGTH_STEP = 20
DEFAULT_STRIDE_FRAC = 0.5
DEFAULT_BINDER_LENGTH = 90
DEFAULT_N_FINAL_DESIGNS = 5

# Low-complexity (Shannon entropy) filter — carried over from
# idp_bindcraft/preprocess_targets.py: repetitive/polyX windows (polyQ, polyS,
# tandem repeats) lack positional identity that BindCraft can hallucinate against.
DEFAULT_LC_ENTROPY_THRESHOLD = 3.0

# Cap on selected crops PER CONTIGUOUS GAP SEGMENT (kept by highest score), applied
# after local-peak selection. None = unlimited (every local score peak survives, per
# the original "overlap allowed, no masking" design). Per-segment (not per-protein)
# so a multi-segment protein's weaker disordered region isn't crowded out by a
# stronger one. Set via --max-crops-per-segment to bound total crop/job volume for a
# full-proteome run.
DEFAULT_MAX_CROPS_PER_SEGMENT = None

# Cap policy: the last N residues at both crop termini must ALL be hydrophilic
# (not in HYDROPHOBIC_AAS). Mirrors idp_bindcraft/preprocess_targets.py's rationale
# (there, N=1): hydrophobic termini get buried by BindCraft's helicity optimization,
# causing exposed-hydrophobic-patch artifacts. N=1 checks only the single terminal
# residue, which doesn't guard against a hydrophobic residue immediately next to it;
# --cap-length lets that buffer be widened.
HYDROPHOBIC_AAS = set('ILVFWMAP')
DEFAULT_CAP_LENGTH = 1

# Boundary-burial check: reject a candidate crop whose start/end sits right next
# to a long run of buried residues within the SAME eligible segment — a local
# safety net for the case where the accepted-domain footprint doesn't fully
# capture a folded stretch (e.g. a domain step 05 rejected for shape/length/
# membrane reasons but which is still a real fold). Deliberately NOT a full fix
# for that case (many rejected PAE-clusters are genuinely low-pLDDT/disordered
# and must stay eligible) — only catches a candidate that would otherwise slice
# straight into a long uninterrupted buried run. RSA-based (raw, unsmoothed),
# independent of domain_stats.csv.
DEFAULT_BOUNDARY_BURIED_RSA = 0.2
DEFAULT_BOUNDARY_BURIED_RUN = 5

# AF-Bind constants (Piovesan et al. 2022)
DEFAULT_AF_BIND_RSA_THRESHOLD = 0.581
DEFAULT_AF_BIND_PLDDT_LOW = 25.0
# Upper pLDDT bound of the AF-Bind MoRF range. idp_bindcraft/preprocess_targets.py
# raised this from 68.8 to 80.0 ("Fix1b: captures near-structured MoRF contacts");
# this script previously still used the stale 68.8 (never actually run in
# production) — brought in line here, overridable via --af-bind-plddt-high.
DEFAULT_AF_BIND_PLDDT_HIGH = 80.0
AF_BIND_RSA_SMOOTH_WIN = 5


# ── argument parsing ──────────────────────────────────────────────────────────
def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stats-csv", type=str, default=str(DEFAULT_DOMAIN_STATS_PATH),
                        help="domain_stats.csv from step 04 (used for full_pdb_path/chain_id/membrane info per protein).")
    parser.add_argument("--filtered-domain-stats-csv", type=str, default=str(DEFAULT_FILTERED_DOMAIN_STATS_PATH),
                        help="filtered_domain_stats.csv from step 05 structured mode — the accepted-domain "
                             "footprint that IDP crops must avoid.")
    parser.add_argument("--crops-dir", type=str, default=str(DEFAULT_CROPS_DIR),
                        help="Output directory for cropped PDB files.")
    parser.add_argument("--settings-dir", type=str, default=str(DEFAULT_SETTINGS_DIR),
                        help="Output directory for BindCraft target settings JSONs.")
    parser.add_argument("--design-root", type=str, required=True,
                        help="Root directory for BindCraft design output (design_path in settings JSON).")
    parser.add_argument("--bindcraft-advanced-settings", type=str, default="",
                        help="Path to BindCraft advanced settings JSON (e.g. strongiface+mv). "
                             "Recorded in crops stats for downstream use.")
    parser.add_argument("--idp-crops-stats-path", type=str, default=str(DEFAULT_IDP_CROPS_STATS_PATH),
                        help="Output CSV path for IDP crop stats (input for step 05 --idp).")
    parser.add_argument("--bindcraft-jobs-file", type=str, default=str(DEFAULT_BINDCRAFT_JOBS_FILE),
                        help="Output text file listing settings JSON paths (one per line).")
    parser.add_argument("--dssp-binary-path", type=str, default=str(DEFAULT_DSSP_BINARY_PATH),
                        help="Path to the DSSP binary.")
    parser.add_argument("--min-crop-length", type=int, default=DEFAULT_MIN_CROP_LENGTH,
                        help="Shortest dynamic IDP crop, after hydrophilic-cap snapping (default: 30).")
    parser.add_argument("--max-crop-length", type=int, default=DEFAULT_MAX_CROP_LENGTH,
                        help="Longest dynamic IDP crop, after hydrophilic-cap snapping (default: 300).")
    parser.add_argument("--crop-length-step", type=int, default=DEFAULT_CROP_LENGTH_STEP,
                        help="Step between candidate crop lengths scanned per gap segment (default: 20).")
    parser.add_argument("--stride-frac", type=float, default=DEFAULT_STRIDE_FRAC,
                        help="Stride as a fraction of each candidate length (default: 0.5).")
    parser.add_argument("--af-bind-plddt-high", type=float, default=DEFAULT_AF_BIND_PLDDT_HIGH,
                        help=f"Upper pLDDT bound of the AF-Bind MoRF range (default: {DEFAULT_AF_BIND_PLDDT_HIGH}). "
                             "Set to 0 or negative to disable the upper bound entirely (only "
                             "pLDDT > AF_BIND_PLDDT_LOW and the RSA gate apply) — recovers more true "
                             "disorder (DisProt-validated: ~87%->93%) but also lets through much more "
                             "confirmed-structured signal (~39%->69%), so mostly relevant if the higher "
                             "recall on genuine disorder outweighs that added false-positive risk.")
    parser.add_argument("--rsa-threshold", type=float, default=DEFAULT_AF_BIND_RSA_THRESHOLD,
                        help=f"RSA (solvent accessibility) threshold of the AF-Bind gate (default: "
                             f"{DEFAULT_AF_BIND_RSA_THRESHOLD}). Higher = stricter (more exposed-only), and "
                             "discriminates confirmed disorder from confirmed-structured-but-low-confidence "
                             "regions much better than pLDDT alone (DisProt-validated: raising to 0.70 costs "
                             "only ~3pp true-disorder recall but cuts confirmed-structured retention by ~10pp "
                             "vs. the 0.581 default).")
    parser.add_argument("--af-bind-plddt-low", type=float, default=DEFAULT_AF_BIND_PLDDT_LOW,
                        help=f"Lower pLDDT bound of the AF-Bind MoRF range (default: {DEFAULT_AF_BIND_PLDDT_LOW}). "
                             "A residue must have pLDDT above this (and below --af-bind-plddt-high, unless "
                             "disabled) to count at all. Raising this excludes the very-lowest-confidence "
                             "residues (fully unstructured tails/linkers) from scoring; not yet validated "
                             "against DisProt the way plddt_high/rsa_threshold were.")
    parser.add_argument("--lc-entropy-threshold", type=float, default=DEFAULT_LC_ENTROPY_THRESHOLD,
                        help="Shannon-entropy threshold below which a crop is flagged low-complexity/repetitive "
                             f"(default: {DEFAULT_LC_ENTROPY_THRESHOLD}).")
    parser.add_argument("--no-lc-filter", action="store_true", default=False,
                        help="Disable the low-complexity filter (allow repetitive/polyX crops).")
    parser.add_argument("--cap-length", type=int, default=DEFAULT_CAP_LENGTH,
                        help="Number of consecutive residues at EACH terminus that must all be hydrophilic. "
                             "Tried first (the strictest/preferred value) for every gap segment "
                             f"(default: {DEFAULT_CAP_LENGTH}, i.e. just the single terminal residue).")
    parser.add_argument("--cap-length-floor", type=int, default=None,
                        help="If a gap segment yields zero crops at --cap-length, retry it at "
                             "--cap-length - 1, then -2, ... down to this floor (inclusive) before giving "
                             "up on that segment. Default: same as --cap-length (no fallback, matches the "
                             "original fixed-cap-length behavior).")
    parser.add_argument("--cap-max-hydrophobic", type=int, default=0,
                        help="Number of hydrophobic residues tolerated ANYWHERE within the --cap-length-"
                             "residue cap window at each terminus (default: 0, i.e. every residue in the "
                             "cap must be hydrophilic — the original all-or-nothing behavior). E.g. "
                             "--cap-length 10 --cap-max-hydrophobic 1 accepts a 10-residue cap with a "
                             "single hydrophobic residue anywhere in it (including the middle), rather "
                             "than requiring an unbroken 10-residue hydrophilic streak.")
    parser.add_argument("--boundary-buried-rsa", type=float, default=DEFAULT_BOUNDARY_BURIED_RSA,
                        help="RSA below which a residue counts as 'buried' for the boundary-cut check "
                             f"(default: {DEFAULT_BOUNDARY_BURIED_RSA}).")
    parser.add_argument("--boundary-buried-run", type=int, default=DEFAULT_BOUNDARY_BURIED_RUN,
                        help="Reject a candidate crop if this many consecutive buried residues "
                             "(RSA < --boundary-buried-rsa) immediately follow its end, or immediately "
                             "precede its start, within the same eligible segment — signals the crop is "
                             "slicing into a well-folded region that wasn't excluded as an accepted "
                             f"structured domain (default: {DEFAULT_BOUNDARY_BURIED_RUN}). Set to 0 to disable.")
    parser.add_argument("--max-crops-per-segment", type=int, default=DEFAULT_MAX_CROPS_PER_SEGMENT,
                        help="Cap on selected crops per contiguous gap segment, kept by highest score, "
                             "applied after local-peak selection (default: unlimited). Applied per-segment "
                             "(not per-protein) so every distinct disordered region of a multi-segment "
                             "protein keeps its own top crops, rather than one high-scoring segment crowding "
                             "out the others. Use this to control total crop/job volume for a full-proteome run "
                             "without changing the underlying overlap policy.")
    parser.add_argument("--keep-nested-crops", action="store_true", default=False,
                        help="By default, a crop that is fully contained within an already-selected "
                             "higher-scoring crop of the same gap segment is dropped (it would likely "
                             "produce a near-duplicate binder for extra compute). Partial overlaps between "
                             "crops are still always allowed — this only removes strict containment. Pass "
                             "this flag to keep nested crops too.")
    parser.add_argument("--binder-length", type=int, default=DEFAULT_BINDER_LENGTH,
                        help="Binder length for BindCraft settings JSON.")
    parser.add_argument("--n-final-designs", type=int, default=DEFAULT_N_FINAL_DESIGNS,
                        help="number_of_final_designs in BindCraft settings JSON.")
    parser.add_argument("--local-data-base-dir", type=str, default=str(DEFAULT_LOCAL_DATA_BASE_DIR),
                        help="Base directory for local paths (for local→remote mapping).")
    parser.add_argument("--remote-data-base-dir", type=str, default=str(DEFAULT_REMOTE_DATA_BASE_DIR),
                        help="Base directory for remote paths.")
    parser.add_argument("--max-workers", type=int, default=None,
                        help="Max parallel workers (default: SLURM, then CPU count, else 8).")
    parser.add_argument("--debug", action=argparse.BooleanOptionalAction, default=False,
                        help="Process only the first 10 proteins.")
    parser.add_argument("--uniprot-list-file", type=str, default=None,
                        help="Path to a text file with one UniProt accession per line — restricts "
                             "processing to exactly these proteins (e.g. a targeted second pass with "
                             "different thresholds for proteins that got 0 crops in a prior run). "
                             "Default: process every protein in --stats-csv.")
    return parser.parse_args()


# ── helpers ───────────────────────────────────────────────────────────────────
def _slurm_cpu_count() -> int | None:
    for key in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE"):
        val = os.environ.get(key)
        if not val:
            continue
        try:
            return max(1, int(str(val).split(",")[0]))
        except ValueError:
            continue
    return None


def _resolve_max_workers(requested: int | None) -> int:
    if requested and requested > 0:
        return requested
    slurm = _slurm_cpu_count()
    if slurm:
        return slurm
    cpu = os.cpu_count()
    return cpu if cpu else 8


def _parse_chopping_ranges(chopping) -> list[tuple[int, int]]:
    """"100-250" or "100-150_200-250" -> [(100,250)] / [(100,150),(200,250)]."""
    ranges = []
    for part in str(chopping).split("_"):
        if "-" not in part:
            continue
        a, b = part.split("-")
        ranges.append((int(a), int(b)))
    return ranges


def _parse_membrane_ranges(raw) -> list[tuple[int, int]]:
    """membrane_embedded_residue_ranges is a JSON string, e.g. "[[100, 130]]"."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return []
    try:
        parsed = json.loads(raw) if isinstance(raw, str) else raw
    except (json.JSONDecodeError, TypeError):
        return []
    return [(int(a), int(b)) for a, b in parsed]


def _in_any_range(resnum: int, ranges: list[tuple[int, int]]) -> bool:
    return any(a <= resnum <= b for a, b in ranges)


def _contiguous_segments(resnums: list[int]) -> list[list[int]]:
    """[5,6,7,10,11] -> [[5,6,7],[10,11]]."""
    if not resnums:
        return []
    resnums = sorted(resnums)
    segments = [[resnums[0]]]
    for rn in resnums[1:]:
        if rn == segments[-1][-1] + 1:
            segments[-1].append(rn)
        else:
            segments.append([rn])
    return segments


class _ResidueRangeSelect(Select):
    def __init__(self, chain_id: str, resnums: set):
        self.chain_id = chain_id
        self.resnums = resnums

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, res):
        hetflag, resseq, _ = res.get_id()
        return hetflag == " " and resseq in self.resnums


# ── per-residue scoring (whole chain) ─────────────────────────────────────────
def _score_all_residues(pdb_path: str, chain_id: str, dssp_bin: str) -> list[dict]:
    """Score every standard residue in the chain: resnum, aa, plddt, rsa.
    (Whole-chain DSSP pass — RSA needs the full structural context; gap-vs-
    domain masking and AF-Bind smoothing happen afterwards, per segment.)"""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    model = struct[0]
    chain = next((c for c in model.get_chains() if c.id == chain_id), None)
    if chain is None:
        return []

    plddt_map = {}
    aa_map = {}
    for res in chain:
        hetflag, resseq, _ = res.get_id()
        if hetflag != " " or "CA" not in res:
            continue
        plddt_map[resseq] = res["CA"].get_bfactor()
        aa_map[resseq] = seq1(res.get_resname())

    dssp = DSSP(model, pdb_path, dssp=dssp_bin)

    records = []
    for (dc, res_id) in dssp.keys():
        if dc != chain_id:
            continue
        hetflag, resseq, _ = res_id
        if hetflag != " ":
            continue
        entry = dssp[(dc, res_id)]
        rsa = entry[3]
        rsa = float(rsa) if rsa is not None else 0.0
        records.append({
            "resnum": resseq,
            "aa": aa_map.get(resseq, "X"),
            "plddt": plddt_map.get(resseq, 50.0),
            "rsa": rsa,
        })

    records.sort(key=lambda r: r["resnum"])
    return records


def _smoothed_rsa(seg_records: list[dict]) -> np.ndarray:
    """AF_BIND_RSA_SMOOTH_WIN-residue (default 5) rolling average of raw RSA,
    computed WITHIN this segment only (never blend across an excluded
    domain/membrane run) — a single incidentally-buried residue surrounded by
    exposed neighbors doesn't pull the local value down on its own."""
    n = len(seg_records)
    rsa_arr = np.array([r["rsa"] for r in seg_records])
    w = min(AF_BIND_RSA_SMOOTH_WIN, n)
    return np.convolve(rsa_arr, np.ones(w) / w, mode="same")


def _compute_af_bind_w(seg_records: list[dict], plddt_low: float, plddt_high: float | None, rsa_threshold: float) -> None:
    """In-place: smooth RSA WITHIN this segment only (never blend across an
    excluded domain/membrane run), then set af_bind_w per record.

    plddt_high=None disables the upper pLDDT bound entirely — a residue only
    needs pLDDT > plddt_low and to clear the RSA gate. Score normalization
    then divides by 100 (max possible pLDDT) instead of plddt_high, since
    there's no longer a boundary value to normalize against."""
    rsa_smooth = _smoothed_rsa(seg_records)
    norm_high = plddt_high if plddt_high is not None else 100.0

    for i, rec in enumerate(seg_records):
        p = rec["plddt"]
        rs = float(rsa_smooth[i])
        in_morf = (plddt_low < p) if plddt_high is None else (plddt_low < p < plddt_high)
        if in_morf and rs > rsa_threshold * 0.7:
            rec["af_bind_w"] = (p / norm_high) * (rs / rsa_threshold)
        else:
            rec["af_bind_w"] = 0.0


def _adjoins_buried_run(seg_records: list[dict], boundary_idx: int, direction: int,
                        rsa_threshold: float, run_length: int) -> bool:
    """True iff >= run_length consecutive residues with rsa < rsa_threshold start
    immediately adjacent to boundary_idx, scanning outward in `direction` (+1 =
    residues after the crop's end index; -1 = residues before its start index).
    Stops as soon as an exposed residue is hit or the segment edge is reached —
    a real domain reads as an uninterrupted buried run, not scattered buried
    residues, so this doesn't fire on ordinary disordered regions with a few
    incidentally-buried residues."""
    if run_length <= 0:
        return False
    n = len(seg_records)
    count = 0
    i = boundary_idx + direction
    while 0 <= i < n and count < run_length:
        if seg_records[i]["rsa"] < rsa_threshold:
            count += 1
            i += direction
        else:
            break
    return count >= run_length


# ── dynamic, hydrophilic-capped window generation ─────────────────────────────
def _hydrophilic_run_starts(is_hydrophilic: list[bool], cap_length: int, max_hydrophobic: int = 0) -> list[bool]:
    """run_ok[i] = True iff is_hydrophilic[i : i + cap_length] contains at most
    max_hydrophobic hydrophobic residues (anywhere in that span, not just at
    the very ends), a sliding-window, not a consecutive-run streak.
    max_hydrophobic=0 (default) reproduces the original "all cap_length
    residues must be hydrophilic" behavior exactly."""
    n = len(is_hydrophilic)
    run_ok = [False] * n
    if cap_length > n:
        return run_ok
    hydrophobic = [0 if h else 1 for h in is_hydrophilic]
    window_sum = sum(hydrophobic[0:cap_length])
    run_ok[0] = window_sum <= max_hydrophobic
    for i in range(1, n - cap_length + 1):
        window_sum += hydrophobic[i + cap_length - 1] - hydrophobic[i - 1]
        run_ok[i] = window_sum <= max_hydrophobic
    return run_ok


def _generate_dynamic_windows(seg_records: list[dict], config: "IdpCropConfig", cap_length: int) -> list[dict]:
    """Scan candidate lengths [min_crop_length, max_crop_length] (step
    crop_length_step) across this gap segment; snap each candidate's start
    forward / end backward so that the first and last cap_length residues
    contain at most config.cap_max_hydrophobic hydrophobic residues (0 =
    ALL must be hydrophilic); keep it if the snapped length is still in
    range and it isn't low-complexity. cap_length is passed explicitly (rather than read
    from config) so a caller can retry a segment at several cap lengths."""
    n = len(seg_records)
    if n < config.min_crop_length:
        return []

    af_bind_w = np.array([r["af_bind_w"] for r in seg_records])
    plddt_arr = np.array([r["plddt"] for r in seg_records])
    is_hydrophilic = [r["aa"] not in HYDROPHOBIC_AAS for r in seg_records]
    run_ok = _hydrophilic_run_starts(is_hydrophilic, cap_length, config.cap_max_hydrophobic)

    windows = []
    seen: set[tuple[int, int]] = set()
    for wlen in range(config.min_crop_length, config.max_crop_length + 1, config.crop_length_step):
        if wlen > n:
            break
        stride = max(1, int(wlen * config.stride_frac))
        for start in range(0, n - wlen + 1, stride):
            end = start + wlen - 1

            # ALL s >= start whose cap window [s, s+cap_length-1] is within the
            # hydrophobic tolerance AND whose own first residue (s itself — the
            # actual crop terminus) is hydrophilic (never counts against the budget).
            # Trying every valid s (not just the nearest one to `start`) matters:
            # a much better-scoring anchor can sit deeper inside the raw window —
            # stopping at the first match would silently miss it.
            valid_s = [s for s in range(start, end - cap_length + 2)
                       if run_ok[s] and is_hydrophilic[s]]

            for s in valid_s:
                # largest e (= e_start + cap_length - 1) whose cap window [e_start, e]
                # is within tolerance AND whose own last residue (e itself — the
                # actual crop terminus) is hydrophilic.
                e_start = end - cap_length + 1
                e = e_start + cap_length - 1
                while e_start >= s and not (run_ok[e_start] and is_hydrophilic[e]):
                    e_start -= 1
                    e = e_start + cap_length - 1
                if e_start < s or s >= e:
                    continue

                length = e - s + 1
                if not (config.min_crop_length <= length <= config.max_crop_length):
                    continue

                # Reject crops that slice into a well-folded-but-unaccepted domain:
                # a long buried run immediately outside either terminus, still within
                # this same eligible segment, means the sliding window stopped in the
                # middle of a real fold rather than at a natural exposed/loopy edge.
                if (_adjoins_buried_run(seg_records, s, -1, config.boundary_buried_rsa, config.boundary_buried_run)
                        or _adjoins_buried_run(seg_records, e, +1, config.boundary_buried_rsa, config.boundary_buried_run)):
                    continue

                key = (s, e)
                if key in seen:
                    continue
                seen.add(key)

                w_plddt = plddt_arr[s:e + 1]
                mean_plddt = float(w_plddt.mean())
                score_af = float(af_bind_w[s:e + 1].mean())
                boundary_penalty = float(w_plddt.std()) / 100.0

                window_seq = [seg_records[i]["aa"] for i in range(s, e + 1)]
                aa_counts: dict[str, int] = {}
                for aa in window_seq:
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1
                freqs = np.array(list(aa_counts.values()), dtype=float) / length
                seq_entropy = float(-np.sum(freqs * np.log2(freqs + 1e-12)))
                low_complexity = seq_entropy < config.lc_entropy_threshold

                if config.exclude_low_complexity and low_complexity:
                    continue

                windows.append({
                    "resnum_start": seg_records[s]["resnum"],
                    "resnum_end": seg_records[e]["resnum"],
                    "window_size": length,
                    "mean_plddt": round(mean_plddt, 2),
                    "mean_af_bind_w": round(score_af, 4),
                    "boundary_std": round(boundary_penalty, 4),
                    "score": round(score_af - 0.5 * boundary_penalty, 4),
                    "seq_entropy": round(seq_entropy, 3),
                    "low_complexity": low_complexity,
                })
    return windows


def _select_local_peaks(candidates: list[dict], max_crops_per_segment: int | None = None,
                        drop_nested: bool = True) -> list[dict]:
    """Sort by score descending, keep everything except exact-duplicate
    (resnum_start, resnum_end) pairs — no exclusion masking, partial overlaps
    allowed. If drop_nested is True (default), also skip a candidate that is
    fully contained within an already-selected (necessarily higher- or equal-
    scoring, since we walk in score order) crop of the same gap segment — a
    strict subset of an already-chosen crop would likely yield a near-duplicate
    binder for extra compute, unlike a staggered partial overlap. If
    max_crops_per_segment is set, additionally cap how many crops survive
    within each contiguous gap segment (still ranked by score) — bounds total
    volume without letting one high-scoring segment crowd out a weaker one
    elsewhere in the same protein. Requires "gap_segment_start"/"gap_segment_end"
    to already be set on each candidate (done by the caller, _process_protein)."""
    positive = [w for w in candidates if w["score"] > 0.0]
    positive.sort(key=lambda w: -w["score"])

    selected = []
    seen: set[tuple[int, int]] = set()
    per_segment_count: dict[tuple[int, int], int] = {}
    per_segment_selected: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for w in positive:
        key = (w["resnum_start"], w["resnum_end"])
        if key in seen:
            continue
        seg_key = (w["gap_segment_start"], w["gap_segment_end"])

        if drop_nested:
            already = per_segment_selected.get(seg_key, [])
            # containment in EITHER direction: candidate inside an already-selected
            # crop, or candidate containing one (score order means whichever of the
            # pair was processed first — necessarily higher- or equal-scoring — wins).
            if any((a <= key[0] and key[1] <= b) or (key[0] <= a and b <= key[1])
                   for a, b in already):
                continue

        if max_crops_per_segment is not None and per_segment_count.get(seg_key, 0) >= max_crops_per_segment:
            continue

        seen.add(key)
        per_segment_count[seg_key] = per_segment_count.get(seg_key, 0) + 1
        per_segment_selected.setdefault(seg_key, []).append(key)
        selected.append(w)
    return selected


# ── PDB writing ───────────────────────────────────────────────────────────────
def _write_crop_pdb(full_pdb: str, chain_id: str,
                    resnum_start: int, resnum_end: int, out_path: str) -> None:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", full_pdb)
    resnums = set(range(resnum_start, resnum_end + 1))
    io = PDBIO()
    io.set_structure(struct)
    io.save(out_path, _ResidueRangeSelect(chain_id, resnums))


def _write_settings_json(crop_name: str, crop_pdb: str, design_root: str,
                         binder_length: int, n_final: int, out_path: str) -> None:
    settings = {
        "design_path": os.path.join(design_root, crop_name, ""),
        "binder_name": crop_name,
        "starting_pdb": str(Path(crop_pdb).resolve()),
        "chains": "A",
        "target_hotspot_residues": None,
        "lengths": [binder_length, binder_length],
        "number_of_final_designs": n_final,
    }
    with open(out_path, "w") as f:
        json.dump(settings, f, indent=4)


# ── per-protein processing ─────────────────────────────────────────────────────
@dataclass(frozen=True)
class IdpCropConfig:
    crops_dir: str
    settings_dir: str
    design_root: str
    advanced_settings: str
    dssp_binary_path: str
    min_crop_length: int
    max_crop_length: int
    crop_length_step: int
    stride_frac: float
    cap_length: int
    cap_length_floor: int
    cap_max_hydrophobic: int
    af_bind_plddt_low: float
    af_bind_plddt_high: float | None
    rsa_threshold: float
    boundary_buried_rsa: float
    boundary_buried_run: int
    lc_entropy_threshold: float
    exclude_low_complexity: bool
    max_crops_per_segment: int | None
    drop_nested_crops: bool
    binder_length: int
    n_final_designs: int
    local_base: str
    remote_base: str

    def to_remote(self, path: str) -> str:
        if self.local_base and self.remote_base:
            return path.replace(self.local_base, self.remote_base)
        return path


def _process_protein(protein_row: dict, config: "IdpCropConfig") -> list[dict]:
    """Process one protein in three steps: (1) X = full sequence minus accepted
    domains and membrane ranges; (2) Y = X filtered down to individually
    accessible residues (raw RSA > config.rsa_threshold); (3) Y is split into
    contiguous gap segments and divided into crops via dynamic hydrophilic-
    capped windows + local-peak selection. Returns a list of crop dicts
    (possibly empty)."""
    uniprot = protein_row["uniprot"]
    full_pdb = protein_row.get("full_pdb_path_local", "")
    if not full_pdb or not os.path.exists(full_pdb):
        return []

    chain_id = str(protein_row.get("chain_id", "A"))
    domain_ranges = protein_row.get("_domain_ranges", [])
    membrane_ranges = _parse_membrane_ranges(protein_row.get("membrane_embedded_residue_ranges"))

    residues = _score_all_residues(full_pdb, chain_id, config.dssp_binary_path)
    if not residues:
        return []
    rec_by_resnum = {r["resnum"]: r for r in residues}

    # Step 1: X = accepted domains + membrane ranges removed.
    x_resnums = [
        r["resnum"] for r in residues
        if not _in_any_range(r["resnum"], domain_ranges) and not _in_any_range(r["resnum"], membrane_ranges)
    ]

    # Step 2: Y = X filtered to accessible residues
    # per-residue RSA > rsa_threshold. A single buried residue splits a segment in two
    eligible_resnums = [rn for rn in x_resnums if rec_by_resnum[rn]["rsa"] > config.rsa_threshold]

    # Step 3: Y split into contiguous gap segments -> crops.
    segments = _contiguous_segments(eligible_resnums)

    all_candidates = []
    for seg_resnums in segments:
        seg_records = [rec_by_resnum[rn] for rn in seg_resnums]
        _compute_af_bind_w(seg_records, config.af_bind_plddt_low, config.af_bind_plddt_high, config.rsa_threshold)

        # Try the strictest (preferred) cap_length first; only if a segment yields
        # NO (score > 0) crop, retry it one residue shorter, down to
        # cap_length_floor. A geometrically-valid-but-zero-score window doesn't
        # count as success — _select_local_peaks would drop it anyway, and
        # stopping the fallback there would wrongly forfeit a lower cap_length
        # that might have found real AF-Bind signal. A segment that finds a
        # useful crop at cap_length=10 never sees 9, 8, ... — the fallback is
        # per-segment, not a blend of all levels.
        candidates: list[dict] = []
        cap_length_used = None
        for cap_len in range(config.cap_length, config.cap_length_floor - 1, -1):
            trial = _generate_dynamic_windows(seg_records, config, cap_length=cap_len)
            if any(w["score"] > 0.0 for w in trial):
                candidates = trial
                cap_length_used = cap_len
                break

        for c in candidates:
            c["gap_segment_start"] = seg_resnums[0]
            c["gap_segment_end"] = seg_resnums[-1]
            c["cap_length_used"] = cap_length_used
        all_candidates.extend(candidates)

    selected = _select_local_peaks(all_candidates, max_crops_per_segment=config.max_crops_per_segment,
                                   drop_nested=config.drop_nested_crops)
    if not selected:
        return []

    pae_json = protein_row.get("json_path_local", "") or ""
    full_pdb_remote = protein_row.get("full_pdb_path", config.to_remote(str(Path(full_pdb).resolve())))

    crops = []
    for crop in selected:
        crop_name = f"{uniprot}_gapcrop{crop['resnum_start']:04d}_{crop['resnum_end']:04d}"
        crop_pdb_local = str(Path(os.path.join(config.crops_dir, f"{crop_name}.pdb")).resolve())
        settings_json_local = str(Path(os.path.join(config.settings_dir, f"{crop_name}_idp.json")).resolve())

        if not os.path.exists(crop_pdb_local):
            _write_crop_pdb(full_pdb, chain_id,
                            crop["resnum_start"], crop["resnum_end"], crop_pdb_local)
        if not os.path.exists(settings_json_local):
            _write_settings_json(crop_name, crop_pdb_local, config.design_root,
                                 config.binder_length, config.n_final_designs, settings_json_local)

        crops.append({
            # Core columns matching domain_stats.csv schema (for step 05 compatibility)
            "uniprot": uniprot,
            "domain_name": crop_name,
            "domain_pdb_path_local": crop_pdb_local,
            "domain_pdb_path": config.to_remote(crop_pdb_local),
            "full_pdb_path_local": str(Path(full_pdb).resolve()),
            "full_pdb_path": full_pdb_remote,
            "json_path_local": pae_json,
            "json_path": config.to_remote(pae_json) if pae_json else "",
            "domain_length": crop["window_size"],
            "chain_id": chain_id,
            # IDP-specific crop columns
            "resnum_start": crop["resnum_start"],
            "resnum_end": crop["resnum_end"],
            "mean_plddt": crop["mean_plddt"],
            "mean_af_bind_w": crop["mean_af_bind_w"],
            "window_score": crop["score"],
            "seq_entropy": crop["seq_entropy"],
            "low_complexity": crop["low_complexity"],
            "bindcraft_settings_json": settings_json_local,
            "bindcraft_advanced_settings": config.advanced_settings,
            # Traceability: which contiguous gap segment this crop came from,
            # and which cap_length that segment actually needed (== --cap-length
            # unless the segment had to fall back toward --cap-length-floor).
            "gap_segment_start": crop["gap_segment_start"],
            "gap_segment_end": crop["gap_segment_end"],
            "cap_length_used": crop["cap_length_used"],
        })

    return crops


# ── main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    args = _parse_args()
    print("Parsed arguments:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    stats_csv = Path(args.stats_csv)
    if not stats_csv.exists():
        raise FileNotFoundError(f"domain_stats.csv not found: {stats_csv}")

    filtered_domain_stats_csv = Path(args.filtered_domain_stats_csv)
    if not filtered_domain_stats_csv.exists():
        raise FileNotFoundError(
            f"filtered_domain_stats.csv not found: {filtered_domain_stats_csv}. "
            "Run 05_filter_domains_and_create_sweep_index.py in structured mode first."
        )

    df = pd.read_csv(stats_csv)
    print(f"\nLoaded {len(df)} domain rows from {stats_csv}")

    # Protein-level metadata: one row per uniprot (full_pdb_path/chain_id/membrane
    # info is identical across all domain rows of the same protein).
    protein_df = df.drop_duplicates(subset="uniprot", keep="first")[
        ["uniprot", "full_pdb_path_local", "full_pdb_path", "chain_id",
         "json_path_local", "membrane_embedded_residue_ranges"]
    ].copy()
    print(f"Proteins in universe: {len(protein_df)}")

    # Accepted structured-domain footprint, from step 05 structured-mode output.
    # Deliberately NOT every domain_stats.csv cluster: many rejected clusters are
    # genuinely low-pLDDT/disordered and must stay eligible for IDP cropping —
    # only what step 05 actually accepted as good structured-binder material is
    # carved out here.
    filtered_df = pd.read_csv(filtered_domain_stats_csv)
    footprint: dict[str, list[tuple[int, int]]] = {}
    for uni, grp in filtered_df.groupby("uniprot"):
        ranges = []
        for chopping in grp["chopping"]:
            ranges.extend(_parse_chopping_ranges(chopping))
        footprint[uni] = ranges
    print(f"Accepted-domain footprint loaded for {len(footprint)} proteins "
          f"(of {len(protein_df)} total — the rest have no accepted domain "
          f"and are fully eligible for gap-crop search over their whole sequence)")

    if args.uniprot_list_file:
        with open(args.uniprot_list_file) as f:
            allowed = {line.strip() for line in f if line.strip()}
        before = len(protein_df)
        protein_df = protein_df[protein_df["uniprot"].isin(allowed)]
        print(f"--uniprot-list-file: restricted {before} -> {len(protein_df)} proteins "
              f"({len(allowed)} requested, {len(allowed) - len(protein_df)} not found in --stats-csv)")

    if args.debug:
        protein_df = protein_df.head(10)
        print(f"[debug] Limiting to {len(protein_df)} proteins.")

    crops_dir = Path(args.crops_dir)
    settings_dir = Path(args.settings_dir)
    crops_dir.mkdir(parents=True, exist_ok=True)
    settings_dir.mkdir(parents=True, exist_ok=True)

    local_base = str(Path(args.local_data_base_dir).resolve())
    remote_base = str(Path(args.remote_data_base_dir).resolve())

    config = IdpCropConfig(
        crops_dir=str(crops_dir),
        settings_dir=str(settings_dir),
        design_root=args.design_root,
        advanced_settings=args.bindcraft_advanced_settings,
        dssp_binary_path=args.dssp_binary_path,
        min_crop_length=args.min_crop_length,
        max_crop_length=args.max_crop_length,
        crop_length_step=args.crop_length_step,
        stride_frac=args.stride_frac,
        cap_length=args.cap_length,
        cap_length_floor=args.cap_length_floor if args.cap_length_floor is not None else args.cap_length,
        cap_max_hydrophobic=args.cap_max_hydrophobic,
        af_bind_plddt_low=args.af_bind_plddt_low,
        af_bind_plddt_high=args.af_bind_plddt_high if args.af_bind_plddt_high > 0 else None,
        rsa_threshold=args.rsa_threshold,
        boundary_buried_rsa=args.boundary_buried_rsa,
        boundary_buried_run=args.boundary_buried_run,
        lc_entropy_threshold=args.lc_entropy_threshold,
        exclude_low_complexity=not args.no_lc_filter,
        max_crops_per_segment=args.max_crops_per_segment,
        drop_nested_crops=not args.keep_nested_crops,
        binder_length=args.binder_length,
        n_final_designs=args.n_final_designs,
        local_base=local_base,
        remote_base=remote_base,
    )

    protein_records = protein_df.to_dict(orient="records")
    for row in protein_records:
        row["_domain_ranges"] = footprint.get(row["uniprot"], [])

    max_workers = _resolve_max_workers(args.max_workers)
    print(f"Processing {len(protein_records)} proteins with up to {max_workers} workers...")

    process_func = partial(_process_protein, config=config)
    all_crops: list[dict] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        for crop_list in tqdm(
            ex.map(process_func, protein_records),
            total=len(protein_records),
            desc="Generating IDP crops",
        ):
            all_crops.extend(crop_list)

    print(f"\nTotal crops generated: {len(all_crops)}")

    # Write idp_crops_stats.csv
    crops_stats_path = Path(args.idp_crops_stats_path)
    crops_stats_path.parent.mkdir(parents=True, exist_ok=True)
    crops_df = pd.DataFrame(all_crops)
    crops_df.to_csv(crops_stats_path, index=False)
    print(f"Wrote IDP crop stats: {crops_stats_path}")

    # Write bindcraft jobs file (one settings JSON path per line)
    jobs_file = Path(args.bindcraft_jobs_file)
    jobs_file.parent.mkdir(parents=True, exist_ok=True)
    with jobs_file.open("w") as f:
        for crop in all_crops:
            f.write(crop["bindcraft_settings_json"] + "\n")
    print(f"Wrote BindCraft jobs file: {jobs_file}  ({len(all_crops)} entries)")

    # Summary
    proteins_with_crops = len({c["uniprot"] for c in all_crops})
    print(f"\nSummary:")
    print(f"  Proteins processed:     {len(protein_records)}")
    print(f"  Proteins with >=1 crop: {proteins_with_crops}")
    print(f"  Total crops:            {len(all_crops)}")
    if all_crops:
        af_vals = [c["mean_af_bind_w"] for c in all_crops]
        print(f"  AF-Bind w mean/max:     {np.mean(af_vals):.3f} / {np.max(af_vals):.3f}")


if __name__ == "__main__":
    main()
