"""Assemble PSSM, chemistry, accessibility, and optional regional context arrays."""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
                      

import argparse
import json
import re
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Optional, Sequence, Tuple
import numpy as np
import pandas as pd
from scipy.stats import gamma

                 
from pssm_conservation_idppi import find_pssm_file, load_pssm_with_extras, parse_region_string
from chemical_window_idppi_fast import (
    KD_UVERSKY_0_1,
    PPII_HILSER,
    Segment,
    build_nardini_lookup,
    center_flank_masks,
    features_for_idr,
    idr_terminal_flags,
    load_categories,
    load_regexes,
    load_regions as load_chemical_regions,
    local_patch_strength_for_all_windows,
    regex_features_for_idr,
)
from raw_feature_pipeline import (
    IDR_KEY_COLUMNS,
    NARDINI_FEATURE_NAMES,
    build_idr_manifest,
    default_esm_python,
    default_nardini_python,
    prepare_derived_features,
    stable_frame_hash,
)
from esm_region_pipeline import REGION_KEY_COLUMNS, build_esm_region_manifest
from build_contextual_esm2 import prepare_contextual_esm2
from build_motifae_sparse_windows import build_motifae_sparse_windows

PROJECT_DIR = Path(os.environ.get("PRIORSLIM_DATA_ROOT", Path(__file__).resolve().parents[2] / "data"))
WINDOWS_FILE: Optional[Path] = None
SPLIT_FILES = tuple(PROJECT_DIR / f"new_data/{split}.tsv" for split in ("train", "val", "test"))
PSSM_DIR = PROJECT_DIR / "sample_pssm50"
FRAGMENT_FILE = PROJECT_DIR / "regions.xlsx"
CACHE_DIR = PROJECT_DIR / "feature_cache"
ESM_REGION_METADATA_FILE = CACHE_DIR / "esm2_region_metadata.tsv"
NARDINI_CACHE_FILE = CACHE_DIR / "nardini_features.tsv"
GUSSIDR_CACHE_FILE = CACHE_DIR / "gussidr_ensemble_features.tsv"
GUSSIDR_DIR = PROJECT_DIR / "GussIDR_zsq"
CHEMICAL_CATEGORIES_JSON = Path(__file__).resolve().parent / "chemical_categories_template.json"
MOTIF_REGEX_JSON = Path(__file__).resolve().parent / "motif_regex_template.json"
OUT_PREFIX = PROJECT_DIR / "train/combined"

WINDOW_LEN = 29
LAMBDA_DECAY = 10.0               
LAMBDA0 = 20.0                    
DEFAULT_V = 0.5
EPS = 1e-8
PATCH = 5
TOP_PATCH_K = 3
CENTER_LEN = 15
ALPHA = 0.5
MAD_FLOOR = 0.05
CLIP_REL = 5.0
NULL_SHUFFLE = 20000
SEED = 13
SKIP_BAD_WINDOWS = True
ALLOW_MISSING_IDR = False
PSSM_TAU = 1.0
PSSM_RESIDUE_FEATURE_NAMES = (
    "pssm_information",
    "log1p_real_match_weight",
    "distribution_specificity",
    "observed_specificity",
    "consensus_conservation",
    "consensus_percentile_IDR",
)
PSSM_WINDOW_FEATURE_NAMES = ("local_conservation_peak",)
ORDERED_COLUMN_CANDIDATES = (
    "Final_Non_Disordered_Intervals",
    "Final_Ordered_Intervals",
    "Ordered_Regions",
)
DISORDERED_COLUMN_CANDIDATES = (
    "Final_Disordered_Intervals",
    "Disorder_Regions",
)
VARIANT_RE = re.compile(r"^(?P<base>.+)_\d+_\d+_[VnMFC]\d+$")

                                                                       
CONTEXT_DROP_INDICES = (1343, 1345, 1347, 1349, 1351, 1353)
ENSEMBLE_DROP_INDICES = tuple(index - (1280 + 2 + 36) for index in CONTEXT_DROP_INDICES)
ENSEMBLE_KEEP_INDICES = tuple(
    index for index in range(38) if index not in ENSEMBLE_DROP_INDICES
)

def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)

def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)

def load_ensemble_frame(path: str | Path) -> pd.DataFrame:
    pass
    frame = pd.read_csv(path, sep="\t")
    expected_columns = 3 + 38
    named = ["uniprot", "idr_start", "idr_end"] + [
        f"feat_{index}" for index in range(1, 39)
    ]
    if frame.columns.tolist() != named:
        frame = pd.read_csv(path, sep="\t", header=None)
    if frame.shape[1] != expected_columns:
        raise ValueError(
            f"Ensemble table must have {expected_columns} columns "
            f"(uniprot, idr_start, idr_end, 38 features); got {frame.shape[1]}"
        )
    frame.columns = named
    frame["uniprot"] = frame["uniprot"].astype(str).str.strip()
    frame["idr_start"] = pd.to_numeric(frame["idr_start"], errors="raise").astype(int)
    frame["idr_end"] = pd.to_numeric(frame["idr_end"], errors="raise").astype(int)
    for column in [f"feat_{index}" for index in range(1, 39)]:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    return frame

def first_existing(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    column_set = set(columns)
    for cand in candidates:
        if cand in column_set:
            return cand
    return None

def base_accession(uniprot: str) -> str:
    match = VARIANT_RE.match(str(uniprot))
    return match.group("base") if match else str(uniprot)

def normalize_windows(df: pd.DataFrame, default_label: Optional[str] = None,
                      start_col: Optional[str] = None, end_col: Optional[str] = None) -> pd.DataFrame:
    uid_col = first_existing(df.columns, ["uniprot", "UniProt_ID", "Protein_Accession", "Protein_Acc", "Primary_Acc", "accession"])
    if uid_col is None:
        raise ValueError(f"Cannot find UniProt column. Columns: {df.columns.tolist()}")
    start_col = start_col or first_existing(df.columns, ["window_start", "win_start", "seg_start", "start"])
    end_col = end_col or first_existing(df.columns, ["window_end", "win_end", "seg_end", "end"])
    if start_col is None or end_col is None:
        raise ValueError("Cannot find window coordinate columns.")
    label_col = first_existing(df.columns, ["label", "class", "y"])
    window_id_col = first_existing(df.columns, ["window_id", "id", "segment_id"])
    idr_start_col = first_existing(df.columns, ["idr_start", "IDR_start", "disorder_start"])
    idr_end_col = first_existing(df.columns, ["idr_end", "IDR_end", "disorder_end"])
    out = pd.DataFrame({
        "uniprot": df[uid_col].astype(str).str.strip(),
        "window_start": pd.to_numeric(df[start_col], errors="raise").astype(int),
        "window_end": pd.to_numeric(df[end_col], errors="raise").astype(int),
    })
    out["label"] = df[label_col].astype(str) if label_col else (default_label or "")
    out["window_id"] = df[window_id_col].astype(str) if window_id_col else [f"window_{i:07d}" for i in range(len(df))]
    out["source_row"] = np.arange(len(df), dtype=int)
    if idr_start_col and idr_end_col:
        out["idr_start"] = pd.to_numeric(df[idr_start_col], errors="coerce")
        out["idr_end"] = pd.to_numeric(df[idr_end_col], errors="coerce")
    else:
        out["idr_start"] = np.nan
        out["idr_end"] = np.nan

    for col in ["ELMIdentifier", "slim_sequence", "slim_start", "slim_end", "sequence", "window_center", "idr_index",
                "protein_sequence", "window_seq", "fragment_id", "interaction_class"]:
        if col in df.columns and col not in out.columns:
            out[col] = df[col]
    return out

def load_windows() -> pd.DataFrame:
    if WINDOWS_FILE is not None:
        return normalize_windows(read_table(WINDOWS_FILE))
    frames = []
    for split, path in zip(("train", "val", "test"), SPLIT_FILES):
        frame = normalize_windows(read_table(path))
        frame["split"] = split
        frame["split_row"] = np.arange(len(frame), dtype=int)
        frame["window_id"] = split + ":" + frame["split_row"].astype(str)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)

def label_to_y(label: object) -> float:
    text = str(label).strip().lower()
    if text in {"positive", "pos", "1", "true", "slim"}:
        return 1.0
    if text in {"negative", "neg", "0", "false", "other"}:
        return 0.0
    try:
        value = float(text)
    except (TypeError, ValueError):
        return np.nan
    return value if 0.0 <= value <= 1.0 else np.nan

def percentile_against_reference(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    pass
    values = np.asarray(values, dtype=float)
    reference = np.sort(np.asarray(reference, dtype=float)[np.isfinite(reference)])
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if len(reference) == 0:
        return result
    if len(reference) == 1:
        result[finite] = 1.0
        return result
    left = np.searchsorted(reference, values[finite], side="left")
    right = np.searchsorted(reference, values[finite], side="right")
    average_rank = (left + right - 1.0) / 2.0
    result[finite] = np.clip(average_rank / (len(reference) - 1.0), 0.0, 1.0)
    return result

def normalized_specificity(probabilities: np.ndarray) -> np.ndarray:
    entropy = -np.nansum(probabilities * np.log(probabilities + EPS), axis=1)
    return np.clip(1.0 - entropy / math.log(20.0), 0.0, 1.0)

def build_base_pssm_features(pssm_path: Path) -> pd.DataFrame:
    pssm, sequence, extras = load_pssm_with_extras(pssm_path)
    if not sequence:
        raise ValueError("ASCII PSSM sequence could not be parsed")
    observed = np.asarray(extras.get("observed_pct"), dtype=float)
    information = np.asarray(extras.get("information"), dtype=float)
    real_weight = np.asarray(extras.get("real_match_weight"), dtype=float)
    length = min(len(sequence), pssm.shape[0], len(information), len(real_weight))
    if observed.ndim != 2 or observed.shape[1] != 20:
        raise ValueError(f"observed_pct must be Lx20, got {observed.shape}")
    length = min(length, observed.shape[0])
    sequence = sequence[:length]
    pssm = np.asarray(pssm[:length], dtype=float)
    observed = observed[:length]
    information = information[:length]
    real_weight = real_weight[:length]

    row_mean = np.nanmean(pssm, axis=1, keepdims=True)
    row_std = np.nanstd(pssm, axis=1, keepdims=True)
    zscore = (pssm - row_mean) / (row_std + EPS)
    scaled = zscore / PSSM_TAU
    scaled -= np.nanmax(scaled, axis=1, keepdims=True)
    exp_score = np.exp(scaled)
    distribution = exp_score / (np.nansum(exp_score, axis=1, keepdims=True) + EPS)

    observed_sum = np.nansum(observed, axis=1, keepdims=True)
    observed_probability = np.divide(
        observed,
        observed_sum,
        out=np.full_like(observed, np.nan, dtype=float),
        where=observed_sum > 0,
    )
    result = pd.DataFrame(
        {
            "position": np.arange(1, length + 1, dtype=int),
            "aa": list(sequence),
            "pssm_information": information,
            "real_match_weight": real_weight,
            "log1p_real_match_weight": np.log1p(np.maximum(real_weight, 0.0)),
            "distribution_specificity": normalized_specificity(distribution),
            "observed_specificity": normalized_specificity(observed_probability),
        }
    )
    return result

def add_consensus_for_idr(
    protein_df: pd.DataFrame,
    idr_start: int,
    idr_end: int,
) -> pd.DataFrame:
    pass
    length = len(protein_df)
    start = max(1, int(idr_start))
    end = min(length, int(idr_end))
    if start > end:
        raise ValueError(f"IDR {idr_start}-{idr_end} is outside protein length {length}")
    index = np.arange(start - 1, end)
    info_all = protein_df["pssm_information"].to_numpy(float)
    dist_all = protein_df["distribution_specificity"].to_numpy(float)
    obs_all = protein_df["observed_specificity"].to_numpy(float)
    info_pct = percentile_against_reference(info_all, info_all[index])
    dist_pct = percentile_against_reference(dist_all, dist_all[index])
    obs_pct = percentile_against_reference(obs_all, obs_all[index])
    model_consensus = np.sqrt(np.clip(info_pct * dist_pct, 0.0, 1.0))
    real_weight = protein_df["real_match_weight"].to_numpy(float)
    support_gate = np.clip(real_weight / (1.0 + real_weight + EPS), 0.0, 1.0)
    observed_consensus = np.sqrt(np.clip(model_consensus * obs_pct, 0.0, 1.0))
    consensus_segment = (
        (1.0 - support_gate) * model_consensus
        + support_gate * observed_consensus
    )
    consensus_segment = np.where(
        np.isfinite(observed_consensus), consensus_segment, model_consensus
    )
    consensus_pct = percentile_against_reference(consensus_segment, consensus_segment[index])
    return pd.DataFrame(
        {
            "consensus_conservation": consensus_segment,
            "consensus_percentile_IDR": consensus_pct,
        }
    )

def local_conservation_peak(consensus_window: np.ndarray) -> float:
    pass
    values = np.asarray(consensus_window, dtype=float)
    if values.shape != (WINDOW_LEN,):
        raise ValueError(f"Expected {WINDOW_LEN} consensus values, got {values.shape}")
    center = values[11:18]
    flanks = np.concatenate([values[:7], values[22:]])
    center = center[np.isfinite(center)]
    flanks = flanks[np.isfinite(flanks)]
    if not len(center) or not len(flanks):
        return np.nan
    return float(np.mean(center) - np.mean(flanks))

def build_pssm_tensor(
    windows: pd.DataFrame,
    collect_residue_rows: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pass
    protein_cache: Dict[str, pd.DataFrame] = {}
    pssm_feature_cache: Dict[str, pd.DataFrame] = {}
    consensus_source_cache: Dict[Tuple[str, int, int], pd.DataFrame] = {}
    pssm_path_by_uid: Dict[str, str] = {}
    residue_tensors = []
    window_tensors = []
    valid_masks = []
    metadata_rows = []
    residue_rows = []
    error_rows = []

    for _, window in windows.iterrows():
        uid = str(window["uniprot"]).strip()
        start = int(window["window_start"])
        end = int(window["window_end"])
        window_id = str(window["window_id"])
        try:
            if end - start + 1 != WINDOW_LEN:
                raise ValueError(f"Invalid window length: {start}-{end}")
            if uid not in protein_cache:

                                                                          

                pssm_uid = base_accession(uid)
                pssm_path = find_pssm_file(PSSM_DIR, pssm_uid)
                if pssm_path is None:
                    raise FileNotFoundError(
                        f"No original-protein PSSM file for {uid} (base={pssm_uid})"
                    )
                pssm_key = str(pssm_path.resolve())
                if pssm_key not in pssm_feature_cache:
                    pssm_feature_cache[pssm_key] = build_base_pssm_features(pssm_path)
                protein_cache[uid] = pssm_feature_cache[pssm_key]
                pssm_path_by_uid[uid] = pssm_key
            protein_df = protein_cache[uid]
            idr_start = window.get("idr_start", np.nan)
            idr_end = window.get("idr_end", np.nan)
            if not np.isfinite(idr_start) or not np.isfinite(idr_end):
                if not ALLOW_MISSING_IDR:
                    raise ValueError("Missing IDR bounds")
                idr_start, idr_end = start, end
            key = (int(idr_start), int(idr_end))
            consensus_key = (pssm_path_by_uid[uid], key[0], key[1])
            if consensus_key not in consensus_source_cache:
                consensus_source_cache[consensus_key] = add_consensus_for_idr(protein_df, *key)
            enriched = protein_df.join(consensus_source_cache[consensus_key])
            residue_tensor = np.zeros((WINDOW_LEN, len(PSSM_RESIDUE_FEATURE_NAMES)), dtype=np.float32)
            consensus_window = np.full(WINDOW_LEN, np.nan, dtype=float)
            valid_mask = np.zeros(WINDOW_LEN, dtype=bool)
            window_records = [] if collect_residue_rows else None
            for offset, position in enumerate(range(start, end + 1)):
                if 1 <= position <= len(enriched):
                    row = enriched.iloc[position - 1]
                    residue_tensor[offset] = row[
                        list(PSSM_RESIDUE_FEATURE_NAMES)
                    ].to_numpy(np.float32)
                    consensus_window[offset] = float(row["consensus_conservation"])
                    valid_mask[offset] = True
                    if window_records is not None:
                        window_records.append(
                            (offset, position, str(row["aa"]), residue_tensor[offset])
                        )
            peak = local_conservation_peak(consensus_window)
            if not np.isfinite(residue_tensor[valid_mask]).all() or not np.isfinite(peak):
                raise ValueError("Non-finite PSSM feature after consensus calculation")
            residue_tensors.append(residue_tensor)
            window_tensors.append([peak])
            valid_masks.append(valid_mask)

            meta = window.to_dict()
            meta.update(
                {
                    "idr_start_used": key[0],
                    "idr_end_used": key[1],
                    "pssm_path": pssm_path_by_uid[uid],
                    "local_conservation_peak": peak,
                    "status": "ok",
                }
            )
            metadata_rows.append(meta)
            if collect_residue_rows:
                for offset, position, aa, feature_values in window_records or []:
                    residue_rows.append(
                        {
                            "window_id": window_id,
                            "uniprot": uid,
                            "window_offset": offset + 1,
                            "protein_position": position,
                            "aa": aa,
                            "valid_mask": 1,
                            **{
                                name: float(feature_values[index])
                                for index, name in enumerate(PSSM_RESIDUE_FEATURE_NAMES)
                            },
                        }
                    )
        except Exception as exc:
            error_rows.append(
                {
                    "window_id": window_id,
                    "uniprot": uid,
                    "window_start": start,
                    "window_end": end,
                    "module": "pssm",
                    "error": str(exc),
                }
            )
            if not SKIP_BAD_WINDOWS:
                raise
    if not residue_tensors:
        raise RuntimeError("No PSSM windows processed")
    return (
        np.stack(residue_tensors).astype(np.float32),
        np.asarray(window_tensors, dtype=np.float32),
        np.stack(valid_masks),
        pd.DataFrame(metadata_rows),
        pd.DataFrame(residue_rows),
        pd.DataFrame(error_rows),
    )

def chemical_feature_names(categories: Dict[str, str], include_motif: bool) -> list[str]:
    names = []
    for category in categories:
        for metric in ("density", "X", "center_enrich"):
            base = f"{metric}_{category}"
            names.extend([base, f"{base}_percentile_IDR"])
    names.extend(["NCPR", "FCR", "RK_logratio", "ED_logratio", "hydropathy_sum", "PPII_sum"])
    if include_motif:
        names.extend(["motif_min_abs_offset", "motif_center_enrichment", "motif_count_log"])
    names.append("uversky_score")
    return names

def padded_protein_window(sequence: str, start: int, end: int) -> str:
    return "".join(
        sequence[position - 1] if 1 <= position <= len(sequence) else "X"
        for position in range(start, end + 1)
    )

def single_chemical_window_features(
    sequence: str,
    categories: Dict[str, str],
    lookup: Dict[int, Tuple[float, float]],
    regexes,
    is_n_terminal_idr: bool,
    is_c_terminal_idr: bool,
) -> Dict[str, float]:
    pass
    if len(sequence) != WINDOW_LEN:
        raise ValueError(f"Chemical window must be {WINDOW_LEN} aa, got {len(sequence)}")
    output: Dict[str, float] = {}
    center_mask, flank_mask = center_flank_masks(WINDOW_LEN, CENTER_LEN)
    for name, letters in categories.items():
        allowed = set(letters)
        binary = np.fromiter((aa in allowed for aa in sequence), dtype=np.int8, count=WINDOW_LEN)
        count = int(binary.sum())
        patch_strength = float(
            local_patch_strength_for_all_windows(binary, WINDOW_LEN, PATCH, TOP_PATCH_K)[0]
        )
        null_mean, null_sd = lookup[count]
        output[f"density_{name}"] = count / float(WINDOW_LEN)
        output[f"X_{name}"] = (
            (patch_strength - null_mean) / (null_sd + EPS) if null_sd > EPS else 0.0
        )
        output[f"center_enrich_{name}"] = float(binary @ center_mask - binary @ flank_mask)

    counts = {aa: sequence.count(aa) for aa in "RKHDE"}
    positive = counts["R"] + counts["K"] + counts["H"]
    negative = counts["D"] + counts["E"]
    output["NCPR"] = (positive - negative) / float(WINDOW_LEN)
    output["FCR"] = (positive + negative) / float(WINDOW_LEN)
    output["RK_logratio"] = math.log((counts["R"] + ALPHA) / (counts["K"] + ALPHA))
    output["ED_logratio"] = math.log((counts["E"] + ALPHA) / (counts["D"] + ALPHA))
    output["hydropathy_sum"] = float(sum(KD_UVERSKY_0_1.get(aa, 0.0) for aa in sequence))
    output["PPII_sum"] = float(sum(PPII_HILSER.get(aa, 0.0) for aa in sequence))
    motif = regex_features_for_idr(
        sequence,
        WINDOW_LEN,
        regexes,
        is_n_terminal_idr,
        is_c_terminal_idr,
        CENTER_LEN,
    )
    output.update({name: float(values[0]) for name, values in motif.items()})
    output["uversky_score"] = output["NCPR"] - output["hydropathy_sum"] / WINDOW_LEN
    return output

def build_chemical_tensors(windows: pd.DataFrame) -> Tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    pass

    categories = load_categories(CHEMICAL_CATEGORIES_JSON)
    regexes = load_regexes(MOTIF_REGEX_JSON) if MOTIF_REGEX_JSON else {}
    idrs, ordered_by_pid = load_chemical_regions(FRAGMENT_FILE)
    for uid in windows["uniprot"].astype(str).unique():
        base_uid = base_accession(uid)
        if uid not in ordered_by_pid and base_uid in ordered_by_pid:
            ordered_by_pid[uid] = ordered_by_pid[base_uid]
    idrs_by_uid = {}
    for idr in idrs:
        idrs_by_uid.setdefault(idr.protein_id, []).append(idr)

                                                                          
    seqs = {}
    if "protein_sequence" in windows.columns:
        for row in windows.itertuples():
            uid = str(row.uniprot)
            seq = str(row.protein_sequence).strip().upper()
            if uid and seq and seq.lower()!="nan":
                if uid not in seqs or len(seq) > len(seqs[uid]):
                    seqs[uid] = seq

                                                                          
    for uid in windows["uniprot"].astype(str).unique():
        if uid in seqs:
            continue
        pssm_path = find_pssm_file(PSSM_DIR, base_accession(uid))
        if pssm_path is None:
            continue
        _pssm, pssm_sequence, _extras = load_pssm_with_extras(pssm_path)
        if pssm_sequence:

                                                                         

                           
            chars = list(pssm_sequence)
            uid_rows = windows[windows["uniprot"].astype(str).eq(uid)]
            if base_accession(uid) != uid and "sequence" in uid_rows.columns:
                for variant_row in uid_rows.itertuples(index=False):
                    window_sequence = str(variant_row.sequence)
                    for offset, amino_acid in enumerate(window_sequence):
                        position = int(variant_row.window_start) + offset
                        if 1 <= position <= len(chars) and amino_acid != "X":
                            chars[position - 1] = amino_acid
            seqs[uid] = "".join(chars)

    lookup_df = build_nardini_lookup(WINDOW_LEN, PATCH, TOP_PATCH_K, NULL_SHUFFLE, SEED)
    lookup = {int(r.n): (float(r.mu_null), float(r.sd_null)) for r in lookup_df.itertuples(index=False)}

    names = chemical_feature_names(categories, include_motif=bool(regexes))

    resolved_rows = []
    errors = []
    needed_idrs = {}
    for _, row in windows.iterrows():
        uid = str(row["uniprot"])
        ws = int(row["window_start"]); we = int(row["window_end"])
        if we - ws + 1 != WINDOW_LEN:
            if SKIP_BAD_WINDOWS:
                errors.append({"window_id": row["window_id"], "uniprot": uid, "reason": f"length {we-ws+1}"})
                continue
            else:
                raise ValueError("Bad window length")
        if pd.notna(row.get("idr_start")) and pd.notna(row.get("idr_end")):
            idr = Segment(uid, int(row["idr_start"]), int(row["idr_end"]), "all_idr", "input")
        else:
            idr = None
            for idr_candidate in idrs_by_uid.get(uid, []):
                if idr_candidate.start <= ws and we <= idr_candidate.end:
                    idr = idr_candidate
                    break
            if idr is None:
                errors.append({"window_id": row["window_id"], "uniprot": uid, "reason": "no_containing_idr"})
                continue
        if uid not in seqs:
            errors.append({"window_id": row["window_id"], "uniprot": uid, "reason": "sequence_not_found"})
            continue
        if idr.end > len(seqs[uid]):
            errors.append({"window_id": row["window_id"], "uniprot": uid, "reason": "idr_end_exceeds_sequence"})
            continue
        new_row = row.copy()
        new_row["idr_start"] = idr.start
        new_row["idr_end"] = idr.end
        resolved_rows.append(new_row)
        needed_idrs[(uid, idr.start, idr.end)] = idr
    if not resolved_rows:
        return np.empty((0, len(names)), dtype=float), pd.DataFrame(), pd.DataFrame(errors)
    resolved = pd.DataFrame(resolved_rows).reset_index(drop=True)

    feature_frames = []
    for (uid, st, en), idr in needed_idrs.items():
        df = features_for_idr(
            uid, idr, seqs[uid], categories, lookup,
            SimpleNamespace(
                window=WINDOW_LEN, window_len=WINDOW_LEN, patch=PATCH, top_patch_k=TOP_PATCH_K,
                center_len=CENTER_LEN, alpha=ALPHA, eps=EPS, mad_floor=MAD_FLOOR, clip_rel=CLIP_REL
            ),
            regexes=regexes,
            ordered_by_pid=ordered_by_pid,
        )
        if df is not None:
            feature_frames.append(df)
    if not feature_frames:
        return np.empty((0, len(names)), dtype=float), pd.DataFrame(), pd.DataFrame(errors)
    idr_features = pd.concat(feature_frames, ignore_index=True)

    merge_cols = ["protein_id", "idr_start", "idr_end", "window_start", "window_end", "window_seq"] + [n for n in names if n != "uversky_score"]
    merged = resolved.merge(
        idr_features[merge_cols],
        left_on=["uniprot","idr_start","idr_end","window_start","window_end"],
        right_on=["protein_id","idr_start","idr_end","window_start","window_end"],
        how="left"
    )

                                                                            

                                                     
    missing_feature = merged[names[0]].isna()
    if missing_feature.any():
        background_by_idr = {
            key: group
            for key, group in idr_features.groupby(
                ["protein_id", "idr_start", "idr_end"], sort=False
            )
        }
        relative_bases = [
            f"{metric}_{category}"
            for category in categories
            for metric in ("density", "X", "center_enrich")
        ]
        for row_index in merged.index[missing_feature]:
            row = merged.loc[row_index]
            uid = str(row["uniprot"])
            idr = Segment(uid, int(row["idr_start"]), int(row["idr_end"]), "all_idr", "input")
            sequence = padded_protein_window(
                seqs[uid], int(row["window_start"]), int(row["window_end"])
            )
            terminal = idr_terminal_flags(idr, ordered_by_pid.get(uid, []))
            direct = single_chemical_window_features(
                sequence,
                categories,
                lookup,
                regexes,
                terminal[0],
                terminal[1],
            )
            background = background_by_idr[(uid, idr.start, idr.end)]
            for name, value in direct.items():
                if name in merged.columns:
                    merged.at[row_index, name] = value
            for base_name in relative_bases:
                values = pd.to_numeric(background[base_name], errors="coerce").to_numpy(float)
                values = np.sort(values[np.isfinite(values)])
                percentile = (
                    np.searchsorted(values, direct[base_name], side="right") / len(values)
                    if len(values)
                    else np.nan
                )
                merged.at[row_index, f"{base_name}_percentile_IDR"] = percentile

    merged["uversky_score"] = merged["NCPR"] - merged["hydropathy_sum"] / WINDOW_LEN

    tensor = merged[names].to_numpy(dtype=float)

    meta_list = []
    for i, row in merged.iterrows():
        meta = {
            "tensor_index": i,
            "window_id": row["window_id"],
            "uniprot": row["uniprot"],
            "window_start": int(row["window_start"]),
            "window_end": int(row["window_end"]),
            "idr_start": int(row["idr_start"]),
            "idr_end": int(row["idr_end"]),
            "label": row["label"],
        }
        meta_list.append(meta)
    metadata = pd.DataFrame(meta_list)
    return tensor, metadata, pd.DataFrame(errors)

def build_context_tensors(windows: pd.DataFrame, return_structured: bool = False):
    pass

    ordered_map = {}
    if Path(FRAGMENT_FILE).exists():
        df = pd.read_excel(FRAGMENT_FILE)
        uid_col = first_existing(df.columns, ["Protein_Accession", "UniProt_ID", "uniprot"])
        ord_col = first_existing(df.columns, ORDERED_COLUMN_CANDIDATES)
        print(f"Selected ordered column: {ord_col}")
        if uid_col is not None and ord_col is not None:
            for _, row in df.iterrows():
                uid = str(row[uid_col]).strip()
                if uid and uid.lower()!="nan":
                    intervals = parse_region_string(row[ord_col])
                    if intervals: ordered_map[uid] = intervals

    pattern_map = {}
    if Path(NARDINI_CACHE_FILE).exists():
        pat_df = pd.read_csv(NARDINI_CACHE_FILE, sep='\t')
        for _, row in pat_df.iterrows():
            uid = str(row['uniprot'])
            st = int(row['idr_start']); en = int(row['idr_end'])
            feat = row[list(NARDINI_FEATURE_NAMES)].values.astype(np.float32)
            pattern_map[(uid, st, en)] = feat
    ensemble_map = {}
    if Path(GUSSIDR_CACHE_FILE).exists():
        ens_df = load_ensemble_frame(GUSSIDR_CACHE_FILE)
        for _, row in ens_df.iterrows():
            uid = str(row['uniprot'])
            st = int(row['idr_start']); en = int(row['idr_end'])
            feat = row[[f"feat_{i}" for i in range(1, 39)]].values.astype(np.float32)
            ensemble_map[(uid, st, en)] = feat
    if not Path(ESM_REGION_METADATA_FILE).is_file():
        raise FileNotFoundError(f"ESM region metadata missing: {ESM_REGION_METADATA_FILE}")
    esm_metadata = pd.read_csv(ESM_REGION_METADATA_FILE, sep="\t")
    required_esm_columns = {
        "uniprot", "region_type", "original_start", "original_end",
        "crop_start", "crop_end", "crop_length", "embedding_path",
    }
    missing_esm_columns = required_esm_columns - set(esm_metadata.columns)
    if missing_esm_columns:
        raise ValueError(
            f"ESM region metadata missing columns: {sorted(missing_esm_columns)}"
        )
    esm_metadata_by_key = {}
    for esm_row in esm_metadata.itertuples(index=False):
        key = (
            str(esm_row.uniprot), str(esm_row.region_type),
            int(esm_row.original_start), int(esm_row.original_end),
        )
        if key in esm_metadata_by_key:
            raise ValueError(f"Duplicate ESM biological region metadata: {key}")
        esm_metadata_by_key[key] = esm_row._asdict()
    esm_array_cache = {}
    raw_features = []
    left_domain_features = []
    right_domain_features = []
    pattern_features = []
    ensemble_features = []
    ensemble_full_features = []
    metadata = []
    errors = []
    for idx, row in windows.iterrows():
        uid = str(row['uniprot'])
        idr_s = int(row['idr_start']); idr_e = int(row['idr_end'])
        ws = int(row['window_start']); we = int(row['window_end'])
        window_id = row.get('window_id', f'window_{idx:07d}')

        ordered = ordered_map.get(uid, ordered_map.get(base_accession(uid), []))
        left_domain = None; right_domain = None
        for st, en in ordered:
            if en < idr_s:
                if left_domain is None or en > left_domain[1]:
                    left_domain = (st, en)
            if st > idr_e:
                if right_domain is None or st < right_domain[0]:
                    right_domain = (st, en)

        def get_emb(dom, side):
            if dom is None:
                return None, None
            st, en = dom
            region_type = f"{side}_domain"
            key = (uid, region_type, st, en)
            if key not in esm_metadata_by_key:
                key = (base_accession(uid), region_type, st, en)
            if key not in esm_metadata_by_key:
                raise KeyError(f"ESM metadata missing for {uid} {region_type} {st}-{en}")
            region = esm_metadata_by_key[key]
            esm_path = str(region["embedding_path"])
            if esm_path not in esm_array_cache:
                array = np.load(esm_path)
                expected = (int(region["crop_length"]), 1280)
                if array.shape != expected:
                    raise ValueError(
                        f"Invalid ESM shape for {uid} {region_type}: "
                        f"{array.shape}, expected {expected}"
                    )
                esm_array_cache[esm_path] = np.asarray(array, dtype=np.float32)
            emb = esm_array_cache[esm_path]
            crop_start = int(region["crop_start"])
            crop_end = int(region["crop_end"])
            positions = np.arange(crop_start, crop_end + 1)
            if side == "left":
                dist = idr_s - positions
            elif side == "right":
                dist = positions - idr_e
            else:
                raise ValueError("side must be 'left' or 'right'")
            if np.any(dist <= 0):
                raise ValueError(f"{region_type} crop overlaps IDR {idr_s}-{idr_e}")
            w = np.exp(-dist / LAMBDA_DECAY)
            w = w / (w.sum() + 1e-8)
            return np.sum(emb * w[:, np.newaxis], axis=0), region
        left_emb = None; right_emb = None
        left_region = None; right_region = None
        try:
            if left_domain:
                left_emb, left_region = get_emb(left_domain, "left")
            if right_domain:
                right_emb, right_region = get_emb(right_domain, "right")
        except Exception as e:
            errors.append({"window_id": window_id, "uniprot": uid, "error": str(e)})
            continue

        if left_emb is not None and right_emb is not None:
            domain_emb = (left_emb + right_emb) / 2.0
            has_left, has_right = 1, 1
        elif left_emb is not None:
            domain_emb = left_emb; has_left, has_right = 1, 0
        elif right_emb is not None:
            domain_emb = right_emb; has_left, has_right = 0, 1
        else:
            domain_emb = np.zeros(1280, dtype=np.float32); has_left, has_right = 0, 0
        onehot = np.array([has_left, has_right], dtype=np.float32)

        key = (uid, idr_s, idr_e)

                                                                           

        pattern_key = key
        if pattern_key not in pattern_map:
            errors.append({"window_id": window_id, "uniprot": uid, "error": "Pattern missing"})
            continue
        ensemble_key = key
        if ensemble_key not in ensemble_map:
            errors.append({"window_id": window_id, "uniprot": uid, "error": "Ensemble missing"})
            continue
        pattern_feat = pattern_map[pattern_key]
        ensemble_feat = ensemble_map[ensemble_key]
        raw = np.concatenate([domain_emb, onehot, pattern_feat, ensemble_feat])
        raw_features.append(raw)
        left_domain_features.append(
            left_emb if left_emb is not None else np.zeros(1280, dtype=np.float32)
        )
        right_domain_features.append(
            right_emb if right_emb is not None else np.zeros(1280, dtype=np.float32)
        )
        pattern_features.append(pattern_feat)
        ensemble_features.append(ensemble_feat[list(ENSEMBLE_KEEP_INDICES)])
        ensemble_full_features.append(ensemble_feat)
        metadata.append({
            "window_id": window_id, "uniprot": uid,
            "window_start": ws, "window_end": we,
            "idr_start": idr_s, "idr_end": idr_e, "label": row["label"],
            "left_domain_start": left_domain[0] if left_domain else np.nan,
            "left_domain_end": left_domain[1] if left_domain else np.nan,
            "left_esm_crop_start": left_region["crop_start"] if left_region else np.nan,
            "left_esm_crop_end": left_region["crop_end"] if left_region else np.nan,
            "right_domain_start": right_domain[0] if right_domain else np.nan,
            "right_domain_end": right_domain[1] if right_domain else np.nan,
            "right_esm_crop_start": right_region["crop_start"] if right_region else np.nan,
            "right_esm_crop_end": right_region["crop_end"] if right_region else np.nan,
        })
    if not raw_features:
        base = (
            np.empty((0, 1350), dtype=np.float32),
            pd.DataFrame(),
            pd.DataFrame(errors),
        )
        if not return_structured:
            return base
        return (*base, *(
            np.empty((0, width), dtype=np.float32)
            for width in (1280, 1280, 36, 32, 38)
        ))
    X = np.stack(raw_features, axis=0)
    keep = np.ones(X.shape[1], dtype=bool)
    keep[list(CONTEXT_DROP_INDICES)] = False
    X = X[:, keep]
    if X.shape[1] != 1350:
        raise RuntimeError(f"Expected 1,350 context features, got {X.shape}")
    base = (X, pd.DataFrame(metadata), pd.DataFrame(errors))
    if not return_structured:
        return base
    return (
        *base,
        np.stack(left_domain_features).astype(np.float32),
        np.stack(right_domain_features).astype(np.float32),
        np.stack(pattern_features).astype(np.float32),
        np.stack(ensemble_features).astype(np.float32),
        np.stack(ensemble_full_features).astype(np.float32),
    )

                                                

@dataclass
class Interval:
    start: int; end: int
    def __post_init__(self):
        if self.start > self.end: raise ValueError
    @property
    def length(self): return self.end - self.start + 1

@dataclass
class EnsembleSummary:
    sample_id: str; as_mean: float; ree_median: float; rg_median: float; v: float; spr: float

def load_ensemble_tsv(path):
    df = load_ensemble_frame(path)
    feat_cols = [f'feat_{i}' for i in range(1,39)]
    out = {}
    for _, row in df.iterrows():
        uid = str(row['uniprot'])
        st = int(row['idr_start']); en = int(row['idr_end'])
        vals = np.array([row[c] for c in feat_cols], dtype=float)

        def gamma_median(w, conc, rate):
            valid = (w>0)&(conc>0)&(rate>0)
            if not np.any(valid): return np.nan
            w = w[valid]/np.sum(w[valid]); c=conc[valid]; r=rate[valid]
            def cdf(x): return np.sum(w * gamma.cdf(x, c, scale=1/r))
            low=0.0; high=max(np.sum(w*c/r)*10, 1.0)
            for _ in range(20):
                if cdf(high)<0.5: high*=2
                else: break
            for _ in range(100):
                mid=(low+high)/2
                if cdf(mid)<0.5: low=mid
                else: high=mid
            return (low+high)/2
        as_vals = vals[0:12]
        weights = as_vals[0:4]; weights = weights/np.sum(weights) if np.sum(weights)>0 else np.ones(4)/4
        alpha = as_vals[4:8]; beta=as_vals[8:12]
        as_mean = np.sum(weights * alpha/(alpha+beta+1e-8))
        ree_vals = vals[12:24]; rg_vals = vals[24:36]
        w_ree = ree_vals[0:4]; w_ree = w_ree/np.sum(w_ree) if np.sum(w_ree)>0 else np.ones(4)/4
        ree_median = gamma_median(w_ree, ree_vals[4:8], ree_vals[8:12])
        w_rg = rg_vals[0:4]; w_rg = w_rg/np.sum(w_rg) if np.sum(w_rg)>0 else np.ones(4)/4
        rg_median = gamma_median(w_rg, rg_vals[4:8], rg_vals[8:12])
        v = vals[36]; spr=vals[37]
        out[(uid, st, en)] = EnsembleSummary(sample_id=f"{uid}({st}-{en})", as_mean=as_mean,
                                             ree_median=ree_median, rg_median=rg_median, v=v, spr=spr)
    return out

def build_accessibility_tensors(
    windows: pd.DataFrame,
    rg_norm_ref: Optional[float] = None,
    return_reference: bool = False,
):

    region_db = {}
    if FRAGMENT_FILE.exists():
        df = pd.read_excel(FRAGMENT_FILE)
        uid_col = first_existing(df.columns, ["UniProt_ID","Protein_Accession","uniprot"])
        disorder_col = first_existing(df.columns, DISORDERED_COLUMN_CANDIDATES)
        ordered_col = first_existing(df.columns, ORDERED_COLUMN_CANDIDATES)
        print(f"Selected IDR column: {disorder_col}")
        print(f"Selected ordered column: {ordered_col}")
        if uid_col and disorder_col:
            for _, row in df.iterrows():
                uid = str(row[uid_col]).strip()
                if not uid or uid.lower()=="nan": continue
                dis = []
                for part in re.split(r'[,;，；\s]+', str(row[disorder_col])):
                    if '-' in part:
                        st, en = map(int, part.split('-'))
                        dis.append(Interval(min(st,en), max(st,en)))
                ord_lst = []
                if ordered_col:
                    for part in re.split(r'[,;，；\s]+', str(row[ordered_col])):
                        if '-' in part:
                            st, en = map(int, part.split('-'))
                            ord_lst.append(Interval(min(st,en), max(st,en)))
                region_db[uid] = {"disorder": dis, "ordered": ord_lst}
    ensembles = load_ensemble_tsv(GUSSIDR_CACHE_FILE)

    errors = []
    prepared = []
    for idx, row in windows.iterrows():
        uid = str(row["uniprot"])
        ws = int(row["window_start"]); we = int(row["window_end"])
        if we - ws + 1 != WINDOW_LEN:
            if SKIP_BAD_WINDOWS:
                errors.append({"source_row": idx, "window_id": row["window_id"], "uniprot": uid, "reason": "bad length"})
                continue
            else: raise ValueError("bad length")

        if pd.notna(row.get("idr_start")) and pd.notna(row.get("idr_end")):
            idr = Interval(int(row["idr_start"]), int(row["idr_end"]))
        else:
            idr = None
            region_record = region_db.get(uid, region_db.get(base_accession(uid), {}))
            for idr_candidate in region_record.get("disorder", []):
                if idr_candidate.start <= ws and we <= idr_candidate.end:
                    idr = idr_candidate
                    break
            if idr is None:
                errors.append({"source_row": idx, "window_id": row["window_id"], "uniprot": uid, "reason": "no_containing_idr"})
                continue

                                                                               

                                                                      

        window_center = (ws + we) / 2.0
        if not (idr.start <= window_center <= idr.end):
            errors.append({
                "source_row": idx,
                "window_id": row["window_id"],
                "uniprot": uid,
                "reason": "window_center_outside_idr",
            })
            continue

                                                                           
        summary = ensembles.get((uid, idr.start, idr.end))
        if summary is None:
            errors.append({"source_row": idx, "window_id": row["window_id"], "uniprot": uid, "reason": "ensemble_not_found"})
            continue
        rg_norm = summary.rg_median / (idr.length ** (summary.v if np.isfinite(summary.v) and summary.v>0 else DEFAULT_V))
        ordered = region_db.get(uid, region_db.get(base_accession(uid), {})).get("ordered", [])
        prepared.append({"row": row, "idr": idr, "ordered": ordered, "summary": summary, "rg_norm": rg_norm})
    if not prepared:
        base = (np.empty((0,7), dtype=float), pd.DataFrame(), pd.DataFrame(errors))
        return (*base, np.nan) if return_reference else base
    if rg_norm_ref is None:
        if "split" not in windows.columns:
            raise ValueError(
                "Accessibility rg_norm_ref must be supplied for inference/single-table "
                "extraction; only a combined train/val/test build may fit it"
            )
        training_idrs = {}
        for item in prepared:
            row = item["row"]
            if str(row.get("split", "")) != "train":
                continue
            idr = item["idr"]
            key = (str(row["uniprot"]), idr.start, idr.end)
            if np.isfinite(item["rg_norm"]):
                training_idrs[key] = float(item["rg_norm"])
        if not training_idrs:
            raise ValueError("No finite unique training-IDR rg_norm values were available")
        rg_norm_ref = float(np.median(list(training_idrs.values())))
    elif not np.isfinite(rg_norm_ref) or rg_norm_ref <= 0:
        raise ValueError(f"Invalid accessibility rg_norm_ref: {rg_norm_ref}")

    tensor_rows=[]; meta_rows=[]
    for item in prepared:
        row = item["row"]; idr=item["idr"]; ordered=item["ordered"]; summary=item["summary"]; rg_norm=item["rg_norm"]
        ws=int(row["window_start"]); we=int(row["window_end"])
        c=(ws+we)/2.0
        L=idr.length
        if L<=1: p_i=0.5; d_edge=0.0
        else:
            p_i=(c-idr.start)/(L-1)
            d_edge=min(c-idr.start, idr.end-c)/L

        left_dom=None; right_dom=None
        for dom in ordered:
            if dom.end < idr.start:
                if left_dom is None or dom.end > left_dom.end:
                    left_dom = dom
            if dom.start > idr.end:
                if right_dom is None or dom.start < right_dom.start:
                    right_dom = dom
        m_l = 1.0 if left_dom else 0.0
        m_r = 1.0 if right_dom else 0.0
        n_l = max(0.0, c-left_dom.end) if left_dom else np.nan
        n_r = max(0.0, right_dom.start-c) if right_dom else np.nan
        if np.isfinite(rg_norm) and np.isfinite(rg_norm_ref) and rg_norm>0:
            g_i = rg_norm_ref / (rg_norm + EPS)
        else: g_i = np.nan
        def occ(m, dist):
            if not m or not np.isfinite(g_i) or g_i<=0: return 0.0
            return math.exp(-dist / (LAMBDA0 * g_i + EPS))
        o_l = occ(m_l, n_l); o_r = occ(m_r, n_r); o_cis = 1.0 - (1.0-o_l)*(1.0-o_r)
        X = np.array([p_i, d_edge, m_l, m_r, o_l, o_r, o_cis], dtype=float)
        tensor_rows.append(X)
        meta = {"window_id": row["window_id"], "uniprot": row["uniprot"], "window_start": ws, "window_end": we,
                "idr_start": idr.start, "idr_end": idr.end, "label": row["label"]}
        meta_rows.append(meta)
    base = (np.stack(tensor_rows, axis=0), pd.DataFrame(meta_rows), pd.DataFrame(errors))
    return (*base, float(rg_norm_ref)) if return_reference else base

def context_feature_names() -> list[str]:
    raw = (
        [f"domain_esm_{index:04d}" for index in range(1280)]
        + ["has_left_ordered_domain", "has_right_ordered_domain"]
        + list(NARDINI_FEATURE_NAMES)
        + [f"ensemble_feat_{index}" for index in range(1, 39)]
    )
    return [name for index, name in enumerate(raw) if index not in CONTEXT_DROP_INDICES]

ACCESSIBILITY_FEATURE_NAMES = (
    "p_I",
    "d_edge_norm",
    "m_L",
    "m_R",
    "occlusion_L",
    "occlusion_R",
    "occlusion_cis",
)

def validate_required_inputs(
    windows: pd.DataFrame,
    esm_python: Optional[Path],
    nardini_python: Optional[Path],
) -> pd.DataFrame:
    required_files = (
        FRAGMENT_FILE,
        CHEMICAL_CATEGORIES_JSON,
        MOTIF_REGEX_JSON,
    )
    issues = [f"missing file: {path}" for path in required_files if not Path(path).is_file()]
    window_inputs = (WINDOWS_FILE,) if WINDOWS_FILE is not None else SPLIT_FILES
    issues.extend(f"missing file: {path}" for path in window_inputs if not Path(path).is_file())
    if not PSSM_DIR.is_dir():
        issues.append(f"missing directory: {PSSM_DIR}/")
    if esm_python is None or not esm_python.is_file():
        issues.append(
            "fair-esm Python was not found; pass --esm-python or set ESM_PYTHON"
        )
    if nardini_python is None or not nardini_python.is_file():
        issues.append(
            "NARDINI Python was not found; pass --nardini-python or set NARDINI_PYTHON"
        )
    for name in ("As", "ree", "rg", "v", "spr"):
        for filename in (f"model_best_{name}.pth", f"best_hyperparameters_{name}.json"):
            path = GUSSIDR_DIR / name / filename
            if not path.is_file():
                issues.append(f"missing GussIDR asset: {path}")
    if PSSM_DIR.is_dir():
        unique_uids = windows["uniprot"].astype(str).map(base_accession).unique()
        missing_pssm = [uid for uid in unique_uids if find_pssm_file(PSSM_DIR, uid) is None]
        if missing_pssm:
            issues.append(
                f"PSSM coverage: {len(unique_uids) - len(missing_pssm)}/{len(unique_uids)} proteins"
            )
    if issues:
        raise FileNotFoundError("Feature preflight failed:\n- " + "\n- ".join(issues))

                                                                               
    idrs = build_idr_manifest(windows, PSSM_DIR)
    build_esm_region_manifest(idrs, FRAGMENT_FILE)
    return idrs

def validate_reused_derived_caches(idrs: pd.DataFrame) -> None:
    pass
    issues = []
    expected_idr_keys = set(
        idrs[IDR_KEY_COLUMNS].itertuples(index=False, name=None)
    )
    idr_sequence_hash = stable_frame_hash(
        idrs, [*IDR_KEY_COLUMNS, "sequence"]
    )
    for cache_path, manifest_name in (
        (NARDINI_CACHE_FILE, "nardini_features.manifest.json"),
        (GUSSIDR_CACHE_FILE, "gussidr_ensemble_features.manifest.json"),
    ):
        if not cache_path.is_file():
            issues.append(f"missing cache: {cache_path}")
            continue
        frame = pd.read_csv(cache_path, sep="\t")
        if frame.duplicated(IDR_KEY_COLUMNS).any():
            issues.append(f"duplicate IDR keys: {cache_path}")
        observed_keys = set(
            frame[IDR_KEY_COLUMNS].itertuples(index=False, name=None)
        )
        if observed_keys != expected_idr_keys:
            issues.append(
                f"IDR key coverage mismatch: {cache_path} "
                f"(missing={len(expected_idr_keys - observed_keys)}, "
                f"extra={len(observed_keys - expected_idr_keys)})"
            )
        manifest_path = CACHE_DIR / manifest_name
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            issues.append(f"missing/invalid provenance: {manifest_path}")
        else:
            if manifest.get("idr_sequence_sha256") != idr_sequence_hash:
                issues.append(f"IDR sequence signature mismatch: {cache_path}")

    expected_regions = build_esm_region_manifest(idrs, FRAGMENT_FILE)
    expected_region_keys = set(
        expected_regions[REGION_KEY_COLUMNS].itertuples(index=False, name=None)
    )
    if not ESM_REGION_METADATA_FILE.is_file():
        issues.append(f"missing cache: {ESM_REGION_METADATA_FILE}")
    else:
        esm_metadata = pd.read_csv(ESM_REGION_METADATA_FILE, sep="\t")
        observed_region_keys = set(
            esm_metadata[REGION_KEY_COLUMNS].itertuples(index=False, name=None)
        )
        if (
            esm_metadata.duplicated(REGION_KEY_COLUMNS).any()
            or observed_region_keys != expected_region_keys
        ):
            issues.append(
                f"ESM region key coverage mismatch: {ESM_REGION_METADATA_FILE} "
                f"(missing={len(expected_region_keys - observed_region_keys)}, "
                f"extra={len(observed_region_keys - expected_region_keys)})"
            )
        elif "embedding_path" not in esm_metadata:
            issues.append(f"embedding_path missing: {ESM_REGION_METADATA_FILE}")
        else:
            missing_embeddings = sum(
                not Path(path).is_file() for path in esm_metadata["embedding_path"]
            )
            if missing_embeddings:
                issues.append(f"missing ESM embedding files: {missing_embeddings}")

    esm_manifest_path = CACHE_DIR / "esm2_region_cache.manifest.json"
    expected_region_hash = stable_frame_hash(
        expected_regions,
        [
            *REGION_KEY_COLUMNS, "crop_start", "crop_end", "sequence",
            "embedding_file",
        ],
    )
    try:
        esm_manifest = json.loads(esm_manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        issues.append(f"missing/invalid provenance: {esm_manifest_path}")
    else:
        if esm_manifest.get("region_sequence_sha256") != expected_region_hash:
            issues.append("ESM region sequence signature mismatch")

    if issues:
        raise RuntimeError(
            "Existing derived caches are not valid for the current altered-IDR "
            "dataset. Re-run without --reuse-derived-cache so missing variant "
            "ESM/NARDINI/GussIDR features are generated:\n- "
            + "\n- ".join(issues)
        )

def aligned_block(
    tensor: np.ndarray,
    metadata: pd.DataFrame,
    window_ids: Sequence[str],
) -> np.ndarray:
    if len(tensor) != len(metadata):
        raise ValueError(f"Tensor rows {len(tensor)} != metadata rows {len(metadata)}")
    if metadata["window_id"].duplicated().any():
        raise ValueError("Module metadata contains duplicate window_id values")
    index_by_id = {str(window_id): index for index, window_id in enumerate(metadata["window_id"])}
    return tensor[[index_by_id[str(window_id)] for window_id in window_ids]]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build all SLiM model features")
    parser.add_argument(
        "--windows",
        type=Path,
        default=None,
        help="Optional single window table; default concatenates train/val/test",
    )
    parser.add_argument("--pssm-dir", type=Path, default=PSSM_DIR)
    parser.add_argument("--fragment", type=Path, default=FRAGMENT_FILE)
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    parser.add_argument("--gussidr-dir", type=Path, default=GUSSIDR_DIR)
    parser.add_argument("--esm-python", type=Path, default=None)
    parser.add_argument("--esm-batch-size", type=int, default=1)
    parser.add_argument(
        "--esm-device", choices=("auto", "cpu", "mps", "cuda"), default="auto"
    )
    parser.add_argument("--nardini-python", type=Path, default=None)
    parser.add_argument("--gussidr-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--nardini-num-scrambles", type=int, default=1000)
    parser.add_argument("--nardini-random-seed", type=int, default=1)
    parser.add_argument("--nardini-workers", type=int, default=1)
    parser.add_argument("--gussidr-batch-size", type=int, default=4)
    parser.add_argument(
        "--gussidr-device", choices=("auto", "cpu", "mps", "cuda"), default="auto"
    )
    cache_mode = parser.add_mutually_exclusive_group()
    cache_mode.add_argument("--force-derived-cache", action="store_true")
    cache_mode.add_argument("--reuse-derived-cache", action="store_true")
    parser.add_argument(
        "--accessibility-rg-norm-ref",
        type=float,
        default=None,
        help=(
            "Fixed training-derived accessibility normalization reference for "
            "single-table inference. Required with --windows."
        ),
    )
    parser.add_argument("--prepare-cache-only", action="store_true")
    parser.add_argument("--out-prefix", type=Path, default=OUT_PREFIX)
    parser.add_argument("--limit", type=int, default=None, help="Optional diagnostic row limit")
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument(
        "--save-pssm-residue-table",
        action="store_true",
        help="Write the per-residue feature table",
    )
    parser.add_argument(
        "--prepare-contextual-esm2",
        action="store_true",
        help=(
            "Build contextual ESM2 chunks and window-to-chunk metadata "
            "using the selected crop mode."
        ),
    )
    parser.add_argument(
        "--prepare-motifae-sparse",
        action="store_true",
        help="Build per-window sparse MotifAE activations from contextual ESM2.",
    )
    parser.add_argument(
        "--v4-feature-only",
        action="store_true",
        help="Stop after optional contextual ESM2/MotifAE preparation.",
    )
    parser.add_argument(
        "--contextual-esm2-plan-only",
        action="store_true",
        help="Write crop/window mapping metadata without running ESM2.",
    )
    parser.add_argument(
        "--contextual-esm2-device",
        choices=("auto", "cpu", "mps", "cuda"),
        default="auto",
    )
    parser.add_argument("--contextual-esm2-batch-size", type=int, default=1)
    parser.add_argument(
        "--contextual-esm2-crop-mode", choices=("window", "idr_coverage"), default="idr_coverage",
        help="Use IDR coverage (default) or greedy grouping of scoring windows",
    )
    parser.add_argument("--contextual-esm2-context-margin", type=int, default=128)
    parser.add_argument("--motifae-checkpoint", type=Path, default=None)
    parser.add_argument(
        "--motifae-device", choices=("auto", "cpu", "mps", "cuda"), default="auto"
    )
    parser.add_argument("--motifae-chunk-residues", type=int, default=64)
    parser.add_argument("--force-v4-feature-cache", action="store_true")
    return parser.parse_args()

def main():
    global WINDOWS_FILE, PSSM_DIR, FRAGMENT_FILE, NARDINI_CACHE_FILE
    global GUSSIDR_CACHE_FILE, ESM_REGION_METADATA_FILE, CACHE_DIR, GUSSIDR_DIR, OUT_PREFIX
    args = parse_args()
    WINDOWS_FILE = args.windows.resolve() if args.windows is not None else None
    PSSM_DIR = args.pssm_dir.resolve()
    FRAGMENT_FILE = args.fragment.resolve()
    CACHE_DIR = args.cache_dir.resolve()
    GUSSIDR_DIR = args.gussidr_dir.resolve()
    NARDINI_CACHE_FILE = CACHE_DIR / "nardini_features.tsv"
    GUSSIDR_CACHE_FILE = CACHE_DIR / "gussidr_ensemble_features.tsv"
    ESM_REGION_METADATA_FILE = CACHE_DIR / "esm2_region_metadata.tsv"
    OUT_PREFIX = args.out_prefix.resolve()

    print("Loading windows...")
    windows = load_windows()
    if args.limit is not None:
        windows = windows.iloc[: args.limit].copy()
    print(f"Total windows: {len(windows)}")

    if args.contextual_esm2_crop_mode == "idr_coverage":
        contextual_metadata = CACHE_DIR / "contextual_esm2_idr_coverage_metadata.tsv"
        contextual_chunk_metadata = CACHE_DIR / "contextual_esm2_idr_coverage_chunks.tsv"
        contextual_dir = CACHE_DIR / "contextual_esm2_idr_coverage"
        motifae_sparse_dir = CACHE_DIR / "motifae_sparse_windows_idr_coverage"
    else:
        contextual_metadata = CACHE_DIR / "contextual_esm2_sequence_aware_metadata.tsv"
        contextual_chunk_metadata = CACHE_DIR / "contextual_esm2_sequence_aware_chunks.tsv"
        contextual_dir = CACHE_DIR / "contextual_esm2"
        motifae_sparse_dir = CACHE_DIR / "motifae_sparse_windows_sequence_aware"
    if args.prepare_contextual_esm2:
        print(f"Preparing contextual ESM2 cache mode={args.contextual_esm2_crop_mode}...")
        prepare_contextual_esm2(
            windows=windows,
            pssm_dir=PSSM_DIR,
            output_dir=contextual_dir,
            metadata_path=contextual_metadata,
            chunk_metadata_path=contextual_chunk_metadata,
            device=args.contextual_esm2_device,
            batch_size=args.contextual_esm2_batch_size,
            force=args.force_v4_feature_cache,
            plan_only=args.contextual_esm2_plan_only,
            crop_mode=args.contextual_esm2_crop_mode,
            regions_path=FRAGMENT_FILE,
            context_margin=args.contextual_esm2_context_margin,
        )
    if args.prepare_motifae_sparse:
        if args.contextual_esm2_plan_only:
            raise ValueError("MotifAE features cannot be built from a plan-only ESM cache")
        if args.motifae_checkpoint is None:
            raise ValueError("--prepare-motifae-sparse requires --motifae-checkpoint")
        print("Preparing frozen MotifAE sparse window features...")
        build_motifae_sparse_windows(
            context_metadata_path=contextual_metadata,
            chunk_metadata_path=contextual_chunk_metadata,
            checkpoint=args.motifae_checkpoint.resolve(),
            output_dir=motifae_sparse_dir,
            device_name=args.motifae_device,
            chunk_residues=args.motifae_chunk_residues,
            force=args.force_v4_feature_cache,
        )
    if args.v4_feature_only:
        if not (args.prepare_contextual_esm2 or args.prepare_motifae_sparse):
            raise ValueError(
                "--v4-feature-only requires --prepare-contextual-esm2 and/or "
                "--prepare-motifae-sparse"
            )
        print("[OK] Offline ESM2/MotifAE feature preparation completed")
        return

    esm_python = args.esm_python.resolve() if args.esm_python else default_esm_python()
    nardini_python = args.nardini_python.resolve() if args.nardini_python else default_nardini_python()
    idr_manifest = validate_required_inputs(windows, esm_python, nardini_python)

    if args.reuse_derived_cache:
        validate_reused_derived_caches(idr_manifest)
        print("Reusing sequence-validated derived caches.")

    if args.preflight_only:
        print("[OK] Primary inputs, model assets, and runtime dependencies are present")
        return

    if not args.reuse_derived_cache:
        print("Preparing reproducible IDR-level derived caches...")

        derived = prepare_derived_features(
            windows=windows,
            pssm_dir=PSSM_DIR,
            fragment_file=FRAGMENT_FILE,
            cache_dir=CACHE_DIR,
            gussidr_dir=GUSSIDR_DIR,
            esm_python=esm_python,
            nardini_python=nardini_python,
            gussidr_python=args.gussidr_python.resolve(),
            esm_batch_size=args.esm_batch_size,
            esm_device=args.esm_device,
            nardini_num_scrambles=args.nardini_num_scrambles,
            nardini_random_seed=args.nardini_random_seed,
            nardini_workers=args.nardini_workers,
            gussidr_batch_size=args.gussidr_batch_size,
            gussidr_device=args.gussidr_device,
            force=args.force_derived_cache,
            idr_manifest_frame=idr_manifest,
        )
        NARDINI_CACHE_FILE = derived.nardini_features
        GUSSIDR_CACHE_FILE = derived.gussidr_features
        ESM_REGION_METADATA_FILE = derived.esm_region_metadata
        print(f"ESM region metadata: {ESM_REGION_METADATA_FILE}")
        print(f"NARDINI cache: {NARDINI_CACHE_FILE}")
        print(f"GussIDR cache: {GUSSIDR_CACHE_FILE}")

    if args.prepare_cache_only:
        if args.reuse_derived_cache:
            print(f"[OK] Existing derived feature caches are ready in: {CACHE_DIR}")
        else:
            print(f"[OK] Derived feature caches prepared; provenance: {derived.provenance}")
        return

    print("Building PSSM tensor...")
    pssm_X, pssm_window_X, pssm_valid_mask, pssm_meta, pssm_residue, pssm_err = build_pssm_tensor(
        windows, collect_residue_rows=args.save_pssm_residue_table
    )
    print(f"PSSM residue shape: {pssm_X.shape}; window shape: {pssm_window_X.shape}")

    print("Building Chemical tensor...")
    chem_X, chem_meta, chem_err = build_chemical_tensors(windows)
    print(f"Chemical shape: {chem_X.shape}")

    print("Building Context tensor...")
    (
        ctx_X,
        ctx_meta,
        ctx_err,
        left_domain_X,
        right_domain_X,
        idr_pattern_X,
        idr_ensemble_X,
        idr_ensemble_full_X,
    ) = build_context_tensors(windows, return_structured=True)
    print(f"Context shape: {ctx_X.shape}")

    print("Building Accessibility tensor...")
    acc_X, acc_meta, acc_err, accessibility_rg_norm_ref = build_accessibility_tensors(
        windows,
        rg_norm_ref=args.accessibility_rg_norm_ref,
        return_reference=True,
    )
    print(f"Accessibility shape: {acc_X.shape}")

    common = set(pssm_meta["window_id"].astype(str))
    for module_meta in (chem_meta, ctx_meta, acc_meta):
        common &= set(module_meta["window_id"].astype(str))
    ordered_ids = [str(value) for value in windows["window_id"] if str(value) in common]
    if not ordered_ids:
        raise RuntimeError("No windows were successfully extracted by every module")
    dropped = len(windows) - len(ordered_ids)
    if dropped:
        print(f"Warning: retaining {len(ordered_ids)} common windows; dropping {dropped}")

    pssm_X = aligned_block(pssm_X, pssm_meta, ordered_ids)
    pssm_window_X = aligned_block(pssm_window_X, pssm_meta, ordered_ids)
    pssm_valid_mask = aligned_block(pssm_valid_mask, pssm_meta, ordered_ids)
    chem_X = aligned_block(chem_X, chem_meta, ordered_ids)
    ctx_X = aligned_block(ctx_X, ctx_meta, ordered_ids)
    left_domain_X = aligned_block(left_domain_X, ctx_meta, ordered_ids)
    right_domain_X = aligned_block(right_domain_X, ctx_meta, ordered_ids)
    idr_pattern_X = aligned_block(idr_pattern_X, ctx_meta, ordered_ids)
    idr_ensemble_X = aligned_block(idr_ensemble_X, ctx_meta, ordered_ids)
    idr_ensemble_full_X = aligned_block(idr_ensemble_full_X, ctx_meta, ordered_ids)
    acc_X = aligned_block(acc_X, acc_meta, ordered_ids)
    pssm_flat = pssm_X.reshape(len(ordered_ids), -1)
    combined_X = np.concatenate([pssm_flat, pssm_window_X, chem_X, ctx_X, acc_X], axis=1)
    print(f"Combined X shape: {combined_X.shape}")

    meta = windows.set_index("window_id").loc[ordered_ids].reset_index()
    peak_by_id = pssm_meta.set_index("window_id")["local_conservation_peak"]
    meta["local_conservation_peak"] = [peak_by_id.loc[window_id] for window_id in ordered_ids]
    meta.insert(0, "tensor_index", np.arange(len(meta), dtype=int))
    y = np.array([label_to_y(value) for value in meta["label"]], dtype=np.float32)
    if not np.isfinite(y).all():
        raise ValueError(f"Could not parse {int((~np.isfinite(y)).sum())} labels")

    all_errors = pd.concat([pssm_err, chem_err, ctx_err, acc_err], ignore_index=True)

    out_prefix = Path(OUT_PREFIX)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    categories = load_categories(CHEMICAL_CATEGORIES_JSON)
    chem_names = chemical_feature_names(categories, include_motif=bool(load_regexes(MOTIF_REGEX_JSON)))
    context_names = context_feature_names()
    if chem_X.shape[1] != len(chem_names) or ctx_X.shape[1] != len(context_names):
        raise RuntimeError("Feature-name manifest does not match extracted tensor dimensions")

    slices = {}
    offset = 0
    for name, width in (
        ("pssm_residue_flat", pssm_flat.shape[1]),
        ("pssm_window", pssm_window_X.shape[1]),
        ("chemical", chem_X.shape[1]),
        ("context", ctx_X.shape[1]),
        ("accessibility", acc_X.shape[1]),
    ):
        slices[name] = [offset, offset + width]
        offset += width
    manifest = {
        "combined_shape": [len(ordered_ids), int(combined_X.shape[1])],
        "combined_order": list(slices),
        "combined_slices_half_open": slices,
        "pssm_residue_shape": [29, len(PSSM_RESIDUE_FEATURE_NAMES)],
        "pssm_residue_features": list(PSSM_RESIDUE_FEATURE_NAMES),
        "pssm_valid_mask": "separate (N,29) uint8 array; not a biological feature channel",
        "pssm_window_features": list(PSSM_WINDOW_FEATURE_NAMES),
        "chemical_features": chem_names,
        "context_features": context_names,
        "structured_context": {
            "left_domain_esm": [1280],
            "right_domain_esm": [1280],
            "idr_pattern": [36],
            "idr_ensemble": [32],
            "idr_ensemble_full": [38],
            "idr_ensemble_full_features": [f"feat_{index}" for index in range(1, 39)],
            "ensemble_dropped_original_zero_based_indices": list(ENSEMBLE_DROP_INDICES),
        },
        "accessibility_features": list(ACCESSIBILITY_FEATURE_NAMES),
        "accessibility_rg_norm_ref": accessibility_rg_norm_ref,
        "accessibility_rg_norm_reference_scope": (
            "supplied training reference"
            if args.accessibility_rg_norm_ref is not None
            else "median of unique training IDRs"
        ),
    }

    np.savez_compressed(
        out_prefix.with_name(out_prefix.name + "_tensor.npz"),
        X=combined_X.astype(np.float32),
        X_pssm_residue=pssm_X.astype(np.float32),
        X_pssm_window=pssm_window_X.astype(np.float32),
        X_pssm_valid_mask=pssm_valid_mask.astype(np.uint8),
        X_chemical=chem_X.astype(np.float32),
        X_context=ctx_X.astype(np.float32),
        X_accessibility=acc_X.astype(np.float32),
        pssm_residue=pssm_X.astype(np.float32),
        pssm_valid_mask=pssm_valid_mask.astype(np.uint8),
        pssm_window=pssm_window_X.astype(np.float32),
        chemical=chem_X.astype(np.float32),
        left_domain_esm=left_domain_X.astype(np.float32),
        right_domain_esm=right_domain_X.astype(np.float32),
        idr_pattern=idr_pattern_X.astype(np.float32),
        idr_ensemble=idr_ensemble_X.astype(np.float32),
        idr_ensemble_full=idr_ensemble_full_X.astype(np.float32),
        accessibility=acc_X.astype(np.float32),
        accessibility_rg_norm_ref=np.asarray(accessibility_rg_norm_ref, dtype=np.float32),
        y=y,
        window_ids=np.asarray(ordered_ids, dtype=str),
        uniprots=np.asarray(meta['uniprot'].astype(str).tolist(), dtype=str),
        labels=np.asarray(meta['label'].astype(str).tolist(), dtype=str),
        pssm_residue_feature_names=np.asarray(PSSM_RESIDUE_FEATURE_NAMES, dtype=str),
        pssm_window_feature_names=np.asarray(PSSM_WINDOW_FEATURE_NAMES, dtype=str),
    )
    write_table(meta, out_prefix.with_name(out_prefix.name + "_metadata.tsv"))
    write_table(all_errors, out_prefix.with_name(out_prefix.name + "_errors.tsv"))
    if args.save_pssm_residue_table:
        write_table(pssm_residue, out_prefix.with_name(out_prefix.name + "_pssm_residue_features.tsv.gz"))
    out_prefix.with_name(out_prefix.name + "_feature_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    print(f"[OK] Combined tensor saved to {out_prefix.with_name(out_prefix.name + '_tensor.npz')}")
    print(f"[OK] Metadata saved to {out_prefix.with_name(out_prefix.name + '_metadata.tsv')}")
    print(f"[OK] Manifest saved to {out_prefix.with_name(out_prefix.name + '_feature_manifest.json')}")
    print(f"[OK] Errors saved to {out_prefix.with_name(out_prefix.name + '_errors.tsv')}")
    if len(all_errors):
        print(f"Total errors: {len(all_errors)}")

if __name__ == "__main__":
    main()
