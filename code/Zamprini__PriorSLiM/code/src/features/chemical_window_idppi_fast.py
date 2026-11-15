"""Calculate local chemical and regular-expression features."""
                      

from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
KD_UVERSKY_0_1 = {
    "I": 1.000000,
    "V": 0.966667,
    "L": 0.922222,
    "F": 0.811111,
    "C": 0.777778,
    "M": 0.711111,
    "A": 0.700000,
    "G": 0.455556,
    "T": 0.422222,
    "S": 0.411111,
    "W": 0.400000,
    "Y": 0.355556,
    "P": 0.322222,
    "H": 0.144444,
    "E": 0.111111,
    "Q": 0.111111,
    "D": 0.111111,
    "N": 0.111111,
    "K": 0.066667,
    "R": 0.000000,
}
PPII_HILSER = {
    "I": 0.39,
    "V": 0.39,
    "L": 0.24,
    "F": 0.17,
    "C": 0.25,
    "M": 0.36,
    "A": 0.37,
    "G": 0.13,
    "T": 0.32,
    "S": 0.24,
    "W": 0.25,
    "Y": 0.25,
    "P": 1.00,
    "H": 0.20,
    "E": 0.42,
    "Q": 0.53,
    "D": 0.30,
    "N": 0.27,
    "K": 0.56,
    "R": 0.38,
}
DEFAULT_CATEGORIES = {
    "positive": "RKH",
    "negative": "DE",
    "charged": "RKHDE",
    "aromatic": "FYW",
    "proline": "P",
    "glycine": "G",
    "polar": "STNQCY",
    "hydrophobic": "AILMFWVY",
    "aliphatic": "AILV",
    "small": "AGSTCPDNV",
    "phospho_acceptor": "STY",
}

@dataclass(frozen=True)
class Segment:
    protein_id: str
    start: int
    end: int
    label: str = ""
    source: str = ""
    @property
    def length(self) -> int:
        return self.end - self.start + 1

def parse_args():
    p = argparse.ArgumentParser(description="Fast chemical-window features for IDPPI/IDR SLiM analysis")
    p.add_argument("--idppi", required=True)
    p.add_argument("--regions", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--window", type=int, default=29)
    p.add_argument("--patch", type=int, default=5)
    p.add_argument("--top_patch_k", type=int, default=3)
    p.add_argument("--center_len", type=int, default=15)
    p.add_argument("--alpha", type=float, default=0.5)
    p.add_argument("--eps", type=float, default=1e-8)
    p.add_argument("--mad_floor", type=float, default=0.05)
    p.add_argument("--clip_rel", type=float, default=5.0)
    p.add_argument("--null_shuffle", type=int, default=20000)
    p.add_argument("--seed", type=int, default=13)
    p.add_argument("--categories_json", default=None)
    p.add_argument("--motif_regex_json", default=None)
    p.add_argument("--positive_state_regex", default="remains_disordered")
    p.add_argument("--negative_scope", choices=["all_nonoverlap_idr", "same_protein_nonoverlap_idr"], default="all_nonoverlap_idr")
    p.add_argument("--label_window_by", choices=["center", "overlap"], default="center")
    p.add_argument("--min_overlap_frac", type=float, default=0.5)
    p.add_argument("--make_plots", action="store_true")
    p.add_argument("--max_plot_features", type=int, default=20)
    return p.parse_args()

def parse_fasta(path: str) -> Dict[str, str]:
    seqs, cur = {}, None
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                token = line[1:].split()[0]
                if "|" in token:
                    parts = token.split("|")
                    if len(parts) > 1 and parts[1]:
                        token = parts[1]
                cur = token
                seqs[cur] = []
            else:
                if cur is None:
                    raise ValueError("FASTA sequence before header")
                seqs[cur].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}

def parse_region_string(x: Any) -> List[Tuple[int, int]]:
    if pd.isna(x):
        return []
    out = []
    for part in re.split(r"[,;，；]\s*", str(x)):
        m = re.search(r"(\d+)\s*[-–—]\s*(\d+)", part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            out.append((min(a, b), max(a, b)))
    return out

def load_regions(path: str) -> Tuple[List[Segment], Dict[str, List[Segment]]]:
    df = pd.read_excel(path)
    uid_col = next(
        (name for name in ("Protein_Accession", "UniProt_ID", "uniprot") if name in df.columns),
        None,
    )
    disorder_col = next(
        (
            name
            for name in ("Final_Disordered_Intervals", "Disorder_Regions")
            if name in df.columns
        ),
        None,
    )
    ordered_col = next(
        (
            name
            for name in (
                "Final_Non_Disordered_Intervals",
                "Final_Ordered_Intervals",
                "Ordered_Regions",
            )
            if name in df.columns
        ),
        None,
    )
    if uid_col is None or disorder_col is None:
        raise ValueError(
            "Region file needs a UniProt identifier and disorder intervals. "
            f"Found {list(df.columns)}"
        )
    idrs = []
    ordered_by_pid: Dict[str, List[Segment]] = defaultdict(list)
    for _, r in df.iterrows():
        pid = str(r[uid_col]).strip()
        for s, e in parse_region_string(r[disorder_col]):
            idrs.append(Segment(pid, s, e, "all_idr", "regions"))
        if ordered_col is not None:
            for s, e in parse_region_string(r[ordered_col]):
                ordered_by_pid[pid].append(Segment(pid, s, e, "ordered", "regions"))
    return idrs, ordered_by_pid

def load_idrs(path: str) -> List[Segment]:
    idrs, _ordered_by_pid = load_regions(path)
    return idrs

def parse_protein_segment(x: Any) -> Optional[Tuple[str, int, int]]:
    if pd.isna(x):
        return None
    m = re.search(r"^\s*([^（(\s]+)\s*[（(]\s*(\d+)\s*[-–—]\s*(\d+)\s*[）)]", str(x))
    if not m:
        return None
    pid = m.group(1).strip()
    s, e = int(m.group(2)), int(m.group(3))
    return pid, min(s, e), max(s, e)

def load_positive_segments(path: str, state_regex: str) -> List[Segment]:
    df = pd.read_csv(path, sep="\t")
    rx = re.compile(state_regex, re.I)
    segs = []
    for idx, r in df.iterrows():
        for side in ("1", "2"):
            if str(r.get(f"F{side}", "")).strip().upper() == "IDR" and rx.search(str(r.get(f"P{side}_IDR_complex_state", ""))):
                p = parse_protein_segment(r.get(f"P{side}"))
                if p:
                    segs.append(Segment(p[0], p[1], p[2], "positive", f"idppi_row_{idx}_P{side}"))
    seen, uniq = set(), []
    for s in segs:
        key = (s.protein_id, s.start, s.end)
        if key not in seen:
            seen.add(key); uniq.append(s)
    return uniq

def overlap_len(a: int, b: int, c: int, d: int) -> int:
    return max(0, min(b, d) - max(a, c) + 1)

def overlaps(a: Segment, b: Segment) -> bool:
    return a.protein_id == b.protein_id and overlap_len(a.start, a.end, b.start, b.end) > 0

def parent_idr(seg: Segment, by_pid: Dict[str, List[Segment]]) -> Optional[Segment]:
    cands = by_pid.get(seg.protein_id, [])
    contain = [x for x in cands if x.start <= seg.start and x.end >= seg.end]
    if contain:
        return min(contain, key=lambda x: x.length)
    best, best_ov = None, 0
    for x in cands:
        ov = overlap_len(seg.start, seg.end, x.start, x.end)
        if ov > best_ov:
            best, best_ov = x, ov
    return best if best_ov > 0 else None

def load_categories(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return DEFAULT_CATEGORIES.copy()
    obj = json.load(open(path, encoding="utf-8"))
    out = {}
    for k, v in obj.items():
        letters = v if isinstance(v, str) else "".join(v)
        letters = "".join([aa for aa in str(letters).upper() if aa in AA20])
        if letters:
            out[str(k)] = letters
    if not out:
        raise ValueError("No valid chemical categories found")
    return out

def load_regexes(path: Optional[str]) -> Dict[str, List[re.Pattern]]:
    if not path:
        return {}
    obj = json.load(open(path, encoding="utf-8"))
    out = {}
    for fam, pats in obj.items():
        if isinstance(pats, str):
            pats = [pats]
        out[str(fam)] = [re.compile(str(p)) for p in pats]
    return out

def has_unescaped_anchor(pattern: str, anchor: str) -> bool:
    pass

    escaped = False
    in_class = False
    for ch in pattern:
        if escaped:
            escaped = False
            continue
        if ch == "\\":
            escaped = True
            continue
        if ch == "[":
            in_class = True
            continue
        if ch == "]" and in_class:
            in_class = False
            continue
        if ch == anchor and not in_class:
            return True
    return False

def idr_terminal_flags(idr: Segment, ordered: List[Segment]) -> Tuple[bool, bool]:
    pass

    has_left_ordered = any(dom.end < idr.start for dom in ordered)
    has_right_ordered = any(dom.start > idr.end for dom in ordered)
    return (not has_left_ordered), (not has_right_ordered)

def sliding_sum(arr: np.ndarray, k: int) -> np.ndarray:
    arr = arr.astype(float, copy=False)
    cs = np.concatenate([[0.0], np.cumsum(arr)])
    return cs[k:] - cs[:-k]

def sliding_windows_2d_1d(arr: np.ndarray, width: int) -> np.ndarray:

    try:
        return np.lib.stride_tricks.sliding_window_view(arr, width)
    except AttributeError:
        n = len(arr) - width + 1
        return np.vstack([arr[i:i+width] for i in range(n)])

def topk_mean_axis1(mat: np.ndarray, k: int) -> np.ndarray:
    k = max(1, min(k, mat.shape[1]))
    part = np.partition(mat, mat.shape[1] - k, axis=1)[:, -k:]
    return part.mean(axis=1)

def local_patch_strength_for_all_windows(binary_idr: np.ndarray, K: int, patch: int, top_patch_k: int) -> np.ndarray:

    patch_counts = sliding_sum(binary_idr, patch) / float(patch)                    
    n_patches_inside_window = K - patch + 1
    patch_matrix = sliding_windows_2d_1d(patch_counts, n_patches_inside_window)                   
    return topk_mean_axis1(patch_matrix, top_patch_k)

@lru_cache(maxsize=8)
def build_nardini_lookup(K: int, patch: int, top_patch_k: int, B: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    recs = []
    n_patches = K - patch + 1
    for n in range(K + 1):
        vals = np.empty(B, dtype=float)
        for b in range(B):
            arr = np.zeros(K, dtype=np.int8)
            if n > 0:
                arr[rng.choice(K, size=n, replace=False)] = 1
            patch_counts = sliding_sum(arr, patch) / float(patch)
            vals[b] = np.sort(patch_counts)[-min(top_patch_k, n_patches):].mean()
        recs.append({"n": n, "mu_null": vals.mean(), "sd_null": vals.std(ddof=1) if B > 1 else 0.0})
    return pd.DataFrame(recs)

def center_flank_masks(K: int, center_len: int) -> Tuple[np.ndarray, np.ndarray]:
    if center_len % 2 == 0:
        raise ValueError("center_len should be odd")
    if center_len >= K:
        raise ValueError("center_len must be < window")
    center = K // 2
    half = center_len // 2
    cm = np.zeros(K, dtype=float)
    cm[center-half:center+half+1] = 1.0 / center_len
    fm = np.ones(K, dtype=float)
    fm[center-half:center+half+1] = 0.0
    fm = fm / fm.sum()
    return cm, fm

def relative_within_idr(df: pd.DataFrame, cols: List[str], mad_floor: float, eps: float, clip_rel: float) -> pd.DataFrame:
    for col in cols:
        vals = df[col].to_numpy(dtype=float)
        ok = np.isfinite(vals)
        if ok.sum() == 0:
            df[f"{col}_rel_IDR"] = np.nan
            df[f"{col}_percentile_IDR"] = np.nan
            continue
        med = float(np.nanmedian(vals))
        mad = float(np.nanmedian(np.abs(vals[ok] - med)))
        scale = max(1.4826 * mad, mad_floor)
        rel = (vals - med) / (scale + eps)
        if clip_rel and clip_rel > 0:
            rel = np.clip(rel, -clip_rel, clip_rel)
        sorted_vals = np.sort(vals[ok])
        pct = np.full_like(vals, np.nan, dtype=float)
        pct[ok] = np.searchsorted(sorted_vals, vals[ok], side="right") / ok.sum()
        df[f"{col}_rel_IDR"] = rel
        df[f"{col}_percentile_IDR"] = pct
    return df

def motif_core_mask(pat: re.Pattern, matched: str, cache: Dict[Tuple[str, str], np.ndarray]) -> np.ndarray:
    pass

    key = (pat.pattern, matched)
    if key in cache:
        return cache[key]
    out = np.zeros(len(matched), dtype=bool)
    chars = list(matched)
    for i, original in enumerate(chars):
        allowed = 0
        for aa in AA20:
            chars[i] = aa
            if pat.fullmatch("".join(chars)):
                allowed += 1
        chars[i] = original
        out[i] = allowed < len(AA20)
    cache[key] = out
    return out

def regex_features_for_idr(
    seq: str,
    K: int,
    regexes: Dict[str, List[re.Pattern]],
    is_n_terminal_idr: bool,
    is_c_terminal_idr: bool,
    center_len: int = 15,
) -> Dict[str, np.ndarray]:
    T = len(seq) - K + 1
    out: Dict[str, np.ndarray] = {}
    if not regexes:
        out["motif_min_abs_offset"] = np.full(T, np.nan)
        out["motif_center_enrichment"] = np.zeros(T, dtype=float)
        out["motif_count_log"] = np.zeros(T, dtype=float)
        return out
    center = (K + 1) / 2.0
    denom = (K - 1) / 2.0
    total = np.zeros(T, dtype=np.int16)
    best_abs = np.full(T, np.inf)
    core_covered = np.zeros(len(seq), dtype=np.int8)
    core_cache: Dict[Tuple[str, str], np.ndarray] = {}
    for _fam, pats in regexes.items():
        for pat in pats:
            needs_n_terminal = has_unescaped_anchor(pat.pattern, "^")
            needs_c_terminal = has_unescaped_anchor(pat.pattern, "$")
            if needs_n_terminal and not is_n_terminal_idr:
                continue
            if needs_c_terminal and not is_c_terminal_idr:
                continue
            look = re.compile(r"(?=(" + pat.pattern + r"))")
            for m in look.finditer(seq):
                s0 = m.start(1)
                matched = m.group(1)
                if matched:
                    e0 = s0 + len(matched) - 1
                    lo = max(0, e0 - K + 1)
                    hi = min(T - 1, s0)
                    if lo > hi:
                        continue
                    core_mask = motif_core_mask(pat, matched, core_cache)
                    core_idx = np.flatnonzero(core_mask)
                    if core_idx.size:
                        core_covered[s0 + core_idx] = 1
                        motif_center = float(np.mean(s0 + core_idx + 1))
                    else:
                        motif_center = ((s0 + 1) + (e0 + 1)) / 2.0
                    starts = np.arange(lo, hi + 1)
                    c_match_in_window = motif_center - starts
                    signed = (c_match_in_window - center) / denom
                    ab = np.abs(signed)
                    total[starts] += 1
                    update = ab < best_abs[starts]
                    idx = starts[update]
                    best_abs[idx] = ab[update]
    cmask, fmask = center_flank_masks(K, center_len=center_len)
    center_cov = np.convolve(core_covered.astype(float), cmask[::-1], mode="valid")
    flank_cov = np.convolve(core_covered.astype(float), fmask[::-1], mode="valid")
    best_abs[~np.isfinite(best_abs)] = 1.0
    out["motif_min_abs_offset"] = best_abs
    out["motif_center_enrichment"] = center_cov - flank_cov
    out["motif_count_log"] = np.log1p(total.astype(float))
    return out

def features_for_idr(
    pid: str,
    idr: Segment,
    full_seq: str,
    categories: Dict[str, str],
    lookup: Dict[int, Tuple[float, float]],
    args,
    regexes,
    ordered_by_pid: Dict[str, List[Segment]],
) -> Optional[pd.DataFrame]:
    K = args.window
    seq = full_seq[idr.start-1:idr.end].upper()
    L = len(seq)
    if L < K:
        return None
    T = L - K + 1
    starts = np.arange(idr.start, idr.start + T, dtype=np.int32)
    data: Dict[str, Any] = {
        "protein_id": np.array([pid] * T, dtype=object),
        "idr_start": np.full(T, idr.start, dtype=np.int32),
        "idr_end": np.full(T, idr.end, dtype=np.int32),
        "idr_len": np.full(T, idr.length, dtype=np.int32),
        "window_start": starts,
        "window_end": starts + K - 1,
        "window_center": starts + K // 2,
        "window_offset_in_IDR": np.arange(1, T + 1, dtype=np.int32),
    }

    data["window_seq"] = [seq[i:i+K] for i in range(T)]

    cmask, fmask = center_flank_masks(K, args.center_len)
    rel_cols = []

    for name, letters in categories.items():
        aset = set(letters)
        arr = np.fromiter((1 if aa in aset else 0 for aa in seq), dtype=np.int8, count=L)
        count_w = sliding_sum(arr, K)
        density = count_w / float(K)
        A = local_patch_strength_for_all_windows(arr, K, args.patch, args.top_patch_k)
        mu = np.array([lookup[int(n)][0] for n in count_w], dtype=float)
        sd = np.array([lookup[int(n)][1] for n in count_w], dtype=float)
        X = np.where(sd > args.eps, (A - mu) / (sd + args.eps), 0.0)
        center_d = np.convolve(arr.astype(float), cmask[::-1], mode="valid")
        flank_d = np.convolve(arr.astype(float), fmask[::-1], mode="valid")
        ce = center_d - flank_d
        data[f"density_{name}"] = density
        data[f"count_{name}"] = count_w.astype(np.int16)
        data[f"A_top{args.top_patch_k}_patch_{name}"] = A
        data[f"X_{name}"] = X
        data[f"center_enrich_{name}"] = ce
        rel_cols.extend([f"density_{name}", f"X_{name}", f"center_enrich_{name}"])

    def aa_count(aa: str) -> np.ndarray:
        arr = np.fromiter((1 if x == aa else 0 for x in seq), dtype=np.int8, count=L)
        return sliding_sum(arr, K)
    nR, nK, nH, nD, nE = [aa_count(a) for a in "RKHDE"]
    pos = nR + nK + nH
    neg = nD + nE
    data["NCPR"] = (pos - neg) / float(K)
    data["FCR"] = (pos + neg) / float(K)
    data["RK_logratio"] = np.log((nR + args.alpha) / (nK + args.alpha))
    data["ED_logratio"] = np.log((nE + args.alpha) / (nD + args.alpha))
    data["n_R"], data["n_K"], data["n_H"], data["n_D"], data["n_E"] = nR, nK, nH, nD, nE
    rel_cols.extend(["NCPR", "FCR", "RK_logratio", "ED_logratio"])

                                                                             

                                                                           

                                           
    kd_arr = np.fromiter((KD_UVERSKY_0_1.get(aa, 0.0) for aa in seq), dtype=float, count=L)
    ppii_arr = np.fromiter((PPII_HILSER.get(aa, 0.0) for aa in seq), dtype=float, count=L)
    data["hydropathy_sum"] = sliding_sum(kd_arr, K)
    data["PPII_sum"] = sliding_sum(ppii_arr, K)

    is_n_terminal_idr, is_c_terminal_idr = idr_terminal_flags(idr, ordered_by_pid.get(pid, []))
    data.update(regex_features_for_idr(seq, K, regexes, is_n_terminal_idr, is_c_terminal_idr, args.center_len))

    df = pd.DataFrame(data)
    df = relative_within_idr(df, rel_cols, args.mad_floor, args.eps, args.clip_rel)
    return df

def label_windows(windows: pd.DataFrame, pos_by_pid: Dict[str, List[Segment]], mode: str, min_overlap_frac: float) -> pd.Series:
    labels = np.array(["other"] * len(windows), dtype=object)
    for pid, idx in windows.groupby("protein_id", sort=False).groups.items():
        segs = pos_by_pid.get(pid, [])
        if not segs:
            continue
        sub = windows.loc[idx]
        positive = np.zeros(len(sub), dtype=bool)
        centers = sub["window_center"].to_numpy()
        starts = sub["window_start"].to_numpy()
        ends = sub["window_end"].to_numpy()
        K = ends - starts + 1
        for seg in segs:
            if mode == "center":
                positive |= (centers >= seg.start) & (centers <= seg.end)
            else:
                ov = np.maximum(0, np.minimum(ends, seg.end) - np.maximum(starts, seg.start) + 1)
                positive |= (ov / K) >= min_overlap_frac
        labels[np.asarray(idx)[positive]] = "positive"
    return pd.Series(labels, index=windows.index)

def aggregate_to_segments(windows: pd.DataFrame, segments: List[Segment], all_idrs_by_pid: Dict[str, List[Segment]], numeric_cols: List[str], label: str) -> pd.DataFrame:
    recs = []
    by_pid = {p: sub for p, sub in windows.groupby("protein_id", sort=False)}
    for seg in segments:
        par = parent_idr(seg, all_idrs_by_pid)
        if par is None or seg.protein_id not in by_pid:
            continue
        sub = by_pid[seg.protein_id]
        m = (sub["idr_start"] == par.start) & (sub["idr_end"] == par.end) & (sub["window_center"] >= seg.start) & (sub["window_center"] <= seg.end)
        w = sub.loc[m]
        if w.empty:
            sp = sub.loc[(sub["idr_start"] == par.start) & (sub["idr_end"] == par.end)].copy()
            if sp.empty:
                continue
            center = (seg.start + seg.end) / 2.0
            w = sp.iloc[[int(np.argmin(np.abs(sp["window_center"].to_numpy() - center)))]]
        r = {
            "protein_id": seg.protein_id,
            "segment_start": seg.start,
            "segment_end": seg.end,
            "segment_len": seg.length,
            "parent_idr_start": par.start,
            "parent_idr_end": par.end,
            "group": label,
            "n_windows": len(w),
        }
        for col in numeric_cols:
            vals = pd.to_numeric(w[col], errors="coerce").to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                continue
            r[f"{col}_mean"] = vals.mean()
            r[f"{col}_median"] = np.median(vals)
            r[f"{col}_max"] = vals.max()
            for k in (3, 5):
                kk = min(k, len(vals))
                r[f"{col}_top{k}_mean"] = np.sort(vals)[-kk:].mean()
            r[f"{col}_peak_count_gt1"] = int((vals > 1).sum())
            r[f"{col}_peak_count_gt2"] = int((vals > 2).sum())
            if "percentile_IDR" in col:
                r[f"{col}_peak_count_gt0.8"] = int((vals > 0.8).sum())
                r[f"{col}_peak_count_gt0.9"] = int((vals > 0.9).sum())
                r[f"{col}_peak_count_gt0.95"] = int((vals > 0.95).sum())
        recs.append(r)
    return pd.DataFrame(recs)

def mann_auc(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if len(x) < 1 or len(y) < 1:
        return np.nan, np.nan, np.nan
    try:
        from scipy.stats import mannwhitneyu
        res = mannwhitneyu(x, y, alternative="two-sided")
        U, p = float(res.statistic), float(res.pvalue)
    except Exception:
        vals = np.r_[x, y]
        ranks = pd.Series(vals).rank(method="average").to_numpy()
        U = ranks[:len(x)].sum() - len(x) * (len(x) + 1) / 2.0
        p = np.nan
    auc = U / (len(x) * len(y))
    return p, auc, 2 * auc - 1

def stats(seg_df: pd.DataFrame) -> pd.DataFrame:
    pos = seg_df[seg_df["group"] == "positive"]
    neg = seg_df[seg_df["group"] == "other_idr"]
    recs = []
    for col in seg_df.columns:
        if col in {"protein_id", "group"} or not pd.api.types.is_numeric_dtype(seg_df[col]):
            continue
        xp = pd.to_numeric(pos[col], errors="coerce").dropna().to_numpy(float)
        xn = pd.to_numeric(neg[col], errors="coerce").dropna().to_numpy(float)
        if len(xp) < 2 or len(xn) < 2:
            continue
        p, auc, cliff = mann_auc(xp, xn)
        recs.append({
            "feature": col,
            "n_positive": len(xp), "n_other": len(xn),
            "positive_mean": xp.mean(), "other_mean": xn.mean(),
            "positive_median": np.median(xp), "other_median": np.median(xn),
            "delta_mean_pos_minus_other": xp.mean() - xn.mean(),
            "mannwhitney_p": p,
            "AUC_pos_gt_other": auc,
            "cliffs_delta_pos_vs_other": cliff,
        })
    out = pd.DataFrame(recs)
    if not out.empty:
        out = out.sort_values(["mannwhitney_p", "AUC_pos_gt_other"], ascending=[True, False])
    return out

def make_plots(seg_df: pd.DataFrame, stats_df: pd.DataFrame, out_dir: Path, max_features: int):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    pdir = out_dir / "plots"; pdir.mkdir(exist_ok=True)
    for feat in stats_df.head(max_features)["feature"]:
        a = pd.to_numeric(seg_df.loc[seg_df.group == "positive", feat], errors="coerce").dropna().to_numpy()
        b = pd.to_numeric(seg_df.loc[seg_df.group == "other_idr", feat], errors="coerce").dropna().to_numpy()
        if len(a) == 0 or len(b) == 0:
            continue
        plt.figure(figsize=(5,4))
        plt.boxplot([a, b], labels=["positive", "other_idr"], showfliers=False)
        plt.title(feat); plt.tight_layout()
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", feat)[:160]
        plt.savefig(pdir / f"box_{safe}.png", dpi=200); plt.close()

def main():
    args = parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    cats = load_categories(args.categories_json)
    regexes = load_regexes(args.motif_regex_json)
    json.dump(cats, open(out_dir / "chemical_categories_used.json", "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    seqs = parse_fasta(args.fasta)
    idrs, ordered_by_pid = load_regions(args.regions)
    positives = load_positive_segments(args.idppi, args.positive_state_regex)
    idrs_by_pid = defaultdict(list)
    for x in idrs: idrs_by_pid[x.protein_id].append(x)
    pos_by_pid = defaultdict(list)
    for x in positives: pos_by_pid[x.protein_id].append(x)

    if args.negative_scope == "same_protein_nonoverlap_idr":
        pids = set(pos_by_pid)
        idr_candidates = [x for x in idrs if x.protein_id in pids]
    else:
        idr_candidates = idrs
    negatives = [x for x in idr_candidates if not any(overlaps(x, p) for p in pos_by_pid.get(x.protein_id, []))]

    lookup_df = build_nardini_lookup(args.window, args.patch, args.top_patch_k, args.null_shuffle, args.seed)
    lookup_df.to_csv(out_dir / "nardini_X_null_lookup.tsv", sep="\t", index=False)
    lookup = {int(r.n): (float(r.mu_null), float(r.sd_null)) for r in lookup_df.itertuples(index=False)}

    missing = []
    frames = []
    for i, idr in enumerate(idrs, 1):
        seq = seqs.get(idr.protein_id)
        if seq is None:
            missing.append({"protein_id": idr.protein_id, "start": idr.start, "end": idr.end, "reason": "missing_fasta"}); continue
        if idr.end > len(seq):
            missing.append({"protein_id": idr.protein_id, "start": idr.start, "end": idr.end, "seq_len": len(seq), "reason": "idr_end_exceeds_sequence"}); continue
        if idr.length < args.window:
            missing.append({"protein_id": idr.protein_id, "start": idr.start, "end": idr.end, "idr_len": idr.length, "reason": f"shorter_than_window_{args.window}"}); continue
        df = features_for_idr(idr.protein_id, idr, seq, cats, lookup, args, regexes, ordered_by_pid)
        if df is not None:
            frames.append(df)
        if i % 1000 == 0:
            print(f"[INFO] processed IDR segments: {i}/{len(idrs)}", file=sys.stderr)

    if not frames:
        pd.DataFrame(missing).to_csv(out_dir / "missing_inputs.tsv", sep="\t", index=False)
        raise SystemExit("No windows generated")
    windows = pd.concat(frames, ignore_index=True)
    windows["window_label"] = label_windows(windows, pos_by_pid, args.label_window_by, args.min_overlap_frac)
    windows.to_csv(out_dir / "chemical_window_level.tsv", sep="\t", index=False)

    exclude = {"idr_start", "idr_end", "idr_len", "window_start", "window_end", "window_center", "window_offset_in_IDR"}
    numeric_cols = [c for c in windows.columns if c not in exclude and pd.api.types.is_numeric_dtype(windows[c])]
    pos_seg = aggregate_to_segments(windows, positives, idrs_by_pid, numeric_cols, "positive")
    neg_seg = aggregate_to_segments(windows, negatives, idrs_by_pid, numeric_cols, "other_idr")
    seg_df = pd.concat([pos_seg, neg_seg], ignore_index=True)
    seg_df.to_csv(out_dir / "chemical_segment_level.tsv", sep="\t", index=False)

    st = stats(seg_df)
    st.to_csv(out_dir / "chemical_group_stats.tsv", sep="\t", index=False)
    pd.DataFrame(missing).to_csv(out_dir / "missing_inputs.tsv", sep="\t", index=False)

    summary = {
        "n_sequences": len(seqs), "n_idr_segments": len(idrs),
        "n_positive_segments": len(positives), "n_other_segments": len(negatives),
        "n_windows": int(len(windows)),
        "n_positive_windows": int((windows.window_label == "positive").sum()),
        "window": args.window, "patch": args.patch, "top_patch_k": args.top_patch_k,
        "center_len": args.center_len, "categories": cats,
        "motif_regex_enabled": bool(regexes), "n_motif_families": len(regexes),
    }
    json.dump(summary, open(out_dir / "run_summary.json", "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    if args.make_plots:
        make_plots(seg_df, st, out_dir, args.max_plot_features)
    print("[DONE]", file=sys.stderr)
    print(f"chemical_window_level.tsv: {len(windows)} rows")
    print(f"chemical_segment_level.tsv: {len(seg_df)} rows")
    if not st.empty:
        print(st[["feature", "mannwhitney_p", "AUC_pos_gt_other", "positive_mean", "other_mean"]].head(15).to_string(index=False))

if __name__ == "__main__":
    main()
