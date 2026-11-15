"""Extract residue- and window-level conservation features from ASCII PSSMs."""
                      

from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_ORDER)}

                                                                              

                                                                              
PSIBLAST_ASCII_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
PSIBLAST_TO_INTERNAL = [PSIBLAST_ASCII_AA_ORDER.index(aa) for aa in AA_ORDER]
EPS = 1e-8

        

def parse_fasta(path: str | Path) -> Dict[str, str]:
    pass
    seqs: Dict[str, str] = {}
    cur_id: Optional[str] = None
    cur_seq: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq).upper()
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(re.sub(r"\s+", "", line))
        if cur_id is not None:
            seqs[cur_id] = "".join(cur_seq).upper()
    return seqs

def parse_region_string(x) -> List[Tuple[int, int]]:
    pass
    if pd.isna(x):
        return []
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return []
    out: List[Tuple[int, int]] = []
    for part in re.split(r"[,;，；]\s*", s):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^(\d+)\s*[-–—]\s*(\d+)$", part)
        if not m:
            continue
        a, b = int(m.group(1)), int(m.group(2))
        if a > b:
            a, b = b, a
        out.append((a, b))
    return out

def parse_fragment_label(x: str) -> Optional[Tuple[str, int, int]]:
    pass
    if pd.isna(x):
        return None
    s = str(x).strip()

    m = re.match(r"^\s*([^（(\s]+)\s*[（(]\s*(\d+)\s*[-–—]\s*(\d+)\s*[）)]\s*$", s)
    if not m:
        return None
    uid = m.group(1)
    a, b = int(m.group(2)), int(m.group(3))
    if a > b:
        a, b = b, a
    return uid, a, b

def ranges_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    return max(a[0], b[0]) <= min(a[1], b[1])

def range_len(r: Tuple[int, int]) -> int:
    return r[1] - r[0] + 1

def clamp_range(start: int, end: int, L: int) -> Optional[Tuple[int, int]]:
    start = max(1, int(start))
    end = min(L, int(end))
    if start > end:
        return None
    return start, end

        

def load_all_idrs(region_xlsx: str | Path) -> pd.DataFrame:
    pass
    df = pd.read_excel(region_xlsx)
    required = {"UniProt_ID", "Disorder_Regions"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"region xlsx 缺少列: {missing}. 当前列: {df.columns.tolist()}")

    rows = []
    for _, row in df.iterrows():
        uid = str(row["UniProt_ID"]).strip()
        plen = row.get("Protein_Length", np.nan)
        for idx, (s, e) in enumerate(parse_region_string(row["Disorder_Regions"]), start=1):
            rows.append({
                "uniprot": uid,
                "idr_index": idx,
                "idr_start": s,
                "idr_end": e,
                "idr_len": e - s + 1,
                "protein_length": plen,
            })
    return pd.DataFrame(rows)

def extract_positive_segments(idppi_tsv: str | Path) -> pd.DataFrame:
    pass
    df = pd.read_csv(idppi_tsv, sep="\t")
    required = {"P1", "P2", "F1", "F2", "P1_IDR_complex_state", "P2_IDR_complex_state"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"IDPPI tsv 缺少列: {missing}. 当前列: {df.columns.tolist()}")

    rows = []
    for ridx, row in df.iterrows():
        for side in ["1", "2"]:
            frag = parse_fragment_label(row[f"P{side}"])
            if frag is None:
                continue
            uid, s, e = frag
            ftype = str(row.get(f"F{side}", "")).strip().lower()
            state = str(row.get(f"P{side}_IDR_complex_state", "")).strip().lower()
            is_idr = (ftype == "idr") or ("remains_disordered" in state)
            if not is_idr:
                continue
            rows.append({
                "uniprot": uid,
                "seg_start": s,
                "seg_end": e,
                "seg_len": e - s + 1,
                "source_row": ridx,
                "source_side": f"P{side}",
                "source_type": row.get(f"F{side}"),
                "complex_state": row.get(f"P{side}_IDR_complex_state"),
                "label": "positive",
            })
    pos = pd.DataFrame(rows)
    if pos.empty:
        return pos

    agg = {
        "source_row": lambda x: ";".join(map(str, sorted(set(x))))[:1000],
        "source_side": lambda x: ";".join(sorted(set(map(str, x))))[:200],
        "source_type": lambda x: ";".join(sorted(set(map(str, x))))[:200],
        "complex_state": lambda x: ";".join(sorted(set(map(str, x))))[:200],
        "label": "first",
        "seg_len": "first",
    }
    pos = pos.groupby(["uniprot", "seg_start", "seg_end"], as_index=False).agg(agg)
    pos["seg_len"] = pos["seg_end"] - pos["seg_start"] + 1
    return pos

def build_negative_segments(
    all_idrs: pd.DataFrame,
    positive: pd.DataFrame,
    negative_scope: str = "all",
) -> pd.DataFrame:
    pass
    if positive.empty:
        pos_by_uid: Dict[str, List[Tuple[int, int]]] = {}
        pos_uids = set()
    else:
        pos_by_uid = {
            uid: list(zip(g["seg_start"].astype(int), g["seg_end"].astype(int)))
            for uid, g in positive.groupby("uniprot")
        }
        pos_uids = set(pos_by_uid)

    rows = []
    for _, r in all_idrs.iterrows():
        uid = r["uniprot"]
        if negative_scope == "same_proteins" and uid not in pos_uids:
            continue
        idr_range = (int(r["idr_start"]), int(r["idr_end"]))
        overlaps_pos = any(ranges_overlap(idr_range, pr) for pr in pos_by_uid.get(uid, []))
        if overlaps_pos:
            continue
        rows.append({
            "uniprot": uid,
            "seg_start": idr_range[0],
            "seg_end": idr_range[1],
            "seg_len": range_len(idr_range),
            "source_row": "from_complex_7002.xlsx",
            "source_side": "Disorder_Regions",
            "source_type": "IDR",
            "complex_state": "other_idr",
            "label": "negative",
        })
    return pd.DataFrame(rows)

def assign_parent_idr(segments: pd.DataFrame, all_idrs: pd.DataFrame) -> pd.DataFrame:
    pass
    idrs_by_uid = {uid: g for uid, g in all_idrs.groupby("uniprot")}
    out = []
    for _, seg in segments.iterrows():
        uid = seg["uniprot"]
        s, e = int(seg["seg_start"]), int(seg["seg_end"])
        best = None
        best_ov = 0
        for _, idr in idrs_by_uid.get(uid, pd.DataFrame()).iterrows():
            is_, ie = int(idr["idr_start"]), int(idr["idr_end"])
            ov = max(0, min(e, ie) - max(s, is_) + 1)
            if ov > best_ov:
                best_ov = ov
                best = idr
        row = seg.to_dict()
        if best is not None and best_ov > 0:
            row["idr_start"] = int(best["idr_start"])
            row["idr_end"] = int(best["idr_end"])
            row["idr_len"] = int(best["idr_len"])
            row["idr_overlap_len"] = best_ov
        else:

            row["idr_start"] = s
            row["idr_end"] = e
            row["idr_len"] = e - s + 1
            row["idr_overlap_len"] = 0
        out.append(row)
    return pd.DataFrame(out)

         

def find_pssm_file(pssm_dir: str | Path, uniprot: str) -> Optional[Path]:
    p = Path(pssm_dir)
    candidates = []
    for ext in [".npy", ".npz", ".pssm", ".txt", ".tsv", ".csv"]:
        candidates.append(p / f"{uniprot}{ext}")
        candidates.append(p / f"{uniprot.lower()}{ext}")
        candidates.append(p / f"{uniprot.upper()}{ext}")
    for c in candidates:
        if c.exists():
            return c

    hits = list(p.glob(f"{uniprot}*")) + list(p.glob(f"{uniprot.lower()}*")) + list(p.glob(f"{uniprot.upper()}*"))
    hits = [h for h in hits if h.suffix.lower() in {".npy", ".npz", ".pssm", ".txt", ".tsv", ".csv"}]
    return hits[0] if hits else None

def _parse_psiblast_ascii_full(path: Path) -> Optional[Tuple[np.ndarray, str, Dict[str, np.ndarray]]]:
    rows = []
    aas = []
    observed_pct_rows = []
    information = []
    real_match_weight = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 22:
                continue
            if not parts[0].isdigit():
                continue
            aa = parts[1].upper()

                                                                            

                                                    
            if len(aa) != 1 or not aa.isalpha():
                continue
            vals = []
            ok = True
            for x in parts[2:22]:
                try:
                    vals.append(float(x))
                except ValueError:
                    ok = False
                    break
            if ok and len(vals) == 20:
                rows.append([vals[index] for index in PSIBLAST_TO_INTERNAL])
                aas.append(aa)
                if len(parts) >= 42:
                    try:
                        observed_raw = [float(x) for x in parts[22:42]]
                        observed_pct_rows.append(
                            [observed_raw[index] for index in PSIBLAST_TO_INTERNAL]
                        )
                    except ValueError:
                        observed_pct_rows.append([np.nan] * 20)
                else:
                    observed_pct_rows.append([np.nan] * 20)
                try:
                    information.append(float(parts[42]) if len(parts) > 42 else np.nan)
                except ValueError:
                    information.append(np.nan)
                try:
                    real_match_weight.append(float(parts[43]) if len(parts) > 43 else np.nan)
                except ValueError:
                    real_match_weight.append(np.nan)
    if not rows:
        return None
    extras = {
        "observed_pct": np.asarray(observed_pct_rows, dtype=float),
        "information": np.asarray(information, dtype=float),
        "real_match_weight": np.asarray(real_match_weight, dtype=float),
    }
    return np.asarray(rows, dtype=float), "".join(aas), extras

def _parse_psiblast_ascii(path: Path) -> Optional[Tuple[np.ndarray, str]]:
    parsed = _parse_psiblast_ascii_full(path)
    if parsed is None:
        return None
    pssm, seq, _extras = parsed
    return pssm, seq

def _parse_table_pssm(path: Path) -> Optional[np.ndarray]:

    try:
        df = pd.read_csv(path, sep=None, engine="python", comment="#")
    except Exception:
        try:
            df = pd.read_csv(path, sep="\t", comment="#")
        except Exception:
            return None

    upper_cols = {str(c).strip().upper(): c for c in df.columns}
    if all(aa in upper_cols for aa in AA_ORDER):
        mat = df[[upper_cols[aa] for aa in AA_ORDER]].to_numpy(dtype=float)
        return mat

    numeric = df.apply(pd.to_numeric, errors="coerce")

    valid_cols = [c for c in numeric.columns if numeric[c].notna().sum() >= max(1, len(numeric) * 0.8)]
    if len(valid_cols) == 20:
        return numeric[valid_cols].to_numpy(dtype=float)
    if len(valid_cols) > 20:

        return numeric[valid_cols[-20:]].to_numpy(dtype=float)
    return None

def load_pssm(path: Path) -> Tuple[np.ndarray, Optional[str]]:
    pass
    suffix = path.suffix.lower()
    if suffix == ".npy":
        arr = np.load(path)
        return np.asarray(arr, dtype=float), None
    if suffix == ".npz":
        data = np.load(path)
        key = "pssm" if "pssm" in data.files else data.files[0]
        return np.asarray(data[key], dtype=float), None

    parsed = _parse_psiblast_ascii(path)
    if parsed is not None:
        return parsed

    mat = _parse_table_pssm(path)
    if mat is not None:
        return mat, None

    raise ValueError(f"无法解析 PSSM 文件: {path}")

def load_pssm_with_extras(path: Path) -> Tuple[np.ndarray, Optional[str], Dict[str, np.ndarray]]:
    pass
    suffix = path.suffix.lower()
    if suffix in {".pssm", ".txt"}:
        parsed = _parse_psiblast_ascii_full(path)
        if parsed is not None:
            return parsed
    pssm, seq = load_pssm(path)
    extras = {
        "observed_pct": np.full_like(pssm, np.nan, dtype=float),
        "information": np.full(pssm.shape[0], np.nan, dtype=float),
        "real_match_weight": np.full(pssm.shape[0], np.nan, dtype=float),
    }
    return pssm, seq, extras

       

def softmax_rows(z: np.ndarray, tau: float = 1.0) -> np.ndarray:
    x = z / tau
    x = x - np.nanmax(x, axis=1, keepdims=True)
    expx = np.exp(x)
    return expx / (np.nansum(expx, axis=1, keepdims=True) + EPS)

def percentile_ranks(x: np.ndarray) -> np.ndarray:
    pass
    n = len(x)
    if n <= 1:
        return np.ones(n)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(n, dtype=float)
    return ranks / (n - 1)

def mad(x: np.ndarray) -> float:
    med = np.nanmedian(x)
    return float(np.nanmedian(np.abs(x - med)))

def compute_residue_conservation(
    pssm: np.ndarray,
    seq: str,
    idr_ranges: Sequence[Tuple[int, int]],
    tau: float = 1.0,
) -> pd.DataFrame:
    pass

    pssm = np.asarray(pssm, dtype=float)
    if pssm.ndim != 2 or pssm.shape[1] != 20:
        raise ValueError(f"PSSM 必须是 L×20，当前 shape={pssm.shape}")
    L = min(len(seq), pssm.shape[0])
    if L == 0:
        raise ValueError("空序列或空 PSSM")
    seq = seq[:L]
    pssm = pssm[:L]

    row_mean = np.nanmean(pssm, axis=1, keepdims=True)
    row_std = np.nanstd(pssm, axis=1, keepdims=True)
    z = (pssm - row_mean) / (row_std + EPS)
    q = softmax_rows(z, tau=tau)

    aa_idx = np.array([AA_TO_IDX.get(a, -1) for a in seq], dtype=int)
    valid_aa = aa_idx >= 0

    q_cur = np.zeros(L, dtype=float)
    q_cur[valid_aa] = q[np.arange(L)[valid_aa], aa_idx[valid_aa]]

    current_pref = np.maximum(0.0, (q_cur - 1.0 / 20.0) / (1.0 - 1.0 / 20.0))

    q_other = q.copy()
    q_other[np.arange(L)[valid_aa], aa_idx[valid_aa]] = -np.inf
    q_2nd = np.nanmax(q_other, axis=1)
    q_2nd[~np.isfinite(q_2nd)] = 0.0
    top1_dominance = np.maximum(0.0, (q_cur - q_2nd) / (q_cur + EPS))

    entropy = -np.nansum(q * np.log(q + EPS), axis=1)
    distribution_specificity = 1.0 - entropy / math.log(20.0)
    distribution_specificity = np.clip(distribution_specificity, 0.0, 1.0)

    residue_cons = np.sqrt(
        np.maximum(
            current_pref * (0.7 * top1_dominance + 0.3 * distribution_specificity),
            0.0,
        )
    )

    residue_rel_idr = np.full(L, np.nan, dtype=float)
    residue_percentile_idr = np.full(L, np.nan, dtype=float)
    idr_id_for_pos = np.full(L, -1, dtype=int)

    for idr_idx, (s, e) in enumerate(idr_ranges, start=1):
        cr = clamp_range(s, e, L)
        if cr is None:
            continue
        s2, e2 = cr
        idx = np.arange(s2 - 1, e2)
        vals = distribution_specificity[idx]
        med = np.nanmedian(vals)
        scale = 1.4826 * mad(vals) + EPS
        residue_rel_idr[idx] = (vals - med) / scale
        residue_percentile_idr[idx] = percentile_ranks(vals)
        idr_id_for_pos[idx] = idr_idx

    out = pd.DataFrame({
        "position": np.arange(1, L + 1),
        "aa": list(seq),
        "current_pref": current_pref,
        "top1_dominance": top1_dominance,
        "distribution_specificity": distribution_specificity,
        "residue_cons": residue_cons,
        "residue_rel_IDR": residue_rel_idr,
        "residue_percentile_IDR": residue_percentile_idr,
        "parent_idr_index": idr_id_for_pos,
    })
    return out

def topk_mean(x: np.ndarray, k: int) -> float:
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.nan
    k = min(k, len(x))
    return float(np.mean(np.sort(x)[-k:]))

def summarize_segment(res_df: pd.DataFrame, start: int, end: int) -> Dict[str, float]:
    sub = res_df[(res_df["position"] >= start) & (res_df["position"] <= end)]
    if sub.empty:
        return {}
    x = sub["residue_cons"].to_numpy(float)
    rel = sub["residue_rel_IDR"].to_numpy(float)
    pct = sub["residue_percentile_IDR"].to_numpy(float)
    return {
        "n_residue_scored": int(len(sub)),
        "residue_cons_mean": float(np.nanmean(x)),
        "residue_cons_median": float(np.nanmedian(x)),
        "residue_cons_max": float(np.nanmax(x)),
        "residue_cons_top3_mean": topk_mean(x, 3),
        "residue_cons_top5_mean": topk_mean(x, 5),
        "residue_rel_IDR_mean": float(np.nanmean(rel)),
        "residue_rel_IDR_median": float(np.nanmedian(rel)),
        "residue_rel_IDR_max": float(np.nanmax(rel)),
        "residue_rel_IDR_top3_mean": topk_mean(rel, 3),
        "residue_rel_IDR_top5_mean": topk_mean(rel, 5),
        "residue_rel_IDR_peak_count_gt1": int(np.nansum(rel > 1.0)),
        "residue_rel_IDR_peak_count_gt2": int(np.nansum(rel > 2.0)),
        "residue_percentile_IDR_mean": float(np.nanmean(pct)),
        "residue_percentile_IDR_top3_mean": topk_mean(pct, 3),
        "current_pref_mean": float(np.nanmean(sub["current_pref"])),
        "top1_dominance_mean": float(np.nanmean(sub["top1_dominance"])),
        "distribution_specificity_mean": float(np.nanmean(sub["distribution_specificity"])),
    }

      

def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan

    gt = 0
    lt = 0
    for xi in x:
        gt += np.sum(xi > y)
        lt += np.sum(xi < y)
    return float((gt - lt) / (len(x) * len(y)))

def safe_mannwhitney(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    try:
        from scipy.stats import mannwhitneyu
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]
        if len(x) == 0 or len(y) == 0:
            return np.nan, np.nan
        res = mannwhitneyu(x, y, alternative="two-sided")
        return float(res.statistic), float(res.pvalue)
    except Exception:
        return np.nan, np.nan

def compare_groups(seg_df: pd.DataFrame, metrics: Sequence[str]) -> pd.DataFrame:
    rows = []
    pos = seg_df[seg_df["label"] == "positive"]
    neg = seg_df[seg_df["label"] == "negative"]
    for m in metrics:
        if m not in seg_df.columns:
            continue
        x = pos[m].to_numpy(float)
        y = neg[m].to_numpy(float)
        U, p = safe_mannwhitney(x, y)
        rows.append({
            "metric": m,
            "n_positive": int(np.sum(~np.isnan(x))),
            "n_negative": int(np.sum(~np.isnan(y))),
            "positive_mean": float(np.nanmean(x)) if np.sum(~np.isnan(x)) else np.nan,
            "negative_mean": float(np.nanmean(y)) if np.sum(~np.isnan(y)) else np.nan,
            "positive_median": float(np.nanmedian(x)) if np.sum(~np.isnan(x)) else np.nan,
            "negative_median": float(np.nanmedian(y)) if np.sum(~np.isnan(y)) else np.nan,
            "delta_mean_pos_minus_neg": (float(np.nanmean(x)) - float(np.nanmean(y))) if np.sum(~np.isnan(x)) and np.sum(~np.isnan(y)) else np.nan,
            "mannwhitney_U": U,
            "mannwhitney_p": p,
            "cliffs_delta_pos_vs_neg": cliffs_delta(x, y),
        })
    out = pd.DataFrame(rows)
    if "mannwhitney_p" in out:
        out = out.sort_values(["mannwhitney_p", "metric"], na_position="last")
    return out

def make_basic_plots(seg_df: pd.DataFrame, out_dir: Path, metrics: Sequence[str]) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        print("[WARN] matplotlib 不可用，跳过绘图。", file=sys.stderr)
        return
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for m in metrics:
        if m not in seg_df.columns:
            continue
        data = [
            seg_df.loc[seg_df["label"] == "positive", m].dropna().to_numpy(float),
            seg_df.loc[seg_df["label"] == "negative", m].dropna().to_numpy(float),
        ]
        if len(data[0]) == 0 or len(data[1]) == 0:
            continue
        plt.figure(figsize=(4, 4), dpi=160)
        plt.boxplot(data, labels=["positive", "other_IDR"], showfliers=False)
        plt.ylabel(m)
        plt.title(f"{m}: positive vs other IDR")
        plt.tight_layout()
        plt.savefig(plot_dir / f"boxplot_{m}.png")
        plt.close()

     

def main() -> None:
    parser = argparse.ArgumentParser(description="基于 PSSM 的 IDR 保守性分析")
    parser.add_argument("--idppi", required=True, help="IDPPI_2_remains_disordered.tsv")
    parser.add_argument("--regions", required=True, help="from_complex_7002.xlsx")
    parser.add_argument("--fasta", required=True, help="from_complex_7002.fa")
    parser.add_argument("--pssm_dir", required=True, help="PSSM 矩阵目录")
    parser.add_argument("--out_dir", required=True, help="输出目录")
    parser.add_argument("--tau", type=float, default=1.0, help="PSSM softmax 温度参数")
    parser.add_argument(
        "--negative_scope",
        choices=["all", "same_proteins"],
        default="all",
        help="negative IDR 取所有其他 IDR，还是只取 positive 蛋白中的其他 IDR",
    )
    parser.add_argument("--min_idr_len", type=int, default=1, help="negative IDR 最小长度过滤")
    parser.add_argument("--make_plots", action="store_true", help="输出基础箱线图")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/7] 读取 FASTA / IDR 区域 / IDPPI positive segments ...")
    seqs = parse_fasta(args.fasta)
    all_idrs = load_all_idrs(args.regions)
    all_idrs = all_idrs[all_idrs["idr_len"] >= args.min_idr_len].copy()
    positives = extract_positive_segments(args.idppi)
    negatives = build_negative_segments(all_idrs, positives, negative_scope=args.negative_scope)

    segments = pd.concat([positives, negatives], ignore_index=True)
    segments = assign_parent_idr(segments, all_idrs)
    segments["segment_id"] = [f"seg_{i:07d}" for i in range(len(segments))]

    print(f"  positive segments: {sum(segments['label']=='positive')}")
    print(f"  negative segments: {sum(segments['label']=='negative')}")

    needed_uids = sorted(set(segments["uniprot"]))
    print(f"[2/7] 需要处理 UniProt 数: {len(needed_uids)}")

    pssm_cache: Dict[str, pd.DataFrame] = {}
    missing = []
    errors = []

    idrs_by_uid = {
        uid: list(zip(g["idr_start"].astype(int), g["idr_end"].astype(int)))
        for uid, g in all_idrs.groupby("uniprot")
    }

    print("[3/7] 读取 PSSM 并计算 residue-level conservation ...")
    for idx, uid in enumerate(needed_uids, start=1):
        seq = seqs.get(uid)
        if seq is None:
            missing.append({"uniprot": uid, "reason": "missing_fasta"})
            continue
        pssm_file = find_pssm_file(args.pssm_dir, uid)
        if pssm_file is None:
            missing.append({"uniprot": uid, "reason": "missing_pssm"})
            continue
        try:
            pssm, pssm_seq = load_pssm(pssm_file)
            if pssm.shape[0] != len(seq):

                errors.append({
                    "uniprot": uid,
                    "reason": "length_mismatch_truncated",
                    "fasta_len": len(seq),
                    "pssm_len": pssm.shape[0],
                    "pssm_file": str(pssm_file),
                })
            res_df = compute_residue_conservation(pssm, seq, idrs_by_uid.get(uid, []), tau=args.tau)
            res_df.insert(0, "uniprot", uid)
            pssm_cache[uid] = res_df
        except Exception as e:
            errors.append({"uniprot": uid, "reason": f"pssm_error: {e}", "pssm_file": str(pssm_file)})
        if idx % 200 == 0:
            print(f"  processed {idx}/{len(needed_uids)} proteins")

    pd.DataFrame(missing).to_csv(out_dir / "missing_inputs.tsv", sep="\t", index=False)
    pd.DataFrame(errors).to_csv(out_dir / "pssm_warnings.tsv", sep="\t", index=False)

    print("[4/7] 汇总 segment-level features ...")
    seg_rows = []
    residue_rows = []
    for _, seg in segments.iterrows():
        uid = seg["uniprot"]
        res_df = pssm_cache.get(uid)
        base = seg.to_dict()
        if res_df is None:
            base["status"] = "no_pssm_or_fasta"
            seg_rows.append(base)
            continue
        L = len(res_df)
        cr = clamp_range(int(seg["seg_start"]), int(seg["seg_end"]), L)
        if cr is None:
            base["status"] = "invalid_range"
            seg_rows.append(base)
            continue
        s, e = cr
        summary = summarize_segment(res_df, s, e)
        base.update(summary)
        base["status"] = "ok"
        seg_rows.append(base)

        sub = res_df[(res_df["position"] >= s) & (res_df["position"] <= e)].copy()
        sub.insert(1, "segment_id", seg["segment_id"])
        sub.insert(2, "label", seg["label"])
        sub.insert(3, "seg_start", s)
        sub.insert(4, "seg_end", e)
        residue_rows.append(sub)

    seg_df = pd.DataFrame(seg_rows)
    residue_df = pd.concat(residue_rows, ignore_index=True) if residue_rows else pd.DataFrame()

    seg_df.to_csv(out_dir / "conservation_segment_level.tsv", sep="\t", index=False)
    if not residue_df.empty:
        residue_df.to_csv(out_dir / "conservation_residue_level.tsv", sep="\t", index=False)

    print("[5/7] 统计 positive vs other IDR 差异 ...")
    metrics = [
        "residue_cons_mean", "residue_cons_median", "residue_cons_max",
        "residue_cons_top3_mean", "residue_cons_top5_mean",
        "residue_rel_IDR_mean", "residue_rel_IDR_median", "residue_rel_IDR_max",
        "residue_rel_IDR_top3_mean", "residue_rel_IDR_top5_mean",
        "residue_rel_IDR_peak_count_gt1", "residue_rel_IDR_peak_count_gt2",
        "residue_percentile_IDR_mean", "residue_percentile_IDR_top3_mean",
        "current_pref_mean", "top1_dominance_mean", "distribution_specificity_mean",
    ]
    ok_seg = seg_df[seg_df.get("status", "") == "ok"].copy()
    stats = compare_groups(ok_seg, metrics)
    stats.to_csv(out_dir / "conservation_group_stats.tsv", sep="\t", index=False)

    if args.make_plots:
        print("[6/7] 输出基础箱线图 ...")
        make_basic_plots(ok_seg, out_dir, [
            "residue_cons_mean", "residue_cons_top5_mean",
            "residue_rel_IDR_mean", "residue_rel_IDR_top5_mean",
            "residue_rel_IDR_peak_count_gt2",
        ])
    else:
        print("[6/7] 跳过绘图。若需要图，添加 --make_plots")

    print("[7/7] 完成。输出文件：")
    for name in [
        "conservation_segment_level.tsv",
        "conservation_residue_level.tsv",
        "conservation_group_stats.tsv",
        "missing_inputs.tsv",
        "pssm_warnings.tsv",
    ]:
        f = out_dir / name
        if f.exists():
            print(f"  {f}")

if __name__ == "__main__":
    main()
