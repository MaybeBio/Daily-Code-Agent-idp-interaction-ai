"""
Region-aware feature evaluator — GPU-accelerated.

Core algorithm: sort-based sweep (O(N log N) per feature-chunk) fully
executed on GPU. All N residues and F_c features in a chunk are sorted,
cumsum'd, and searchsorted simultaneously using PyTorch primitives.

For N=1.8M, F_c=256:
  CPU argsort per feature: ~300 ms × 10240 features = 51 min per triplet
  GPU sort (F_c=256 chunk): ~50 ms per chunk × 40 chunks = 2 min per triplet
"""
import logging
from typing import Optional

import numpy as np
import torch

logger = logging.getLogger(__name__)


class RegionAwareEvaluator:
    """
    Sort-based, GPU-accelerated evaluator for SAE feature quality.

    Processes features in chunks on the GPU. For each chunk:
      1. Load (N, F_c) activations → GPU float32.
      2. torch.sort (radix sort) along N dim → sorted acts + sort order.
      3. Gather labels by sort order → sorted labels.
      4. cumsum → cumulative TP/FP → precision at every rank.
      5. searchsorted vs sorted region_maxes → region recall at every rank.
      6. Compute F1, take argmax per feature.

    Memory per chunk (F_c=256, N=1.8M):
      Acts    : 256 × 1.8M × 4 = 1.8 GB
      Indices : 256 × 1.8M × 8 = 3.7 GB  (int64 from torch.sort)
      Labels  : 256 × 1.8M × 4 = 1.8 GB
      Cumsums : 2 × above    = 3.6 GB
      Totals  : ~11 GB
    """

    def __init__(self, chunk_size: int = 256):
        self.chunk_size = chunk_size

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def evaluate_all_features(
        self,
        all_activations: np.ndarray,           # (N, F) mmap or ndarray, any dtype
        labels: np.ndarray,                     # (N,) float32
        regions: list[tuple[int, int]],
        triplet: str,
        region_maxes_matrix: Optional[np.ndarray] = None,  # (R, F) float32
        device: str = "cpu",
        val_residue_mask: Optional[np.ndarray] = None,     # (N,) bool, True = validation
        val_region_mask: Optional[np.ndarray] = None,      # (R,) bool, True = validation
    ) -> list[dict]:
        """
        Evaluate all F features against one triplet label set.

        Args:
            all_activations: (N, F) activation matrix. Mmap OK.
            labels: (N,) binary labels for this triplet.
            regions: (start, end) region tuples in mmap coordinates.
            triplet: Name for logging.
            region_maxes_matrix: (R, F) precomputed per-region max activations.
                If None, computed here from all_activations (row scan).
            device: "cuda" or "cpu". GPU strongly recommended.
            val_residue_mask: (N,) bool. If provided (together with
                ``val_region_mask``), the per-feature activation threshold is
                selected on validation residues/regions only, and the reported
                precision / region_recall / F1 / average_precision are computed
                on the held-out test residues/regions. This removes the
                threshold-selection data leakage of the legacy (no-split) path.
                When None, the legacy full-corpus behaviour is used.
            val_region_mask: (R,) bool aligned with ``regions`` /
                ``region_maxes_matrix`` rows. Required when ``val_residue_mask``
                is given.
        """
        from tqdm import tqdm

        n_residues, n_features = all_activations.shape
        n_regions = len(regions)
        C = self.chunk_size
        labels_f32 = labels.astype(np.float32)

        # ── Precompute region_maxes if not provided ────────────────────────
        if region_maxes_matrix is None and n_regions > 0:
            logger.info(
                f"  Computing region maxes for {n_regions} regions × {n_features} features..."
            )
            region_maxes_matrix = np.zeros((n_regions, n_features), dtype=np.float32)
            for r, (start, end) in enumerate(
                tqdm(regions, desc="Region maxes", leave=False, mininterval=30)
            ):
                if end > start:
                    region_maxes_matrix[r] = (
                        all_activations[start:end].astype(np.float32).max(axis=0)
                    )

        # ── Split-aware path: select threshold on val, report on test ──────
        if val_residue_mask is not None:
            if val_region_mask is None:
                raise ValueError(
                    "val_region_mask must be provided when val_residue_mask is set."
                )
            return self._evaluate_split(
                all_activations=all_activations,
                labels_f32=labels_f32,
                region_maxes_matrix=region_maxes_matrix,
                val_residue_mask=val_residue_mask.astype(bool),
                val_region_mask=val_region_mask.astype(bool),
                triplet=triplet,
                device=device,
            )

        # ── Move labels and sorted region maxes to GPU ─────────────────────
        labels_gpu = torch.from_numpy(labels_f32).to(device)              # (N,)
        n_pos = float(labels_gpu.sum().item())                            # scalar

        # sorted_rm_all: (F, R) ascending — one searchsorted per chunk column
        if region_maxes_matrix is not None and n_regions > 0:
            # Transpose to (F, R) so we can slice F_c rows per chunk
            sorted_rm_all = np.sort(region_maxes_matrix, axis=0).T       # (F, R)
        else:
            sorted_rm_all = None

        results = []
        n_chunks = (n_features + C - 1) // C

        for chunk_idx in tqdm(range(n_chunks), desc=f"Eval {triplet}",
                              leave=False, mininterval=30):
            cs = chunk_idx * C
            ce = min(cs + C, n_features)
            F_c = ce - cs

            # ── Load chunk → GPU ──────────────────────────────────────────
            # C-order mmap: reading adjacent columns is efficient per row
            acts_np = all_activations[:, cs:ce].astype(np.float32)       # (N, F_c) CPU
            acts_gpu = torch.from_numpy(acts_np).to(device)               # (N, F_c) GPU

            # Skip fully dead chunks
            feat_max = acts_gpu.max(dim=0).values                         # (F_c,)
            active = feat_max > 1e-6                                       # (F_c,)
            if not active.any():
                continue

            # ── Sort descending along N (GPU radix sort) ─────────────────
            sorted_acts, sort_order = torch.sort(acts_gpu, dim=0, descending=True)
            # sorted_acts : (N, F_c) descending
            # sort_order  : (N, F_c) int64 — which original row is at each rank

            # ── Gather labels by sort order ───────────────────────────────
            sorted_labels = labels_gpu[sort_order]                         # (N, F_c)

            # ── Cumulative TP/FP → precision ─────────────────────────────
            cum_tp = torch.cumsum(sorted_labels, dim=0)                    # (N, F_c)
            cum_fp = torch.cumsum(1.0 - sorted_labels, dim=0)             # (N, F_c)
            precision = cum_tp / (cum_tp + cum_fp + 1e-10)                 # (N, F_c)

            # ── Region recall via searchsorted ────────────────────────────
            if sorted_rm_all is not None:
                # sorted_rm_chunk: (F_c, R) ascending on GPU
                sorted_rm_chunk = torch.from_numpy(
                    sorted_rm_all[cs:ce].copy()
                ).to(device)                                               # (F_c, R)

                # sorted_acts.T: (F_c, N) — threshold values per feature
                sorted_acts_T = sorted_acts.T.contiguous()                # (F_c, N)

                # For each (feature, rank), n_regions_below_threshold:
                n_below = torch.searchsorted(
                    sorted_rm_chunk, sorted_acts_T                        # sorted seq, queries
                )                                                          # (F_c, N)
                n_hit = n_regions - n_below                               # (F_c, N)
                recall_T = n_hit.float() / max(1, n_regions)              # (F_c, N)

                precision_T = precision.T.contiguous()                    # (F_c, N)
            else:
                precision_T = precision.T.contiguous()                    # (F_c, N)
                recall_T = torch.zeros_like(precision_T)

            # ── F1 and best threshold per feature ────────────────────────
            f1_T = (
                2.0 * precision_T * recall_T
                / (precision_T + recall_T + 1e-10)
            )                                                              # (F_c, N)

            # Mask positions where feature doesn't actually fire
            sorted_acts_T_mask = sorted_acts.T.contiguous() > 1e-6       # (F_c, N)
            f1_T = f1_T * sorted_acts_T_mask

            best_t_idx = f1_T.argmax(dim=1)                               # (F_c,)

            # ── Average Precision from GPU tensors (no CPU sort) ─────────
            # AP = (1/n_pos) × Σ_k [sorted_labels[k,f] × precision[k,f]]
            # This is mathematically equivalent to sklearn's AP but runs on GPU.
            if n_pos > 0:
                ap_chunk = (sorted_labels * precision).sum(dim=0) / n_pos  # (F_c,)
            else:
                ap_chunk = torch.zeros(F_c, device=device)

            # ── Gather results back to CPU ────────────────────────────────
            best_t_np   = best_t_idx.cpu().numpy()
            f1_np       = f1_T.cpu().numpy()
            prec_np     = precision_T.cpu().numpy()
            rec_np      = recall_T.cpu().numpy()
            acts_t_np   = sorted_acts.T.contiguous().cpu().numpy()
            active_np   = active.cpu().numpy()
            ap_np       = ap_chunk.cpu().numpy()

            for local_f in range(F_c):
                if not active_np[local_f]:
                    continue

                bt       = best_t_np[local_f]
                best_f1  = float(f1_np[local_f, bt])
                if best_f1 < 1e-6:
                    continue

                global_f = cs + local_f

                results.append({
                    "feature_idx":       global_f,
                    "triplet":           triplet,
                    "f1":                best_f1,
                    "precision":         float(prec_np[local_f, bt]),
                    "region_recall":     float(rec_np[local_f, bt]),
                    "best_threshold":    float(acts_t_np[local_f, bt]),
                    "average_precision": float(ap_np[local_f]),
                })

        logger.info(f"  → {len(results)} features with F1 > 0 for '{triplet}'")
        return results

    # ------------------------------------------------------------------
    # Split-aware evaluation (threshold on val, metrics on test)
    # ------------------------------------------------------------------

    def _evaluate_split(
        self,
        all_activations: np.ndarray,            # (N, F)
        labels_f32: np.ndarray,                 # (N,) float32
        region_maxes_matrix: Optional[np.ndarray],  # (R, F) float32, aligned to regions
        val_residue_mask: np.ndarray,           # (N,) bool, True = val
        val_region_mask: np.ndarray,            # (R,) bool, True = val
        triplet: str,
        device: str = "cpu",
    ) -> list[dict]:
        """
        For each feature: pick the activation threshold that maximises F1 on the
        *validation* residues/regions, then report precision / region_recall /
        F1 / average_precision computed on the *test* residues/regions at that
        fixed threshold. No information from the test split influences the
        threshold choice, so the reported numbers are leakage-free.
        """
        from tqdm import tqdm

        n_residues, n_features = all_activations.shape
        C = self.chunk_size
        eps = 1e-10

        # Residue index sets (global, feature-independent)
        val_idx_np = np.nonzero(val_residue_mask)[0]
        test_idx_np = np.nonzero(~val_residue_mask)[0]
        val_idx = torch.from_numpy(val_idx_np).long().to(device)
        test_idx = torch.from_numpy(test_idx_np).long().to(device)

        labels_gpu = torch.from_numpy(labels_f32).to(device)            # (N,)
        val_labels = labels_gpu[val_idx]                                # (N_val,)
        test_labels = labels_gpu[test_idx]                             # (N_test,)
        test_labels_bool = test_labels > 0.5
        n_test_pos = float(test_labels.sum().item())

        # Region maxes split into val / test, sorted ascending per feature → (F, R_*)
        has_regions = (
            region_maxes_matrix is not None and region_maxes_matrix.shape[0] > 0
        )
        if has_regions:
            val_rm = region_maxes_matrix[val_region_mask]               # (R_val, F)
            test_rm = region_maxes_matrix[~val_region_mask]            # (R_test, F)
            n_val_regions = val_rm.shape[0]
            n_test_regions = test_rm.shape[0]
            val_rm_sorted_all = np.sort(val_rm, axis=0).T if n_val_regions else None    # (F, R_val)
            # test region maxes kept unsorted: we count >= threshold directly
            test_rm_all = test_rm.T if n_test_regions else None         # (F, R_test)
        else:
            n_val_regions = n_test_regions = 0
            val_rm_sorted_all = test_rm_all = None

        if n_val_regions == 0 or n_test_regions == 0:
            logger.warning(
                f"  '{triplet}': insufficient regions after split "
                f"(val={n_val_regions}, test={n_test_regions}); skipping triplet."
            )
            return []

        results = []
        n_chunks = (n_features + C - 1) // C

        for chunk_idx in tqdm(range(n_chunks), desc=f"Eval(split) {triplet}",
                              leave=False, mininterval=30):
            cs = chunk_idx * C
            ce = min(cs + C, n_features)
            F_c = ce - cs

            acts_np = all_activations[:, cs:ce].astype(np.float32)      # (N, F_c)
            acts_gpu = torch.from_numpy(acts_np).to(device)            # (N, F_c)

            active = acts_gpu.max(dim=0).values > 1e-6                  # (F_c,)
            if not active.any():
                continue

            val_acts = acts_gpu[val_idx]                               # (N_val, F_c)
            test_acts = acts_gpu[test_idx]                            # (N_test, F_c)

            # ── Step A: best threshold on validation ──────────────────────
            sorted_val_acts, val_order = torch.sort(val_acts, dim=0, descending=True)
            sorted_val_labels = val_labels[val_order]                  # (N_val, F_c)

            cum_tp = torch.cumsum(sorted_val_labels, dim=0)
            cum_fp = torch.cumsum(1.0 - sorted_val_labels, dim=0)
            val_prec = cum_tp / (cum_tp + cum_fp + eps)                # (N_val, F_c)

            val_rm_chunk = torch.from_numpy(
                val_rm_sorted_all[cs:ce].copy()
            ).to(device)                                               # (F_c, R_val)
            sva_T = sorted_val_acts.T.contiguous()                     # (F_c, N_val)
            n_below = torch.searchsorted(val_rm_chunk, sva_T)         # (F_c, N_val)
            val_recall_T = (n_val_regions - n_below).float() / max(1, n_val_regions)

            val_f1_T = 2.0 * val_prec.T * val_recall_T / (val_prec.T + val_recall_T + eps)
            val_f1_T = val_f1_T * (sva_T > 1e-6)                       # mask non-firing
            best_idx = val_f1_T.argmax(dim=1)                          # (F_c,)
            thresholds = sva_T.gather(1, best_idx[:, None]).squeeze(1)  # (F_c,) activation value

            # ── Step B: evaluate on test at fixed threshold ───────────────
            thr = thresholds[None, :]                                  # (1, F_c)
            test_pred = test_acts >= thr                              # (N_test, F_c)
            tp = (test_pred & test_labels_bool[:, None]).sum(dim=0).float()
            fp = (test_pred & ~test_labels_bool[:, None]).sum(dim=0).float()
            test_prec = tp / (tp + fp + eps)                          # (F_c,)

            test_rm_chunk = torch.from_numpy(
                test_rm_all[cs:ce].copy()
            ).to(device)                                              # (F_c, R_test)
            test_hits = (test_rm_chunk >= thresholds[:, None]).sum(dim=1).float()
            test_recall = test_hits / max(1, n_test_regions)          # (F_c,)

            test_f1 = 2.0 * test_prec * test_recall / (test_prec + test_recall + eps)

            # Average precision on the test split (GPU, sklearn-equivalent)
            sorted_test_acts, test_order = torch.sort(test_acts, dim=0, descending=True)
            sorted_test_labels = test_labels[test_order]
            t_cum_tp = torch.cumsum(sorted_test_labels, dim=0)
            t_cum_fp = torch.cumsum(1.0 - sorted_test_labels, dim=0)
            t_prec = t_cum_tp / (t_cum_tp + t_cum_fp + eps)
            if n_test_pos > 0:
                ap = (sorted_test_labels * t_prec).sum(dim=0) / n_test_pos
            else:
                ap = torch.zeros(F_c, device=device)

            # ── Gather → CPU ──────────────────────────────────────────────
            active_np = active.cpu().numpy()
            thr_np = thresholds.cpu().numpy()
            f1_np = test_f1.cpu().numpy()
            prec_np = test_prec.cpu().numpy()
            rec_np = test_recall.cpu().numpy()
            ap_np = ap.cpu().numpy()

            for local_f in range(F_c):
                if not active_np[local_f]:
                    continue
                # Reject features that don't fire in the VALIDATION split: their
                # threshold collapses to ~0, which degenerates to "predict all
                # positive" (precision=base-rate, recall=1) — a phantom F1, not a
                # real detector. A meaningful threshold must be clearly positive.
                if thr_np[local_f] <= 1e-6:
                    continue
                best_f1 = float(f1_np[local_f])
                if best_f1 < 1e-6:
                    continue
                results.append({
                    "feature_idx":       cs + local_f,
                    "triplet":           triplet,
                    "f1":                best_f1,
                    "precision":         float(prec_np[local_f]),
                    "region_recall":     float(rec_np[local_f]),
                    "best_threshold":    float(thr_np[local_f]),
                    "average_precision": float(ap_np[local_f]),
                })

        logger.info(
            f"  → {len(results)} features with test-F1 > 0 for '{triplet}' "
            f"(val/test split enforced)"
        )
        return results

    # Compatibility wrapper for single-feature callers

    def evaluate_feature(
        self,
        acts_col: np.ndarray,
        labels: np.ndarray,
        regions: list[tuple[int, int]],
        region_maxes: Optional[np.ndarray] = None,
        n_regions: int = 0,
    ) -> Optional[dict]:
        """Single-feature wrapper — CPU only, for compatibility."""
        max_act = float(acts_col.max())
        if max_act < 1e-6:
            return None

        order        = np.argsort(-acts_col, kind="stable")
        sorted_acts  = acts_col[order]
        sorted_labels = labels.astype(np.float32)[order]

        cum_tp    = np.cumsum(sorted_labels)
        cum_fp    = np.cumsum(1.0 - sorted_labels)
        precision = cum_tp / (cum_tp + cum_fp + 1e-10)

        if n_regions > 0 and region_maxes is not None:
            sorted_rm    = np.sort(region_maxes)
            n_hit        = n_regions - np.searchsorted(sorted_rm, sorted_acts, side="left")
            region_recall = n_hit / n_regions
        else:
            region_recall = np.zeros(len(acts_col))

        f1       = 2.0 * precision * region_recall / (precision + region_recall + 1e-10)
        valid    = sorted_acts > 1e-6
        if not valid.any():
            return None
        best_idx = int(np.argmax(f1 * valid))
        best_f1  = float(f1[best_idx])
        if best_f1 < 1e-6:
            return None

        try:
            ap = float(average_precision_score(labels, acts_col))
        except ValueError:
            ap = 0.0

        return {
            "f1":               best_f1,
            "precision":        float(precision[best_idx]),
            "region_recall":    float(region_recall[best_idx]),
            "best_threshold":   float(sorted_acts[best_idx]),
            "average_precision": ap,
        }
