"""
Atlas generation: builds the final output tables and flags missing annotations.

The Feature Atlas is the primary deliverable of IDPInterp — a table mapping
each SAE latent feature to its best-matching MobiDB biological concept,
with F1 scores, top-activating proteins, and discovery candidates.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class AtlasBuilder:
    """
    Builds the Disorder Mechanistic Atlas from per-triplet evaluation results.
    """

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def build_atlas(
        self,
        triplet_results: dict[str, pd.DataFrame],
    ) -> pd.DataFrame:
        """
        Combine per-triplet evaluation results into a unified atlas.
        
        For each feature, find the triplet with the highest F1 score.
        
        Args:
            triplet_results: Dict mapping triplet name to evaluation DataFrame
                            (output of RegionAwareEvaluator.evaluate_all_features).
        
        Returns:
            Atlas DataFrame with one row per feature, showing the best-matching
            triplet and its scores.
        """
        all_results = []
        for triplet, df in triplet_results.items():
            if df is not None and len(df) > 0:
                all_results.append(df)

        if not all_results:
            logger.warning("No evaluation results to build atlas from.")
            return pd.DataFrame()

        combined = pd.concat(all_results, ignore_index=True)

        # For each feature, select the triplet with the highest F1
        idx = combined.groupby("feature_idx")["f1"].idxmax()
        atlas = combined.loc[idx].reset_index(drop=True)
        atlas = atlas.sort_values("f1", ascending=False).reset_index(drop=True)

        logger.info(
            f"Atlas built: {len(atlas)} features, "
            f"{(atlas['f1'] > 0.5).sum()} with F1 > 0.5, "
            f"{(atlas['f1'] > 0.3).sum()} with F1 > 0.3"
        )

        return atlas

    def flag_missing_annotations(
        self,
        atlas: pd.DataFrame,
        all_activations: np.ndarray,
        all_labels: np.ndarray,
        metadata: pd.DataFrame,
        precision_threshold: float = 0.9,
    ) -> pd.DataFrame:
        """
        Flag candidate missing annotations.
        
        A feature is flagged if:
        1. It has high precision (>90%) on the validation set.
        2. It fires strongly on residues that are NOT labeled.
        
        Args:
            atlas: The feature atlas DataFrame.
            all_activations: Shape (n_residues, d_latent).
            all_labels: Shape (n_residues,), binary.
            metadata: Metadata index with protein boundaries.
            precision_threshold: Minimum precision to consider.
        
        Returns:
            DataFrame of candidate missing annotations.
        """
        high_precision = atlas[atlas["precision"] >= precision_threshold]
        candidates = []

        for _, row in high_precision.iterrows():
            feat_idx = row["feature_idx"]
            threshold = row["best_threshold"]
            feat_acts = all_activations[:, feat_idx]
            if hasattr(feat_acts, "toarray"):
                feat_acts = feat_acts.toarray().flatten()
            feat_acts = feat_acts.astype(np.float32)

            # Find residues where feature fires but no label exists
            fires = feat_acts >= threshold
            unlabeled = ~all_labels.astype(bool)
            novel = fires & unlabeled

            if not np.any(novel):
                continue

            # Map back to proteins (vectorized via searchsorted on sorted starts)
            novel_indices = np.where(novel)[0]
            if len(novel_indices) == 0:
                continue

            # metadata is sorted by mmap_start_idx (it was built in order)
            starts = metadata["mmap_start_idx"].values
            ends   = metadata["mmap_end_idx"].values

            # For each novel residue, find which protein it belongs to
            prot_idx = np.searchsorted(starts, novel_indices, side="right") - 1
            prot_idx = np.clip(prot_idx, 0, len(metadata) - 1)

            # Validate: residue must actually be within the found protein's range
            valid = novel_indices < ends[prot_idx]

            for ni, pi, ok in zip(novel_indices, prot_idx, valid):
                if not ok:
                    continue
                prot = metadata.iloc[pi]
                candidates.append({
                    "feature_idx": feat_idx,
                    "uniprot_id":  prot["uniprot_id"],
                    "residue_idx": ni - prot["mmap_start_idx"],
                    "activation":  float(feat_acts[ni]),
                    "precision":   row["precision"],
                    "triplet":     row["triplet"],
                })

        candidates_df = pd.DataFrame(candidates)

        if len(candidates_df) > 0:
            # Aggregate to protein-level spans
            candidates_df = (
                candidates_df.groupby(["feature_idx", "uniprot_id", "triplet"])
                .agg(
                    residue_start=("residue_idx", "min"),
                    residue_end=("residue_idx", "max"),
                    max_activation=("activation", "max"),
                    precision=("precision", "first"),
                    n_residues=("residue_idx", "count"),
                )
                .reset_index()
                .sort_values("max_activation", ascending=False)
            )

        logger.info(f"Flagged {len(candidates_df)} candidate missing annotations.")
        return candidates_df

    def save(self, atlas: pd.DataFrame, missing: pd.DataFrame | None = None) -> None:
        """Save atlas and missing annotations as Parquet files."""
        atlas_path = self.output_dir / "feature_atlas.parquet"
        atlas.to_parquet(atlas_path, engine="pyarrow", index=False)
        logger.info(f"Saved feature atlas: {atlas_path}")

        # Also save as CSV for easy human reading
        csv_path = self.output_dir / "feature_atlas.csv"
        atlas.to_csv(csv_path, index=False)
        logger.info(f"Saved feature atlas CSV: {csv_path}")

        if missing is not None and len(missing) > 0:
            missing_path = self.output_dir / "missing_annotations.parquet"
            missing.to_parquet(missing_path, engine="pyarrow", index=False)
            logger.info(f"Saved missing annotations: {missing_path}")
