import gzip
import json
import os
import numpy as np
from typing import List, Dict, Optional, Tuple

class PAEService:
    """
    Service for analyzing AlphaFold 3 Predicted Aligned Error (PAE) matrices.
    
    Identifies 'Rigid Blocks' (domains/regions with low internal error) and 
    calculates inter-domain stability.
    """

    def __init__(self, static_base_path: str = "backend/static/data"):
        self.static_base_path = static_base_path

    def load_pae_matrix(self, relative_path: str) -> Optional[np.ndarray]:
        """
        Loads the PAE matrix from a JSON file.
        Handles both AFDB format (nested list) and standard AF3 format.
        """
        # Resolve path robustly
        from app.services.storage import STATIC_DATA_DIR
        
        # Standardize URL/Path logic
        if relative_path.startswith("/static/data/"):
            clean_path = relative_path.replace("/static/data/", "", 1).lstrip('/')
            full_path = STATIC_DATA_DIR / clean_path
        elif relative_path.startswith("static/data/"):
            clean_path = relative_path.replace("static/data/", "", 1).lstrip('/')
            full_path = STATIC_DATA_DIR / clean_path
        elif relative_path.startswith("/static/"):
            clean_path = relative_path.replace("/static/", "", 1).lstrip('/')
            full_path = STATIC_DATA_DIR.parent / clean_path
        else:
            # Fallback to direct joins
            clean_path = relative_path.lstrip('/')
            full_path = os.path.join(self.static_base_path, clean_path)
            
        full_path = str(full_path)

        # Transparent .gz fallback: if the plain .json is gone but a sibling
        # <file>.gz exists, read that instead. Enables the storage-cleanup
        # gzip pass (Phase: pre-server-transfer) to drop the originals.
        opener = open
        if not os.path.exists(full_path):
            if os.path.exists(full_path + ".gz"):
                full_path = full_path + ".gz"
                opener = lambda p, _mode: gzip.open(p, "rt")
            else:
                return None
        elif full_path.endswith(".gz"):
            opener = lambda p, _mode: gzip.open(p, "rt")
        elif os.path.exists(full_path + ".gz") and os.path.getsize(full_path) < 4096:
            # 9,685 ccRCC dirs carry a ~250-byte `_confidences.json` stub (it holds
            # the summary keys, no PAE) next to the real `.json.gz`; the stub won
            # and the isoform got no PAE mean/blocks. Prefer the gz in that case.
            full_path = full_path + ".gz"
            opener = lambda p, _mode: gzip.open(p, "rt")

        try:
            with opener(full_path, 'r') as f:
                data = json.load(f)
            
            pae_data = None
            if isinstance(data, list) and len(data) > 0 and 'predicted_aligned_error' in data[0]:
                pae_data = data[0]['predicted_aligned_error']
            elif 'predicted_aligned_error' in data:
                pae_data = data['predicted_aligned_error']
            elif 'pae' in data:
                pae_data = data['pae']
            elif isinstance(data, list) and isinstance(data[0], list):
                pae_data = data
                
            if pae_data:
                return np.array(pae_data)
        except Exception as e:
            print(f"Error loading PAE matrix: {e}")
            return None
        return None

    # Adjacent-merge thresholds tightened in audit-6.
    # The previous values (gap ≤ 5 residues, combined_mean < threshold + 2)
    # over-fragmented contiguous rigid regions: live data showed proteins
    # decomposed into 17 adjacent blocks where the right answer was "one
    # rigid region 302-1170".  The tighter values below collapse the common
    # case of 1-3 residue inter-block gaps with near-identical PAE.
    ADJACENT_MERGE_MAX_GAP = 3       # was 5
    ADJACENT_MERGE_TOLERANCE = 1.0   # was 2.0 (additive over threshold)

    def find_rigid_blocks(self, matrix: np.ndarray, threshold: float = 6.0) -> List[Dict]:
        """
        Identifies 'Rigid Blocks' using a Seed-and-Scale approach.

          1. Find every residue where the local 10×10 window has mean PAE
             below `threshold`.  Mark it as a seed.
          2. From each seed, extend the right edge to the largest j such
             that mean(matrix[i:j, i:j]) < threshold.
          3. Merge adjacent blocks whose inter-block gap is ≤ 3 residues AND
             whose combined-region mean stays under threshold + 1.0
             (audit-6 tightening — previously gap ≤ 5 and tolerance 2.0).

        Edge bleed (the right-edge extension keeps the block as long as the
        BLOCK mean is under threshold, not each residue individually) is a
        known limitation; documented in the structural-scaffolding doc page.
        """
        n = matrix.shape[0]
        blocks = []

        # Summed-area table: the mean of any diagonal block matrix[a:b, a:b] is
        # then O(1) instead of an O((b-a)^2) np.mean. The original per-j loop
        # called np.mean ~16k times per isoform (57 of 68 s in the profile of
        # the structural-metrics backfill); the semantics below are identical.
        S = np.zeros((n + 1, n + 1), dtype=np.float64)
        S[1:, 1:] = np.cumsum(np.cumsum(np.asarray(matrix, dtype=np.float64), axis=0), axis=1)

        def block_mean(a: int, b: int) -> float:
            return float((S[b, b] - S[a, b] - S[b, a] + S[a, a]) / ((b - a) ** 2))

        # 1. Start scanning residues
        i = 0
        while i < n - 20:
            # Look for a start where the immediate window is stable
            if block_mean(i, i + 10) < threshold:
                # We found a potential start. Now find the largest j (n ≥ j > i+20)
                # whose block mean is under the threshold — vectorised over all j.
                js = np.arange(i + 21, n + 1)
                means = (S[js, js] - S[i, js] - S[js, i] + S[i, i]) / ((js - i) ** 2)
                hits = js[means < threshold]
                best_j = int(hits.max()) if hits.size else i + 20

                if best_j > i + 20:
                    blocks.append({
                        "start": i + 1,
                        "end": best_j,
                        "mean_pae": block_mean(i, best_j)
                    })
                    i = best_j
                    continue
            i += 1

        # 2. Merge adjacent blocks with the audit-6 tightened thresholds.
        if len(blocks) > 1:
            merged = []
            current = blocks[0]
            for next_b in blocks[1:]:
                if next_b['start'] <= current['end'] + self.ADJACENT_MERGE_MAX_GAP:
                    s, e = current['start'] - 1, next_b['end']
                    combined_mean = block_mean(s, e)
                    if combined_mean < threshold + self.ADJACENT_MERGE_TOLERANCE:
                        current = {
                            "start": current['start'],
                            "end": next_b['end'],
                            "mean_pae": combined_mean
                        }
                    else:
                        merged.append(current)
                        current = next_b
                else:
                    merged.append(current)
                    current = next_b
            merged.append(current)
            return merged

        return blocks

    def get_peptide_stability(self, matrix: np.ndarray, start: int, end: int) -> float:
        """
        Calculates the mean PAE within a specific residue window.
        Lower is better (more stable).
        """
        # Adjust 1-based to 0-based
        s0 = max(0, start - 1)
        e0 = min(matrix.shape[0], end)
        
        if s0 >= e0: return 99.0
        
        sub = matrix[s0:e0, s0:e0]
        return float(np.mean(sub))

    def get_inter_domain_pae(self, matrix: np.ndarray, domain1: Tuple[int, int], domain2: Tuple[int, int]) -> float:
        """
        Calculates the mean PAE between two domains.
        Used to identify if two domains are floating relative to each other (high PAE)
        or rigidly docked (low PAE).
        """
        s1, e1 = domain1[0]-1, domain1[1]
        s2, e2 = domain2[0]-1, domain2[1]
        
        if s1 >= e1 or s2 >= e2: return 99.0
        
        sub = matrix[s1:e1, s2:e2]
        return float(np.mean(sub))
