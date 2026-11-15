import numpy as np


class RandomBaseline:
    """A baseline that makes random predictions."""
    
    @staticmethod
    def evaluate(labels: np.ndarray, regions: list[tuple[int, int]]) -> dict:
        """Evaluate a random predictor."""
        preds = np.random.rand(len(labels))
        
        # We can just threshold at 0.5 since it's random
        binary = (preds > 0.5).astype(np.float32)
        
        tp = np.sum(binary * labels)
        fp = np.sum(binary * (1 - labels))
        
        precision = tp / (tp + fp + 1e-10)
        
        # Region recall
        region_hits = 0
        for start, end in regions:
            if np.any(binary[start:end]):
                region_hits += 1
                
        region_recall = region_hits / max(1, len(regions))
        
        f1 = 2 * (precision * region_recall) / (precision + region_recall + 1e-10)
        
        return {
            "f1": float(f1),
            "precision": float(precision),
            "region_recall": float(region_recall),
        }
