import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def create_splits(
    df: pd.DataFrame,
    split_mode: str = "eval_only",
    test_fraction: float = 0.5,
    random_seed: int = 42,
) -> pd.DataFrame:
    """
    Split the dataset based on the split_mode.
    
    Args:
        split_mode: One of ["eval_only", "train_only", "train_val_test"]
    """
    np.random.seed(random_seed)
    df = df.copy()
    
    if split_mode == "train_only":
        df["split"] = "train"
    elif split_mode == "eval_only":
        mask = np.random.rand(len(df)) < test_fraction
        df["split"] = np.where(mask, "test", "val")
    else:
        # 80/10/10 split
        rand = np.random.rand(len(df))
        df["split"] = np.where(
            rand < 0.8, "train",
            np.where(rand < 0.9, "val", "test")
        )
        
    return df
