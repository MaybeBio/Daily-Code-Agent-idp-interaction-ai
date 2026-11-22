"""
Inference for EmbeddingMultiTaskRegressor: pass a PLM embedding (or a CSV of
sequence+embedding rows), get back the 5 predicted properties.

Usage
-----
From a CSV where col 0 = sequence, cols 1..end = embedding dims (no property
columns needed for inference):
    python predict_from_embeddings.py --csv new_embeddings.csv --out predictions.csv

As a library:
    from predict_from_embeddings import load_model, predict_embeddings
    model, embedding_dim, device = load_model("best_model.pt")
    results = predict_embeddings(model, device, sequences, embeddings)  # embeddings: [N, D] array
"""

import argparse
import json

import pandas as pd
import numpy as np
import torch

from model_2 import EmbeddingMultiTaskRegressor

CHECKPOINT_PATH_DEFAULT = "outputs/best_model.pt"
PROPERTY_NAMES = EmbeddingMultiTaskRegressor.PROPERTY_NAMES


def load_model(checkpoint_path=CHECKPOINT_PATH_DEFAULT, device=None):
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    checkpoint = torch.load(checkpoint_path, map_location=device)
    embedding_dim = checkpoint["embedding_dim"]

    model = EmbeddingMultiTaskRegressor(embedding_dim=embedding_dim).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    print(
        f"Loaded model from {checkpoint_path} (epoch {checkpoint.get('epoch', '?')}, "
        f"val_loss={checkpoint.get('val_loss', float('nan')):.4f}, embedding_dim={embedding_dim})"
    )
    return model, embedding_dim, device


@torch.no_grad()
def predict_embeddings(model, device, sequences, embeddings, batch_size=64):
    """
    sequences: list of N sequence strings (used only for labeling output rows)
    embeddings: array-like [N, embedding_dim]

    Returns: list of dicts {"sequence": ..., "asphericity": ..., "scaled_re": ...,
                            "scaled_rg": ..., "scaling_exponent": ...}
    """
    embeddings = np.asarray(embeddings, dtype=np.float32)
    results = []

    for start in range(0, len(sequences), batch_size):
        seq_chunk = sequences[start:start + batch_size]
        emb_chunk = torch.tensor(embeddings[start:start + batch_size]).to(device)

        outputs = model(emb_chunk)
        final_pred = outputs["final_pred"].cpu()

        for seq, pred in zip(seq_chunk, final_pred):
            entry = {"sequence": seq}
            for i, name in enumerate(PROPERTY_NAMES):
                entry[name] = float(pred[i])
            results.append(entry)

    return results


def _load_csv(path):
    """col 0 = sequence, cols 1..end = embedding dims."""
    df = pd.read_csv(path, header=0)
    sequences = df.iloc[:, 0].astype(str).tolist()
    embeddings = df.iloc[:, 1:].to_numpy(dtype=np.float32)
    return sequences, embeddings


def main():
    parser = argparse.ArgumentParser(description="Predict protein properties from PLM embeddings.")
    parser.add_argument("--checkpoint", default=CHECKPOINT_PATH_DEFAULT, help="Path to model checkpoint (.pt)")
    parser.add_argument("--csv", required=True, help="CSV: col0=sequence, cols1..end=embedding dims")
    parser.add_argument("--out", help="Optional path to save predictions as CSV")
    args = parser.parse_args()

    model, embedding_dim, device = load_model(args.checkpoint)
    sequences, embeddings = _load_csv(args.csv)

    if embeddings.shape[1] != embedding_dim:
        raise ValueError(
            f"CSV has {embeddings.shape[1]} embedding columns but model expects {embedding_dim}."
        )

    results = predict_embeddings(model, device, sequences, embeddings)

    for r in results:
        print(json.dumps(r, indent=2))

    if args.out:
        pd.DataFrame(results).to_csv(args.out, index=False)
        print(f"Saved predictions to {args.out}")


if __name__ == "__main__":
    main()
