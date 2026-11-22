"""
Training pipeline for EmbeddingMultiTaskRegressor.

Expects a CSV at DATA_CSV_PATH with a header row and columns in this order:
    col 0:      sequence
    col 1:      asphericity
    col 2:      prefactor
    col 3:      RE  (scaled end-to-end distance)
    col 4:      RG  (scaled radius of gyration)
    col 5:      scaling exponent
    col 6..end: PLM embedding dimensions (e.g. emb_0, emb_1, ..., emb_{D-1})

Column NAMES are ignored except to size the embedding — only column POSITION
is used, so the header can be anything.

Outputs (written to OUTPUT_DIR):
    best_model.pt   - checkpoint with lowest validation loss
    last_model.pt   - checkpoint after final epoch
    losses.csv      - per-epoch train/val total/final/aux losses
    loss_curve.png  - plot of train vs val loss over epochs
"""

import os
import csv

import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader, random_split
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model_2 import EmbeddingMultiTaskRegressor, compute_loss

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA_CSV_PATH = "../data/train_with_prop_emb.csv"
OUTPUT_DIR = "sequence_500_multi_prop_regression"

VAL_FRACTION = 0.2
BATCH_SIZE = 512
NUM_EPOCHS = 2500
LEARNING_RATE = 1e-4
WEIGHT_DECAY = 1e-5
AUX_LOSS_WEIGHT = 0.1
SEED = 42

os.makedirs(OUTPUT_DIR, exist_ok=True)
torch.manual_seed(SEED)


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

# ==========================================================
# Dataset
# ==========================================================

class EmbeddingPropertyDataset(Dataset):
    """
    Dataset returning:
        embedding : FloatTensor [embedding_dim]
        labels    : FloatTensor [5]
    """

    PROPERTY_NAMES = ["asphericity","prefactor","scaled_re","scaled_rg","scaling_exponent"]

    def __init__(self, records):
        """
        records : list of dictionaries
        """
        self.records = records

    def __len__(self):
        return len(self.records)

    def __getitem__(self, idx):
        rec = self.records[idx]

        embedding = torch.tensor(rec["embedding"],dtype=torch.float32)

        labels = torch.tensor([rec[name] for name in self.PROPERTY_NAMES],dtype=torch.float32)

        return {"embedding": embedding,"labels": labels}


# ==========================================================
# Collate Function
# ==========================================================

def collate_fn(batch):
    embeddings = torch.stack([b["embedding"] for b in batch])
    labels = torch.stack([b["labels"] for b in batch])

    return {
        "embedding": embeddings,
        "labels": labels,
    }



# ==========================================================
# Load records for training
# ==========================================================

def load_records(csv_path):
    """
    Loads merged CSV into a list of dictionaries.

    CSV format:

    Sequence,
    Asphericity,
    Prefactor,
    Scaled_RE,
    Scaled_RG,
    Scaling_Exponent,
    emb_0,
    emb_1,
    ...
    """

    df = pd.read_csv(csv_path)

    if df.shape[1] < 7:
        raise ValueError(
            "CSV should contain:\n"
            "Sequence, Asphericity, Prefactor, "
            "Scaled_RE, Scaled_RG, Scaling_Exponent, "
            "embedding columns...")

    embedding_dim = df.shape[1] - 6

    records = []

    for _, row in df.iterrows():

        records.append(
            {
                "sequence": str(row["sequence"]).strip(),
                "asphericity": float(row["Asphericity"]),
                "prefactor": float(row["Prefactor"]),
                "scaled_re": float(row["Scaled_RE"]),
                "scaled_rg": float(row["Scaled_RG"]),
                "scaling_exponent": float(row["Scaling_exponent"]),
                "embedding": row.iloc[6:].to_numpy(dtype=np.float32),
            }
        )

    return records, embedding_dim


# ---------------------------------------------------------------------------
# Train / eval loop
# ---------------------------------------------------------------------------

def run_epoch(model, loader, device, optimizer=None, aux_weight=AUX_LOSS_WEIGHT):
    is_train = optimizer is not None
    model.train() if is_train else model.eval()

    total_loss_sum = final_loss_sum = aux_loss_sum = 0.0
    n_batches = 0

    context = torch.enable_grad() if is_train else torch.no_grad()
    with context:
        for batch in loader:
            embedding = batch["embedding"].to(device)
            labels = batch["labels"].to(device)

            if is_train:
                optimizer.zero_grad()

            outputs = model(embedding)
            total_loss, final_loss, aux_loss = compute_loss(outputs, labels, aux_weight)

            if is_train:
                total_loss.backward()
                optimizer.step()

            total_loss_sum += total_loss.item()
            final_loss_sum += final_loss.item()
            aux_loss_sum += aux_loss.item()
            n_batches += 1

    return (
        total_loss_sum / n_batches,
        final_loss_sum / n_batches,
        aux_loss_sum / n_batches,
    )


def main():
    device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    records, embedding_dim = load_records(DATA_CSV_PATH)
    print(f"Loaded {len(records)} sequences from {DATA_CSV_PATH}, embedding_dim={embedding_dim}")

    full_dataset = EmbeddingPropertyDataset(records)

    val_size = max(1, int(len(full_dataset) * VAL_FRACTION))
    train_size = len(full_dataset) - val_size
    train_ds, val_ds = random_split(
        full_dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(SEED),
    )
    print(f"Train size: {train_size}, Val size: {val_size}")

    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True, collate_fn=collate_fn)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate_fn)

    model = EmbeddingMultiTaskRegressor(embedding_dim=embedding_dim).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)

    history = []
    best_val_loss = float("inf")

    for epoch in range(1, NUM_EPOCHS + 1):
        train_total, train_final, train_aux = run_epoch(model, train_loader, device, optimizer=optimizer)
        val_total, val_final, val_aux = run_epoch(model, val_loader, device, optimizer=None)

        history.append({
            "epoch": epoch,
            "train_total_loss": train_total, "train_final_loss": train_final, "train_aux_loss": train_aux,
            "val_total_loss": val_total, "val_final_loss": val_final, "val_aux_loss": val_aux,
        })

        print(
            f"Epoch {epoch:3d}/{NUM_EPOCHS} | "
            f"train_loss={train_total:.4f} (final={train_final:.4f}, aux={train_aux:.4f}) | "
            f"val_loss={val_total:.4f} (final={val_final:.4f}, aux={val_aux:.4f})" , flush=True
        )

        if val_total < best_val_loss:
            best_val_loss = val_total
            torch.save(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "val_loss": val_total,
                    "embedding_dim": embedding_dim,
                },
                os.path.join(OUTPUT_DIR, "best_model.pt"),
            )

    torch.save(
        {
            "epoch": NUM_EPOCHS,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "val_loss": val_total,
            "embedding_dim": embedding_dim,
        },
        os.path.join(OUTPUT_DIR, "last_model.pt"),
    )

    losses_csv_path = os.path.join(OUTPUT_DIR, "losses.csv")
    with open(losses_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)
    print(f"Saved loss history to {losses_csv_path}")

    epochs = [h["epoch"] for h in history]
    train_losses = [h["train_total_loss"] for h in history]
    val_losses = [h["val_total_loss"] for h in history]

    plt.figure(figsize=(8, 5))
    plt.plot(epochs, train_losses, label="Train Loss")
    plt.plot(epochs, val_losses, label="Val Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Total Loss (final + aux)")
    plt.title("Training vs Validation Loss")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plot_path = os.path.join(OUTPUT_DIR, "loss_curve.png")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    print(f"Saved loss curve to {plot_path}")

    print(f"Best val loss: {best_val_loss:.4f}")
    print(f"Best model saved to {os.path.join(OUTPUT_DIR, 'best_model.pt')}")


if __name__ == "__main__":
    main()
