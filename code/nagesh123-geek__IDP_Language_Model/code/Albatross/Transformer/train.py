"""
Training pipeline for ProteinMultiTaskRegressor.

Expects a CSV at DATA_CSV_PATH with a header row and columns in this order:
    col 0: sequence
    col 1: asphericity
    col 2: RE   (scaled end-to-end distance)
    col 3: RG   (scaled radius of gyration)
    col 4: scaling exponent

Column NAMES in the header are ignored; only column POSITION is used

Outputs ( written to OUTPUT_DIR):
    vocab.json        - token -> id mapping used for encoding sequences
    best_model.pt      - model checkpoint with lowest validation loss
    last_model.pt      - model checkpoint after final epoch
    losses.csv          - per-epoch train/val total/final/aux losses
    loss_curve.png      - plot of train vs val loss over epochs
"""

import os
import json
import csv

import pandas as pd
import torch
from torch.utils.data import DataLoader, random_split
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from protein_transformer_regression import (
    ProteinPropertyDataset,
    ProteinMultiTaskRegressor,
    collate_fn,
    compute_loss,
    VOCAB,
    TOKEN_TO_IDX,
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATA_CSV_PATH = "../data/train_all_prop_minmax.csv"
OUTPUT_DIR = "with_minmax"

VAL_FRACTION = 0.2
BATCH_SIZE = 256
NUM_EPOCHS = 2500
LEARNING_RATE = 1e-4
WEIGHT_DECAY = 1e-5
AUX_LOSS_WEIGHT = 0.1
SEED = 42

os.makedirs(OUTPUT_DIR, exist_ok=True)
torch.manual_seed(SEED)


def load_records(csv_path):
    """Read data.csv by column POSITION (header names are ignored)."""
    df = pd.read_csv(csv_path, header=0)
    if df.shape[1] < 5:
        raise ValueError(
            f"Expected at least 5 columns (sequence, asphericity, RE, RG, "
            f"scaling_exponent) but got {df.shape[1]} in {csv_path}"
        )

    records = []
    for _, row in df.iterrows():
        seq = str(row.iloc[0]).strip().upper()
        records.append(
            {
                "Sequence": seq,
                "Asphericity": float(row.iloc[1]),
                "Prefactor":float(row.iloc[2]),
                "Scaled_RE": float(row.iloc[3]),
                "Scaled_RG": float(row.iloc[4]),
                "Scaling_Exponent": float(row.iloc[5]),
            }
        )
    return records


def save_vocab(path):
    with open(path, "w") as f:
        json.dump(
            {"vocab": VOCAB, "token_to_idx": TOKEN_TO_IDX},
            f,
            indent=2,
        )


def run_epoch(model, loader, device, optimizer=None, aux_weight=AUX_LOSS_WEIGHT):
    """One pass over the data. If optimizer is provided, trains; else evaluates."""
    is_train = optimizer is not None
    model.train() if is_train else model.eval()

    total_loss_sum = 0.0
    final_loss_sum = 0.0
    aux_loss_sum = 0.0
    n_batches = 0

    context = torch.enable_grad() if is_train else torch.no_grad()
    with context:
        for batch in loader:
            input_ids = batch["input_ids"].to(device)
            attention_mask = batch["attention_mask"].to(device)
            labels = batch["labels"].to(device)

            if is_train:
                optimizer.zero_grad()

            outputs = model(input_ids, attention_mask)
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
    device = torch.device("cuda:2" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # --- Data ---
    records = load_records(DATA_CSV_PATH)
    print(f"Loaded {len(records)} sequences from {DATA_CSV_PATH}")

    full_dataset = ProteinPropertyDataset(records)
    max_len = full_dataset.max_len
    print(f"Max sequence length (incl. CLS): {max_len}")

    val_size = max(1, int(len(full_dataset) * VAL_FRACTION))
    train_size = len(full_dataset) - val_size
    train_ds, val_ds = random_split(
        full_dataset,
        [train_size, val_size],
        generator=torch.Generator().manual_seed(SEED),
    )
    print(f"Train size: {train_size}, Val size: {val_size}")

    train_loader = DataLoader(
        train_ds, batch_size=BATCH_SIZE, shuffle=True, collate_fn=collate_fn
    )
    val_loader = DataLoader(
        val_ds, batch_size=BATCH_SIZE, shuffle=False, collate_fn=collate_fn
    )

    # --- Model ---
    model = ProteinMultiTaskRegressor(max_len=max_len).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )

    # --- Save vocab up front (doesn't depend on training) ---
    save_vocab(os.path.join(OUTPUT_DIR, "vocab.json"))

    # --- Training loop ---
    history = []
    best_val_loss = float("inf")

    for epoch in range(1, NUM_EPOCHS + 1):
        train_total, train_final, train_aux = run_epoch(
            model, train_loader, device, optimizer=optimizer
        )
        val_total, val_final, val_aux = run_epoch(
            model, val_loader, device, optimizer=None
        )

        history.append(
            {
                "epoch": epoch,
                "train_total_loss": train_total,
                "train_final_loss": train_final,
                "train_aux_loss": train_aux,
                "val_total_loss": val_total,
                "val_final_loss": val_final,
                "val_aux_loss": val_aux,
            }
        )

        print(
            f"Epoch {epoch:3d}/{NUM_EPOCHS} | "
            f"train_loss={train_total:.4f} (final={train_final:.4f}, aux={train_aux:.4f}) | "
            f"val_loss={val_total:.4f} (final={val_final:.4f}, aux={val_aux:.4f})" , flush=True
        )

        # Save best model by validation total loss
        if val_total < best_val_loss:
            best_val_loss = val_total
            torch.save(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "val_loss": val_total,
                    "max_len": max_len,
                },
                os.path.join(OUTPUT_DIR, "best_model.pt"),
            )

    # Save final-epoch model too
    torch.save(
        {
            "epoch": NUM_EPOCHS,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "val_loss": val_total,
            "max_len": max_len,
        },
        os.path.join(OUTPUT_DIR, "last_model.pt"),
    )

    # --- Save loss history to CSV ---
    losses_csv_path = os.path.join(OUTPUT_DIR, "losses.csv")
    with open(losses_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)
    print(f"Saved loss history to {losses_csv_path}")

    # --- Plot loss curve ---
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
