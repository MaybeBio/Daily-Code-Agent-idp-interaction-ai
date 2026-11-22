"""
Training script for the BiLSTM multi-task protein property regressor.

Expected input: a CSV file with columns:
    Sequence, Asphericity, Prefactor, Scaled_RE, Scaled_RG, Scaling_Exponent

Usage:
    python train.py \
        --data_csv path/to/train_data.csv \
        --out_dir  checkpoints \
        --epochs 50 \
        --batch_size 32 \
        --lr 1e-3 \
        --val_frac 0.1

Saves:
    <out_dir>/best_model.pt   -- state dict + config + vocab, best val loss
    <out_dir>/last_model.pt   -- state dict + config + vocab, final epoch
"""

import argparse
import csv
import os
import random

import torch
from torch.utils.data import DataLoader, random_split

from lstm import (
    ProteinMultiTaskRegressor,
    ProteinPropertyDataset,
    collate_fn,
    compute_loss,
    VOCAB_SIZE,
    PAD_IDX,
)


def set_seed(seed: int):
    random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def load_records_from_csv(path):
    property_names = ProteinPropertyDataset.PROPERTY_NAMES
    records = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rec = {"Sequence": row["Sequence"].strip()}
            for name in property_names:
                rec[name] = float(row[name])
            records.append(rec)
    return records


def run_epoch(model, loader, device, optimizer=None, aux_weight=0.1):
    is_train = optimizer is not None
    model.train() if is_train else model.eval()

    total_loss, total_final, total_aux, n_batches = 0.0, 0.0, 0.0, 0

    with torch.set_grad_enabled(is_train):
        for batch in loader:
            input_ids = batch["input_ids"].to(device)
            attention_mask = batch["attention_mask"].to(device)
            lengths = batch["lengths"]  # kept on CPU for pack_padded_sequence
            labels = batch["labels"].to(device)

            outputs = model(input_ids, attention_mask, lengths)
            loss, final_loss, aux_loss = compute_loss(outputs, labels, aux_weight)

            if is_train:
                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()

            total_loss += loss.item()
            total_final += final_loss.item()
            total_aux += aux_loss.item() if isinstance(aux_loss, torch.Tensor) else aux_loss
            n_batches += 1

    return total_loss / n_batches, total_final / n_batches, total_aux / n_batches


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_csv", type=str, required=True)
    parser.add_argument("--out_dir", type=str, default="checkpoints")
    parser.add_argument("--epochs", type=int, default=3000)
    parser.add_argument("--batch_size", type=int, default=512)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--val_frac", type=float, default=0.2)
    parser.add_argument("--aux_weight", type=float, default=0.1)
    parser.add_argument("--max_len", type=int, default=None,
                         help="Override max sequence length (incl. CLS). "
                              "Defaults to longest sequence in the data + 1.")
    parser.add_argument("--hidden_size", type=int, default=256)
    parser.add_argument("--num_layers", type=int, default=4)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--patience", type=int, default=1000,
                         help="Early-stopping patience on val loss (epochs).")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    set_seed(args.seed)
    os.makedirs(args.out_dir, exist_ok=True)
    device = torch.device("cuda:7" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # ---------------- Data ----------------
    records = load_records_from_csv(args.data_csv)
    print(f"Loaded {len(records)} records from {args.data_csv}")

    dataset = ProteinPropertyDataset(records, max_len=args.max_len)
    max_len = dataset.max_len
    print(f"max_len (incl. CLS) = {max_len}")

    n_val = max(1, int(len(dataset) * args.val_frac))
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(
        dataset, [n_train, n_val],
        generator=torch.Generator().manual_seed(args.seed),
    )
    print(f"Train: {n_train} | Val: {n_val}")

    train_loader = DataLoader(
        train_set, batch_size=args.batch_size, shuffle=True, collate_fn=collate_fn
    )
    val_loader = DataLoader(
        val_set, batch_size=args.batch_size, shuffle=False, collate_fn=collate_fn
    )

    # ---------------- Model ----------------
    model = ProteinMultiTaskRegressor(
        vocab_size=VOCAB_SIZE,
        hidden_size=args.hidden_size,
        num_layers=args.num_layers,
        dropout=args.dropout,
        pad_idx=PAD_IDX,
    ).to(device)

    print(f"Using one-hot token encoding; LSTM input_size = {model.encoder.input_dim}")

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Model parameters: {n_params:,}")

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.lr, weight_decay=args.weight_decay
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=3
    )

    # ---------------- Train loop ----------------
    best_val_loss = float("inf")
    epochs_no_improve = 0

    config = {
        "vocab_size": VOCAB_SIZE,
        "hidden_size": args.hidden_size,
        "num_layers": args.num_layers,
        "dropout": args.dropout,
        "pad_idx": PAD_IDX,
        "max_len": max_len,
    }

    for epoch in range(1, args.epochs + 1):
        train_loss, train_final, train_aux = run_epoch(
            model, train_loader, device, optimizer, args.aux_weight
        )
        val_loss, val_final, val_aux = run_epoch(
            model, val_loader, device, optimizer=None, aux_weight=args.aux_weight
        )
        scheduler.step(val_loss)
        current_lr = optimizer.param_groups[0]["lr"]

        print(
            f"Epoch {epoch:3d}/{args.epochs} | "
            f"train_loss {train_loss:.4f} (final {train_final:.4f}, aux {train_aux:.4f}) | "
            f"val_loss {val_loss:.4f} (final {val_final:.4f}, aux {val_aux:.4f}) | "
            f"lr {current_lr:.2e}" , flush=True
        )

        # Save "last" checkpoint every epoch
        torch.save(
            {"model_state": model.state_dict(), "config": config},
            os.path.join(args.out_dir, "last_model.pt"),
        )

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            epochs_no_improve = 0
            torch.save(
                {"model_state": model.state_dict(), "config": config},
                os.path.join(args.out_dir, "best_model.pt"),
            )
            print(f"  -> new best val_loss {best_val_loss:.4f}, checkpoint saved")
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= args.patience:
                print(f"Early stopping: no val improvement for {args.patience} epochs.")
                break

    print(f"Training complete. Best val_loss = {best_val_loss:.4f}")
    print(f"Checkpoints saved to: {args.out_dir}")


if __name__ == "__main__":
    main()
