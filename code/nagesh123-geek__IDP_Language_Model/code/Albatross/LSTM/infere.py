"""
Inference script for the BiLSTM multi-task protein property regressor.

Usage (predict on a CSV of sequences):
    python infer.py \
        --checkpoint checkpoints/best_model.pt \
        --input_csv  path/to/sequences.csv \
        --output_csv predictions.csv

    The input CSV must have a 'Sequence' column (any other columns are
    ignored/passed through).

Usage (predict on sequences given directly on the command line):
    python infer.py \
        --checkpoint checkpoints/best_model.pt \
        --sequences MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR \
                    GQPLYIQAEVDDTVEDNLKLYQKLKATSSLL
"""

import argparse
import csv

import torch

from lstm import (ProteinMultiTaskRegressor,encode_sequence)


def load_model(checkpoint_path, device):
    ckpt = torch.load(checkpoint_path, map_location=device)
    config = ckpt["config"]

    model = ProteinMultiTaskRegressor(
        vocab_size=config["vocab_size"],
        hidden_size=config["hidden_size"],
        num_layers=config["num_layers"],
        dropout=config["dropout"],
        pad_idx=config["pad_idx"],
    ).to(device)
    model.load_state_dict(ckpt["model_state"])
    model.eval()

    return model, config["max_len"]


@torch.no_grad()
def predict_batch(model, sequences, max_len, device, batch_size=32):
    """
    sequences: list[str]
    Returns: list[dict[property_name -> float]]
    """
    property_names = ProteinMultiTaskRegressor.PROPERTY_NAMES
    results = []

    for start in range(0, len(sequences), batch_size):
        chunk = sequences[start:start + batch_size]

        input_ids_list, attn_list, len_list = [], [], []
        for seq in chunk:
            ids, attn, length = encode_sequence(seq, max_len)
            input_ids_list.append(ids)
            attn_list.append(attn)
            len_list.append(length)

        input_ids = torch.stack(input_ids_list).to(device)
        attention_mask = torch.stack(attn_list).to(device)
        lengths = torch.tensor(len_list, dtype=torch.long)  # CPU, for pack_padded_sequence

        outputs = model(input_ids, attention_mask, lengths)
        final_pred = outputs["final_pred"].cpu().numpy()  # [B, 5]

        for row in final_pred:
            results.append({name: float(val) for name, val in zip(property_names, row)})

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=str, required=True)
    parser.add_argument("--input_csv", type=str, default=None,
                         help="CSV file with a 'Sequence' column.")
    parser.add_argument("--sequences", type=str, nargs="+", default=None,
                         help="One or more sequences given directly on the CLI.")
    parser.add_argument("--output_csv", type=str, default="predictions.csv")
    parser.add_argument("--batch_size", type=int, default=32)
    args = parser.parse_args()

    if not args.input_csv and not args.sequences:
        raise ValueError("Provide either --input_csv or --sequences.")

    device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    model, max_len = load_model(args.checkpoint, device)
    print(f"Loaded model (max_len={max_len}) from {args.checkpoint}")

    # ---------------- Gather sequences ----------------
    rows = []
    if args.input_csv:
        with open(args.input_csv, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows.append(dict(row))
    else:
        rows = [{"Sequence": s} for s in args.sequences]

    sequences = [r["Sequence"].strip() for r in rows]

    # Warn (but still run) if any sequence would be truncated to max_len
    n_over = sum(1 for s in sequences if len(s) + 1 > max_len)
    if n_over:
        print(
            f"Warning: {n_over} sequence(s) exceed the model's training max_len "
            f"({max_len} incl. CLS) and will be truncated."
        )

    # ---------------- Predict ----------------
    predictions = predict_batch(model, sequences, max_len, device, args.batch_size)

    property_names = ProteinMultiTaskRegressor.PROPERTY_NAMES
    for row, pred in zip(rows, predictions):
        row.update(pred)

    # ---------------- Write output ----------------
    fieldnames = list(rows[0].keys())
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} predictions to {args.output_csv}")

    # Also print a small preview to stdout
    preview_n = min(5, len(rows))
    print(f"\nPreview (first {preview_n} rows):")
    for row in rows[:preview_n]:
        seq_preview = row["Sequence"][:30] + ("..." if len(row["Sequence"]) > 30 else "")
        props = ", ".join(f"{name}={row[name]:.4f}" for name in property_names)
        print(f"  {seq_preview}: {props}")


if __name__ == "__main__":
    main()
