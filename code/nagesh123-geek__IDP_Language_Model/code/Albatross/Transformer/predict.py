"""
Inference script for ProteinMultiTaskRegressor.

Usage 
--------------
Single sequence, from the command line:
    python predict.py --sequence "MKTAYIAKQRQISFVKSHFSRQ"

Multiple sequences from a text file (one sequence per line):
    python predict.py --file sequences.txt

Multiple sequences from a CSV (first column = sequence):
    python predict.py --csv new_sequences.csv --out predictions.csv

As a library, from your own code:
    from predict import load_model, predict_sequences
    model, max_len, device = load_model("best_model.pt")
    results = predict_sequences(model, max_len, device, ["MKTAYIAK..."])
"""

import argparse
import csv
import json

import torch

from protein_transformer_regression import (
    ProteinMultiTaskRegressor,
    encode_sequence,
)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"     # 20 standard residues
EXTRA_TOKENS = "XBZUO"                    # non-standard / ambiguous residues
PAD_TOKEN = "<PAD>"
CLS_TOKEN = "<CLS>"

VOCAB = [PAD_TOKEN, CLS_TOKEN] + list(AMINO_ACIDS) + list(EXTRA_TOKENS)
TOKEN_TO_IDX = {tok: i for i, tok in enumerate(VOCAB)}
PAD_IDX = TOKEN_TO_IDX[PAD_TOKEN]
CLS_IDX = TOKEN_TO_IDX[CLS_TOKEN]
UNK_IDX = TOKEN_TO_IDX["X"]
VOCAB_SIZE = len(VOCAB)



CHECKPOINT_PATH_DEFAULT = "seq_to_property_albatross_model/best_model.pt"
PROPERTY_NAMES = ProteinMultiTaskRegressor.PROPERTY_NAMES  # order of outputs


def load_model(checkpoint_path=CHECKPOINT_PATH_DEFAULT, device=None):
    """Load a trained model from a checkpoint saved by train.py."""
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    checkpoint = torch.load(checkpoint_path, map_location=device)
    max_len = checkpoint["max_len"]

    model = ProteinMultiTaskRegressor(max_len=max_len).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    print(f"Loaded model from {checkpoint_path} (epoch {checkpoint.get('epoch', '?')}, "
          f"val_loss={checkpoint.get('val_loss', float('nan')):.4f}, max_len={max_len})")

    return model, max_len, device


@torch.no_grad()
def predict_sequences(model, max_len, device, sequences, batch_size=32):
    """
    sequences: list of raw protein sequence strings (no CLS/PAD needed, handled internally)

    Returns: list of dicts, one per sequence, e.g.
        {"sequence": "...", "asphericity": 0.42, "scaled_re": 1.05,
         "scaled_rg": 0.88, "scaling_exponent": 0.55}
    """
    results = []

    for start in range(0, len(sequences), batch_size):
        chunk = sequences[start:start + batch_size]

        input_ids_list = []
        attn_list = []
        for seq in chunk:
            seq = seq.strip().upper()
            if len(seq) + 1 > max_len:
                raise ValueError(
                    f"Sequence of length {len(seq)} + CLS exceeds model's max_len={max_len}. "
                    f"Sequence: {seq[:30]}..."
                )
            ids, attn = encode_sequence(seq, max_len)
            input_ids_list.append(ids)
            attn_list.append(attn)

        input_ids = torch.stack(input_ids_list).to(device)
        attention_mask = torch.stack(attn_list).to(device)

        outputs = model(input_ids, attention_mask)
        final_pred = outputs["final_pred"].cpu()  # [B, 4], order = PROPERTY_NAMES

        for seq, pred in zip(chunk, final_pred):
            entry = {"sequence": seq}
            for i, name in enumerate(PROPERTY_NAMES):
                entry[name] = float(pred[i])
            results.append(entry)

    return results


def _read_sequences_from_file(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def _read_sequences_from_csv(path):
    seqs = []
    with open(path) as f:
        reader = csv.reader(f)
        rows = list(reader)
    # skip header if the first cell of the first row doesn't look like a sequence
    start_idx = 1 if rows and not rows[0][0].strip().isalpha() is False and rows[0][0].strip().upper() != rows[0][0].strip() else 0
    # Simpler/robust rule: try to detect header by checking if row0 col0 contains only
    # standard amino-acid letters; if not, treat as header.
    valid_aa = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
    first_cell = rows[0][0].strip().upper() if rows else ""
    has_header = not (first_cell and all(c in valid_aa for c in first_cell))
    data_rows = rows[1:] if has_header else rows
    for row in data_rows:
        if row and row[0].strip():
            seqs.append(row[0].strip())
    return seqs


def main():
    parser = argparse.ArgumentParser(description="Predict protein properties from sequence(s).")
    parser.add_argument("--checkpoint", default=CHECKPOINT_PATH_DEFAULT, help="Path to model checkpoint (.pt)")
    parser.add_argument("--sequence", help="A single protein sequence")
    parser.add_argument("--file", help="Text file with one sequence per line")
    parser.add_argument("--csv", help="CSV file with sequences in the first column")
    parser.add_argument("--out", help="Optional path to save predictions as CSV")
    args = parser.parse_args()

    if not (args.sequence or args.file or args.csv):
        parser.error("Provide one of --sequence, --file, or --csv")

    sequences = []
    if args.sequence:
        sequences.append(args.sequence)
    if args.file:
        sequences.extend(_read_sequences_from_file(args.file))
    if args.csv:
        sequences.extend(_read_sequences_from_csv(args.csv))

    model, max_len, device = load_model(args.checkpoint)
    results = predict_sequences(model, max_len, device, sequences)

    for r in results:
        print(json.dumps(r, indent=2))

    if args.out:
        with open(args.out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["sequence"] + PROPERTY_NAMES)
            writer.writeheader()
            writer.writerows(results)
        print(f"Saved predictions to {args.out}")


if __name__ == "__main__":
    main()
