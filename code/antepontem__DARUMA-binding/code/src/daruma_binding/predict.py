"""
This script takes a FASTA file containing one or more amino acid sequences and predicts binding regions using DARUMA-binding.

Execution command:
    Run on CPU:
        `python predict.py -i input.fasta -o output.txt -d cpu`
    Run on GPU:
        `python predict.py -i input.fasta -o output.txt -d 0`

Arguments:
    -i, --input: Path to the input FASTA file (e.g., protein1.fasta)
    -o, --output: Path to the output prediction file (e.g., protein1.out)
    -d, --device: Device to use (e.g., cpu, 0, 1, ... for GPU)

Output:
    The prediction probabilities will be written to the output file specified in the arguments.
"""

import sys,os
import argparse
import torch
import numpy as np

dir = os.path.dirname(__file__)
sys.path.append(dir)

from backend.models import daruma_binding
from backend.functions import read_fasta, SequenceFeaturizer, ResultWriter


def parse_args():
    parser = argparse.ArgumentParser(description="Run prediction using ProtT5 embedding vectors.")
    parser.add_argument("-i", "--input", required=False, help="Input FASTA file, e.g. protein1.fasta")
    parser.add_argument("-o", "--output", required=True, help="Output prediction file, e.g. protein1.out")
    parser.add_argument("-d", "--device", default="0", help="Use device, e.g. cpu, 0, 1, ... for GPU")
    return parser.parse_args()

args = parse_args()
input_fasta_path = args.input
output_file_path = args.output
target_device = args.device


if target_device == "cpu":
    device = torch.device("cpu")
elif torch.cuda.is_available() and target_device != "cpu":
    torch.cuda.set_device(int(target_device))
    device = torch.device(f"cuda:{target_device}")
else:
    device = torch.device("cpu")
print(f"Using {device} device.")


# Load input data
input_datas = read_fasta(input_fasta_path)

# Initialize the featurizer
featurizer = SequenceFeaturizer(feature_type="ProtT5", device=device)

# Load the model
model = daruma_binding(device,featurizer=featurizer)

result_writer = ResultWriter(output_file_path, model.threshold)


for ac, seq in input_datas.items():
    print(f"Processing {ac} ...")

    result = model.predict_from_seq(seq)

    result_writer.append_write(ac, result, model.threshold)

