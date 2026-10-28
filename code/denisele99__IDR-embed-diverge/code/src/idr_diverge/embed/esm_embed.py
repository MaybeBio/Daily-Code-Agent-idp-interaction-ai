import glob
import math
import multiprocessing
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from Bio import SeqIO


#Sections of code based on implementation from:
#https://github.com/facebookresearch/esm 
#Modified by Denise Le


# ---------------------------
# Global configuration
# ---------------------------

OUTPUT_DIR = Path("/home/moseslab/denise/embeddings/esm1b/")
ESM_MODEL_NAME = "esm1b_t33_650M_UR50S"
ESM_REPR_LAYER = 33 #specific for ESM
ESM_EMBED_DIM = 1280
ESM_MAX_TOKENS = 1022   # max residues per chunk for ESM-1b
BATCH_SIZE = 100        # number of sequences to process before writing to HDF5


# ---------------------------
# Load model once
# ---------------------------

model, alphabet = torch.hub.load("facebookresearch/esm:main","esm1b_t33_650M_UR50S")
batch_converter = alphabet.get_batch_converter()
model = model.eval()



def read_fasta(file:str)-> list: #return list of SeqRecord Objects without the alignments
    all_rec=[]
    for rec in SeqIO.parse(file,"fasta"):
        rec.seq=rec.seq.replace('-','')
        all_rec.append(rec)
    return all_rec

def get_esm_model():
    """
    Prepare the globally loaded ESM model for inference and return it
    together with the batch converter.

    On GPU:move model to CUDA.

    On CPU:enable a few inference optimizations.
    """
    global model, batch_converter

    if torch.cuda.is_available():
        model.to("cuda")
    else:
        torch._C._jit_set_profiling_mode(False)
        torch.set_num_threads(multiprocessing.cpu_count())
        model = torch.jit.freeze(model)
        model = torch.jit.optimize_for_inference(model)

    return model, batch_converter


def embed_sequence_residue_level(sequence: str) -> np.ndarray:
    """
    Generate per-residue embeddings for a single amino acid sequence.

    Returns
    -------
    np.ndarray
        Array of shape (sequence_length, ESM_EMBED_DIM), excluding
        BOS/CLS and EOS tokens.
    """
    model, batch_converter = get_esm_model()

    _, _, tokens = batch_converter([("protein", sequence)])

    if torch.cuda.is_available():
        tokens = tokens.to(device="cuda", non_blocking=True)

    with torch.no_grad():
        results = model(tokens, repr_layers=[ESM_REPR_LAYER])

    residue_embeddings = (
        results["representations"][ESM_REPR_LAYER]
        .to(device="cpu")[0]
        .detach()
        .numpy()
    )

    # Remove BOS/CLS and EOS token embeddings
    residue_embeddings = residue_embeddings[1:-1]

    return residue_embeddings


def embed_long_sequence(sequence: str, max_chunk_length: int = ESM_MAX_TOKENS) -> np.ndarray:
    """
    Embed a sequence at residue level, splitting it into chunks if it is
    longer than the ESM model's maximum input length.

    Returns
    -------
    np.ndarray
        Per-residue embeddings for the full sequence.
    """
    seq_len = len(sequence)

    if seq_len <= max_chunk_length:
        return embed_sequence_residue_level(sequence)

    n_chunks = math.ceil(seq_len / max_chunk_length)
    chunk_size = seq_len / n_chunks

    all_embeddings = []
    for chunk_idx in range(n_chunks):
        start = int(chunk_idx * chunk_size)
        stop = int((chunk_idx + 1) * chunk_size)
        chunk_embeddings = embed_sequence_residue_level(sequence[start:stop])
        all_embeddings.append(chunk_embeddings)

    return np.concatenate(all_embeddings, axis=0)


def load_fasta_sequences(fasta_path: str, is_directory: bool = False):
    """
    Read sequences from either:
    - a single FASTA file
    - all FASTA files inside a directory
    """
    if is_directory:
        all_sequences = []
        fasta_files = glob.glob(f"{fasta_path}/*")
        for fasta_file in fasta_files:
            all_sequences.extend(read_fasta(fasta_file))
        return all_sequences

    return read_fasta(fasta_path)


def resolve_output_path(out_name: str, suffix: str) -> str:
    """
    Build an output path. If out_name already includes a path, use it directly.
    Otherwise, place it inside OUTPUT_DIR.
    """
    if "/" in out_name:
        return f"{out_name}{suffix}"
    return str(OUTPUT_DIR / f"{out_name}{suffix}")


def iter_batches(items, batch_size: int):
    """
    Yield successive batches from a list.
    """
    for start_idx in range(0, len(items), batch_size):
        yield start_idx, items[start_idx:start_idx + batch_size]


def get_sequence_embedding(sequence: str) -> np.ndarray:
    """
    Compute a single per-sequence embedding by averaging per-residue embeddings.
    """
    residue_embeddings = embed_long_sequence(sequence)
    return residue_embeddings.mean(axis=0)


def get_esm_embeddings_csv(fasta_path, out_name, is_directory=False):
    """
    Generate one averaged embedding per sequence and save results to CSV.

    Output format:
        row index = sequence name
        columns   = embedding dimensions
    """
    sequence_records = load_fasta_sequences(fasta_path, is_directory=is_directory)
    output_path = resolve_output_path(out_name, ".csv")

    # Create/clear file
    pd.DataFrame().to_csv(output_path, mode="w", header=False)

    for record in sequence_records:
        sequence_str = str(record.seq)
        sequence_embedding = get_sequence_embedding(sequence_str)

        embed_df = pd.DataFrame(
            sequence_embedding.reshape(1, -1),
            index=[record.name]
        )
        embed_df.to_csv(output_path, mode="a", index=True, header=False)

    return output_path


def get_esm_embeddings_hdf5(fasta_path, out_name, is_directory=False):
    """
    Generate one averaged embedding per sequence and save results to an HDF5 file.

    HDF5 structure:
        /embeddings : float32 array of shape (n_sequences, ESM_EMBED_DIM)
        /ids        : sequence identifiers
    """
    sequence_records = load_fasta_sequences(fasta_path, is_directory=is_directory)
    output_path = resolve_output_path(out_name, ".h5")

    # Initialize output file with extendable datasets
    with h5py.File(output_path, "w") as h5f:
        h5f.create_dataset(
            "embeddings",
            shape=(0, ESM_EMBED_DIM),
            maxshape=(None, ESM_EMBED_DIM),
            dtype="float32",
            chunks=True,
            compression="gzip",
        )
        h5f.create_dataset(
            "ids",
            shape=(0,),
            maxshape=(None,),
            dtype="S100",
            chunks=True,
            compression="gzip",
        )

    total_batches = math.ceil(len(sequence_records) / BATCH_SIZE)

    for batch_num, (start_idx, batch_records) in enumerate(
        iter_batches(sequence_records, BATCH_SIZE),
        start=1
    ):
        print(
            f"Processing batch {batch_num}/{total_batches} "
            f"with {len(batch_records)} sequences..."
        )

        batch_embeddings = []
        batch_ids = []

        for record in batch_records:
            sequence_str = str(record.seq)
            sequence_embedding = get_sequence_embedding(sequence_str)

            batch_embeddings.append(sequence_embedding)
            batch_ids.append(record.name.encode())

        batch_embeddings = np.asarray(batch_embeddings, dtype=np.float32)
        batch_ids = np.asarray(batch_ids, dtype="S100")

        end_idx = start_idx + len(batch_records)

        with h5py.File(output_path, "a") as h5f:
            embedding_ds = h5f["embeddings"]
            ids_ds = h5f["ids"]

            embedding_ds.resize((end_idx, ESM_EMBED_DIM))
            ids_ds.resize((end_idx,))

            embedding_ds[start_idx:end_idx] = batch_embeddings
            ids_ds[start_idx:end_idx] = batch_ids

    return output_path


if __name__ == "__main__":
    fasta_path = sys.argv[1]
    out_name = sys.argv[2]

    # Example:
    # python script.py input.fasta output_name
    get_esm_embeddings_hdf5(
        fasta_path=fasta_path,
        out_name=out_name,
        is_directory=False,
    )