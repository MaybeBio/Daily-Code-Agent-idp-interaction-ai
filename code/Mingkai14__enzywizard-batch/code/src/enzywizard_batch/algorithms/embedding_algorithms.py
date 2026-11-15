from __future__ import annotations

from typing import Dict, List, Any, Optional
import torch
from Bio.Data.IUPACData import protein_letters_1to3
from ..utils.logging_utils import Logger
from ..utils.embedding_utils import load_esm2, postprocess_embedding_report_to_schema
from ..utils.sequence_utils import  normalize_aa_name_to_one_letter

ESM2_EMBEDDING_DIMS = {
    "esm2_t6_8M_UR50D": 320,
    "esm2_t12_35M_UR50D": 480,
    "esm2_t30_150M_UR50D": 640,
}

def generate_embedding(sequence_dict: Dict[str, str],logger: Logger, model_name: str = "esm2_t6_8M_UR50D", device: Optional[str] = None) -> List[Dict[str, Any]] | None:

    # ---------- check input ----------
    if model_name not in ESM2_EMBEDDING_DIMS:
        logger.print("[ERROR] Model name parameter has to be in esm2_t6_8M_UR50D,esm2_t12_35M_UR50D, or esm2_t30_150M_UR50D.")
        return None

    if not sequence_dict:
        logger.print("[ERROR] Sequence_dict is empty.")
        return None

    header = sequence_dict.get("header")
    sequence = sequence_dict.get("sequence")

    if not sequence:
        logger.print("[ERROR] Empty sequence.")
        return None

    sequence = sequence.strip().upper()
    header = header or "protein"

    allowed = set("ACDEFGHIKLMNPQRSTVWY")
    bad = sorted(set(sequence) - allowed)
    if bad:
        logger.print(f"[ERROR] Sequence contains non-standard residues: {bad}")
        return None

    # ---------- device ----------
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    # ---------- load model ----------
    try:
        model, alphabet = load_esm2(model_name)
    except Exception as e:
        logger.print(f"[ERROR] Failed to load ESM2 model: {e}")
        return None

    batch_converter = alphabet.get_batch_converter()
    model = model.to(device)
    model.eval()

    # ---------- tokenize ----------
    data = [(header, sequence)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(device)

    # ---------- layer ----------
    try:
        n_layers = int(model_name.split("_")[1].lstrip("t"))
    except Exception as e:
        logger.print(f"[ERROR] Failed to parse ESM2 layer count from model name {model_name}: {e}")
        n_layers = model.num_layers

    # ---------- inference ----------
    try:
        with torch.no_grad():
            out = model(tokens, repr_layers=[n_layers], return_contacts=False)
            reps = out["representations"][n_layers]  # [1, T, D]
    except Exception as e:
        logger.print(f"[ERROR] ESM2 inference failed: {e}")
        return None

    # ---------- extract ----------
    L = len(sequence)
    per_res = reps[0, 1:1 + L, :].detach().cpu()  # [L, D]
    embedding_dim = per_res.shape[1]
    expected_dim = ESM2_EMBEDDING_DIMS[model_name]
    if embedding_dim != expected_dim:
        logger.print(f"[ERROR] Unexpected embedding dimension for {model_name}: expected {expected_dim}, got {embedding_dim}.")
        return None

    # ---------- format output ----------
    result: List[Dict[str, Any]] = []

    for i in range(L):

        result.append({
            "aa_id": i + 1,
            "aa_name": normalize_aa_name_to_one_letter(protein_letters_1to3[sequence[i]].upper()),
            "embedding": per_res[i].tolist(),
        })

    return result


def generate_embedding_report(embeddings: List[Dict[str, Any]]) -> dict | None:
    raw_report = {
        "output_type": "enzywizard_embedding",
        "embeddings": embeddings,
    }

    return postprocess_embedding_report_to_schema(raw_report)
