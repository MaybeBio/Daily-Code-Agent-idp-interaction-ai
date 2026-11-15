from __future__ import annotations
import torch
import esm
from typing import Tuple, Any, Dict, List
from functools import lru_cache

@lru_cache(maxsize=4)
def load_esm2(model_name: str = "esm2_t6_8M_UR50D") -> Tuple[torch.nn.Module, esm.data.Alphabet]:
    """
    Load ESM2 model + alphabet once (cached).
    """
    loader = getattr(esm.pretrained, model_name, None)
    if loader is None:
        raise ValueError(f"Unknown ESM2 model name: {model_name}")

    model, alphabet = loader()
    model.eval()
    return model, alphabet


def postprocess_embedding_report_to_schema(
    raw_report: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Postprocess the raw EnzyWizard-Embedding report into the official JSON Schema field names.
    """

    sequence_embeddings: List[Dict[str, Any]] = []

    for raw_embedding in raw_report["embeddings"]:
        sequence_embeddings.append(
            {
                "residue_index": raw_embedding["aa_id"],
                "residue_name": raw_embedding["aa_name"],
                "residue_embedding": raw_embedding["embedding"],
            }
        )

    schema_report: Dict[str, Any] = {
        "report_type": raw_report["output_type"],
        "sequence_embeddings": sequence_embeddings,
    }

    return schema_report