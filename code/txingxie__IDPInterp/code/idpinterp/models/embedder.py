import logging
import os
import sys
from pathlib import Path

import torch
from transformers import AutoModel, AutoTokenizer

logger = logging.getLogger(__name__)

# ESM2-35M is the backbone SeqDance / ESMDance are built on (its tokenizer).
SEQDANCE_BACKBONE = "facebook/esm2_t12_35M_UR50D"


class PLMEmbedder:
    """Base wrapper for generating embeddings from Protein Language Models."""
    def __init__(self, model_name: str, hf_path: str, device: str, max_seq_len: int):
        self.model_name = model_name
        self.hf_path = hf_path
        self.device = device
        self.max_seq_len = max_seq_len

        logger.info(f"Loading {hf_path} to {device}...")
        self.tokenizer = AutoTokenizer.from_pretrained(hf_path)
        self.model = AutoModel.from_pretrained(hf_path).to(device).eval()

    def get_hidden_dim(self) -> int:
        return self.model.config.hidden_size

    def encode_batch(self, sequences: list[str], layer: int) -> tuple[torch.Tensor, list[int]]:
        """
        Encode a batch of sequences and extract hidden states from a specific layer.
        
        Returns:
            Flat tensor of embeddings (shape: [total_residues, hidden_dim])
            List of actual sequence lengths (to track where each sequence ends)
        """
        # Truncate sequences to max length
        trunc_seqs = [seq[:self.max_seq_len] for seq in sequences]
        
        inputs = self.tokenizer(
            trunc_seqs, 
            return_tensors="pt", 
            padding=True, 
            truncation=True,
            max_length=self.max_seq_len + 2  # +2 for <cls> and <eos>
        ).to(self.device)

        with torch.no_grad():
            outputs = self.model(
                **inputs, 
                output_hidden_states=True, 
                return_dict=True
            )

        # Layer 0 is embeddings, Layer 1 is first transformer block
        hidden_states = outputs.hidden_states[layer]

        embeddings = []
        lengths = []
        
        for i, seq in enumerate(trunc_seqs):
            seq_len = len(seq)
            # Slice out <cls> at index 0 and <eos>/padding at the end
            emb = hidden_states[i, 1 : seq_len + 1, :].cpu()
            embeddings.append(emb)
            lengths.append(seq_len)

        flat_embeddings = torch.cat(embeddings, dim=0)
        return flat_embeddings, lengths


class ESM2Embedder(PLMEmbedder):
    """Specific implementation for ESM-2 models (already compatible with base)."""
    pass


class SeqDanceEmbedder(PLMEmbedder):
    """
    Embedder for SeqDance / ESMDance (Hou & Shen, PNAS 2026).

    SeqDance is `ESMwrap(nn.Module, PyTorchModelHubMixin)` from the ShenLab/SeqDance
    repo; it wraps a standard transformers ``EsmModel`` (ESM2-35M architecture,
    12 layers, hidden 480) as ``self.esm2`` plus dynamics-prediction heads. For
    interpretability we only need the inner transformer's per-residue hidden
    states, so we pull out ``.esm2`` and reuse the base ``encode_batch`` (which
    calls ``output_hidden_states=True`` and slices off <cls>/<eos>).

    Setup: clone the repo so ``ESMwrap`` is importable, then point the embedder at it:
        git clone https://github.com/ShenLab/SeqDance
        export SEQDANCE_REPO=$PWD/SeqDance
    """

    def __init__(self, model_name: str, hf_path: str, device: str, max_seq_len: int):
        self.model_name = model_name
        self.hf_path = hf_path
        self.device = device
        self.max_seq_len = max_seq_len

        repo = os.environ.get("SEQDANCE_REPO", "./SeqDance")
        repo_path = Path(repo).expanduser().resolve()
        if not (repo_path / "model" / "model.py").exists():
            raise RuntimeError(
                f"SeqDance repo not found at {repo_path}. Clone it and set "
                f"SEQDANCE_REPO:\n  git clone https://github.com/ShenLab/SeqDance\n"
                f"  export SEQDANCE_REPO=$PWD/SeqDance"
            )
        # SeqDance import is fiddly. <repo>/model/ has NO __init__.py, so it is a
        # *namespace* package: `from model.model import ESMwrap` resolves only when
        # the repo ROOT (and not <repo>/model) is on sys.path -- putting
        # <repo>/model on the path makes model.py shadow the namespace package and
        # you get "'model' is not a package". But model/model.py also does a bare
        # `from config import config`, where config.py lives at <repo>/model/config.py.
        # Solution: keep only the repo root on the path, and pre-load config.py
        # directly into sys.modules['config'] (it's a standalone dict module, no
        # imports of its own) so the bare import resolves without touching the path.
        import importlib.util
        repo_str = str(repo_path)
        while repo_str in sys.path:
            sys.path.remove(repo_str)
        sys.path.insert(0, repo_str)
        # Drop any stale entries cached by an earlier failed import in this kernel.
        for stale in ("model", "model.model", "config"):
            sys.modules.pop(stale, None)
        cfg_path = repo_path / "model" / "config.py"
        spec = importlib.util.spec_from_file_location("config", cfg_path)
        cfg_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cfg_mod)
        sys.modules["config"] = cfg_mod   # satisfies `from config import config`
        from model.model import ESMwrap  # namespace package via repo root on path

        logger.info(f"Loading SeqDance wrapper {hf_path} (repo: {repo_path})...")
        wrapper = ESMwrap.from_pretrained(hf_path)
        # Inner transformer (standard EsmModel) carries the representations.
        self.model = wrapper.esm2.to(device).eval()
        self.tokenizer = AutoTokenizer.from_pretrained(SEQDANCE_BACKBONE)
        logger.info(
            f"SeqDance inner EsmModel ready: hidden={self.model.config.hidden_size}, "
            f"layers={self.model.config.num_hidden_layers}"
        )


def get_embedder(model_name: str, hf_path: str, device: str, max_seq_len: int) -> PLMEmbedder:
    if "esm2" in model_name:
        return ESM2Embedder(model_name, hf_path, device, max_seq_len)
    elif "seqdance" in model_name or "esmdance" in model_name:
        return SeqDanceEmbedder(model_name, hf_path, device, max_seq_len)
    else:
        # Fallback to base
        return PLMEmbedder(model_name, hf_path, device, max_seq_len)
