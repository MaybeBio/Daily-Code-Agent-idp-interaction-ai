"""Generate GussIDR ensemble predictions from IDR sequences."""
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
                      

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.nn.utils.rnn import pad_sequence

def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def load_state(path: Path, device):
    try:
        return torch.load(path, map_location=device, weights_only=True)
    except TypeError:
        return torch.load(path, map_location=device)

def load_models(root: Path, device):
    modules = {
        name: load_module(root / name / f"predict_{name}.py", f"gussidr_{name}")
        for name in ("As", "ree", "rg", "v", "spr")
    }
    models = {}
    for name, module in modules.items():
        hparams_path = root / name / f"best_hyperparameters_{name}.json"
        weight_path = root / name / f"model_best_{name}.pth"
        hparams = json.loads(hparams_path.read_text(encoding="utf-8"))
        feedforward = hparams.get("d_model", 256) * hparams.get(
            "dim_feedforward_factor", 2
        )
        common = dict(
            input_dim=hparams.get("input_dim", 1280),
            emb_dim=hparams.get("emb_dim", 128),
            d_model=hparams.get("d_model", 256),
            nhead=hparams.get("nhead", 4),
            transformer_layers=hparams.get("transformer_layers", 2),
            dim_feedforward=feedforward,
            mlp_hidden_dim=hparams.get("mlp_hidden_dim", 128),
            dropout_rate=hparams.get("dropout_rate", 0.15),
        )
        if name in {"As", "ree", "rg"}:
            model = module.TransformerMDNModel(
                **common, num_mixtures=hparams.get("num_mixtures", 4)
            )
        else:
            model = module.TransformerMoEModel(
                **common, num_experts=hparams.get("num_experts", 8)
            )
        model.load_state_dict(load_state(weight_path, device))
        model.to(device).eval()
        models[name] = model
    return modules, models

def load_idr_embedding(row, metadata_by_key):
    uid = str(row.uniprot)
    start, end = int(row.idr_start), int(row.idr_end)
    key = (uid, "idr", start, end)
    if key not in metadata_by_key:
        raise KeyError(f"No IDR ESM metadata for {uid}:{start}-{end}")
    metadata = metadata_by_key[key]
    embedding = np.array(np.load(metadata["embedding_path"]), dtype=np.float32, copy=True)
    expected = (int(metadata["crop_length"]), 1280)
    if embedding.shape != expected:
        raise ValueError(f"Expected IDR ESM shape {expected}, got {embedding.shape}")
    full_idr_sequence = str(row.sequence)
    offset = int(metadata["crop_start"]) - start
    crop_sequence = full_idr_sequence[offset : offset + int(metadata["crop_length"])]
    if len(crop_sequence) != embedding.shape[0]:
        raise ValueError(f"IDR crop sequence/ESM mismatch for {uid}:{start}-{end}")
    return embedding, crop_sequence

def physical_features(sequence: str, aa_phy):
    length = len(sequence)
    result = np.zeros((length, 3), dtype=np.float32)
    for index, aa in enumerate(sequence):
        charge, hydropathy = aa_phy.get(aa, (0.0, 0.0))
        result[index] = [charge, hydropathy / 5.0, index / max(length - 1, 1)]
    return result

def make_batch(items, aa_phy, device):
    embeddings = [torch.from_numpy(item["embedding"]) for item in items]
    physical = [torch.from_numpy(physical_features(item["sequence"], aa_phy)) for item in items]
    lengths = torch.tensor([len(item["sequence"]) for item in items], dtype=torch.float32)
    padded_embedding = pad_sequence(embeddings, batch_first=True, padding_value=0.0)
    padded_physical = pad_sequence(physical, batch_first=True, padding_value=0.0)
    mask = torch.arange(padded_embedding.shape[1])[None, :] >= lengths[:, None]
    return (
        padded_embedding.to(device),
        padded_physical.to(device),
        lengths.to(device),
        mask.to(device),
    )

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--idr-manifest", required=True, type=Path)
    parser.add_argument("--esm-region-metadata", required=True, type=Path)
    parser.add_argument("--gussidr-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--device", choices=("auto", "cpu", "mps", "cuda"), default="auto")
    parser.add_argument("--v-mean", type=float, default=0.5331)
    parser.add_argument("--v-std", type=float, default=0.0377)
    parser.add_argument("--spr-mean", type=float, default=10.1507)
    parser.add_argument("--spr-std", type=float, default=0.1471)
    return parser.parse_args()

def main():
    args = parse_args()
    if args.device == "auto":
        device = torch.device(
            "cuda" if torch.cuda.is_available() else "mps"
            if torch.backends.mps.is_available() else "cpu"
        )
    else:
        device = torch.device(args.device)
    modules, models = load_models(args.gussidr_dir, device)
    frame = pd.read_csv(args.idr_manifest, sep="\t")
    esm_metadata = pd.read_csv(args.esm_region_metadata, sep="\t")
    metadata_by_key = {
        (
            str(row.uniprot), str(row.region_type), int(row.original_start),
            int(row.original_end),
        ): row._asdict()
        for row in esm_metadata.itertuples(index=False)
    }
    items = []
    for original_index, row in enumerate(frame.itertuples(index=False)):
        embedding, crop_sequence = load_idr_embedding(row, metadata_by_key)
        items.append(
            {
                "original_index": original_index,
                "uniprot": str(row.uniprot),
                "idr_start": int(row.idr_start),
                "idr_end": int(row.idr_end),
                "sequence": crop_sequence,
                "embedding": embedding,
            }
        )
    items.sort(key=lambda item: len(item["sequence"]))
    rows = []
    aa_phy = modules["As"].AA_PHY
    with torch.no_grad():
        for offset in range(0, len(items), args.batch_size):
            batch_items = items[offset : offset + args.batch_size]
            embedding, physical, lengths, mask = make_batch(batch_items, aa_phy, device)
            _, as_params, _ = models["As"](embedding, mask, physical, lengths)
            _, ree_params, _ = models["ree"](embedding, mask, physical, lengths)
            _, rg_params, _ = models["rg"](embedding, mask, physical, lengths)
            v = models["v"](embedding, mask, physical, lengths).squeeze(-1)
            spr = models["spr"](embedding, mask, physical, lengths).squeeze(-1)
            v = v * args.v_std + args.v_mean
            spr = spr * args.spr_std + args.spr_mean
            parameter_batches = [
                [part.detach().cpu().numpy() for part in params]
                for params in (as_params, ree_params, rg_params)
            ]
            v = v.detach().cpu().numpy()
            spr = spr.detach().cpu().numpy()
            for index, item in enumerate(batch_items):
                vector = []
                for params in parameter_batches:
                    vector.extend(np.concatenate([part[index] for part in params]).tolist())
                vector.extend([float(v[index]), float(spr[index])])
                if len(vector) != 38 or not np.isfinite(vector).all():
                    raise ValueError(f"Invalid GussIDR vector for {item['uniprot']}: {len(vector)}")
                rows.append(
                    [
                        item["original_index"], item["uniprot"], item["idr_start"],
                        item["idr_end"], *vector,
                    ]
                )
    output = pd.DataFrame(
        rows,
        columns=["original_index", "uniprot", "idr_start", "idr_end"]
        + [f"feat_{index}" for index in range(1, 39)],
    ).sort_values("original_index").drop(columns="original_index")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(
        f"[OK] GussIDR cache rows={len(output)} device={device} "
        f"embedding_source=biological_idr_crop output={args.output}"
    )

if __name__ == "__main__":
    main()
