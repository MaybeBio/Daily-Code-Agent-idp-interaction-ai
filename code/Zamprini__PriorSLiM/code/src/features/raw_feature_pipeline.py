"""Build IDR manifests and derived sequence-feature resources."""
                      

from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from pssm_conservation_idppi import find_pssm_file, load_pssm_with_extras
from esm_region_pipeline import (
    ESM_EMBEDDING_DIM,
    ESM_MODEL_NAME,
    ESM_REPR_LAYER,
    REGION_KEY_COLUMNS,
    base_accession,
    build_esm_region_manifest,
)

IDR_KEY_COLUMNS = ["uniprot", "idr_start", "idr_end"]
NARDINI_GROUPS = ("STNQCH", "ILMV", "RK", "ED", "FWY", "A", "P", "G")
NARDINI_FEATURE_NAMES = tuple(
    f"nardini_{left}__{right}"
    for left_index, left in enumerate(NARDINI_GROUPS)
    for right_index, right in enumerate(NARDINI_GROUPS)
    if left_index <= right_index
)
GUSSIDR_FEATURE_NAMES = tuple(f"feat_{index}" for index in range(1, 39))

@dataclass(frozen=True)
class DerivedFeaturePaths:
    idr_manifest: Path
    esm_region_manifest: Path
    esm_region_metadata: Path
    nardini_features: Path
    gussidr_features: Path
    provenance: Path

def stable_frame_hash(frame: pd.DataFrame, columns: Sequence[str]) -> str:
    data = frame[list(columns)].to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(data).hexdigest()

def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def json_hash(value) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()

def build_idr_manifest(windows: pd.DataFrame, pssm_dir: Path) -> pd.DataFrame:
    required = {"uniprot", "idr_start", "idr_end"}
    missing = required - set(windows.columns)
    if missing:
        raise ValueError(f"Window table is missing IDR columns: {sorted(missing)}")
    keys = windows[IDR_KEY_COLUMNS].drop_duplicates().copy()
    if keys[["idr_start", "idr_end"]].isna().any().any():
        raise ValueError("Every window must have finite IDR bounds")
    keys["idr_start"] = keys["idr_start"].astype(int)
    keys["idr_end"] = keys["idr_end"].astype(int)
    protein_sequences = {}
    pssm_paths = {}
    base_uids = keys["uniprot"].astype(str).map(base_accession).unique()
    for uid in base_uids:
        pssm_path = find_pssm_file(pssm_dir, uid)
        if pssm_path is None:
            raise FileNotFoundError(f"No PSSM file for {uid}")
        _matrix, sequence, _extras = load_pssm_with_extras(pssm_path)
        if not sequence:
            raise ValueError(f"PSSM does not contain a sequence: {pssm_path}")
        protein_sequences[uid] = str(sequence).upper()
        pssm_paths[uid] = str(Path(pssm_path).resolve())
    records = []
    for row in keys.itertuples(index=False):
        uid, start, end = str(row.uniprot), int(row.idr_start), int(row.idr_end)
        base_uid = base_accession(uid)
        protein_sequence = protein_sequences[base_uid]
        if start < 1 or end < start or end > len(protein_sequence):
            raise ValueError(
                f"Invalid IDR {uid}:{start}-{end} for PSSM sequence length "
                f"{len(protein_sequence)}"
            )

                                                                             
        chars = list(protein_sequence)
        if uid != base_uid:
            variant_rows = windows[
                windows["uniprot"].astype(str).eq(uid)
                & pd.to_numeric(windows["idr_start"], errors="coerce").eq(start)
                & pd.to_numeric(windows["idr_end"], errors="coerce").eq(end)
            ]
            if "sequence" not in variant_rows.columns or variant_rows.empty:
                raise ValueError(f"Variant {uid}:{start}-{end} has no window sequence")
            assigned = {}
            for variant_row in variant_rows.itertuples(index=False):
                window_sequence = str(variant_row.sequence).upper()
                window_start = int(variant_row.window_start)
                for offset, amino_acid in enumerate(window_sequence):
                    position = window_start + offset
                    if amino_acid == "X" or not (1 <= position <= len(chars)):
                        continue
                    prior = assigned.get(position)
                    if prior is not None and prior != amino_acid:
                        raise ValueError(
                            f"Conflicting variant residues for {uid} at position {position}"
                        )
                    assigned[position] = amino_acid
                    chars[position - 1] = amino_acid
        sequence = "".join(chars[start - 1 : end])
        records.append(
            {
                "uniprot": uid,
                "idr_start": start,
                "idr_end": end,
                "sequence": sequence,
                "sequence_length": len(sequence),
                "protein_length": len(protein_sequence),
                "pssm_path": pssm_paths[base_uid],
            }
        )
    return pd.DataFrame(records).sort_values(IDR_KEY_COLUMNS).reset_index(drop=True)

def default_nardini_python() -> Path | None:
    override = os.environ.get("NARDINI_PYTHON")
    if override:
        candidate = Path(override).expanduser().resolve()
        if not candidate.is_file():
            raise FileNotFoundError(f"NARDINI_PYTHON does not exist: {candidate}")
        candidates = [candidate]
    else:
        current = Path(sys.executable).resolve()
        candidates = [current]
        if len(current.parents) > 2 and current.parents[2].name == "envs":
            candidates = [current.parents[2] / name / "bin/python"
                          for name in ("priorslim-nardini", "nardini")] + candidates
    for candidate in candidates:
        if not candidate.is_file():
            continue
        check = subprocess.run(
            [str(candidate), "-c", "from nardini.score_and_plot import calculate_zscore; import Bio"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if check.returncode == 0:
            return candidate
    return None

def default_esm_python() -> Path | None:
    override = os.environ.get("ESM_PYTHON")
    if override and Path(override).is_file():
        return Path(override).resolve()
    current = Path(sys.executable).resolve()
    candidates = [current]
    if len(current.parents) > 2 and current.parents[2].name == "envs":
        candidates.append(current.parents[2] / "IDPPI/bin/python")
    for candidate in candidates:
        if not candidate.is_file():
            continue
        check = subprocess.run(
            [
                str(candidate), "-c",
                "import esm, torch; assert hasattr(esm.pretrained, 'esm2_t33_650M_UR50D')",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            env={**os.environ, "KMP_DUPLICATE_LIB_OK": "TRUE"},
        )
        if check.returncode == 0:
            return candidate
    return None

def cache_is_current(metadata_path: Path, expected: dict) -> bool:
    if not metadata_path.is_file():
        return False
    try:
        observed = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return observed.get("cache_signature") == json_hash(expected)

def cache_settings_compatible(
    metadata_path: Path, expected: dict, ignored_keys: set[str],
) -> bool:
    pass
    if not metadata_path.is_file():
        return False
    try:
        observed = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return all(
        observed.get(key) == value
        for key, value in expected.items()
        if key not in ignored_keys
    )

def reusable_base_rows(
    idrs: pd.DataFrame, output: Path, feature_names: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pass
    key_columns = list(IDR_KEY_COLUMNS)
    expected_columns = [*key_columns, *feature_names]
    cached = pd.read_csv(output, sep="\t")
    if cached.columns.tolist() != expected_columns:
        return pd.DataFrame(columns=expected_columns), idrs.copy()
    cached = cached[
        cached["uniprot"].astype(str).map(base_accession)
        .eq(cached["uniprot"].astype(str))
    ].drop_duplicates(key_columns, keep="first")
    wanted = idrs[key_columns]
    reusable = wanted.merge(cached, on=key_columns, how="inner", validate="one_to_one")
    pending = idrs.merge(
        reusable[key_columns].assign(_reuse=1), on=key_columns, how="left"
    )
    pending = pending[pending["_reuse"].isna()].drop(columns="_reuse")
    return reusable, pending

def combine_cached_and_generated(
    idrs: pd.DataFrame, reusable: pd.DataFrame, generated: pd.DataFrame,
    feature_names: Sequence[str],
) -> pd.DataFrame:
    expected_columns = [*IDR_KEY_COLUMNS, *feature_names]
    combined = pd.concat([reusable, generated], ignore_index=True)
    if combined.duplicated(IDR_KEY_COLUMNS).any():
        raise RuntimeError("Incremental derived cache contains duplicate IDR keys")
    ordered = idrs[IDR_KEY_COLUMNS].merge(
        combined, on=IDR_KEY_COLUMNS, how="left", validate="one_to_one"
    )
    if ordered.columns.tolist() != expected_columns or ordered[list(feature_names)].isna().any().any():
        raise RuntimeError("Incremental derived cache is incomplete")
    return ordered

def run_worker(command: list[str], env: dict | None = None) -> None:
    subprocess.run(command, env=env, check=True)

def ensure_esm_region_cache(
    regions: pd.DataFrame,
    cache_dir: Path,
    worker_script: Path,
    python_path: Path,
    batch_size: int,
    device: str,
    force: bool,
) -> Path:
    region_manifest = cache_dir / "esm2_region_manifest.tsv"
    metadata = cache_dir / "esm2_region_metadata.tsv"
    cache_manifest = cache_dir / "esm2_region_cache.manifest.json"
    embedding_dir = cache_dir / "esm2_regions"
    regions.to_csv(region_manifest, sep="\t", index=False)
    settings = {
        "schema": 1,
        "generator": str(worker_script.resolve()),
        "generator_sha256": file_hash(worker_script),
        "region_sequence_sha256": stable_frame_hash(
            regions,
            [
                *REGION_KEY_COLUMNS, "crop_start", "crop_end", "sequence",
                "embedding_file",
            ],
        ),
        "esm_model": ESM_MODEL_NAME,
        "esm_repr_layer": ESM_REPR_LAYER,
        "embedding_dim": ESM_EMBEDDING_DIM,
    }
    if not force and metadata.is_file() and cache_is_current(cache_manifest, settings):
        try:
            observed = pd.read_csv(metadata, sep="\t")
            expected_keys = set(
                regions[REGION_KEY_COLUMNS].itertuples(index=False, name=None)
            )
            observed_keys = set(
                observed[REGION_KEY_COLUMNS].itertuples(index=False, name=None)
            )
            complete = len(observed) == len(regions) and observed_keys == expected_keys
            if complete:
                for row in observed.itertuples(index=False):
                    array = np.load(Path(row.embedding_path), mmap_mode="r")
                    if array.shape != (int(row.crop_length), ESM_EMBEDDING_DIM):
                        complete = False
                        break
            if complete:
                return metadata
        except (KeyError, OSError, ValueError):

                                                                   
            pass
    command = [
        str(python_path), str(worker_script), "--region-manifest", str(region_manifest),
        "--embedding-dir", str(embedding_dir), "--metadata-output", str(metadata),
        "--batch-size", str(batch_size), "--device", device,
    ]
    if force:
        command.append("--force")
    environment = os.environ.copy()
    environment.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    environment.setdefault("OMP_NUM_THREADS", "1")
    run_worker(command, env=environment)
    generated = pd.read_csv(metadata, sep="\t")
    if len(generated) != len(regions):
        raise RuntimeError("ESM worker returned incomplete region metadata")
    cache_manifest.write_text(
        json.dumps({**settings, "cache_signature": json_hash(settings)}, indent=2) + "\n",
        encoding="utf-8",
    )
    return metadata

def ensure_nardini_cache(
    idrs: pd.DataFrame,
    cache_dir: Path,
    worker_script: Path,
    python_path: Path,
    num_scrambles: int,
    random_seed: int,
    workers: int,
    force: bool,
) -> Path:
    output = cache_dir / "nardini_features.tsv"
    metadata = cache_dir / "nardini_features.manifest.json"
    settings = {
        "schema": 1,
        "generator": str(worker_script.resolve()),
        "generator_sha256": file_hash(worker_script),
        "idr_sequence_sha256": stable_frame_hash(
            idrs, [*IDR_KEY_COLUMNS, "sequence"]
        ),
        "num_scrambles": int(num_scrambles),
        "random_seed": int(random_seed),
        "feature_names": list(NARDINI_FEATURE_NAMES),
    }
    if not force and output.is_file() and cache_is_current(metadata, settings):
        return output
    reusable = pd.DataFrame(columns=[*IDR_KEY_COLUMNS, *NARDINI_FEATURE_NAMES])
    pending_idrs = idrs
    if (
        not force and output.is_file()
        and cache_settings_compatible(metadata, settings, {"idr_sequence_sha256"})
    ):
        reusable, pending_idrs = reusable_base_rows(
            idrs, output, NARDINI_FEATURE_NAMES
        )
    pending = cache_dir / "_pending_nardini_idrs.tsv"
    temporary = cache_dir / "_pending_nardini_features.tsv"
    pending_idrs[[*IDR_KEY_COLUMNS, "sequence"]].to_csv(pending, sep="\t", index=False)
    runtime_cache = cache_dir / "runtime_cache"
    runtime_cache.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment.setdefault("MPLCONFIGDIR", str(runtime_cache / "matplotlib"))
    environment.setdefault("XDG_CACHE_HOME", str(runtime_cache))
    if len(pending_idrs):
        run_worker(
            [
                str(python_path), str(worker_script), "--idr-manifest", str(pending),
                "--output", str(temporary), "--num-scrambles", str(num_scrambles),
                "--random-seed", str(random_seed), "--workers", str(workers),
            ],
            env=environment,
        )
        generated = pd.read_csv(temporary, sep="\t")
    else:
        generated = pd.DataFrame(columns=[*IDR_KEY_COLUMNS, *NARDINI_FEATURE_NAMES])
    expected_columns = [*IDR_KEY_COLUMNS, *NARDINI_FEATURE_NAMES]
    if generated.columns.tolist() != expected_columns or len(generated) != len(pending_idrs):
        raise RuntimeError("NARDINI worker returned an incompatible cache")
    combined = combine_cached_and_generated(
        idrs, reusable, generated, NARDINI_FEATURE_NAMES
    )
    combined.to_csv(output, sep="\t", index=False)
    temporary.unlink(missing_ok=True)
    pending.unlink(missing_ok=True)
    metadata.write_text(
        json.dumps({**settings, "cache_signature": json_hash(settings)}, indent=2) + "\n",
        encoding="utf-8",
    )
    return output

def region_embedding_signature(metadata: pd.DataFrame):
    records = []
    for row in metadata.itertuples(index=False):
        path = Path(row.embedding_path)
        if not path.is_file():
            raise FileNotFoundError(f"Missing ESM region embedding: {path}")
        stat = path.stat()
        records.append(
            (
                str(row.uniprot), str(row.region_type), int(row.original_start),
                int(row.original_end), int(row.crop_start), int(row.crop_end),
                stat.st_size, stat.st_mtime_ns,
            )
        )
    return json_hash(records)

def ensure_gussidr_cache(
    idrs: pd.DataFrame,
    cache_dir: Path,
    worker_script: Path,
    gussidr_dir: Path,
    esm_region_metadata: Path,
    python_path: Path,
    batch_size: int,
    device: str,
    force: bool,
) -> Path:
    output = cache_dir / "gussidr_ensemble_features.tsv"
    metadata = cache_dir / "gussidr_ensemble_features.manifest.json"
    model_assets = []
    for name in ("As", "ree", "rg", "v", "spr"):
        for filename in (f"model_best_{name}.pth", f"best_hyperparameters_{name}.json"):
            path = gussidr_dir / name / filename
            if not path.is_file():
                raise FileNotFoundError(f"Missing GussIDR model asset: {path}")
            model_assets.append((str(path.resolve()), file_hash(path)))
    esm_metadata = pd.read_csv(esm_region_metadata, sep="\t")
    idr_esm_metadata = esm_metadata.loc[esm_metadata["region_type"] == "idr"].copy()
    settings = {
        "schema": 1,
        "generator": str(worker_script.resolve()),
        "generator_sha256": file_hash(worker_script),
        "idr_sequence_sha256": stable_frame_hash(
            idrs, [*IDR_KEY_COLUMNS, "sequence"]
        ),
        "esm_source": "biological_idr_crop",
        "esm_source_sha256": region_embedding_signature(idr_esm_metadata),
        "model_assets": model_assets,
        "feature_names": list(GUSSIDR_FEATURE_NAMES),
    }
    if not force and output.is_file() and cache_is_current(metadata, settings):
        return output
    reusable = pd.DataFrame(columns=[*IDR_KEY_COLUMNS, *GUSSIDR_FEATURE_NAMES])
    pending_idrs = idrs
    if (
        not force and output.is_file()
        and cache_settings_compatible(
            metadata, settings, {"idr_sequence_sha256", "esm_source_sha256"}
        )
    ):
        reusable, pending_idrs = reusable_base_rows(
            idrs, output, GUSSIDR_FEATURE_NAMES
        )
    pending = cache_dir / "_pending_gussidr_idrs.tsv"
    temporary = cache_dir / "_pending_gussidr_features.tsv"
    pending_idrs[[*IDR_KEY_COLUMNS, "sequence"]].to_csv(pending, sep="\t", index=False)
    command = [
        str(python_path), str(worker_script), "--idr-manifest", str(pending),
        "--esm-region-metadata", str(esm_region_metadata),
        "--gussidr-dir", str(gussidr_dir),
        "--output", str(temporary), "--batch-size", str(batch_size),
        "--device", device,
    ]
    environment = os.environ.copy()

                                                                               
    environment.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
    environment.setdefault("OMP_NUM_THREADS", "1")
    runtime_cache = cache_dir / "runtime_cache"
    runtime_cache.mkdir(exist_ok=True)
    environment.setdefault("MPLCONFIGDIR", str(runtime_cache / "matplotlib"))
    environment.setdefault("XDG_CACHE_HOME", str(runtime_cache))
    if len(pending_idrs):
        run_worker(command, env=environment)
        generated = pd.read_csv(temporary, sep="\t")
    else:
        generated = pd.DataFrame(columns=[*IDR_KEY_COLUMNS, *GUSSIDR_FEATURE_NAMES])
    expected_columns = [*IDR_KEY_COLUMNS, *GUSSIDR_FEATURE_NAMES]
    if generated.columns.tolist() != expected_columns or len(generated) != len(pending_idrs):
        raise RuntimeError("GussIDR worker returned an incompatible cache")
    combined = combine_cached_and_generated(
        idrs, reusable, generated, GUSSIDR_FEATURE_NAMES
    )
    combined.to_csv(output, sep="\t", index=False)
    temporary.unlink(missing_ok=True)
    pending.unlink(missing_ok=True)
    metadata.write_text(
        json.dumps({**settings, "cache_signature": json_hash(settings)}, indent=2) + "\n",
        encoding="utf-8",
    )
    return output

def prepare_derived_features(
    windows: pd.DataFrame,
    pssm_dir: Path,
    fragment_file: Path,
    cache_dir: Path,
    gussidr_dir: Path,
    esm_python: Path,
    nardini_python: Path,
    gussidr_python: Path,
    esm_batch_size: int = 1,
    esm_device: str = "auto",
    nardini_num_scrambles: int = 1000,
    nardini_random_seed: int = 1,
    nardini_workers: int = 1,
    gussidr_batch_size: int = 4,
    gussidr_device: str = "auto",
    force: bool = False,
    idr_manifest_frame: pd.DataFrame | None = None,
) -> DerivedFeaturePaths:
    cache_dir.mkdir(parents=True, exist_ok=True)
    idrs = (
        idr_manifest_frame.copy()
        if idr_manifest_frame is not None
        else build_idr_manifest(windows, pssm_dir)
    )
    idr_manifest = cache_dir / "idr_manifest.tsv"
    idrs.to_csv(idr_manifest, sep="\t", index=False)
    feature_dir = Path(__file__).resolve().parent
    esm_regions = build_esm_region_manifest(idrs, fragment_file)
    esm_region_metadata = ensure_esm_region_cache(
        esm_regions, cache_dir, feature_dir / "generate_esm2_region_cache.py",
        esm_python, esm_batch_size, esm_device, force,
    )
    esm_region_manifest = cache_dir / "esm2_region_manifest.tsv"
    nardini_features = ensure_nardini_cache(
        idrs, cache_dir, feature_dir / "generate_nardini_cache.py", nardini_python,
        nardini_num_scrambles, nardini_random_seed, nardini_workers, force,
    )
    gussidr_features = ensure_gussidr_cache(
        idrs, cache_dir, feature_dir / "generate_gussidr_cache.py", gussidr_dir,
        esm_region_metadata, gussidr_python, gussidr_batch_size,
        gussidr_device, force,
    )
    provenance = cache_dir / "pipeline_manifest.json"
    payload = {
        "schema": 1,
        "primary_inputs": {
            "windows": "train/val/test or --windows",
            "pssm_dir": str(pssm_dir.resolve()),
            "fragment_file": str(fragment_file.resolve()),
            "esm_model": ESM_MODEL_NAME,
            "esm_repr_layer": ESM_REPR_LAYER,
        },
        "derived_caches": {
            "idr_manifest": str(idr_manifest.resolve()),
            "esm_region_manifest": str(esm_region_manifest.resolve()),
            "esm_region_metadata": str(esm_region_metadata.resolve()),
            "nardini_features": str(nardini_features.resolve()),
            "gussidr_features": str(gussidr_features.resolve()),
        },
        "idr_count": len(idrs),
        "idr_sequence_sha256": stable_frame_hash(idrs, [*IDR_KEY_COLUMNS, "sequence"]),
    }
    provenance.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return DerivedFeaturePaths(
        idr_manifest, esm_region_manifest, esm_region_metadata,
        nardini_features, gussidr_features, provenance,
    )
