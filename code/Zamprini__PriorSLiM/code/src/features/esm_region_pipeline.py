"""Construct sequence and coordinate manifests for regional ESM2 encoding."""
                      

from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import hashlib
import re
from pathlib import Path

import pandas as pd

from chemical_window_idppi_fast import load_regions
from pssm_conservation_idppi import load_pssm_with_extras

ESM_MODEL_NAME = "esm2_t33_650M_UR50D"
ESM_REPR_LAYER = 33
ESM_EMBEDDING_DIM = 1280
ESM_MAX_RESIDUES = 1022
REGION_TYPES = ("idr", "left_domain", "right_domain")
VARIANT_RE = re.compile(r"^(?P<base>.+)_\d+_\d+_[VnMFC]\d+$")
REGION_KEY_COLUMNS = ["uniprot", "region_type", "original_start", "original_end"]

def base_accession(uniprot: str) -> str:
    match = VARIANT_RE.match(str(uniprot))
    return match.group("base") if match else str(uniprot)

def centered_crop(start: int, end: int, max_len: int = ESM_MAX_RESIDUES):
    pass
    length = end - start + 1
    if length <= 0:
        raise ValueError(f"Invalid inclusive interval: {start}-{end}")
    if length <= max_len:
        return start, end
    left_trim = (length - max_len) // 2
    crop_start = start + left_trim
    return crop_start, crop_start + max_len - 1

def crop_adjacent_domain(
    domain_start: int,
    domain_end: int,
    side: str,
    max_len: int = ESM_MAX_RESIDUES,
):
    length = domain_end - domain_start + 1
    if length <= 0:
        raise ValueError(f"Invalid inclusive interval: {domain_start}-{domain_end}")
    if length <= max_len:
        return domain_start, domain_end
    if side == "left":
        return domain_end - max_len + 1, domain_end
    if side == "right":
        return domain_start, domain_start + max_len - 1
    raise ValueError("side must be 'left' or 'right'")

def nearest_adjacent_domains(ordered, idr_start: int, idr_end: int):
    left = None
    right = None
    for domain in ordered:
        start, end = int(domain.start), int(domain.end)
        if end < idr_start and (left is None or end > left[1]):
            left = (start, end)
        if start > idr_end and (right is None or start < right[0]):
            right = (start, end)
    return left, right

def embedding_filename(row: dict) -> str:
    identity = "|".join(
        str(row[name])
        for name in (
            "uniprot", "region_type", "original_start", "original_end",
            "crop_start", "crop_end",
        )
    )

                                                                           

                                                        
    if (
        str(row["region_type"]) == "idr"
        and base_accession(str(row["uniprot"])) != str(row["uniprot"])
    ):
        sequence_digest = hashlib.sha256(
            str(row["sequence"]).encode("utf-8")
        ).hexdigest()
        identity += f"|sequence_sha256={sequence_digest}"
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16]
    safe_uid = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(row["uniprot"]))
    return f"{safe_uid}__{row['region_type']}__{digest}.npy"

def region_record(
    uid: str,
    region_type: str,
    original_start: int,
    original_end: int,
    crop_start: int,
    crop_end: int,
    protein_sequence: str,
):
    original_length = original_end - original_start + 1
    crop_length = crop_end - crop_start + 1
    if crop_length > ESM_MAX_RESIDUES:
        raise ValueError(f"ESM crop exceeds {ESM_MAX_RESIDUES}: {uid}:{crop_start}-{crop_end}")
    if crop_start < original_start or crop_end > original_end:
        raise ValueError("Crop must be contained inside its biological region")
    if crop_start < 1 or crop_end > len(protein_sequence):
        raise ValueError(
            f"Crop {uid}:{crop_start}-{crop_end} exceeds protein length {len(protein_sequence)}"
        )
    sequence = protein_sequence[crop_start - 1 : crop_end]
    if len(sequence) != crop_length:
        raise RuntimeError("1-based inclusive to Python slice conversion failed")
    row = {
        "uniprot": uid,
        "region_type": region_type,
        "original_start": original_start,
        "original_end": original_end,
        "crop_start": crop_start,
        "crop_end": crop_end,
        "original_length": original_length,
        "crop_length": crop_length,
        "was_truncated": int(original_length > ESM_MAX_RESIDUES),
        "was_sequence_clipped": int(
            crop_start > original_start or crop_end < original_end
        ),
        "sequence": sequence,
    }
    row["embedding_file"] = embedding_filename(row)
    return row

def build_esm_region_manifest(idrs: pd.DataFrame, fragment_file: Path) -> pd.DataFrame:
    pass
    required = {
        "uniprot", "idr_start", "idr_end", "protein_length", "pssm_path"
    }
    missing = required - set(idrs.columns)
    if missing:
        raise ValueError(f"IDR manifest is missing columns: {sorted(missing)}")
    _idrs_from_regions, ordered_by_uid = load_regions(str(fragment_file))
    protein_sequences = {}
    for row in idrs.drop_duplicates("pssm_path").itertuples(index=False):
        _matrix, sequence, _extras = load_pssm_with_extras(Path(row.pssm_path))
        sequence = str(sequence).upper()
        if len(sequence) != int(row.protein_length):
            raise ValueError(f"PSSM sequence length changed for {row.uniprot}")
        protein_sequences[base_accession(str(row.uniprot))] = sequence

    records = []
    for row in idrs.itertuples(index=False):
        uid = str(row.uniprot)
        base_uid = base_accession(uid)
        idr_start, idr_end = int(row.idr_start), int(row.idr_end)
        base_sequence = protein_sequences[base_uid]
        chars = list(base_sequence)
        idr_sequence = str(row.sequence).upper()
        if len(idr_sequence) != idr_end - idr_start + 1:
            raise ValueError(f"IDR sequence length mismatch for {uid}:{idr_start}-{idr_end}")
        chars[idr_start - 1 : idr_end] = idr_sequence
        sequence = "".join(chars)
        crop_start, crop_end = centered_crop(idr_start, idr_end)
        records.append(
            region_record(
                uid, "idr", idr_start, idr_end, crop_start, crop_end, sequence
            )
        )
        ordered = ordered_by_uid.get(base_uid, [])
        left, right = nearest_adjacent_domains(ordered, idr_start, idr_end)
        if left is not None:
            available_left = (max(1, left[0]), min(len(sequence), left[1]))
            if available_left[0] <= available_left[1]:
                crop_start, crop_end = crop_adjacent_domain(
                    *available_left, side="left"
                )
                records.append(
                    region_record(
                        base_uid, "left_domain", *left, crop_start, crop_end, base_sequence
                    )
                )
        if right is not None:
            available_right = (max(1, right[0]), min(len(sequence), right[1]))
            if available_right[0] <= available_right[1]:
                crop_start, crop_end = crop_adjacent_domain(
                    *available_right, side="right"
                )
                records.append(
                    region_record(
                        base_uid, "right_domain", *right, crop_start, crop_end, base_sequence
                    )
                )
    frame = pd.DataFrame(records)
    frame = frame.drop_duplicates(REGION_KEY_COLUMNS, keep="first")
    frame = frame.sort_values(REGION_KEY_COLUMNS).reset_index(drop=True)
    if frame.duplicated(REGION_KEY_COLUMNS).any():
        raise RuntimeError("ESM region manifest contains duplicate biological units")
    if (frame["crop_length"] > ESM_MAX_RESIDUES).any():
        raise RuntimeError("ESM region crop length invariant was violated")
    return frame
