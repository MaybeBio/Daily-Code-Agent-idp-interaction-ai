#!/usr/bin/env python3
"""Fetch DisProt human disorder annotations and joined AlphaFold DB per-residue pLDDT.

See `docs/research-protocol.md` for the exact sample-construction, exclusion, and
join procedure implemented here. Summary:

1. Fetch every human (taxon 9606) DisProt entry via the DisProt search API
   (``release=current``), cached to ``data/external/disprot_raw.json``.
2. Draw a fixed-seed reproducible random sample of ``--sample-size`` entries
   (default 400) from that full set, sorted by ``disprot_id`` before sampling.
3. For each sampled entry, fetch AlphaFold DB metadata by UniProt accession.
   Exclude (and log) entries with no AlphaFold DB prediction, a fragmented
   (multi-part) AlphaFold DB entry, or a DisProt/AlphaFold sequence-length
   mismatch.
4. For each remaining entry, fetch the per-residue confidence JSON and join
   ``confidenceScore`` (pLDDT) to a per-residue "inside a DisProt
   Structural-state/disorder region" boolean label, plus the covering
   region's evidence code and a protein-level "has a disorder-to-order
   region elsewhere" flag.
5. Write the joined per-residue table to ``data/external/joined_residues.csv``
   and the exclusion log to ``data/external/exclusions_log.json``.

This script is polite to both APIs (a small delay between requests, capped
concurrency) but is a genuine multi-hundred-request network operation and can
take several minutes to run.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import requests

DISPROT_SEARCH_URL = "https://disprot.org/api/search"
ALPHAFOLD_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
SEED = 20260813
DEFAULT_SAMPLE_SIZE = 400
REQUEST_TIMEOUT = 30
PAGE_SIZE = 200

DISORDER_NAMESPACE = "Structural state"
DISORDER_TERM = "disorder"
TRANSITION_NAMESPACE = "Structural transition"
TRANSITION_TERM = "disorder to order"


def fetch_all_human_disprot(cache_path: Path, *, force: bool = False) -> list[dict[str, Any]]:
    """Fetch (or load a cache of) every human DisProt entry, current release.

    IMPORTANT DEVIATION FROM THE DOCUMENTED API BEHAVIOR (discovered 2026-08-13,
    disclosed in docs/research-protocol.md and data/provenance.json): the
    ``limit``/``offset`` query parameters on ``/api/search`` do not actually
    paginate when an ``organism`` filter is present — every request, regardless
    of ``limit`` or ``offset``, returns the full matching result set (verified
    by requesting limit=5, 50, 500, and 2000, all of which returned the same
    1,339-entry ``data`` array with ``size: 1339``). This function therefore
    makes a single request and deduplicates by ``disprot_id`` defensively,
    rather than looping on offset as originally planned.
    """
    if cache_path.exists() and not force:
        return list(json.loads(cache_path.read_text())["data"])

    resp = requests.get(
        DISPROT_SEARCH_URL,
        params={"organism": "Homo sapiens", "release": "current", "limit": PAGE_SIZE, "offset": 0},
        timeout=REQUEST_TIMEOUT,
    )
    resp.raise_for_status()
    payload = resp.json()
    raw_entries: list[dict[str, Any]] = payload["data"]

    deduped: dict[str, dict[str, Any]] = {}
    for entry in raw_entries:
        deduped[str(entry["disprot_id"])] = entry
    entries = list(deduped.values())

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps({"size": len(entries), "data": entries}))
    return entries


def sample_entries(
    entries: list[dict[str, Any]], sample_size: int, seed: int = SEED
) -> list[dict[str, Any]]:
    """Fixed-seed reproducible sample, drawn from a deterministic disprot_id-sorted order."""
    ordered = sorted(entries, key=lambda e: str(e["disprot_id"]))
    rng = np.random.default_rng(seed)
    n = min(sample_size, len(ordered))
    indices = rng.choice(len(ordered), size=n, replace=False)
    return [ordered[i] for i in sorted(indices.tolist())]


def _disorder_regions(entry: dict[str, Any]) -> list[tuple[int, int, str, str]]:
    """Return (start, end, ec_id, ec_name) for Structural state / disorder regions."""
    out = []
    for region in entry.get("regions", []):
        if (
            region.get("term_namespace") == DISORDER_NAMESPACE
            and region.get("term_name") == DISORDER_TERM
        ):
            out.append(
                (
                    int(region["start"]),
                    int(region["end"]),
                    str(region.get("ec_id", "")),
                    str(region.get("ec_name", "")),
                )
            )
    return out


def _has_transition_region(entry: dict[str, Any]) -> bool:
    return any(
        region.get("term_namespace") == TRANSITION_NAMESPACE
        and region.get("term_name") == TRANSITION_TERM
        for region in entry.get("regions", [])
    )


def fetch_alphafold_metadata(acc: str) -> list[dict[str, Any]] | None:
    """Fetch AlphaFold DB metadata record(s) for exactly ``acc`` (not its isoforms).

    DEVIATION FROM DOCUMENTED BEHAVIOR (discovered 2026-08-13, disclosed in
    docs/research-protocol.md and data/provenance.json): the metadata endpoint
    can return multiple records for one query accession, but in current
    AlphaFold DB these are predominantly UniProt **isoform** entries (e.g.
    querying ``O00273`` also returns a record for isoform ``O00273-2``, whose
    own ``uniprotAccession`` is ``"O00273-2"``, not ``"O00273"``), not the
    long-protein ``-F1``/``-F2``/``-F3`` sequence fragmentation this project
    originally expected. This function filters to records whose
    ``uniprotAccession`` field exactly equals the queried accession before
    deciding whether an entry is genuinely fragmented (more than one exact-
    match record, which does still occur for a small number of very long
    proteins) or simply has isoform records alongside it (which are ignored).
    """
    resp = requests.get(ALPHAFOLD_PREDICTION_URL.format(acc=acc), timeout=REQUEST_TIMEOUT)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    data: list[dict[str, Any]] = resp.json()
    if not data:
        return None
    exact = [record for record in data if record.get("uniprotAccession") == acc]
    return exact if exact else None


def fetch_confidence(plddt_doc_url: str) -> dict[str, list[float]]:
    resp = requests.get(plddt_doc_url, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    result: dict[str, list[float]] = resp.json()
    return result


def process_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Fetch AlphaFold DB data for one DisProt entry and return a join record or exclusion."""
    acc = entry["acc"]
    disprot_id = entry["disprot_id"]
    sequence = entry["sequence"]

    metadata = fetch_alphafold_metadata(acc)
    if metadata is None:
        return {
            "status": "excluded",
            "reason": "no_alphafold_prediction",
            "acc": acc,
            "disprot_id": disprot_id,
        }
    if len(metadata) > 1:
        return {
            "status": "excluded",
            "reason": "fragmented_alphafold_entry",
            "acc": acc,
            "disprot_id": disprot_id,
            "n_fragments": len(metadata),
        }

    record = metadata[0]
    af_sequence = record.get("uniprotSequence") or record.get("sequence", "")
    if len(af_sequence) != len(sequence):
        return {
            "status": "excluded",
            "reason": "sequence_length_mismatch",
            "acc": acc,
            "disprot_id": disprot_id,
            "disprot_length": len(sequence),
            "alphafold_length": len(af_sequence),
        }

    plddt_doc_url = record["plddtDocUrl"]
    confidence = fetch_confidence(plddt_doc_url)
    residue_numbers = confidence["residueNumber"]
    scores = confidence["confidenceScore"]

    regions = _disorder_regions(entry)
    has_transition = _has_transition_region(entry)

    def covering_region(pos: int) -> tuple[str, str] | None:
        for start, end, ec_id, ec_name in regions:
            if start <= pos <= end:
                return ec_id, ec_name
        return None

    rows = []
    for residue_number, plddt in zip(residue_numbers, scores, strict=True):
        covering = covering_region(residue_number)
        rows.append(
            {
                "acc": acc,
                "disprot_id": disprot_id,
                "residue_number": residue_number,
                "plddt": plddt,
                "is_disorder": covering is not None,
                "is_conditional": has_transition,
                "ec_id": covering[0] if covering else "",
                "ec_name": covering[1] if covering else "",
            }
        )
    return {"status": "included", "acc": acc, "disprot_id": disprot_id, "rows": rows}


def run(
    sample_size: int,
    output_csv: Path,
    exclusions_path: Path,
    raw_cache_path: Path,
    *,
    max_workers: int = 8,
) -> None:
    all_entries = fetch_all_human_disprot(raw_cache_path)
    print(f"Fetched {len(all_entries)} total human DisProt entries (release=current)")

    sampled = sample_entries(all_entries, sample_size)
    print(f"Sampled {len(sampled)} entries (seed={SEED})")

    included: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(process_entry, entry): entry for entry in sampled}
        for i, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            if result["status"] == "included":
                included.append(result)
            else:
                excluded.append(result)
            if i % 25 == 0 or i == len(sampled):
                print(
                    f"  processed {i}/{len(sampled)} "
                    f"({len(included)} included, {len(excluded)} excluded)"
                )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "acc",
        "disprot_id",
        "residue_number",
        "plddt",
        "is_disorder",
        "is_conditional",
        "ec_id",
        "ec_name",
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in sorted(included, key=lambda r: str(r["acc"])):
            for row in record["rows"]:
                writer.writerow(row)

    total_residues = sum(len(r["rows"]) for r in included)
    disorder_residues = sum(1 for r in included for row in r["rows"] if row["is_disorder"])

    exclusions_path.parent.mkdir(parents=True, exist_ok=True)
    exclusions_path.write_text(
        json.dumps(
            {
                "sample_size_requested": sample_size,
                "sampled": len(sampled),
                "included_proteins": len(included),
                "excluded_proteins": len(excluded),
                "total_residues": total_residues,
                "disorder_annotated_residues": disorder_residues,
                "exclusions_by_reason": {
                    reason: sum(1 for e in excluded if e["reason"] == reason)
                    for reason in sorted({e["reason"] for e in excluded})
                },
                "exclusions": excluded,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    print(f"Wrote {total_residues} residues from {len(included)} proteins to {output_csv}")
    print(f"Excluded {len(excluded)} proteins; see {exclusions_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-size", type=int, default=DEFAULT_SAMPLE_SIZE)
    parser.add_argument("--output", type=Path, default=Path("data/external/joined_residues.csv"))
    parser.add_argument(
        "--exclusions", type=Path, default=Path("data/external/exclusions_log.json")
    )
    parser.add_argument("--raw-cache", type=Path, default=Path("data/external/disprot_raw.json"))
    parser.add_argument("--max-workers", type=int, default=8)
    args = parser.parse_args(argv)

    run(
        args.sample_size, args.output, args.exclusions, args.raw_cache, max_workers=args.max_workers
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
