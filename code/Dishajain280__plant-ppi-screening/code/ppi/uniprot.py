"""Resolve Arabidopsis gene symbols to UniProt accessions and sequences.

Uses the public UniProt REST API. Every lookup is cached under .cache/ so a
repeated run costs nothing and the pipeline still works offline once populated.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

ARABIDOPSIS_TAXID = 3702
API = "https://rest.uniprot.org"
CACHE_DIR = Path(".cache") / "uniprot"


class ResolutionError(RuntimeError):
    """Raised when a symbol cannot be resolved unambiguously."""


@dataclass(frozen=True)
class Entry:
    gene: str
    accession: str
    protein_name: str
    length: int
    reviewed: bool

    def as_dict(self) -> dict:
        return {
            "gene": self.gene,
            "accession": self.accession,
            "protein_name": self.protein_name,
            "length": self.length,
            "reviewed": self.reviewed,
        }


def _cache_path(gene: str) -> Path:
    safe = "".join(c if c.isalnum() else "_" for c in gene)
    return CACHE_DIR / f"{safe}.json"


def _get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "plant-ppi-screening/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read()


def _query(term: str, gene: str, taxid: int) -> list[Entry]:
    query = f"{term}:{gene} AND organism_id:{taxid} AND reviewed:true"
    url = f"{API}/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": query, "fields": "accession,protein_name,gene_names,length",
         "format": "json", "size": 10}
    )
    payload = json.loads(_get(url))
    return [
        Entry(
            gene=gene,
            accession=item["primaryAccession"],
            protein_name=(
                item.get("proteinDescription", {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value", "")
            ),
            length=item["sequence"]["length"],
            reviewed=True,
        )
        for item in payload.get("results", [])
    ]


def lookup(gene: str, taxid: int = ARABIDOPSIS_TAXID, use_cache: bool = True) -> list[Entry]:
    """All reviewed entries matching ``gene``, searched by exact name then synonym.

    The synonym pass matters in practice: Arabidopsis AGB1 is stored under the
    gene name GB1, so an exact-name search alone misses it.
    """
    cache = _cache_path(gene)
    if use_cache and cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
        cached = [Entry(**item) for item in data.get("results", [])]
        # An empty cache usually means the exact-name search failed, not that the
        # protein is absent, so let the synonym pass run before concluding.
        if cached:
            return cached

    entries: list[Entry] = []
    try:
        for term in ("gene_exact", "gene"):
            for entry in _query(term, gene, taxid):
                if all(entry.accession != seen.accession for seen in entries):
                    entries.append(entry)
            if entries:
                break
    except (urllib.error.URLError, TimeoutError) as exc:
        raise ResolutionError(
            f"{gene}: UniProt unreachable ({exc}). Populate .cache/uniprot/ "
            f"or pass the accession explicitly in candidates.csv."
        ) from exc

    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"results": [e.as_dict() for e in entries]}, indent=2),
                     encoding="utf-8")
    return entries


def resolve(gene: str, taxid: int = ARABIDOPSIS_TAXID, use_cache: bool = True) -> Entry:
    """Return the single reviewed UniProt entry for ``gene``.

    Raises ResolutionError with an actionable message when the symbol is
    missing, ambiguous, or when the network is unavailable and no cache exists.
    """
    entries = lookup(gene, taxid, use_cache)
    if not entries:
        raise ResolutionError(
            f"{gene}: no reviewed Arabidopsis entry. Check the symbol, or put the "
            f"accession directly in candidates.csv."
        )
    # More than one hit means the symbol is not unique; refuse to guess.
    if len(entries) > 1:
        options = ", ".join(f"{e.accession} ({e.protein_name or 'unnamed'})" for e in entries)
        raise ResolutionError(
            f"{gene}: ambiguous - {options}. Put the intended accession in candidates.csv."
        )
    return entries[0]


def fetch_fasta(accession: str, out_dir: Path, use_cache: bool = True) -> Path:
    """Download the canonical FASTA for ``accession`` and return the local path."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{accession}.fasta"
    if use_cache and path.exists():
        return path
    try:
        raw = _get(f"{API}/uniprotkb/{accession}.fasta").decode("utf-8")
    except (urllib.error.URLError, TimeoutError) as exc:
        raise ResolutionError(f"{accession}: could not download sequence ({exc})") from exc
    path.write_text(raw, encoding="utf-8")
    return path


def read_fasta(path: Path) -> str:
    """Return the sequence from a single-record FASTA file."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def slice_sequence(sequence: str, start: int | None, end: int | None) -> str:
    """1-based inclusive slice; None on both ends returns the full sequence."""
    if start is None and end is None:
        return sequence
    first = 1 if start is None else start
    last = len(sequence) if end is None else end
    if first < 1 or last > len(sequence) or first > last:
        raise ValueError(f"range {first}-{last} outside 1-{len(sequence)}")
    return sequence[first - 1:last]
