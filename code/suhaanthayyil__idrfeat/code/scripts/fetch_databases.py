#!/usr/bin/env python
"""Download the annotation databases idrfeat can use and record provenance.

Stdlib only, so it runs without the project installed. Idempotent: a file already present and
non-empty is skipped. A source that fails to download prints a manual instruction instead of
failing the whole run. Versions, dates, URLs, and checksums are written to SOURCES.txt.

Usage:
    python scripts/fetch_databases.py --out data/raw
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

DBPTM = "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment"
UNIPROT = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/"
    "reference_proteomes/Eukaryota/UP000005640"
)
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38"


@dataclass
class Source:
    filename: str
    url: str
    note: str
    fallback: str | None = None


SOURCES = [
    Source(
        "elm_instances.tsv",
        "http://elm.eu.org/instances.tsv?q=*&taxon=Homo+sapiens",
        "ELM curated motif instances, Homo sapiens. Academic, non-commercial license.",
    ),
    Source(
        "elm_classes.tsv",
        "http://elm.eu.org/elms/elms_index.tsv",
        "ELM motif class definitions.",
    ),
    Source(
        "dbPTM_Phosphorylation.gz", f"{DBPTM}/Phosphorylation.gz", "dbPTM phosphorylation sites."
    ),
    Source("dbPTM_Acetylation.gz", f"{DBPTM}/Acetylation.gz", "dbPTM acetylation sites."),
    Source("dbPTM_Ubiquitination.gz", f"{DBPTM}/Ubiquitination.gz", "dbPTM ubiquitination sites."),
    Source("dbPTM_Methylation.gz", f"{DBPTM}/Methylation.gz", "dbPTM methylation sites."),
    Source(
        "phasepro_full.json",
        "https://phasepro.elte.hu/download_full.json",
        "PhaSePro curated phase-separation regions, keyed by UniProt accession.",
    ),
    Source(
        "drllps_LLPS.txt",
        "http://llps.biocuckoo.cn/download/LLPS.txt",
        "DrLLPS LLPS-associated proteins and their roles.",
    ),
    Source(
        "UP000005640_9606.fasta.gz",
        f"{UNIPROT}/UP000005640_9606.fasta.gz",
        "UniProt human reference proteome canonical sequences.",
    ),
    Source(
        "hg38.phyloP100way.bw",
        f"{UCSC}/phyloP100way/hg38.phyloP100way.bw",
        "UCSC phyloP 100-way vertebrate per-base conservation. Optional, large.",
    ),
    Source(
        "hg38.phastCons100way.bw",
        f"{UCSC}/phastCons100way/hg38.phastCons100way.bw",
        "UCSC phastCons 100-way vertebrate conserved-element probability. Optional, large.",
    ),
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path) -> bool:
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "idrfeat-fetch"})
        with urllib.request.urlopen(req, timeout=120) as resp, open(tmp, "wb") as fh:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
        tmp.rename(dest)
        return True
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        print(f"  download failed: {exc}")
        if tmp.exists():
            tmp.unlink()
        return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/raw", help="directory to download into")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    today = dt.date.today().isoformat()
    records = []

    for source in SOURCES:
        dest = out / source.filename
        if dest.exists() and dest.stat().st_size > 0:
            print(f"skip {source.filename} (present)")
        else:
            print(f"fetch {source.filename}")
            if not _download(source.url, dest):
                if source.fallback:
                    print(f"  {source.fallback}")
                else:
                    print(f"  download manually from {source.url} into {out}/")
                continue
        records.append(f"{source.filename}\t{today}\t{source.url}\t{_sha256(dest)}\t{source.note}")

    sources_txt = out / "SOURCES.txt"
    header = "# filename\tfetch_date\turl\tsha256\tnote\n"
    sources_txt.write_text(header + "\n".join(records) + "\n")
    print(f"wrote provenance to {sources_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
