"""
_sources_registry.py — declarative registry of the data sources ingested into
SPLISOFORMS. Single source of truth for both `reset_db.py` (Source-table seeding)
and `ingest_from_metadata.py` (ingest order + per-source settings).

Adding a future dataset:
  1. put structures under  data/sources/<dir>/structures/<md5>/<md5>_model.cif
  2. produce  data/sources/<dir>/metadata.tsv  conforming to the contract in
     docs/unified_metadata_schema.md
  3. append one SourceDef entry below (pick the next `order`)
No ingest code changes are needed.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

# Repo root = .../splisoforms ; this file is at backend/scripts/_sources_registry.py
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SOURCES_DIR = REPO_ROOT / "data" / "sources"


@dataclass(frozen=True)
class SourceDef:
    key: str                 # short slug + data/sources/<dir> name (key == dir)
    name: str                # Source.name (unique) — must match reset_db seed
    source_type: str         # SourceTypeEnum value
    is_reference: bool       # reference layer? (always ingested; dedup anchor)
    order: int               # ingest order — reference first
    struct_url_prefix: str   # /static prefix that serves this source's structures
    cancer_type: str = ""
    tissue_type: str = ""
    sequencing_platform: str = ""
    notes: str = ""

    @property
    def dir(self) -> Path:
        return SOURCES_DIR / self.key

    @property
    def metadata_tsv(self) -> Path:
        return self.dir / "metadata.tsv"

    @property
    def structures_dir(self) -> Path:
        return self.dir / "structures"


SOURCES: List[SourceDef] = [
    SourceDef(
        key="gencode_v47",
        name="GENCODE v47",
        source_type="other",
        is_reference=True,
        order=1,
        struct_url_prefix="/static/data/structures/gencode_v47",
        notes=("GENCODE v47 (= Ensembl 114) protein-coding reference. The dedup "
               "anchor: ingested first so cross-source identical proteins attach "
               "here as memberships."),
    ),
    SourceDef(
        key="ccrcc",
        name="ccRCC PacBio MAS-seq (Karakulak 2025)",
        source_type="pacbio_longread",
        is_reference=False,
        order=2,
        struct_url_prefix="/static/data/structures/ccrcc",
        cancer_type="ccRCC",
        tissue_type="kidney",
        sequencing_platform="PacBio MAS-seq",
        notes="Long-read isoforms from ccRCC patient samples (Normal / ccRCC2..5).",
    ),
    SourceDef(
        key="breast_cancer",
        name="Breast Cancer Nanopore ESPRESSO (Head 2026)",
        source_type="ont_longread",
        is_reference=False,
        order=3,
        struct_url_prefix="/static/data/structures/bc",
        cancer_type="breast",
        tissue_type="breast",
        sequencing_platform="Nanopore",
        notes=("Long-read isoforms from breast tumor / healthy breast / cultured "
               "fibroblast samples."),
    ),
]

SOURCES_BY_KEY = {s.key: s for s in SOURCES}


def ordered_sources() -> List[SourceDef]:
    return sorted(SOURCES, key=lambda s: s.order)
