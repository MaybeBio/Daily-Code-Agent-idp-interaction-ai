"""
PTM site data management — dual-source loader.

Primary source: IsoProDB (Mapping) + UniProtKB Features REST API.
Fallback: Automated UniProt retrieval if dedicated IsoProDB PTM files are unavailable.

Note: PhosphoSitePlus TSV files are no longer used due to licensing restrictions.
"""
import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

import httpx

logger = logging.getLogger(__name__)

DEFAULT_PTM_DIR = Path(__file__).parent.parent.parent / "data" / "ptm"

# MOD_RSD format: "S100-p", "K270-ub", "Y302-p", "K100-ac"
MOD_RSD_PATTERN = re.compile(r"^([A-Z])(\d+)-(\w+)$")

# Map PhosphoSitePlus modification abbreviations to canonical names
PSP_MOD_MAP = {
    "p": "phosphorylation",
    "ub": "ubiquitination",
    "ac": "acetylation",
    "m1": "methylation",
    "m2": "dimethylation",
    "m3": "trimethylation",
    "sm": "sumoylation",
    "ga": "O-GalNAc",
    "gl": "O-GlcNAc",
}

# PTM types we care about for this analysis
RELEVANT_PTM_TYPES = {"phosphorylation", "ubiquitination", "acetylation"}

# Map UniProt modification descriptions to canonical names
UNIPROT_MOD_MAP = {
    "phosphoserine": "phosphorylation",
    "phosphothreonine": "phosphorylation",
    "phosphotyrosine": "phosphorylation",
    "n6-acetyllysine": "acetylation",
    "glycyl lysine isopeptide": "ubiquitination",
}


@dataclass
class PTMSite:
    """A single known PTM site."""
    uniprot_id: str
    residue_position: int   # 1-based
    residue_type: str       # single-letter AA code
    ptm_type: str           # phosphorylation, ubiquitination, acetylation
    source_db: str          # "PhosphoSitePlus" or "UniProt"


def get_ptm_data_dir() -> Path:
    """Return the configured path to PTM data files."""
    env_path = os.environ.get("PTM_DATA_PATH")
    if env_path:
        return Path(env_path)
    return DEFAULT_PTM_DIR


# Process-lifetime cache for PhosphoSitePlus data
_psp_cache: Optional[Dict[str, List[PTMSite]]] = None


def load_phosphositeplus() -> Dict[str, List[PTMSite]]:
    """
    Load PhosphoSitePlus datasets into a dict keyed by UniProt accession.

    Parses Phosphorylation, Ubiquitination, and Acetylation datasets.
    Returns empty dict if files are not found (fallback to UniProt API).
    """
    global _psp_cache
    if _psp_cache is not None:
        return _psp_cache

    import pandas as pd

    data_dir = get_ptm_data_dir()
    result: Dict[str, List[PTMSite]] = {}

    files = {
        "Phosphorylation_site_dataset.gz": "phosphorylation",
        "Ubiquitination_site_dataset.gz": "ubiquitination",
        "Acetylation_site_dataset.gz": "acetylation",
    }

    loaded_any = False
    for filename, expected_type in files.items():
        filepath = data_dir / filename
        if not filepath.exists():
            logger.debug(f"PhosphoSitePlus file not found: {filepath}")
            continue

        try:
            # PhosphoSitePlus files have 3 header lines before the actual data
            df = pd.read_csv(filepath, sep="\t", skiprows=3, low_memory=False)

            # Filter to human entries only
            if "ORGANISM" in df.columns:
                df = df[df["ORGANISM"] == "human"]

            for _, row in df.iterrows():
                acc_id = str(row.get("ACC_ID", "")).strip()
                mod_rsd = str(row.get("MOD_RSD", "")).strip()

                if not acc_id or not mod_rsd:
                    continue

                # Strip isoform suffix (e.g., "P41181-8" → "P41181")
                base_acc = acc_id.split("-")[0]

                match = MOD_RSD_PATTERN.match(mod_rsd)
                if not match:
                    continue

                aa, pos_str, mod_abbr = match.groups()
                ptm_type = PSP_MOD_MAP.get(mod_abbr, mod_abbr)

                if ptm_type not in RELEVANT_PTM_TYPES:
                    continue

                site = PTMSite(
                    uniprot_id=base_acc,
                    residue_position=int(pos_str),
                    residue_type=aa,
                    ptm_type=ptm_type,
                    source_db="PhosphoSitePlus",
                )

                if base_acc not in result:
                    result[base_acc] = []
                result[base_acc].append(site)

            loaded_any = True
            logger.info(f"Loaded {filename}: {len(df)} human entries")
        except Exception as e:
            logger.warning(f"Failed to parse {filename}: {e}")

    if loaded_any:
        _psp_cache = result
        total = sum(len(v) for v in result.values())
        logger.info(f"PhosphoSitePlus loaded: {total} PTM sites across {len(result)} proteins")
    else:
        _psp_cache = {}
        logger.info("No PhosphoSitePlus files found, will use UniProt API fallback")

    return _psp_cache


def fetch_uniprot_ptms(uniprot_id: str) -> List[PTMSite]:
    """
    Fetch PTM annotations from UniProt REST API as a fallback.

    Queries the UniProt entry and extracts "Modified residue" features
    for phosphorylation, ubiquitination, and acetylation.
    """
    try:
        with httpx.Client(timeout=10.0) as client:
            res = client.get(
                f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
                headers={"Accept": "application/json"},
            )
            if res.status_code != 200:
                logger.debug(f"UniProt API returned {res.status_code} for {uniprot_id}")
                return []

            data = res.json()
    except Exception as e:
        logger.warning(f"UniProt API request failed for {uniprot_id}: {e}")
        return []

    sites: List[PTMSite] = []
    features = data.get("features", [])

    for feat in features:
        if feat.get("type") != "Modified residue":
            continue

        location = feat.get("location", {})
        start = location.get("start", {}).get("value")
        end = location.get("end", {}).get("value")

        # PTM sites are single-residue features
        if start is None or start != end:
            continue

        description = feat.get("description", "").lower()

        # Map the UniProt description to our canonical PTM type
        ptm_type = None
        for keyword, mapped_type in UNIPROT_MOD_MAP.items():
            if keyword in description:
                ptm_type = mapped_type
                break

        if ptm_type is None or ptm_type not in RELEVANT_PTM_TYPES:
            continue

        # Try to determine residue type from description
        residue_type = "?"
        if "serine" in description or "phosphoserine" in description:
            residue_type = "S"
        elif "threonine" in description or "phosphothreonine" in description:
            residue_type = "T"
        elif "tyrosine" in description or "phosphotyrosine" in description:
            residue_type = "Y"
        elif "lysine" in description or "acetyllysine" in description:
            residue_type = "K"

        sites.append(PTMSite(
            uniprot_id=uniprot_id,
            residue_position=int(start),
            residue_type=residue_type,
            ptm_type=ptm_type,
            source_db="UniProt",
        ))

    logger.info(f"UniProt API: {len(sites)} PTM sites for {uniprot_id}")
    return sites


def fetch_proteomics_ptms(uniprot_id: str) -> List[PTMSite]:
    """
    Fetch high-throughput PTM annotations (from PRIDE/PTMeXchange) via the Proteins REST API.
    """
    try:
        with httpx.Client(timeout=10.0) as client:
            res = client.get(
                f"https://www.ebi.ac.uk/proteins/api/proteomics/ptm/{uniprot_id}",
                headers={"Accept": "application/json"},
            )
            if res.status_code != 200:
                return []
            
            data = res.json()
            if isinstance(data, dict) and "errorMessage" in data:
                 return []
                 
            features = []
            if isinstance(data, dict):
                features = data.get('features', [])
            elif isinstance(data, list):
                if len(data) > 0 and isinstance(data[0], dict) and 'features' in data[0]:
                    features = data[0].get('features', [])
                else:
                    features = data
                    
    except Exception as e:
        logger.warning(f"Proteomics PTM API request failed for {uniprot_id}: {e}")
        return []

    sites: List[PTMSite] = []

    for feat in features:
        if feat.get("type") != "PROTEOMICS_PTM":
            continue
            
        try:
            begin = int(feat.get("begin", 0))
        except ValueError:
            continue
            
        peptide = feat.get("peptide", "")
        
        for ptm in feat.get("ptms", []):
            ptm_name = ptm.get("name", "").lower()
            pep_pos = ptm.get("position")
            
            if not pep_pos:
                continue
                
            global_pos = begin + pep_pos - 1
            
            # Map the PTM name to canonical type
            canonical_type = None
            if "phospho" in ptm_name:
                canonical_type = "phosphorylation"
            elif "ubiquitin" in ptm_name:
                canonical_type = "ubiquitination"
            elif "acetyl" in ptm_name:
                canonical_type = "acetylation"
                
            if not canonical_type:
                continue
                
            res_type = peptide[pep_pos - 1] if 0 < pep_pos <= len(peptide) else "?"
            
            sites.append(PTMSite(
                uniprot_id=uniprot_id,
                residue_position=global_pos,
                residue_type=res_type,
                ptm_type=canonical_type,
                source_db="ProteomicsDB",
            ))
            
    logger.info(f"Proteomics API: {len(sites)} PTM sites for {uniprot_id}")
    return sites



def get_ptm_sites(uniprot_id: str) -> List[PTMSite]:
    """
    Get PTM sites for a protein.

    Site-source policy:
    1. Uses UniProtKB Features REST API for manually curated sites.
    2. Uses Proteins Proteomics API for high-throughput MS (PRIDE/PTMeXchange) sites.

    Results are cached per accession, in process and on disk
    (backend/data/ptm_cache/ptm_sites_by_accession.json), so a backfill fetches
    each accession once instead of once per isoform (review PTM-4/PTM-7), and a
    re-run reproduces the same site set. A failed/empty fetch is NOT cached, so
    transient errors are retried next time rather than frozen as "no PTMs".
    """
    cache = _sites_cache()
    if uniprot_id in cache:
        return list(cache[uniprot_id])

    curated = fetch_uniprot_ptms(uniprot_id)
    high_throughput = fetch_proteomics_ptms(uniprot_id)

    # Deduplicate by position and type
    seen = set()
    unique_sites = []

    for site in curated + high_throughput:
        key = (site.residue_position, site.ptm_type)
        if key not in seen:
            seen.add(key)
            unique_sites.append(site)

    if unique_sites:
        cache[uniprot_id] = unique_sites
        _maybe_flush_sites_cache()
    return unique_sites


# ── per-accession site cache (process + disk) ───────────────────────────────
SITES_CACHE_PATH = Path(__file__).parent.parent.parent / "data" / "ptm_cache" / "ptm_sites_by_accession.json"
_SITES_CACHE: Optional[Dict[str, List[PTMSite]]] = None
_SITES_CACHE_DIRTY = 0
_SITES_CACHE_FLUSH_EVERY = 100


def _sites_cache() -> Dict[str, List[PTMSite]]:
    global _SITES_CACHE
    if _SITES_CACHE is None:
        _SITES_CACHE = {}
        try:
            if SITES_CACHE_PATH.exists():
                import json
                raw = json.loads(SITES_CACHE_PATH.read_text())
                for acc, rows in raw.items():
                    _SITES_CACHE[acc] = [
                        PTMSite(uniprot_id=acc, residue_position=int(r[0]),
                                residue_type=r[1], ptm_type=r[2], source_db=r[3])
                        for r in rows
                    ]
                logger.info(f"PTM site cache: {len(_SITES_CACHE)} accessions from {SITES_CACHE_PATH}")
        except Exception as e:                      # corrupt cache → ignore, rebuild
            logger.warning(f"PTM site cache unreadable ({e}); starting empty")
            _SITES_CACHE = {}
    return _SITES_CACHE


def _maybe_flush_sites_cache(force: bool = False) -> None:
    global _SITES_CACHE_DIRTY
    _SITES_CACHE_DIRTY += 1
    if not force and _SITES_CACHE_DIRTY < _SITES_CACHE_FLUSH_EVERY:
        return
    flush_sites_cache()


def flush_sites_cache() -> None:
    """Write the per-accession site cache to disk (atomic replace)."""
    global _SITES_CACHE_DIRTY
    if not _SITES_CACHE:
        return
    import json, os, tempfile
    SITES_CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    payload = {acc: [[s.residue_position, s.residue_type, s.ptm_type, s.source_db] for s in sites]
               for acc, sites in _SITES_CACHE.items()}
    fd, tmp = tempfile.mkstemp(dir=str(SITES_CACHE_PATH.parent), suffix=".tmp")
    with os.fdopen(fd, "w") as fh:
        json.dump(payload, fh)
    os.replace(tmp, SITES_CACHE_PATH)
    _SITES_CACHE_DIRTY = 0
