"""
Pfam-A.hmm database download, caching, and loading via pyhmmer.

The HMM database is loaded once into process memory (~500MB for all ~20k profiles)
and cached for the lifetime of the process. Subsequent scans reuse the cached profiles.

Clan mapping
------------
``load_clan_map()`` reads ``Pfam-A.hmm.dat`` (the sidecar metadata file distributed
alongside ``Pfam-A.hmm`` on the EBI FTP server).  Each entry in that file contains::

    #=GF AC   PFxxxxx.version
    #=GF CL   CLxxxx          ← present only when the family belongs to a clan

The function returns a dict mapping Pfam accession *without* version suffix
(e.g. ``"PF00069"``) to its clan identifier (e.g. ``"CL0016"``), or an empty string
when the family has no clan.  The result is cached in ``_clan_map_cache`` for the
lifetime of the process.

Configuration:
    PFAM_HMM_PATH env var overrides the default location at backend/data/pfam/Pfam-A.hmm
    PFAM_DAT_PATH env var overrides the default location at backend/data/pfam/Pfam-A.hmm.dat
"""
import os
import gzip
import shutil
import logging
from pathlib import Path
from typing import Optional

import pyhmmer

logger = logging.getLogger(__name__)

DEFAULT_PFAM_DIR = Path(__file__).parent.parent.parent / "data" / "pfam"
PFAM_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_DAT_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"

# Known kinase clans and accessions for the is_kinase flag
KINASE_CLANS = {"CL0016"}  # Protein kinase superfamily clan
KINASE_ACCESSIONS = {
    "PF00069",   # Pkinase
    "PF07714",   # Pkinase_Tyr (Tyrosine kinase)
    # PF00081 was listed here as "Lipid kinase" — it is Sod_Fe_N (iron/manganese
    # superoxide dismutase, review DOM-7). PI3/PI4 kinases (PF00454) are covered
    # by clan CL0016 above.
}


def get_pfam_hmm_path() -> Path:
    """Return the configured path to Pfam-A.hmm, raising if not found."""
    env_path = os.environ.get("PFAM_HMM_PATH")
    if env_path:
        p = Path(env_path)
        if p.exists():
            return p
        raise FileNotFoundError(f"PFAM_HMM_PATH={env_path} does not exist")

    default_path = DEFAULT_PFAM_DIR / "Pfam-A.hmm"
    if default_path.exists():
        return default_path

    raise FileNotFoundError(
        f"Pfam-A.hmm not found at {default_path}. "
        f"Set PFAM_HMM_PATH env var or run: "
        f"python -c \"from app.analysis.pfam_db import download_pfam_hmm; download_pfam_hmm()\""
    )


def get_pfam_dat_path() -> Optional[Path]:
    """Return the configured path to Pfam-A.hmm.dat, or None if not found.

    Unlike the HMM file, the .dat file is optional — its absence degrades
    clan assignment to None rather than raising an error.
    """
    env_path = os.environ.get("PFAM_DAT_PATH")
    if env_path:
        p = Path(env_path)
        return p if p.exists() else None

    default_path = DEFAULT_PFAM_DIR / "Pfam-A.hmm.dat"
    return default_path if default_path.exists() else None


def download_pfam_hmm(dest_dir: Optional[Path] = None) -> Path:
    """Download Pfam-A.hmm.gz from EBI FTP and decompress it."""
    import httpx

    dest_dir = dest_dir or DEFAULT_PFAM_DIR
    dest_dir.mkdir(parents=True, exist_ok=True)

    gz_path = dest_dir / "Pfam-A.hmm.gz"
    hmm_path = dest_dir / "Pfam-A.hmm"

    if hmm_path.exists():
        logger.info(f"Pfam-A.hmm already exists at {hmm_path}")
        return hmm_path

    logger.info(f"Downloading Pfam-A.hmm.gz from {PFAM_HMM_URL}...")
    with httpx.stream("GET", PFAM_HMM_URL, follow_redirects=True, timeout=300.0) as response:
        response.raise_for_status()
        with open(gz_path, "wb") as f:
            for chunk in response.iter_bytes(chunk_size=8192):
                f.write(chunk)

    logger.info("Decompressing Pfam-A.hmm.gz...")
    with gzip.open(gz_path, "rb") as f_in:
        with open(hmm_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()
    logger.info(f"Pfam-A.hmm ready at {hmm_path}")
    return hmm_path


def download_pfam_dat(dest_dir: Optional[Path] = None) -> Path:
    """Download Pfam-A.hmm.dat.gz from EBI FTP and decompress it.

    The .dat file is the sidecar metadata companion to Pfam-A.hmm.  It
    contains clan assignments and other per-family metadata but is NOT used
    for HMM scanning — only for post-processing (clan lookup).

    Size: ~100 MB compressed, ~400 MB decompressed.
    """
    import httpx

    dest_dir = dest_dir or DEFAULT_PFAM_DIR
    dest_dir.mkdir(parents=True, exist_ok=True)

    gz_path  = dest_dir / "Pfam-A.hmm.dat.gz"
    dat_path = dest_dir / "Pfam-A.hmm.dat"

    if dat_path.exists():
        logger.info(f"Pfam-A.hmm.dat already exists at {dat_path}")
        return dat_path

    logger.info(f"Downloading Pfam-A.hmm.dat.gz from {PFAM_DAT_URL}...")
    with httpx.stream("GET", PFAM_DAT_URL, follow_redirects=True, timeout=600.0) as response:
        response.raise_for_status()
        with open(gz_path, "wb") as f:
            for chunk in response.iter_bytes(chunk_size=65536):
                f.write(chunk)

    logger.info("Decompressing Pfam-A.hmm.dat.gz...")
    with gzip.open(gz_path, "rb") as f_in:
        with open(dat_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()
    logger.info(f"Pfam-A.hmm.dat ready at {dat_path}")
    return dat_path


# ── Clan map (accession-without-version → clan ID) ────────────────────────────

_clan_map_cache: Optional[dict] = None


def load_clan_map() -> dict:
    """Parse Pfam-A.hmm.dat and return a {pfam_accession: clan_id} mapping.

    The returned dict maps each Pfam accession (without version suffix, e.g.
    ``"PF00069"``) to its clan identifier (e.g. ``"CL0016"``), or to an empty
    string ``""`` for families with no clan.

    The .dat file uses HMMER Stockholm-style notation::

        #=GF ID   Pkinase
        #=GF AC   PF00069.29      ← accession with version
        #=GF CL   CL0016          ← present only when family belongs to a clan
        //                        ← record separator

    If the file does not exist locally, an automatic download is attempted.
    If download fails (e.g. offline), an empty dict is returned with a warning
    so that the rest of the pipeline can still function without clan info.

    The result is cached in ``_clan_map_cache`` for the process lifetime.
    """
    global _clan_map_cache
    if _clan_map_cache is not None:
        return _clan_map_cache

    dat_path = get_pfam_dat_path()
    if dat_path is None:
        logger.info(
            "Pfam-A.hmm.dat not found locally — attempting download. "
            "To skip, set PFAM_DAT_PATH to an existing file or an empty string."
        )
        try:
            dat_path = download_pfam_dat()
        except Exception as exc:
            logger.warning(
                f"Could not download Pfam-A.hmm.dat: {exc}. "
                "Clan information will be unavailable (pfam_clan=None)."
            )
            _clan_map_cache = {}
            return _clan_map_cache

    logger.info(f"Parsing clan assignments from {dat_path}...")
    clan_map: dict = {}
    current_acc: Optional[str] = None
    current_clan: Optional[str] = None

    with open(dat_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#=GF AC"):
                # e.g. "#=GF AC   PF00069.29"
                raw_acc = line.split()[-1]
                current_acc = raw_acc.split(".")[0]  # strip .version
                current_clan = None
            elif line.startswith("#=GF CL"):
                # e.g. "#=GF CL   CL0016"
                current_clan = line.split()[-1]
            elif line.startswith("//"):
                # End of record
                if current_acc is not None:
                    clan_map[current_acc] = current_clan or ""
                current_acc = None
                current_clan = None

    logger.info(
        f"Loaded {len(clan_map):,} Pfam family→clan entries "
        f"({sum(1 for v in clan_map.values() if v):,} have a clan assignment)"
    )
    _clan_map_cache = clan_map
    return _clan_map_cache


# Process-lifetime singleton cache for the HMM database
_hmm_db_cache: Optional[list] = None


def load_hmm_database() -> list:
    """Load the entire Pfam HMM database into memory (cached after first call)."""
    global _hmm_db_cache
    if _hmm_db_cache is not None:
        return _hmm_db_cache

    hmm_path = get_pfam_hmm_path()
    logger.info(f"Loading Pfam HMM database from {hmm_path}...")

    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hmm_file:
        _hmm_db_cache = list(hmm_file)

    logger.info(f"Loaded {len(_hmm_db_cache)} Pfam HMM profiles")
    return _hmm_db_cache


def is_kinase_domain(pfam_accession: str, pfam_clan: Optional[str]) -> bool:
    """Check if a Pfam domain is a kinase domain."""
    if pfam_accession in KINASE_ACCESSIONS:
        return True
    if pfam_clan and pfam_clan in KINASE_CLANS:
        return True
    return False
