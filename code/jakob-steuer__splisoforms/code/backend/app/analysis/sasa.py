"""
Solvent Accessible Surface Area (SASA) calculation utilities.

Uses gemmi to parse mmCIF structures and freesasa to calculate per-residue SASA.
FreeSASA's Python API doesn't support mmCIF directly, so structures are converted
to PDB format via gemmi in a temporary file.

AF3/AFDB structures are single-chain predictions, so we use model[0], chain[0].
"""
import logging
import tempfile
from pathlib import Path
from typing import Dict, Optional

import gemmi
import freesasa

from app.services.storage import STATIC_DATA_DIR

logger = logging.getLogger(__name__)


def calculate_residue_sasa(cif_path: Path) -> Dict[int, float]:
    """
    Calculate per-residue total SASA (Å²) from an mmCIF structure file.

    Loads the mmCIF with gemmi, converts to a temporary PDB file, then
    runs FreeSASA to get per-residue solvent accessibility.

    Args:
        cif_path: Path to the mmCIF (.cif) file.

    Returns:
        Dict mapping 1-based residue position to total SASA in Å².
        Empty dict if calculation fails.
    """
    if not cif_path.exists():
        logger.warning(f"Structure file not found: {cif_path}")
        return {}

    try:
        # Parse mmCIF with gemmi
        structure = gemmi.read_structure(str(cif_path))
        if len(structure) == 0:
            logger.warning(f"No models in structure: {cif_path}")
            return {}

        # Write to temporary PDB for FreeSASA
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=True, mode="w") as tmp:
            structure.write_pdb(tmp.name)
            tmp.flush()

            # Calculate SASA with FreeSASA
            fs_structure = freesasa.Structure(tmp.name)
            result = freesasa.calc(fs_structure)

        # Extract per-residue SASA
        residue_areas = result.residueAreas()
        sasa_map: Dict[int, float] = {}

        # AF structures are single-chain — iterate over first (and usually only) chain
        for chain_id, residues in residue_areas.items():
            for res_num_str, area in residues.items():
                try:
                    res_num = int(res_num_str)
                    sasa_map[res_num] = area.total
                except ValueError:
                    # Skip residues with insertion codes or non-numeric IDs
                    continue
            break  # Only process first chain

        logger.debug(f"SASA calculated for {len(sasa_map)} residues from {cif_path.name}")
        return sasa_map

    except Exception as e:
        logger.error(f"SASA calculation failed for {cif_path}: {e}")
        return {}


def calculate_residue_plddt(cif_path: Path) -> Dict[int, float]:
    """
    Extract per-residue pLDDT (CA-atom B-factor) from an mmCIF structure.

    Both AF3 and AFDB store the per-residue confidence score (pLDDT, 0-100)
    in the B-factor column of CA atoms.  Returns {1-based residue position
    → pLDDT in [0, 100]}.  Empty dict if extraction fails.

    pLDDT < 70 is the conventional "low-confidence" cutoff (Tunyasuvunakool
    et al., Nature 2021); residues below this should be treated cautiously
    in SASA / SASA-delta analyses because side-chain placement is unreliable.
    """
    if not cif_path.exists():
        return {}

    try:
        structure = gemmi.read_structure(str(cif_path))
        if len(structure) == 0:
            return {}
        plddt_map: Dict[int, float] = {}
        # AF structures are single-chain; pull pLDDT from the CA atom of every
        # residue in the first chain of the first model.
        model = structure[0]
        if not len(model):
            return {}
        chain = model[0]
        for residue in chain:
            for atom in residue:
                if atom.name == "CA":
                    res_num = residue.seqid.num
                    plddt_map[res_num] = float(atom.b_iso)
                    break
        return plddt_map
    except Exception as e:
        logger.error(f"pLDDT extraction failed for {cif_path}: {e}")
        return {}


def resolve_structure_path(af3_cif_url: Optional[str]) -> Optional[Path]:
    """
    Resolve a Structure.af3_cif_url to an absolute file path on disk.
    
    Standardizes "/static/data/..." or "/static/af3_results/..." URLs
    to their physical locations in the backend filesystem.
    """
    if not af3_cif_url:
        return None

    # Resolve path relative to backend/static
    if af3_cif_url.startswith("/static/data/"):
        clean_path = af3_cif_url.replace("/static/data/", "", 1).lstrip('/')
        filepath = STATIC_DATA_DIR / clean_path
    elif af3_cif_url.startswith("static/data/"):
        clean_path = af3_cif_url.replace("static/data/", "", 1).lstrip('/')
        filepath = STATIC_DATA_DIR / clean_path
    elif af3_cif_url.startswith("/static/af3_results/"):
        clean_path = af3_cif_url.replace("/static/af3_results/", "", 1).lstrip('/')
        filepath = (STATIC_DATA_DIR.parent / "af3_results") / clean_path
    elif af3_cif_url.startswith("/static/"):
        # Generic /static/ mapping
        clean_path = af3_cif_url.replace("/static/", "", 1).lstrip('/')
        filepath = STATIC_DATA_DIR.parent / clean_path
    else:
        # Fallback to filename-only resolution if path is unknown
        filename = af3_cif_url.split("/")[-1]
        filepath = STATIC_DATA_DIR / filename

    if filepath.exists():
        return filepath

    logger.debug(f"Structure path resolution failed for: {af3_cif_url} (Checked {filepath})")
    return None
