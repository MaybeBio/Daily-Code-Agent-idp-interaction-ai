import os
import re
import json
import logging
import numpy as np
from pathlib import Path
from sqlalchemy.orm import Session
from app.models.bio import Isoform, StructuralMetrics, Structure

logger = logging.getLogger(__name__)

def parse_plddt_from_cif(cif_path):
    """
    Portable parser to extract per-residue pLDDT from AlphaFold mmCIF files.
    Skips external dependencies like gemmi.
    """
    plddt_by_res = {} # {seq_id: plddt}
    try:
        with open(cif_path, 'r') as f:
            in_atom_site = False
            col_map = {}
            col_idx = 0
            
            for line in f:
                if line.startswith('loop_'):
                    in_atom_site = False
                    col_map = {}
                    col_idx = 0
                    continue
                
                if line.startswith('_atom_site.'):
                    in_atom_site = True
                    col_map[line.strip()] = col_idx
                    col_idx += 1
                    continue
                
                if in_atom_site and (line.startswith('ATOM') or line.startswith('HETATM')):
                    parts = line.split()
                    try:
                        seq_id = int(parts[col_map['_atom_site.label_seq_id']])
                        plddt = float(parts[col_map['_atom_site.B_iso_or_equiv']])
                        # CIF has an entry for every atom; just take the first one per residue
                        if seq_id not in plddt_by_res:
                            plddt_by_res[seq_id] = plddt
                    except (ValueError, KeyError, IndexError):
                        continue
                        
        if not plddt_by_res:
            return None
            
        # Ensure it's sorted and contains all residues (filling gaps if any, though AFDB shouldn't have any)
        max_res = max(plddt_by_res.keys())
        plddt_series = [plddt_by_res.get(i, 0.0) for i in range(1, max_res + 1)]
        return plddt_series
        
    except Exception as e:
        logger.error(f"Portable CIF parse failed for {cif_path}: {e}")
        return None

def analyze_isoform_portable(db: Session, isoform: Isoform):
    """
    Portable version of structural impact analysis (pLDDT only).
    """
    metrics = isoform.structural_metrics
    if not metrics:
        metrics = StructuralMetrics(isoform_id=isoform.id)
        db.add(metrics)
        db.flush()

    struct = db.query(Structure).filter(Structure.isoform_id == isoform.id).first()
    if not struct:
        return None

    # Resolve local path
    # The backend directory (…/backend) that `/static/...` URLs are relative to.
    # Override with SPLISOFORMS_BACKEND_ROOT when the app is not run from its
    # source tree (e.g. in a container the static tree is mounted elsewhere).
    base_path = Path(
        os.environ.get("SPLISOFORMS_BACKEND_ROOT")
        or Path(__file__).resolve().parents[2]
    )
    cif_path = None
    if struct.af3_cif_url:
        cif_path = base_path / struct.af3_cif_url.lstrip('/')
    
    if not cif_path or not cif_path.exists():
        # Try legacy static/data path
        cif_path = base_path / "static/data" / f"{struct.isoform.transcript_id.replace('AFDB_', '')}.cif"
    
    if not cif_path.exists():
        return None

    plddt_series = parse_plddt_from_cif(cif_path)
    if plddt_series:
        metrics.plddt_series = plddt_series
        metrics.plddt_mean = float(np.mean(plddt_series))
        metrics.plddt_min = float(np.min(plddt_series))
        metrics.plddt_median = float(np.median(plddt_series))
        db.commit()
        return True
    
    return False
