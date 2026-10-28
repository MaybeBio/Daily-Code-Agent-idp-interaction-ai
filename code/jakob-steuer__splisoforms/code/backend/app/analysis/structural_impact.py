"""
Per-residue structural impact analysis using DSSP secondary structure assignment.

Pipeline
--------
  1. Parse AF3 / AFDB mmCIF with gemmi; extract N/CA/C/O backbone coordinates.
  2. Run PyDSSP to assign 3-class secondary structure (H=helix, E=sheet,
     C=coil).  Disordered residues (`-`) collapse to C.
  3. Extract phi/psi backbone torsion angles and per-residue pLDDT from the
     CA B-factor column.
  4. When a canonical structure exists, align novel-vs-canonical via the
     **splice-aware aligner** (same algorithm used by the MSA panel, ported
     in audit-2).  For each novel residue, classify the structural change:
       same                — DSSP unchanged
       coil_to_helix / …   — H/E/C transition between confidently-modelled
                             residues (min pLDDT ≥ 70)
       low_confidence      — either canonical or novel residue has pLDDT
                             < 70; the call is unreliable so we refuse
                             to commit to a transition direction
       novel_residue       — the novel residue has no canonical counterpart
                             (splicing introduced new sequence here)

Why 3-class?
  PyDSSP can emit 8 codes (G=3-10, I=π, T=turn, S=bend, B=isolated bridge);
  we collapse to H/E/C for storage simplicity.  Extended codes (audit-4
  Option B) would land in a separate `dssp_extended` column.

Exon-level aggregation is deferred — exon_number remains NULL.
"""
import logging
import math
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

import numpy as np
import gemmi
from sqlalchemy.orm import Session

from app.models.bio import Isoform, ExonStructuralImpact, Structure, StructuralMetrics
from app.analysis.pfam_domains import find_canonical_isoform, get_canonical_sequence
from app.analysis.sasa import resolve_structure_path

logger = logging.getLogger(__name__)

# pydssp with out_type='c3' (default) returns character strings: 'H', 'E', 'C', or '-'
# '-' indicates a disordered/unresolved residue; map it to 'C' (coil).
# Valid codes to store in DB:
VALID_DSSP_CODES = {"H", "E", "C"}

# pLDDT below this → structural-change call is treated as `low_confidence`
# (Tunyasuvunakool et al., Nature 2021).  Applied to both canonical and
# novel residues — if either is unreliable, the call is unreliable.
LOW_PLDDT_CUTOFF = 70.0


def _normalize_dssp(code) -> str:
    """Normalize a pydssp character code to one of H/E/C."""
    s = str(code).strip().upper()
    return s if s in VALID_DSSP_CODES else "C"


@dataclass
class ResidueStructure:
    """Per-residue structural features extracted from a structure."""
    residue_number: int   # 1-based
    dssp: str             # H, E, or C
    phi: Optional[float]  # degrees in (-180, 180]; None at the N-terminus
    psi: Optional[float]  # degrees in (-180, 180]; None at the C-terminus
    plddt: Optional[float]


def extract_backbone_coords(cif_path) -> Optional[np.ndarray]:
    """
    Extract N, CA, C, O backbone atom coordinates from an mmCIF file.

    PyDSSP expects an array of shape (n_residues, 4, 3) where the 4 atoms
    are [N, CA, C, O] in order.

    Args:
        cif_path: Path to mmCIF file.

    Returns:
        NumPy array of shape (n_residues, 4, 3), or None if extraction fails.
    """
    try:
        structure = gemmi.read_structure(str(cif_path))
        if len(structure) == 0:
            return None

        model = structure[0]
        if len(model) == 0:
            return None

        chain = model[0]
        atom_names = ["N", "CA", "C", "O"]
        coords = []

        for residue in chain:
            # Skip non-polymer residues (water, ligands)
            if not residue.name.isalpha() or len(residue.name) > 3:
                continue

            atom_coords = []
            found_all = True
            for aname in atom_names:
                atom = residue.find_atom(aname, "\0")
                if atom is None:
                    found_all = False
                    break
                atom_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])

            if found_all:
                coords.append(atom_coords)

        if len(coords) == 0:
            return None

        return np.array(coords, dtype=np.float32)

    except Exception as e:
        logger.error(f"Failed to extract backbone coords from {cif_path}: {e}")
        return None


def run_dssp(coords: np.ndarray) -> List[str]:
    """
    Run PyDSSP secondary structure assignment on backbone coordinates.

    Args:
        coords: Array of shape (n_residues, 4, 3) with N/CA/C/O positions.

    Returns:
        List of DSSP codes ('H', 'E', 'C') per residue.
    """
    try:
        import pydssp

        # pydssp.assign with out_type='c3' (default) returns character codes:
        # 'H' (helix), 'E' (sheet), 'C' (coil), '-' (disordered/unresolved).
        # Do NOT cast to int — that raised a silent ValueError causing all-coil fallback.
        if coords.ndim == 3:
            coords_batch = coords[np.newaxis, ...]
        else:
            coords_batch = coords

        assignments = pydssp.assign(coords_batch)

        # Extract first batch if batched
        if assignments.ndim > 1:
            result = assignments[0]
        else:
            result = assignments

        return [_normalize_dssp(a) for a in result]

    except Exception as e:
        logger.error(f"PyDSSP assignment failed: {e}")
        n_residues = coords.shape[0] if coords is not None else 0
        return ["C"] * n_residues




def extract_residue_features(cif_path) -> List[ResidueStructure]:
    """
    Extract per-residue structural features from an mmCIF structure.

    Combines:
    1. DSSP secondary structure assignment (via PyDSSP)
    2. Phi/psi backbone torsion angles (via gemmi)
    3. Per-residue pLDDT confidence (from B-factor column)

    Args:
        cif_path: Path to mmCIF file.

    Returns:
        List of ResidueStructure, one per polymer residue.
    """
    # Get backbone coords and run DSSP
    backbone = extract_backbone_coords(cif_path)
    if backbone is None:
        return []

    dssp_codes = run_dssp(backbone)

    # Also parse with gemmi for phi/psi and pLDDT
    try:
        structure = gemmi.read_structure(str(cif_path))
        model = structure[0]
        chain = model[0]
        polymer = chain.get_polymer()
    except Exception as e:
        logger.error(f"Failed to parse structure for phi/psi: {e}")
        polymer = None

    results = []
    residue_idx = 0

    if polymer is not None:
        n_res = len(polymer)
        for i, residue in enumerate(polymer):
            if residue_idx >= len(dssp_codes):
                break

            # Phi/Psi via gemmi.calculate_phi_psi(prev, res, next) → radians;
            # NaN at chain termini / missing backbone atoms.  Stored in DEGREES.
            # (Review STR-3: the previous code called gemmi.calculate_phi /
            # calculate_psi, which do not exist in gemmi ≥ 0.6, inside a bare
            # `except Exception: pass`, so every phi/psi was NULL.)  No try/except
            # here on purpose — a real failure must surface, not become NULL.
            prev_res = polymer[i - 1] if i > 0 else None
            next_res = polymer[i + 1] if i + 1 < n_res else None
            phi_rad, psi_rad = gemmi.calculate_phi_psi(prev_res, residue, next_res)
            phi = None if math.isnan(phi_rad) else round(math.degrees(phi_rad), 2)
            psi = None if math.isnan(psi_rad) else round(math.degrees(psi_rad), 2)

            # pLDDT from B-factor of CA atom
            plddt = None
            ca = residue.find_atom("CA", "\0")
            if ca is not None:
                plddt = round(ca.b_iso, 2)

            results.append(ResidueStructure(
                residue_number=residue_idx + 1,  # 1-based
                dssp=dssp_codes[residue_idx],
                phi=phi,
                psi=psi,
                plddt=plddt,
            ))
            residue_idx += 1
    else:
        # Fallback: no polymer info, just use DSSP codes
        for idx, code in enumerate(dssp_codes):
            results.append(ResidueStructure(
                residue_number=idx + 1,
                dssp=code,
                phi=None,
                psi=None,
                plddt=None,
            ))

    return results


def classify_structural_change(novel_dssp: str, canonical_dssp: str) -> str:
    """
    Classify the structural change between a novel and canonical residue.

    Args:
        novel_dssp: DSSP code for the novel isoform residue.
        canonical_dssp: DSSP code for the canonical isoform residue.

    Returns:
        Classification string describing the change.
    """
    if novel_dssp == canonical_dssp:
        return "same"

    return f"{_dssp_name(canonical_dssp)}_to_{_dssp_name(novel_dssp)}"


def _dssp_name(code: str) -> str:
    """Human-readable name for DSSP code."""
    return {"H": "helix", "E": "sheet", "C": "coil"}.get(code, "unknown")


def compare_with_canonical(
    novel_features: List[ResidueStructure],
    canonical_features: List[ResidueStructure],
    novel_sequence: str,
    canonical_sequence: str,
) -> Dict[int, Tuple[Optional[str], Optional[str]]]:
    """
    Align novel-vs-canonical via the splice-aware aligner and classify
    per-residue secondary-structure changes.

    For each novel residue position, returns a (canonical_dssp,
    structural_change) tuple where structural_change ∈ {
        same, coil_to_helix, helix_to_coil, coil_to_sheet, sheet_to_coil,
        helix_to_sheet, sheet_to_helix, low_confidence, novel_residue
    }.

    `low_confidence` is emitted when EITHER the canonical or novel residue
    pLDDT is below LOW_PLDDT_CUTOFF (= 70) — DSSP assignments at low
    confidence are unreliable.  `novel_residue` is emitted when the novel
    position has no canonical counterpart (the splice-aware aligner placed
    it in a gap of canonical).
    """
    if not novel_sequence or not canonical_sequence:
        return {}

    # Splice-aware alignment (same algorithm as frontend MSA + PTM module);
    # k-mer-anchored chaining keeps novel exon insertions as proper gaps
    # instead of force-aligning them to unrelated canonical residues.
    from app.analysis.splice_aware_align import map_canonical_to_novel

    canonical_to_novel = map_canonical_to_novel(canonical_sequence, novel_sequence)
    if not canonical_to_novel:
        return {}

    # Invert: novel_pos → canonical_pos (for residues with a canonical match)
    novel_to_canonical: Dict[int, int] = {
        nov: can for can, nov in canonical_to_novel.items() if nov is not None
    }

    # Per-residue lookups
    canonical_features_by_pos = {f.residue_number: f for f in canonical_features}
    novel_features_by_pos     = {f.residue_number: f for f in novel_features}

    result: Dict[int, Tuple[Optional[str], Optional[str]]] = {}

    for nov_feat in novel_features:
        nov_pos = nov_feat.residue_number
        can_pos = novel_to_canonical.get(nov_pos)

        if can_pos is None:
            # No canonical counterpart — splicing introduced this residue.
            result[nov_pos] = (None, "novel_residue")
            continue

        can_feat = canonical_features_by_pos.get(can_pos)
        if can_feat is None:
            # Canonical position aligned but structure-side feature missing
            # (rare; backbone-atom dropout).  Treat as low-confidence.
            result[nov_pos] = (None, "low_confidence")
            continue

        # Confidence gate — min pLDDT across the two sites
        min_plddt = None
        if nov_feat.plddt is not None and can_feat.plddt is not None:
            min_plddt = min(nov_feat.plddt, can_feat.plddt)
        elif nov_feat.plddt is not None:
            min_plddt = nov_feat.plddt
        elif can_feat.plddt is not None:
            min_plddt = can_feat.plddt

        if min_plddt is not None and min_plddt < LOW_PLDDT_CUTOFF:
            result[nov_pos] = (can_feat.dssp, "low_confidence")
            continue

        change = classify_structural_change(nov_feat.dssp, can_feat.dssp)
        result[nov_pos] = (can_feat.dssp, change)

    return result


def analyze_isoform(
    db: Session,
    isoform: Isoform,
    force: bool = False,
) -> List[ExonStructuralImpact]:
    """
    Run per-residue structural impact analysis for a single isoform.

    Steps:
        1. Resolve CIF structure file path
        2. Extract backbone coords, run DSSP, get phi/psi/pLDDT
        3. Find canonical isoform + structure, compare DSSP if available
        4. Store per-residue ExonStructuralImpact records

    Args:
        db: SQLAlchemy session.
        isoform: The isoform to analyze.
        force: If True, delete existing results and re-analyze.

    Returns:
        List of ExonStructuralImpact ORM objects.
    """
    # Check for existing results
    existing = db.query(ExonStructuralImpact).filter(
        ExonStructuralImpact.isoform_id == isoform.id
    ).first()
    if existing and not force:
        logger.info(f"Structural impact already exists for {isoform.transcript_id}, skipping")
        return db.query(ExonStructuralImpact).filter(
            ExonStructuralImpact.isoform_id == isoform.id
        ).order_by(ExonStructuralImpact.residue_number).all()

    if force:
        db.query(ExonStructuralImpact).filter(
            ExonStructuralImpact.isoform_id == isoform.id
        ).delete()
        db.flush()

    # Get structure for this isoform
    struct = db.query(Structure).filter(Structure.isoform_id == isoform.id).first()
    if not struct:
        logger.warning(f"No structure found for {isoform.transcript_id}, skipping")
        return []

    cif_path = resolve_structure_path(struct.af3_cif_url)
    if cif_path is None:
        logger.warning(f"CIF file not found for {isoform.transcript_id}: {struct.af3_cif_url}")
        return []

    # Extract features
    novel_features = extract_residue_features(cif_path)
    if not novel_features:
        logger.warning(f"No residue features extracted for {isoform.transcript_id}")
        return []

    # Canonical comparison
    canonical = find_canonical_isoform(db, isoform.gene_id)
    canonical_id: Optional[int] = None
    comparison: Dict[int, Tuple[Optional[str], Optional[str]]] = {}

    if canonical and canonical.id != isoform.id:
        canonical_id = canonical.id
        canonical_struct = db.query(Structure).filter(
            Structure.isoform_id == canonical.id
        ).first()

        if canonical_struct:
            canonical_cif = resolve_structure_path(canonical_struct.af3_cif_url)
            if canonical_cif is not None:
                canonical_features = extract_residue_features(canonical_cif)
                if canonical_features and isoform.sequence and canonical.sequence:
                    comparison = compare_with_canonical(
                        novel_features, canonical_features,
                        isoform.sequence, canonical.sequence,
                    )

    # Store results
    results: List[ExonStructuralImpact] = []
    for feat in novel_features:
        can_dssp, change = comparison.get(feat.residue_number, (None, None))

        # If this IS the canonical and no comparison target, mark as "canonical"
        if canonical and canonical.id == isoform.id:
            change = None
            can_dssp = None
            canonical_id = None

        record = ExonStructuralImpact(
            isoform_id=isoform.id,
            residue_number=feat.residue_number,
            dssp_assignment=feat.dssp,
            phi_angle=feat.phi,
            psi_angle=feat.psi,
            residue_plddt=feat.plddt,
            exon_number=None,  # populated when exon boundary data is available
            canonical_dssp=can_dssp,
            structural_change=change,
            # Keep the comparison target on every row — including novel residues
            # with no canonical counterpart — so "compared against X, no match"
            # is distinguishable from "no canonical available" (review STR-6).
            canonical_isoform_id=canonical_id,
        )
        results.append(record)

    # Bulk insert for performance (can be hundreds of residues)
    db.bulk_save_objects(results)
    db.flush()

    # Update/Create StructuralMetrics entry for flattened retrieval
    metrics = isoform.structural_metrics
    if not metrics:
        metrics = StructuralMetrics(isoform_id=isoform.id)
        db.add(metrics)
        db.flush()

    plddt_series = [r.residue_plddt for r in results]
    metrics.plddt_series = plddt_series
    
    if plddt_series:
        p_arr = np.array([p for p in plddt_series if p is not None])
        if p_arr.size > 0:
            metrics.plddt_mean = float(np.mean(p_arr))
            metrics.plddt_min = float(np.min(p_arr))
            metrics.plddt_median = float(np.median(p_arr))

    db.commit()

    # Log summary
    if results:
        dssp_counts = {}
        change_counts = {}
        for r in results:
            dssp_counts[r.dssp_assignment] = dssp_counts.get(r.dssp_assignment, 0) + 1
            if r.structural_change:
                change_counts[r.structural_change] = change_counts.get(r.structural_change, 0) + 1

        logger.info(
            f"Structural impact for {isoform.transcript_id}: "
            f"{len(results)} residues, DSSP={dssp_counts}, changes={change_counts}"
        )

    return results
