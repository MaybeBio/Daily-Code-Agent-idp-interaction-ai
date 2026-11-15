from __future__ import annotations
from rdkit import Chem
from openmm.app import Modeller
from openmm.app.element import hydrogen as ELEMENT_H
from openmm import Vec3
from ..utils.substrate_utils import is_valid_mol_3d, is_valid_mol_h
from typing import Tuple, List, Any, Dict, Set
import numpy as np
from ..resources.aa_resources import AA3_STANDARD, HBOND_SIDE_ACCEPTORS, HBOND_SIDE_DONOR_HEAVY, PROTEIN_IONIC_RESIDUES, VDW_RADIUS_A, PROTEIN_PIPI_AROMATIC_RESIDUES, PROTEIN_PIPI_RING_ATOMS, DISULFIDE_RESNAME, DISULFIDE_ATOM_NAME
from ..resources.aa_resources import PROTEIN_PICATION_LYS_ATOM, PROTEIN_PICATION_ARG_CENTER_ATOMS, PROTEIN_PICATION_ARG_PLANE_ATOMS, PROTEIN_PICATION_CATION_RESIDUES
from scipy.spatial import cKDTree
from functools import lru_cache
from rdkit.Chem import ChemicalFeatures
from ..utils.logging_utils import Logger
from rdkit import RDConfig
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

'''
helper
'''




def filter_valid_docked_substrates(
    substrate_name_list: List[str],
    ligand_mol_list: List[Chem.Mol],
    modeller: Modeller,
    logger: Logger,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
) -> Tuple[List[str], List[Chem.Mol]] | None:

    if not isinstance(substrate_name_list, list) or not isinstance(ligand_mol_list, list):
        logger.print("[ERROR] substrate_name_list and ligand_mol_list must be lists.")
        return None

    if len(substrate_name_list) != len(ligand_mol_list):
        logger.print("[ERROR] substrate_name_list and ligand_mol_list must have the same length.")
        return None

    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(docked_heavy_atom_distance_cutoff_A, (int, float)) or float(docked_heavy_atom_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] docked_heavy_atom_distance_cutoff_A must be a positive number.")
        return None

    valid_substrate_name_list: List[str] = []
    valid_ligand_mol_list: List[Chem.Mol] = []

    for substrate_name, lig_mol in zip(substrate_name_list, ligand_mol_list):

        # ---- check mol 3D ----
        if not is_valid_mol_3d(lig_mol, logger):
            logger.print(f"[WARNING] Invalid Mol(3D) for substrate '{substrate_name}'. Skipped. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")
            continue

        # ---- check hydrogen ----
        if not is_valid_mol_h(lig_mol, logger):
            logger.print(f"[WARNING] Substrate '{substrate_name}' does not contain valid explicit hydrogen atoms. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")


        # ---- check docking ----
        if not is_protein_substrate_docked(
            modeller=modeller,
            ligand_mol=lig_mol,
            logger=logger,
            heavy_atom_distance_cutoff_A=float(docked_heavy_atom_distance_cutoff_A),
        ):
            logger.print(f"[WARNING] Substrate '{substrate_name}' is not spatially docked to protein. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")



        valid_substrate_name_list.append(substrate_name)
        valid_ligand_mol_list.append(lig_mol)

    return valid_substrate_name_list, valid_ligand_mol_list

def to_angstrom(v: Vec3) -> np.ndarray:
    return np.array([v.x, v.y, v.z], dtype=float) * 10.0

def safe_int(s: str) -> Any:
    try:
        return int(s)
    except Exception:
        return s

def angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    v1 = a - b
    v2 = c - b
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)

    if n1 == 0.0 or n2 == 0.0:
        return 0.0

    cosang = float(np.dot(v1, v2) / (n1 * n2))
    cosang = max(-1.0, min(1.0, cosang))
    return float(np.degrees(np.arccos(cosang)))


def is_protein_residue(resname: str) -> bool:
    return isinstance(resname, str) and resname in AA3_STANDARD


def get_atom_in_residue(res, atom_name: str):
    for atom in res.atoms():
        if atom.name == atom_name:
            return atom
    return None


def edge_key_pair(a: Any, b: Any) -> Tuple[Any, Any]:
    if isinstance(a, int) and isinstance(b, int):
        return (a, b) if a < b else (b, a)

    sa = str(a)
    sb = str(b)
    return (a, b) if sa < sb else (b, a)


def edge_sort_key(edge: Dict[str, Any]) -> Tuple[int, str, str]:
    return (0, str(edge.get("node1_index")), str(edge.get("node2_index")))

def build_openmm_tables(modeller: Modeller) -> Tuple[Any, List[Vec3], List[Any], Dict[Any, int], np.ndarray, Dict[int, List[int]]]:
    topology = modeller.topology
    positions_nm = list(modeller.positions)
    atoms = list(topology.atoms())
    atom_index: Dict[Any, int] = {atom: i for i, atom in enumerate(atoms)}
    coords_A = np.vstack([to_angstrom(positions_nm[i]) for i in range(len(atoms))])

    bonded: Dict[int, List[int]] = {i: [] for i in range(len(atoms))}
    for bond in topology.bonds():
        i = atom_index[bond[0]]
        j = atom_index[bond[1]]
        bonded[i].append(j)
        bonded[j].append(i)

    return topology, positions_nm, atoms, atom_index, coords_A, bonded


def build_substrate_tables(mol: Chem.Mol) -> Tuple[np.ndarray, Dict[int, List[int]]]:
    conf = mol.GetConformer()
    coords_A = np.array(
        [
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ],
        dtype=float,
    )

    bonded: Dict[int, List[int]] = {i: [] for i in range(mol.GetNumAtoms())}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bonded[i].append(j)
        bonded[j].append(i)

    return coords_A, bonded

@lru_cache(maxsize=1)
def _get_feature_factory():
    fdef = f"{RDConfig.RDDataDir}/BaseFeatures.fdef"
    return ChemicalFeatures.BuildFeatureFactory(fdef)


def unit_vector(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n == 0.0:
        return v
    return v / n


def angle_deg_between_vectors(u: np.ndarray, v: np.ndarray) -> float:
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu == 0.0 or nv == 0.0:
        return 0.0

    cosang = float(np.dot(u, v) / (nu * nv))
    cosang = max(-1.0, min(1.0, cosang))
    return float(np.degrees(np.arccos(cosang)))


def min_angle_0_90(deg: float) -> float:
    deg = abs(float(deg))
    if deg > 90.0:
        deg = 180.0 - deg
    return deg


def fit_plane_normal(points: np.ndarray) -> np.ndarray:
    if not isinstance(points, np.ndarray):
        return np.zeros(3, dtype=float)

    if points.ndim != 2 or points.shape[0] < 3 or points.shape[1] != 3:
        return np.zeros(3, dtype=float)

    center = points.mean(axis=0)
    x = points - center
    cov = x.T @ x
    w, v = np.linalg.eigh(cov)
    normal = v[:, int(np.argmin(w))]
    return unit_vector(normal)

'''
'''

def is_protein_substrate_docked(modeller: Modeller,ligand_mol: Chem.Mol, logger,heavy_atom_distance_cutoff_A: float = 6.5) -> bool:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return False

    if not is_valid_mol_3d(ligand_mol, logger):
        return False

    if not isinstance(heavy_atom_distance_cutoff_A, (int, float)) or float(heavy_atom_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] heavy_atom_distance_cutoff_A must be a positive number.")
        return False

    try:
        topology, positions_nm, atoms, atom_index, coords_A, _ = build_openmm_tables(modeller)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build protein tables for docking check: {e}")
        return False

    try:
        lig_coords_A, _ = build_substrate_tables(ligand_mol)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build ligand tables for docking check: {e}")
        return False

    prot_idx: List[int] = []
    for i, atom in enumerate(atoms):
        if not is_protein_residue(atom.residue.name):
            continue
        if atom.element is not None and atom.element == ELEMENT_H:
            continue
        prot_idx.append(i)

    lig_idx: List[int] = []
    for i, atom in enumerate(ligand_mol.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue
        lig_idx.append(i)

    if len(prot_idx) == 0:
        logger.print("[ERROR] No protein heavy atoms found in modeller.")
        return False

    if len(lig_idx) == 0:
        logger.print("[ERROR] No ligand heavy atoms found in mol.")
        return False

    prot_xyz = coords_A[prot_idx]
    lig_xyz = lig_coords_A[lig_idx]

    try:
        kdtree = cKDTree(lig_xyz)
        for p in prot_xyz:
            cand = kdtree.query_ball_point(p, r=float(heavy_atom_distance_cutoff_A))
            if len(cand) > 0:
                return True
        return False
    except Exception:
        for p in prot_xyz:
            d = np.linalg.norm(lig_xyz - p[None, :], axis=1)
            if float(d.min()) <= float(heavy_atom_distance_cutoff_A):
                return True
        return False

'''
H-bond
'''

def protein_bonded_hydrogens(atom_obj,atom_index: Dict[Any, int],atoms: List[Any],coords_A: np.ndarray,bonded: Dict[int, List[int]],bonded_h_min_distance_A: float,bonded_h_max_distance_A: float) -> List[Tuple[str, np.ndarray]]:
    if atom_obj is None:
        return []

    i = atom_index[atom_obj]
    p_i = coords_A[i]
    out: List[Tuple[str, np.ndarray]] = []

    for j in bonded[i]:
        atom_j = atoms[j]
        if atom_j.element is None or atom_j.element != ELEMENT_H:
            continue

        d = float(np.linalg.norm(coords_A[j] - p_i))
        if float(bonded_h_min_distance_A) <= d <= float(bonded_h_max_distance_A):
            out.append((atom_j.name, coords_A[j]))

    return out


def ligand_bonded_hydrogens(mol: Chem.Mol,coords_A: np.ndarray,bonded: Dict[int, List[int]],atom_idx: int,bonded_h_min_distance_A: float,bonded_h_max_distance_A: float) -> List[Tuple[str, np.ndarray]]:
    out: List[Tuple[str, np.ndarray]] = []
    p_i = coords_A[atom_idx]

    for j in bonded[atom_idx]:
        atom_j = mol.GetAtomWithIdx(j)
        if atom_j.GetAtomicNum() != 1:
            continue

        d = float(np.linalg.norm(coords_A[j] - p_i))
        if float(bonded_h_min_distance_A) <= d <= float(bonded_h_max_distance_A):
            out.append((atom_j.GetSymbol(), coords_A[j]))

    return out


def collect_protein_hbond_sites(topology,positions_nm: List[Vec3],atoms: List[Any],atom_index: Dict[Any, int],coords_A: np.ndarray,bonded: Dict[int, List[int]],bonded_h_min_distance_A: float,bonded_h_max_distance_A: float) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    acceptor_atoms: List[Dict[str, Any]] = []
    donor_entries: List[Dict[str, Any]] = []

    for res in topology.residues():
        resname = res.name
        if not is_protein_residue(resname):
            continue

        resid = safe_int(res.id)

        o_atom = get_atom_in_residue(res, "O")
        if o_atom is not None:
            o_xyz = to_angstrom(positions_nm[atom_index[o_atom]])
            acceptor_atoms.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "atom": "O",
                    "xyz": o_xyz,
                }
            )

        if resname != "PRO":
            n_atom = get_atom_in_residue(res, "N")
            if n_atom is not None:
                h_list = protein_bonded_hydrogens(
                    atom_obj=n_atom,
                    atom_index=atom_index,
                    atoms=atoms,
                    coords_A=coords_A,
                    bonded=bonded,
                    bonded_h_min_distance_A=bonded_h_min_distance_A,
                    bonded_h_max_distance_A=bonded_h_max_distance_A,
                )
                if len(h_list) > 0:
                    n_xyz = to_angstrom(positions_nm[atom_index[n_atom]])
                    donor_entries.append(
                        {
                            "res_id": resid,
                            "res_name": resname,
                            "atom": "N",
                            "xyz": n_xyz,
                            "H": [{"atom": hn, "xyz": hxyz} for hn, hxyz in h_list],
                        }
                    )

        for atom_name in HBOND_SIDE_ACCEPTORS.get(resname, []):
            atom_obj = get_atom_in_residue(res, atom_name)
            if atom_obj is None:
                continue

            a_xyz = to_angstrom(positions_nm[atom_index[atom_obj]])
            acceptor_atoms.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "atom": atom_name,
                    "xyz": a_xyz,
                }
            )

        for atom_name in HBOND_SIDE_DONOR_HEAVY.get(resname, []):
            atom_obj = get_atom_in_residue(res, atom_name)
            if atom_obj is None:
                continue

            h_list = protein_bonded_hydrogens(
                atom_obj=atom_obj,
                atom_index=atom_index,
                atoms=atoms,
                coords_A=coords_A,
                bonded=bonded,
                bonded_h_min_distance_A=bonded_h_min_distance_A,
                bonded_h_max_distance_A=bonded_h_max_distance_A,
            )
            if len(h_list) == 0:
                continue

            d_xyz = to_angstrom(positions_nm[atom_index[atom_obj]])
            donor_entries.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "atom": atom_name,
                    "xyz": d_xyz,
                    "H": [{"atom": hn, "xyz": hxyz} for hn, hxyz in h_list],
                }
            )

    return acceptor_atoms, donor_entries

def collect_substrate_hbond_sites(mol: Chem.Mol,coords_A: np.ndarray,bonded: Dict[int, List[int]],bonded_h_min_distance_A: float,bonded_h_max_distance_A: float) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]] | None:
    try:
        factory = _get_feature_factory()
        feats = factory.GetFeaturesForMol(mol)
    except Exception:
        return None

    acceptor_idx: Set[int] = set()
    donor_idx: Set[int] = set()

    for feat in feats:
        family = feat.GetFamily()
        atom_ids = feat.GetAtomIds()

        if family == "Acceptor":
            acceptor_idx.update(atom_ids)
        elif family == "Donor":
            donor_idx.update(atom_ids)

    acceptor_atoms: List[Dict[str, Any]] = []
    donor_entries: List[Dict[str, Any]] = []

    for i in sorted(acceptor_idx):
        atom = mol.GetAtomWithIdx(i)
        acceptor_atoms.append(
            {
                "atom_idx": i,
                "atom": atom.GetSymbol(),
                "xyz": coords_A[i],
            }
        )

    for i in sorted(donor_idx):
        atom = mol.GetAtomWithIdx(i)

        if atom.GetAtomicNum() == 1:
            continue

        h_list = ligand_bonded_hydrogens(
            mol=mol,
            coords_A=coords_A,
            bonded=bonded,
            atom_idx=i,
            bonded_h_min_distance_A=bonded_h_min_distance_A,
            bonded_h_max_distance_A=bonded_h_max_distance_A,
        )
        if len(h_list) == 0:
            continue

        donor_entries.append(
            {
                "atom_idx": i,
                "atom": atom.GetSymbol(),
                "xyz": coords_A[i],
                "H": [{"atom": hn, "xyz": hxyz} for hn, hxyz in h_list],
            }
        )

    return acceptor_atoms, donor_entries


def find_hbond_hits_from_donors_to_acceptors(donor_entries: List[Dict[str, Any]],acceptor_atoms: List[Dict[str, Any]],da_max_distance_A: float,ha_max_distance_A: float,dha_min_angle_deg: float) -> List[Tuple[Dict[str, Any], Dict[str, Any]]]:
    if len(donor_entries) == 0 or len(acceptor_atoms) == 0:
        return []

    acc_xyz = np.vstack([item["xyz"] for item in acceptor_atoms])
    hits: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []

    use_kdtree = False
    kdtree = None
    try:
        from scipy.spatial import cKDTree  # type: ignore
        kdtree = cKDTree(acc_xyz)
        use_kdtree = True
    except Exception:
        use_kdtree = False

    for donor in donor_entries:
        d_xyz = donor["xyz"]

        if use_kdtree and kdtree is not None:
            cand_idx = kdtree.query_ball_point(d_xyz, r=float(da_max_distance_A))
        else:
            diff = acc_xyz - d_xyz[None, :]
            cand_idx = np.where(np.linalg.norm(diff, axis=1) < float(da_max_distance_A))[0].tolist()

        if len(cand_idx) == 0:
            continue

        for j in cand_idx:
            acceptor = acceptor_atoms[j]
            a_xyz = acceptor["xyz"]

            if float(np.linalg.norm(d_xyz - a_xyz)) >= float(da_max_distance_A):
                continue

            ok = False
            for h in donor["H"]:
                h_xyz = h["xyz"]

                if float(np.linalg.norm(h_xyz - a_xyz)) >= float(ha_max_distance_A):
                    continue

                angle = angle_deg(d_xyz, h_xyz, a_xyz)
                if angle > float(dha_min_angle_deg):
                    ok = True
                    break

            if ok:
                hits.append((donor, acceptor))

    return hits

'''
'''

'''
IONIC-Bond
'''

def is_protonated_his(res) -> bool:
    """
    Same rule as original code:
    treat HIS as cation only when both ND1 and NE2 are protonated
    (HD1 and HE2 both present).
    """
    if res is None or res.name != "HIS":
        return False

    has_hd1 = get_atom_in_residue(res, "HD1") is not None
    has_he2 = get_atom_in_residue(res, "HE2") is not None
    return has_hd1 and has_he2

def collect_protein_ionic_centers(topology,positions_nm: List[Vec3],atom_index: Dict[Any, int]) -> List[Dict[str, Any]]:
    charged_centers: List[Dict[str, Any]] = []

    for res in topology.residues():
        resname = res.name
        if not is_protein_residue(resname):
            continue
        if resname not in PROTEIN_IONIC_RESIDUES:
            continue

        resid = safe_int(res.id)

        if resname == "ASP":
            od1 = get_atom_in_residue(res, "OD1")
            od2 = get_atom_in_residue(res, "OD2")
            if od1 is None or od2 is None:
                continue

            p1 = to_angstrom(positions_nm[atom_index[od1]])
            p2 = to_angstrom(positions_nm[atom_index[od2]])
            xyz = 0.5 * (p1 + p2)

            charged_centers.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "charge": "anion",
                    "xyz": xyz,
                    "center_id": "ASP_CARBOXYLATE",
                }
            )

        elif resname == "GLU":
            oe1 = get_atom_in_residue(res, "OE1")
            oe2 = get_atom_in_residue(res, "OE2")
            if oe1 is None or oe2 is None:
                continue

            p1 = to_angstrom(positions_nm[atom_index[oe1]])
            p2 = to_angstrom(positions_nm[atom_index[oe2]])
            xyz = 0.5 * (p1 + p2)

            charged_centers.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "charge": "anion",
                    "xyz": xyz,
                    "center_id": "GLU_CARBOXYLATE",
                }
            )

        elif resname == "LYS":
            nz = get_atom_in_residue(res, "NZ")
            if nz is None:
                continue

            xyz = to_angstrom(positions_nm[atom_index[nz]])
            charged_centers.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "charge": "cation",
                    "xyz": xyz,
                    "center_id": "LYS_NZ",
                }
            )

        elif resname == "ARG":
            ne = get_atom_in_residue(res, "NE")
            nh1 = get_atom_in_residue(res, "NH1")
            nh2 = get_atom_in_residue(res, "NH2")
            heavy = [a for a in (ne, nh1, nh2) if a is not None]
            if len(heavy) == 0:
                continue

            xyzs = np.vstack([to_angstrom(positions_nm[atom_index[a]]) for a in heavy])
            xyz = xyzs.mean(axis=0)

            charged_centers.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "charge": "cation",
                    "xyz": xyz,
                    "center_id": "ARG_GUANIDINIUM",
                }
            )

        elif resname == "HIS":
            if not is_protonated_his(res):
                continue

            nd1 = get_atom_in_residue(res, "ND1")
            ne2 = get_atom_in_residue(res, "NE2")
            if nd1 is None or ne2 is None:
                continue

            p1 = to_angstrom(positions_nm[atom_index[nd1]])
            p2 = to_angstrom(positions_nm[atom_index[ne2]])
            xyz = 0.5 * (p1 + p2)

            charged_centers.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "charge": "cation",
                    "xyz": xyz,
                    "center_id": "HIS_IMIDAZOLIUM",
                }
            )

    return charged_centers

def collect_substrate_ionic_centers(mol: Chem.Mol,coords_A: np.ndarray) -> List[Dict[str, Any]]:
    charged_centers: List[Dict[str, Any]] = []

    for i, atom in enumerate(mol.GetAtoms()):
        formal_charge = int(atom.GetFormalCharge())

        if formal_charge > 0:
            charged_centers.append(
                {
                    "atom_idx": int(i),
                    "atom": atom.GetSymbol(),
                    "charge": "cation",
                    "xyz": coords_A[i],
                    "center_id": f"ATOM_{i}",
                }
            )
        elif formal_charge < 0:
            charged_centers.append(
                {
                    "atom_idx": int(i),
                    "atom": atom.GetSymbol(),
                    "charge": "anion",
                    "xyz": coords_A[i],
                    "center_id": f"ATOM_{i}",
                }
            )

    return charged_centers


def find_ionic_hits_between_centers(center_entries_1: List[Dict[str, Any]],center_entries_2: List[Dict[str, Any]],distance_cutoff_A: float) -> List[Tuple[Dict[str, Any], Dict[str, Any]]]:
    if len(center_entries_1) == 0 or len(center_entries_2) == 0:
        return []

    xyz2 = np.vstack([item["xyz"] for item in center_entries_2])
    hits: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []

    use_kdtree = False
    kdtree = None
    try:

        kdtree = cKDTree(xyz2)
        use_kdtree = True
    except Exception:
        use_kdtree = False

    for center1 in center_entries_1:
        xyz1 = center1["xyz"]

        if use_kdtree and kdtree is not None:
            cand_idx = kdtree.query_ball_point(xyz1, r=float(distance_cutoff_A))
        else:
            diff = xyz2 - xyz1[None, :]
            cand_idx = np.where(np.linalg.norm(diff, axis=1) <= float(distance_cutoff_A))[0].tolist()

        if len(cand_idx) == 0:
            continue

        for j in cand_idx:
            center2 = center_entries_2[j]
            if float(np.linalg.norm(xyz1 - center2["xyz"])) <= float(distance_cutoff_A):
                hits.append((center1, center2))

    return hits

'''
'''

'''
VDW
'''

def get_openmm_atom_element_symbol(atom_obj, vdw_radius_map: Dict[str, float]) -> str:
    if atom_obj is None:
        return ""

    try:
        if getattr(atom_obj, "element", None) is not None and atom_obj.element is not None:
            sym = atom_obj.element.symbol
            if isinstance(sym, str) and sym.strip():
                return sym.strip().upper()
    except Exception:
        pass

    try:
        atom_name = str(getattr(atom_obj, "name", "")).strip()
        if not atom_name:
            return ""

        head2 = atom_name[:2].upper()
        head1 = atom_name[:1].upper()

        if head2 in vdw_radius_map:
            return head2
        if head1 in vdw_radius_map:
            return head1
        return ""
    except Exception:
        return ""

def get_rdkit_atom_element_symbol(atom_obj) -> str:
    if atom_obj is None:
        return ""

    try:
        sym = atom_obj.GetSymbol()
        if isinstance(sym, str) and sym.strip():
            return sym.strip().upper()
        return ""
    except Exception:
        return ""

def collect_protein_vdw_atoms(topology,atoms: List[Any],atom_index: Dict[Any, int],coords_A: np.ndarray,vdw_radius_map: Dict[str, float]=VDW_RADIUS_A) -> List[Dict[str, Any]]:
    atom_entries: List[Dict[str, Any]] = []

    for res in topology.residues():
        resname = res.name
        if not is_protein_residue(resname):
            continue

        try:
            res_id = int(res.id)
        except Exception:
            res_id = res.id

        for atom_obj in res.atoms():
            elem_symbol = get_openmm_atom_element_symbol(atom_obj, vdw_radius_map)
            if not elem_symbol:
                continue

            radius_A = vdw_radius_map.get(elem_symbol)
            if radius_A is None:
                continue

            idx = atom_index.get(atom_obj)
            if idx is None:
                continue

            atom_entries.append(
                {
                    "res_id": res_id,
                    "res_name": resname,
                    "atom": str(atom_obj.name),
                    "element": elem_symbol,
                    "radius_A": float(radius_A),
                    "xyz": coords_A[idx],
                }
            )

    return atom_entries

def collect_substrate_vdw_atoms(mol: Chem.Mol,coords_A: np.ndarray,vdw_radius_map: Dict[str, float]=VDW_RADIUS_A) -> List[Dict[str, Any]]:
    atom_entries: List[Dict[str, Any]] = []

    for i, atom_obj in enumerate(mol.GetAtoms()):
        elem_symbol = get_rdkit_atom_element_symbol(atom_obj)
        if not elem_symbol:
            continue

        radius_A = vdw_radius_map.get(elem_symbol)
        if radius_A is None:
            continue

        atom_entries.append(
            {
                "atom_idx": int(i),
                "atom": elem_symbol,
                "element": elem_symbol,
                "radius_A": float(radius_A),
                "xyz": coords_A[i],
            }
        )

    return atom_entries

def find_vdw_hits_between_atom_entries(atom_entries_1: List[Dict[str, Any]],atom_entries_2: List[Dict[str, Any]],mu: float,allow_same_entry_list: bool = False) -> List[Tuple[Dict[str, Any], Dict[str, Any]]]:
    if len(atom_entries_1) == 0 or len(atom_entries_2) == 0:
        return []

    xyz2 = np.vstack([item["xyz"] for item in atom_entries_2])
    radius2 = np.array([float(item["radius_A"]) for item in atom_entries_2], dtype=float)

    rmax_1 = max(float(item["radius_A"]) for item in atom_entries_1)
    rmax_2 = float(radius2.max())
    d_max_global = (1.0 + float(mu)) * (float(rmax_1) + float(rmax_2))

    hits: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []

    use_kdtree = False
    kdtree = None
    try:

        kdtree = cKDTree(xyz2)
        use_kdtree = True
    except Exception:
        use_kdtree = False

    for i, item1 in enumerate(atom_entries_1):
        xyz1 = item1["xyz"]
        r1 = float(item1["radius_A"])

        if use_kdtree and kdtree is not None:
            cand_idx = kdtree.query_ball_point(xyz1, r=float(d_max_global))
        else:
            diff = xyz2 - xyz1[None, :]
            cand_idx = np.where(np.linalg.norm(diff, axis=1) <= float(d_max_global))[0].tolist()

        if len(cand_idx) == 0:
            continue

        for j in cand_idx:
            if allow_same_entry_list and atom_entries_1 is atom_entries_2 and j <= i:
                continue

            item2 = atom_entries_2[j]
            rsum = r1 + float(radius2[j])
            low = (1.0 - float(mu)) * rsum
            high = (1.0 + float(mu)) * rsum

            d = float(np.linalg.norm(xyz1 - item2["xyz"]))
            if low <= d <= high:
                hits.append((item1, item2))

    return hits


'''
'''

'''
PIPISTACK
'''

def classify_pipi_geometry(theta_deg: float, delta_deg: float, gamma_deg: float) -> str:
    if float(gamma_deg) > 50.0:
        return "T"
    if float(gamma_deg) < 30.0:
        if float(theta_deg) > 80.0 or float(delta_deg) > 80.0:
            return "S"
        return "P"
    return "I"


def collect_protein_ring_coords(res,positions_nm: List[Any],atom_index: Dict[Any, int],atom_names: List[str]) -> np.ndarray | None:
    pts: List[np.ndarray] = []

    for atom_name in atom_names:
        atom_obj = get_atom_in_residue(res, atom_name)
        if atom_obj is None:
            return None

        idx = atom_index.get(atom_obj)
        if idx is None:
            return None

        pts.append(to_angstrom(positions_nm[idx]))

    if len(pts) == 0:
        return None

    try:
        return np.vstack(pts)
    except Exception:
        return None


def collect_protein_pipi_rings(topology,positions_nm: List[Any],atom_index: Dict[Any, int]) -> List[Dict[str, Any]]:
    ring_entries: List[Dict[str, Any]] = []

    for res in topology.residues():
        resname = res.name

        if not is_protein_residue(resname):
            continue

        if resname not in PROTEIN_PIPI_AROMATIC_RESIDUES:
            continue

        resid = safe_int(res.id)
        ring_def_map = PROTEIN_PIPI_RING_ATOMS.get(resname, {})

        for ring_id, atom_names in ring_def_map.items():
            pts = collect_protein_ring_coords(
                res=res,
                positions_nm=positions_nm,
                atom_index=atom_index,
                atom_names=atom_names,
            )
            if pts is None:
                continue

            center = pts.mean(axis=0)
            normal = fit_plane_normal(pts)

            ring_entries.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "ring_id": ring_id,
                    "center": center,
                    "normal": normal,
                }
            )

    return ring_entries


def collect_substrate_pipi_rings(mol: Chem.Mol,coords_A: np.ndarray) -> List[Dict[str, Any]]:
    ring_entries: List[Dict[str, Any]] = []

    if mol is None:
        return ring_entries

    try:
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
    except Exception:
        return ring_entries

    for ring_idx, atom_ids in enumerate(atom_rings):
        if not isinstance(atom_ids, tuple):
            atom_ids = tuple(atom_ids)

        if len(atom_ids) < 5:
            continue

        try:
            atom_objs = [mol.GetAtomWithIdx(int(i)) for i in atom_ids]
        except Exception:
            continue

        if not all(atom.GetIsAromatic() for atom in atom_objs):
            continue

        try:
            pts = np.vstack([coords_A[int(i)] for i in atom_ids])
        except Exception:
            continue

        center = pts.mean(axis=0)
        normal = fit_plane_normal(pts)

        ring_entries.append(
            {
                "ring_idx": int(ring_idx),
                "atom_ids": tuple(int(i) for i in atom_ids),
                "center": center,
                "normal": normal,
            }
        )

    return ring_entries


def find_pipi_hits_between_ring_entries(ring_entries_1: List[Dict[str, Any]],ring_entries_2: List[Dict[str, Any]],ring_center_distance_cutoff_A: float,allow_same_entry_list: bool = False) -> List[Tuple[Dict[str, Any], Dict[str, Any]]]:
    if len(ring_entries_1) == 0 or len(ring_entries_2) == 0:
        return []

    xyz2 = np.vstack([item["center"] for item in ring_entries_2])
    hits: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []

    use_kdtree = False
    kdtree = None
    try:
        kdtree = cKDTree(xyz2)
        use_kdtree = True
    except Exception:
        use_kdtree = False

    for i, ring1 in enumerate(ring_entries_1):
        center1 = ring1["center"]

        if use_kdtree and kdtree is not None:
            cand_idx = kdtree.query_ball_point(center1, r=float(ring_center_distance_cutoff_A))
        else:
            diff = xyz2 - center1[None, :]
            cand_idx = np.where(np.linalg.norm(diff, axis=1) <= float(ring_center_distance_cutoff_A))[0].tolist()

        if len(cand_idx) == 0:
            continue

        for j in cand_idx:
            if allow_same_entry_list and ring_entries_1 is ring_entries_2 and j <= i:
                continue

            ring2 = ring_entries_2[j]
            dcent = float(np.linalg.norm(center1 - ring2["center"]))
            if dcent > float(ring_center_distance_cutoff_A):
                continue

            hits.append((ring1, ring2))

    return hits

'''
'''

'''
PICATION
'''

def classify_pication_arg_geometry(gamma_deg: float) -> str:
    gamma = float(gamma_deg)

    if (0.0 <= gamma < 30.0) or (150.0 < gamma <= 180.0):
        return "S"

    if 60.0 <= gamma <= 120.0:
        return "T"

    return "I"


def collect_protein_pication_centers(topology,positions_nm: List[Any],atom_index: Dict[Any, int]) -> List[Dict[str, Any]]:
    cation_entries: List[Dict[str, Any]] = []

    for res in topology.residues():
        resname = res.name

        if not is_protein_residue(resname):
            continue

        if resname not in PROTEIN_PICATION_CATION_RESIDUES:
            continue

        resid = safe_int(res.id)

        if resname == "LYS":
            nz = get_atom_in_residue(res, PROTEIN_PICATION_LYS_ATOM)
            if nz is None:
                continue

            center = to_angstrom(positions_nm[atom_index[nz]])
            cation_entries.append(
                {
                    "res_id": resid,
                    "res_name": resname,
                    "center_id": "LYS_NZ",
                    "center": center,
                }
            )

        elif resname == "ARG":
            center_atoms = []
            for atom_name in PROTEIN_PICATION_ARG_CENTER_ATOMS:
                atom_obj = get_atom_in_residue(res, atom_name)
                if atom_obj is not None:
                    center_atoms.append(atom_obj)

            if len(center_atoms) < 3:
                continue

            center_pts = np.vstack(
                [to_angstrom(positions_nm[atom_index[a]]) for a in center_atoms]
            )
            center = center_pts.mean(axis=0)

            plane_pts_list: List[np.ndarray] = []
            for atom_name in PROTEIN_PICATION_ARG_PLANE_ATOMS:
                atom_obj = get_atom_in_residue(res, atom_name)
                if atom_obj is None:
                    continue
                plane_pts_list.append(to_angstrom(positions_nm[atom_index[atom_obj]]))

            guanidinium_normal = None
            if len(plane_pts_list) >= 3:
                try:
                    plane_pts = np.vstack(plane_pts_list)
                    guanidinium_normal = fit_plane_normal(plane_pts)
                except Exception:
                    guanidinium_normal = None

            out_item = {
                "res_id": resid,
                "res_name": resname,
                "center_id": "ARG_GUANIDINIUM",
                "center": center,
            }

            if isinstance(guanidinium_normal, np.ndarray) and guanidinium_normal.shape == (3,):
                out_item["guan_norm"] = guanidinium_normal

            cation_entries.append(out_item)

    return cation_entries


def collect_substrate_pication_centers(mol: Chem.Mol,coords_A: np.ndarray) -> List[Dict[str, Any]]:
    cation_entries: List[Dict[str, Any]] = []

    if mol is None:
        return cation_entries

    for i, atom in enumerate(mol.GetAtoms()):
        if int(atom.GetFormalCharge()) <= 0:
            continue

        cation_entries.append(
            {
                "atom_idx": int(i),
                "atom": atom.GetSymbol(),
                "center_id": f"ATOM_{i}",
                "center": coords_A[i],
            }
        )

    return cation_entries


def find_pication_hits_between_ring_entries_and_cation_entries(
    ring_entries: List[Dict[str, Any]],
    cation_entries: List[Dict[str, Any]],
    ring_cation_distance_cutoff_A: float,
    ring_cation_angle_cutoff_deg: float,
) -> List[Tuple[Dict[str, Any], Dict[str, Any]]]:
    if len(ring_entries) == 0 or len(cation_entries) == 0:
        return []

    center2 = np.vstack([item["center"] for item in cation_entries])
    hits: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []

    use_kdtree = False
    kdtree = None
    try:
        kdtree = cKDTree(center2)
        use_kdtree = True
    except Exception:
        use_kdtree = False

    for ring in ring_entries:
        rc = ring["center"]
        rn = ring["normal"]

        if use_kdtree and kdtree is not None:
            cand_idx = kdtree.query_ball_point(rc, r=float(ring_cation_distance_cutoff_A))
        else:
            diff = center2 - rc[None, :]
            cand_idx = np.where(
                np.linalg.norm(diff, axis=1) <= float(ring_cation_distance_cutoff_A)
            )[0].tolist()

        if len(cand_idx) == 0:
            continue

        for j in cand_idx:
            cat = cation_entries[j]
            cc = cat["center"]

            dcc = float(np.linalg.norm(cc - rc))
            if dcc > float(ring_cation_distance_cutoff_A):
                continue

            v = cc - rc
            alpha = min_angle_0_90(angle_deg_between_vectors(rn, v))
            if alpha > float(ring_cation_angle_cutoff_deg):
                continue

            hits.append((ring, cat))

    return hits

'''
'''

'''
SSBOND
'''

def collect_protein_disulfide_sites(topology,positions_nm: List[Vec3],atom_index: Dict[Any, int]) -> List[Dict[str, Any]]:
    cys_sg_atoms: List[Dict[str, Any]] = []

    for res in topology.residues():
        resname = res.name
        if not is_protein_residue(resname):
            continue

        if resname != DISULFIDE_RESNAME:
            continue

        resid = safe_int(res.id)

        sg_atom = get_atom_in_residue(res, DISULFIDE_ATOM_NAME)
        if sg_atom is None:
            continue

        sg_xyz = to_angstrom(positions_nm[atom_index[sg_atom]])
        cys_sg_atoms.append(
            {
                "res_id": resid,
                "res_name": resname,
                "atom": DISULFIDE_ATOM_NAME,
                "xyz": sg_xyz,
            }
        )

    return cys_sg_atoms


def find_disulfide_bond_hits(cys_sg_atoms: List[Dict[str, Any]],ss_max_distance_A: float,min_residue_index_gap: int) -> List[Dict[str, Any]]:
    if len(cys_sg_atoms) < 2:
        return []

    res_ids = [item["res_id"] for item in cys_sg_atoms]
    sg_xyz = np.vstack([item["xyz"] for item in cys_sg_atoms])

    result_edge_list: List[Dict[str, Any]] = []
    used_sg_atoms_per_pair: Set[Tuple[str, str, str]] = set()

    use_kdtree = False
    kdtree = None
    try:
        from scipy.spatial import cKDTree  # type: ignore
        kdtree = cKDTree(sg_xyz)
        use_kdtree = True
    except Exception:
        use_kdtree = False

    def passes_residue_gap(a: Any, b: Any) -> bool:
        if a == b:
            return False
        if isinstance(a, int) and isinstance(b, int):
            return abs(a - b) >= int(min_residue_index_gap)
        return True

    if use_kdtree and kdtree is not None:
        pair_indices = kdtree.query_pairs(r=float(ss_max_distance_A))

        for i, j in pair_indices:
            res_i = res_ids[i]
            res_j = res_ids[j]

            if not passes_residue_gap(res_i, res_j):
                continue

            pair_a, pair_b = edge_key_pair(res_i, res_j)
            pair_key = (str(pair_a), str(pair_b))

            sg_i_key = (pair_key[0], pair_key[1], str(res_i))
            sg_j_key = (pair_key[0], pair_key[1], str(res_j))

            if sg_i_key in used_sg_atoms_per_pair:
                continue
            if sg_j_key in used_sg_atoms_per_pair:
                continue

            used_sg_atoms_per_pair.add(sg_i_key)
            used_sg_atoms_per_pair.add(sg_j_key)

            result_edge_list.append(
                {
                    "node1_index": pair_a,
                    "node1_type": "amino_acid",
                    "node2_index": pair_b,
                    "node2_type": "amino_acid",
                    "interaction_type": "SSBOND",
                }
            )

    else:
        n = len(cys_sg_atoms)
        for i in range(n):
            for j in range(i + 1, n):
                res_i = res_ids[i]
                res_j = res_ids[j]

                if not passes_residue_gap(res_i, res_j):
                    continue

                d_ss = float(np.linalg.norm(sg_xyz[i] - sg_xyz[j]))
                if d_ss > float(ss_max_distance_A):
                    continue

                pair_a, pair_b = edge_key_pair(res_i, res_j)
                pair_key = (str(pair_a), str(pair_b))

                sg_i_key = (pair_key[0], pair_key[1], str(res_i))
                sg_j_key = (pair_key[0], pair_key[1], str(res_j))

                if sg_i_key in used_sg_atoms_per_pair:
                    continue
                if sg_j_key in used_sg_atoms_per_pair:
                    continue

                used_sg_atoms_per_pair.add(sg_i_key)
                used_sg_atoms_per_pair.add(sg_j_key)

                result_edge_list.append(
                    {
                        "node1_index": pair_a,
                        "node1_type": "amino_acid",
                        "node2_index": pair_b,
                        "node2_type": "amino_acid",
                        "interaction_type": "SSBOND",
                    }
                )

    return result_edge_list

'''
'''


def postprocess_interaction_report_to_schema(raw_report: Dict[str, Any]) -> Dict[str, Any]:
    """
    Postprocess the raw EnzyWizard-Interaction report into the new JSON Schema format.
    """

    interaction_type_field_map: Dict[str, str] = {
        "HBOND": "hydrogen_bond_count",
        "IONIC": "ionic_bond_count",
        "VDW": "van_der_waals_contact_count",
        "PIPISTACK": "pi_pi_stacking_count",
        "PICATION": "pi_cation_interaction_count",
        "SSBOND": "disulfide_bond_count",
    }

    statistics_scope_field_map: Dict[str, str] = {
        "overall": "overall_molecular_interaction_statistics",
        "intra_protein": "intra_enzyme_interaction_statistics",
        "protein_substrate": "enzyme_substrate_interaction_statistics",
    }

    def map_node(raw_node: Dict[str, Any]) -> Dict[str, Any]:
        raw_node_type = raw_node.get("node_type")

        if raw_node_type == "amino_acid":
            return {
                "residue_index": raw_node.get("aa_index"),
                "residue_name": raw_node.get("aa_name"),
                "node_type": "residue",
            }

        return {
            "substrate_index": raw_node.get("substrate_index"),
            "substrate_name": raw_node.get("substrate_name"),
            "node_type": raw_node.get("node_type"),
        }

    def map_interaction_count(raw_count: Dict[str, Any]) -> Dict[str, Any]:
        return {
            schema_field: raw_count.get(raw_interaction_type)
            for raw_interaction_type, schema_field in interaction_type_field_map.items()
        }

    def map_scope_statistics(raw_scope_statistics: Dict[str, Any]) -> Dict[str, Any]:
        return map_interaction_count(raw_scope_statistics)

    molecular_interactions: List[Dict[str, Any]] = []

    for raw_interaction in raw_report.get("interactions", []):
        molecular_interactions.append(
            {
                "molecular_interaction_type": raw_interaction.get("interaction"),
                "source_node": map_node(raw_interaction.get("node1", {})),
                "target_node": map_node(raw_interaction.get("node2", {})),
            }
        )

    raw_statistics = raw_report.get("interactions_statistics", {})

    molecular_interaction_statistics: Dict[str, Any] = {}

    for raw_scope, schema_scope in statistics_scope_field_map.items():
        molecular_interaction_statistics[schema_scope] = map_scope_statistics(
            raw_statistics.get(raw_scope, {})
        )

    schema_report: Dict[str, Any] = {
        "report_type": raw_report.get("output_type"),
        "molecular_interactions": molecular_interactions,
        "molecular_interaction_statistics": molecular_interaction_statistics,
    }

    return schema_report
