from __future__ import annotations

from typing import NamedTuple, Tuple, Dict, List, Optional, Any

import math
from collections import Counter

import numpy as np
import networkx as nx
from scipy.spatial.distance import cdist

from Bio.PDB.Chain import Chain
from Bio.PDB.NeighborSearch import NeighborSearch

from ..utils.logging_utils import Logger
from ..utils.structure_utils import get_residues_by_chain, get_single_chain




class Cluster(NamedTuple):
    area: float
    residues: List[Tuple[str, int, str]]


sel = "protein and chain A and not backbone and noh and resname ILE VAL LEU"
_ATOMIC_RADII = {'C': 1.88}
water_radius = 1.4
sphere_radius_carbon = _ATOMIC_RADII['C'] + water_radius
sphere_points = 610
sphere_area_const = 4.0 * math.pi * (sphere_radius_carbon ** 2) / sphere_points


def extract_ilv_residue_keys(chain: Chain, logger: Logger) -> List[Tuple[str, int, str]] | None:
    residue_list = get_residues_by_chain(chain, logger)
    if residue_list is None:
        logger.print("[ERROR] Failed to get residue list from chain.")
        return None

    ilv_residue_keys: List[Tuple[str, int, str]] = []
    for residue_key, resname, _ in residue_list:
        if resname in ("ILE", "VAL", "LEU"):
            ilv_residue_keys.append(residue_key)

    if len(ilv_residue_keys) == 0:
        return []

    return ilv_residue_keys

def extract_protein_non_h_atoms(chain: Chain, logger: Logger) -> List[Dict[str, Any]] | None:
    residue_list = get_residues_by_chain(chain, logger)
    if residue_list is None:
        logger.print("[ERROR] Failed to get residue list from chain when extracting protein non-hydrogen atoms.")
        return None

    valid_residue_key_set = {residue_key for residue_key, _, _ in residue_list}

    atom_list: List[Dict[str, Any]] = []

    for res in chain.get_residues():
        residue_key = res.id

        if residue_key not in valid_residue_key_set:
            continue

        for atom in res:
            element = str(atom.element).strip().upper()
            if element == "H":
                continue

            atom_list.append({
                "coord": np.asarray(atom.get_coord(), dtype=float),
                "residue_key": residue_key,
                "bio_atom": atom,
            })

    if len(atom_list) == 0:
        logger.print("[ERROR] No protein non-hydrogen atoms were extracted.")
        return None

    return atom_list


def extract_ilv_atoms(chain: Chain, logger: Logger) -> List[Dict[str, Any]] | None:
    residue_list = get_residues_by_chain(chain, logger)
    if residue_list is None:
        logger.print("[ERROR] Failed to get residue list from chain when extracting ILV atoms.")
        return None

    residue_key_set = set()
    for residue_key, resname, _ in residue_list:
        if resname in ("ILE", "VAL", "LEU"):
            residue_key_set.add(residue_key)

    if len(residue_key_set) == 0:
        return []

    atom_list: List[Dict[str, Any]] = []
    backbone_names = {"N", "CA", "C", "O"}

    for res in chain.get_residues():
        residue_key = res.id

        if residue_key not in residue_key_set:
            continue

        for atom in res:
            atom_name = atom.get_name().strip()
            if atom_name in backbone_names:
                continue

            element = str(atom.element).strip().upper()
            if element == "H":
                continue

            atom_list.append({
                "coord": np.asarray(atom.get_coord(), dtype=float),
                "residue_key": residue_key,
                "bio_atom": atom,
            })

    if len(atom_list) == 0:
        return []

    return atom_list


def generate_sphere_points(coords: np.ndarray,n: int = 610,radius: float = 1.88) -> np.ndarray:
    total_radius = radius + water_radius
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)

    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y * y)
        phi = k * inc
        points.append([math.cos(phi) * r, y, math.sin(phi) * r])

    vec = np.asarray(points, dtype=float)
    vec *= total_radius
    vec += coords
    return vec


class ILVAtom:
    radius = _ATOMIC_RADII["C"]

    def __init__(
            self,
            index: int,
            ilv_atom_list: List[Dict[str, Any]],
            protein_atom_list: List[Dict[str, Any]],
            ns: NeighborSearch,
            protein_atom_id_to_index: Dict[int, int],
            logger: Logger,
    ) -> None:
        self.logger = logger
        self.index = index
        self.coords = ilv_atom_list[index]["coord"]
        self.residue_key = ilv_atom_list[index]["residue_key"]
        self.point_coords = generate_sphere_points(
            self.coords,
            sphere_points,
            self.radius,
        )
        self.neighbor_indices = self.get_neighbors(
            ns,
            protein_atom_id_to_index,
            protein_atom_list,
        )

    def get_neighbors(
        self,
        ns: NeighborSearch,
        protein_atom_id_to_index: Dict[int, int],
        protein_atom_list: List[Dict[str, Any]],
    ) -> np.ndarray:
        nearby_atoms = ns.search(self.coords, 6.56, level="A")

        neighbors: List[int] = []
        for bio_atom in nearby_atoms:
            atom_id = id(bio_atom)
            if atom_id not in protein_atom_id_to_index:
                self.logger.print("[ERROR] Neighbor atom not found in protein atom index mapping.")
                raise ValueError("Neighbor atom not found in protein atom index mapping.")

            protein_idx = protein_atom_id_to_index[atom_id]

            if protein_atom_list[protein_idx]["residue_key"] == self.residue_key:
                continue

            neighbors.append(protein_idx)

        return np.asarray(neighbors, dtype=int)


def retrieve_neighbor_positions(atom: ILVAtom,protein_atom_list: List[Dict[str, Any]]) -> Tuple[np.ndarray, Dict[int, Tuple[str, int, str]]]:
    positions = np.asarray(
        [protein_atom_list[i]["coord"] for i in atom.neighbor_indices],
        dtype=float,
    )
    position_index_to_resid = {
        idx: protein_atom_list[i]["residue_key"]
        for idx, i in enumerate(atom.neighbor_indices)
    }
    return positions, position_index_to_resid


def retrieve_indices(matrix: np.ndarray,coords: np.ndarray,neighborpositions: np.ndarray,radius: float = 1.88) -> List[int]:
    sphere_radius = water_radius + radius
    dist_center_atoms = cdist(np.reshape(coords, (1, 3)), neighborpositions)
    ranking = np.argsort(dist_center_atoms)

    valid = matrix <= sphere_radius
    idx2: List[int] = []

    for row in valid:
        if row.any():
            idx2.append(ranking[0][np.where(np.isin(ranking, np.where(row)))[1][0]])

    return idx2


def fill_matrices(atom: ILVAtom,protein_atom_list: List[Dict[str, Any]],resid_matrix: np.ndarray,protein_to_ilv_index: Dict[int, int],atom_to_residposition: Dict[int, int],logger:Logger) -> np.ndarray | None:
    neighbor_positions, _ = retrieve_neighbor_positions(atom, protein_atom_list)
    distances = cdist(atom.point_coords, neighbor_positions)

    column_indices = retrieve_indices(distances, atom.coords, neighbor_positions)
    colpos_occurrences = Counter(column_indices)

    for colpos, occurrences in colpos_occurrences.items():
        protein_neighbor_idx = atom.neighbor_indices[colpos]

        if protein_neighbor_idx not in protein_to_ilv_index:
            continue

        ilv_neighbor_idx = protein_to_ilv_index[protein_neighbor_idx]

        if atom.index not in atom_to_residposition:
            logger.print(f"[ERROR] ILV atom index {atom.index} not found in atom-to-residue-position mapping.")
            return None

        if ilv_neighbor_idx not in atom_to_residposition:
            logger.print(
                f"[ERROR] ILV neighbor atom index {ilv_neighbor_idx} not found in atom-to-residue-position mapping.")
            return None

        area = sphere_area_const * occurrences
        index_i = atom_to_residposition[atom.index]
        index_j = atom_to_residposition[ilv_neighbor_idx]
        resid_matrix[index_i, index_j] += area

    return resid_matrix


def create_graph(resid_matrix: np.ndarray,resid_list: List[Tuple[str, int, str]],cutoff_area: float = 10.0) -> nx.DiGraph:
    g = nx.DiGraph()

    for i, residue_key in enumerate(resid_list):
        g.add_node(i, residue_key=residue_key)

    for row_index, row in enumerate(resid_matrix):
        for column_index, area in enumerate(row):
            if row_index != column_index and area > cutoff_area:
                g.add_edge(row_index, column_index, area=float(area))

    return g


def add_clusters(g: nx.DiGraph) -> List[Cluster]:
    clusters: List[Cluster] = []

    for component in nx.weakly_connected_components(g):
        if len(component) < 2:
            continue

        sub = g.subgraph(component)
        area = float(sum(data["area"] for _, _, data in sub.edges(data=True)))
        residues = [g.nodes[i]["residue_key"] for i in component]

        clusters.append(Cluster(
            area=area,
            residues=residues,
        ))

    return clusters

def postprocess_hydrocluster_report_to_schema(
    raw_report: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Postprocess the raw EnzyWizard-Hydrocluster report into the new JSON Schema format.
    """

    raw_statistics = raw_report.get("hydrophobic_cluster_statistics", {})
    raw_clusters = raw_report.get("hydrophobic_cluster", [])

    hydrophobic_clusters: List[Dict[str, Any]] = []

    for raw_cluster in raw_clusters:
        residues: List[Dict[str, Any]] = []

        for raw_residue in raw_cluster.get("residues", []):
            residues.append(
                {
                    "residue_index": raw_residue.get("aa_id"),
                    "residue_name": raw_residue.get("aa_name"),
                }
            )

        hydrophobic_clusters.append(
            {
                "hydrophobic_cluster_area": raw_cluster.get("area"),
                "residues": residues,
            }
        )

    schema_report: Dict[str, Any] = {
        "report_type": raw_report.get("output_type", "enzywizard_hydrocluster"),
        "hydrophobic_cluster_statistics": {
            "hydrophobic_cluster_count": raw_statistics.get("cluster_num"),
            "max_hydrophobic_cluster_area": raw_statistics.get("max_cluster_area"),
            "total_hydrophobic_cluster_area": raw_statistics.get("total_cluster_area"),
        },
        "hydrophobic_clusters": hydrophobic_clusters,
    }

    return schema_report