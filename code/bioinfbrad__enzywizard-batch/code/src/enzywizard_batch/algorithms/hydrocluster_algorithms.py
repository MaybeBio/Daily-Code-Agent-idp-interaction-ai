from __future__ import annotations
from ..utils.hydrocluster_utils import *

from Bio.PDB.Structure import Structure
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


def compute_hydrophobic_clusters(struct: Structure,logger: Logger,cutoff_area: float = 10.0) -> List[Cluster] | None:

    chain = get_single_chain(struct, logger)
    if chain is None:
        return None

    resid_list = extract_ilv_residue_keys(chain, logger)
    if resid_list is None:
        return None

    ilv_atom_list = extract_ilv_atoms(chain, logger)
    if ilv_atom_list is None:
        return None

    protein_atom_list = extract_protein_non_h_atoms(chain, logger)
    if protein_atom_list is None:
        return None

    if len(protein_atom_list) == 0:
        logger.print("[ERROR] No protein non-hydrogen atoms found in structure.")
        return None

    if len(resid_list) == 0:
        logger.print("[INFO] No ILE, LEU, or VAL residues found; no hydrophobic clusters to report.")
        return []

    if len(ilv_atom_list) == 0:
        logger.print("[INFO] No non-backbone, non-hydrogen ILE/LEU/VAL atoms found; no hydrophobic clusters to report.")
        return []

    dims = len(resid_list)
    indices = np.arange(len(ilv_atom_list), dtype=int)

    atom_to_residposition: Dict[int, int] = {}
    residue_key_to_position = {residue_key: i for i, residue_key in enumerate(resid_list)}

    for idx, atom in enumerate(ilv_atom_list):
        residue_key = atom["residue_key"]
        if residue_key not in residue_key_to_position:
            logger.print(f"[ERROR] Atom residue key {residue_key} not found in ILV residue list.")
            return None
        atom_to_residposition[idx] = residue_key_to_position[residue_key]

    protein_atom_id_to_index: Dict[int, int] = {}
    for protein_idx, atom in enumerate(protein_atom_list):
        protein_atom_id_to_index[id(atom["bio_atom"])] = protein_idx

    if len(protein_atom_id_to_index) == 0:
        logger.print("[ERROR] Protein atom id to index mapping is empty.")
        return None

    ilv_atom_id_to_index: Dict[int, int] = {}
    for ilv_idx, atom in enumerate(ilv_atom_list):
        ilv_atom_id_to_index[id(atom["bio_atom"])] = ilv_idx

    if len(ilv_atom_id_to_index) == 0:
        logger.print("[ERROR] ILV atom id to index mapping is empty.")
        return None

    protein_to_ilv_index: Dict[int, int] = {}
    for protein_idx, atom in enumerate(protein_atom_list):
        atom_id = id(atom["bio_atom"])
        if atom_id in ilv_atom_id_to_index:
            protein_to_ilv_index[protein_idx] = ilv_atom_id_to_index[atom_id]

    if len(protein_to_ilv_index) == 0:
        logger.print("[ERROR] No overlap found between protein atom list and ILV atom list.")
        return None

    ns = NeighborSearch([atom["bio_atom"] for atom in protein_atom_list])

    contacts = np.zeros((dims, dims), dtype=float)

    for index in indices:
        try:
            a = ILVAtom(
                index=index,
                ilv_atom_list=ilv_atom_list,
                protein_atom_list=protein_atom_list,
                ns=ns,
                protein_atom_id_to_index=protein_atom_id_to_index,
                logger=logger,
            )
        except Exception as e:
            logger.print(f"[ERROR] Failed to build ILVAtom for index {index}: {e}")
            return None
        if len(a.neighbor_indices) == 0:
            continue

        updated_contacts = fill_matrices(
            atom=a,
            protein_atom_list=protein_atom_list,
            resid_matrix=contacts,
            protein_to_ilv_index=protein_to_ilv_index,
            atom_to_residposition=atom_to_residposition,
            logger=logger,
        )
        if updated_contacts is None:
            logger.print(f"[ERROR] Failed to fill contact matrix for ILV atom index {index}.")
            return None

        contacts = updated_contacts

    graph = create_graph(contacts, resid_list, cutoff_area=cutoff_area)
    clusters = add_clusters(graph)

    clusters.sort(key=lambda x: x.area, reverse=True)

    return clusters


def calculate_hydrocluster_statistics(clusters: List[Cluster], logger: Logger) -> Dict[str, Any] | None:
    if not isinstance(clusters, list):
        logger.print("[ERROR] clusters must be a list.")
        return None

    for cluster in clusters:
        if not isinstance(cluster, Cluster):
            logger.print("[ERROR] Invalid cluster item in clusters.")
            return None

        if not hasattr(cluster, "area"):
            logger.print("[ERROR] Cluster item missing area attribute.")
            return None

        if not isinstance(cluster.area, (int, float)):
            logger.print("[ERROR] cluster.area must be a number.")
            return None

    cluster_num = len(clusters)
    max_cluster_area = 0.0
    total_cluster_area = 0.0

    if cluster_num > 0:
        max_cluster_area = max(float(cluster.area) for cluster in clusters)
        total_cluster_area = sum(float(cluster.area) for cluster in clusters)

    return {
        "cluster_num": cluster_num,
        "max_cluster_area": max_cluster_area,
        "total_cluster_area": total_cluster_area,
    }

def generate_hydrocluster_report(clusters: List[Cluster],struct: Structure,logger: Logger) -> Dict[str, Any] | None:
    chain = get_single_chain(struct, logger)
    if chain is None:
        logger.print("[ERROR] Failed to get single chain when generating hydrophobic cluster report.")
        return None

    residue_list = get_residues_by_chain(chain, logger)
    if residue_list is None:
        logger.print("[ERROR] Failed to get residue list when generating hydrophobic cluster report.")
        return None

    hydrophobic_cluster_statistics = calculate_hydrocluster_statistics(clusters, logger)
    if hydrophobic_cluster_statistics is None:
        logger.print("[ERROR] Failed to calculate hydrophobic cluster statistics.")
        return None

    residue_key_to_name: Dict[Tuple[str, int, str], str] = {}
    for residue_key, resname, _ in residue_list:
        residue_key_to_name[residue_key] = resname

    hydrophobic_cluster_list: List[Dict[str, Any]] = []

    for cluster in sorted(clusters, key=lambda x: x.area, reverse=True):
        residue_dict_list: List[Dict[str, Any]] = []

        for residue_key in cluster.residues:
            if residue_key not in residue_key_to_name:
                logger.print(f"[ERROR] Residue key {residue_key} not found in structure residue list.")
                return None

            _, resseq, _ = residue_key
            resname = residue_key_to_name[residue_key]

            residue_dict_list.append({
                "aa_id": resseq,
                "aa_name": normalize_aa_name_to_one_letter(resname),
            })

        hydrophobic_cluster_list.append({
            "area": cluster.area,
            "residues": residue_dict_list,
        })

    raw_report = {
        "output_type": "enzywizard_hydrocluster",
        "hydrophobic_cluster_statistics": hydrophobic_cluster_statistics,
        "hydrophobic_cluster": hydrophobic_cluster_list,
    }

    return postprocess_hydrocluster_report_to_schema(
        raw_report=raw_report
    )
