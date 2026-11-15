from __future__ import annotations
from rdkit import Chem
from openmm.app import Modeller

from typing import Tuple, List, Any, Dict, Set
from Bio.PDB.Structure import Structure
from Bio.Data.IUPACData import protein_letters_3to1
from ..utils.structure_utils import get_single_chain, get_residues_by_chain

from ..utils.interaction_utils import build_openmm_tables, collect_protein_hbond_sites, collect_substrate_hbond_sites, find_hbond_hits_from_donors_to_acceptors, edge_key_pair, is_protein_substrate_docked, build_substrate_tables, edge_sort_key
from ..utils.interaction_utils import collect_protein_ionic_centers,collect_substrate_ionic_centers, find_ionic_hits_between_centers
from ..utils.interaction_utils import collect_substrate_vdw_atoms,collect_protein_vdw_atoms,find_vdw_hits_between_atom_entries
from ..utils.interaction_utils import collect_protein_pipi_rings, collect_substrate_pipi_rings, find_pipi_hits_between_ring_entries, min_angle_0_90, angle_deg_between_vectors, classify_pipi_geometry
from ..utils.interaction_utils import collect_protein_pication_centers, classify_pication_arg_geometry, find_pication_hits_between_ring_entries_and_cation_entries, collect_substrate_pication_centers
from ..utils.interaction_utils import collect_protein_disulfide_sites, find_disulfide_bond_hits
from ..utils.interaction_utils import postprocess_interaction_report_to_schema
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


def calculate_hydrogen_bond_network(
    modeller: Modeller,
    ligand_mol_list: List[Chem.Mol],
    logger,
    bonded_h_min_distance_A: float = 0.8,
    bonded_h_max_distance_A: float = 1.3,
    da_max_distance_A: float = 3.9,
    ha_max_distance_A: float = 2.5,
    dha_min_angle_deg: float = 90,
    docked_heavy_atom_distance_cutoff_A: float = 6.5
) -> List[Dict[str, Any]] | None:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(ligand_mol_list, list):
        logger.print("[ERROR] ligand_mol_list must be a list.")
        return None

    if not isinstance(bonded_h_min_distance_A, (int, float)) or float(bonded_h_min_distance_A) <= 0.0:
        logger.print("[ERROR] bonded_h_min_distance_A must be a positive number.")
        return None

    if not isinstance(bonded_h_max_distance_A, (int, float)) or float(bonded_h_max_distance_A) <= 0.0:
        logger.print("[ERROR] bonded_h_max_distance_A must be a positive number.")
        return None

    if float(bonded_h_min_distance_A) >= float(bonded_h_max_distance_A):
        logger.print("[ERROR] bonded_h_min_distance_A must be smaller than bonded_h_max_distance_A.")
        return None

    if not isinstance(da_max_distance_A, (int, float)) or float(da_max_distance_A) <= 0.0:
        logger.print("[ERROR] da_max_distance_A must be a positive number.")
        return None

    if not isinstance(ha_max_distance_A, (int, float)) or float(ha_max_distance_A) <= 0.0:
        logger.print("[ERROR] ha_max_distance_A must be a positive number.")
        return None

    if not isinstance(dha_min_angle_deg, (int, float)) or float(dha_min_angle_deg) < 0.0:
        logger.print("[ERROR] dha_min_angle_deg must be a non-negative number.")
        return None

    if not isinstance(docked_heavy_atom_distance_cutoff_A, (int, float)) or float(docked_heavy_atom_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] docked_heavy_atom_distance_cutoff_A must be a positive number.")
        return None

    try:
        topology, positions_nm, atoms, atom_index, coords_A_prot, bonded_prot = build_openmm_tables(modeller)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM tables from modeller: {e}")
        return None

    try:
        prot_acceptors, prot_donors = collect_protein_hbond_sites(
            topology=topology,
            positions_nm=positions_nm,
            atoms=atoms,
            atom_index=atom_index,
            coords_A=coords_A_prot,
            bonded=bonded_prot,
            bonded_h_min_distance_A=float(bonded_h_min_distance_A),
            bonded_h_max_distance_A=float(bonded_h_max_distance_A),
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect protein hydrogen-bond sites: {e}")
        return None

    result_edge_list: List[Dict[str, Any]] = []

    # -----------------------------
    # part 1) intra-protein H bonds
    # -----------------------------
    if len(prot_acceptors) > 0 and len(prot_donors) > 0:
        used_donor_atoms_per_pair: Set[Tuple[str, str, str, str]] = set()
        used_acceptor_atoms_per_pair: Set[Tuple[str, str, str, str]] = set()

        pp_hits = find_hbond_hits_from_donors_to_acceptors(
            donor_entries=prot_donors,
            acceptor_atoms=prot_acceptors,
            da_max_distance_A=float(da_max_distance_A),
            ha_max_distance_A=float(ha_max_distance_A),
            dha_min_angle_deg=float(dha_min_angle_deg),
        )

        for donor, acceptor in pp_hits:
            d_res = donor["res_id"]
            a_res = acceptor["res_id"]

            if d_res == a_res:
                continue

            if isinstance(d_res, int) and isinstance(a_res, int):
                if abs(d_res - a_res) < 3:
                    continue

            pair_a, pair_b = edge_key_pair(d_res, a_res)
            pair_key = (str(pair_a), str(pair_b))

            donor_key = (pair_key[0], pair_key[1], str(d_res), str(donor["atom"]))
            acceptor_key = (pair_key[0], pair_key[1], str(a_res), str(acceptor["atom"]))

            if donor_key in used_donor_atoms_per_pair:
                continue
            if acceptor_key in used_acceptor_atoms_per_pair:
                continue

            used_donor_atoms_per_pair.add(donor_key)
            used_acceptor_atoms_per_pair.add(acceptor_key)

            result_edge_list.append(
                {
                    "node1_index": pair_a,
                    "node1_type": "amino_acid",
                    "node2_index": pair_b,
                    "node2_type": "amino_acid",
                    "interaction_type": "HBOND",
                }
            )

    # -----------------------------
    # part 2) protein-ligand H bonds
    # -----------------------------
    for ligand_index, lig_mol in enumerate(ligand_mol_list, start=1):
        if lig_mol is None:
            logger.print(f"[WARNING] Skip invalid ligand mol at index {ligand_index}.")
            continue

        if not is_protein_substrate_docked(
            modeller=modeller,
            ligand_mol=lig_mol,
            logger=logger,
            heavy_atom_distance_cutoff_A=float(docked_heavy_atom_distance_cutoff_A),
        ):
            logger.print(f"[WARNING] Substrate at index {ligand_index} is not spatially docked to protein. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")

        try:
            lig_coords_A, bonded_lig = build_substrate_tables(lig_mol)
        except Exception as e:
            logger.print(f"[WARNING] Failed to build ligand tables at index {ligand_index}: {e}")
            continue

        lig_sites = collect_substrate_hbond_sites(
            mol=lig_mol,
            coords_A=lig_coords_A,
            bonded=bonded_lig,
            bonded_h_min_distance_A=float(bonded_h_min_distance_A),
            bonded_h_max_distance_A=float(bonded_h_max_distance_A),
        )
        if lig_sites is None:
            logger.print(f"[WARNING] Failed to collect ligand hydrogen-bond sites at index {ligand_index}.")
            continue

        lig_acceptors, lig_donors = lig_sites
        if (len(lig_acceptors) == 0 and len(lig_donors) == 0) or (len(prot_acceptors) == 0 and len(prot_donors) == 0):
            continue

        used_protein_residue_set: Set[str] = set()

        # protein donor -> ligand acceptor
        if len(prot_donors) > 0 and len(lig_acceptors) > 0:
            hits = find_hbond_hits_from_donors_to_acceptors(
                donor_entries=prot_donors,
                acceptor_atoms=lig_acceptors,
                da_max_distance_A=float(da_max_distance_A),
                ha_max_distance_A=float(ha_max_distance_A),
                dha_min_angle_deg=float(dha_min_angle_deg),
            )
            for donor, _ in hits:
                res_id = donor["res_id"]
                key = f"{str(res_id)}__{ligand_index}"
                if key in used_protein_residue_set:
                    continue
                used_protein_residue_set.add(key)

                result_edge_list.append(
                    {
                        "node1_index": res_id,
                        "node1_type": "amino_acid",
                        "node2_index": ligand_index,
                        "node2_type": "substrate",
                        "interaction_type": "HBOND",
                    }
                )

        # ligand donor -> protein acceptor
        if len(lig_donors) > 0 and len(prot_acceptors) > 0:
            hits = find_hbond_hits_from_donors_to_acceptors(
                donor_entries=lig_donors,
                acceptor_atoms=prot_acceptors,
                da_max_distance_A=float(da_max_distance_A),
                ha_max_distance_A=float(ha_max_distance_A),
                dha_min_angle_deg=float(dha_min_angle_deg),
            )
            for _, acceptor in hits:
                res_id = acceptor["res_id"]
                key = f"{str(res_id)}__{ligand_index}"
                if key in used_protein_residue_set:
                    continue
                used_protein_residue_set.add(key)

                result_edge_list.append(
                    {
                        "node1_index": res_id,
                        "node1_type": "amino_acid",
                        "node2_index": ligand_index,
                        "node2_type": "substrate",
                        "interaction_type": "HBOND",
                    }
                )

    result_edge_list.sort(key=edge_sort_key)
    return result_edge_list

def calculate_ionic_bond_network(
    modeller: Modeller,
    ligand_mol_list: List[Chem.Mol],
    logger,
    ionic_distance_cutoff_A: float = 4.0,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
    min_residue_index_gap: int = 3,
) -> List[Dict[str, Any]] | None:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(ligand_mol_list, list):
        logger.print("[ERROR] ligand_mol_list must be a list.")
        return None

    if not isinstance(ionic_distance_cutoff_A, (int, float)) or float(ionic_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] ionic_distance_cutoff_A must be a positive number.")
        return None

    if not isinstance(docked_heavy_atom_distance_cutoff_A, (int, float)) or float(docked_heavy_atom_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] docked_heavy_atom_distance_cutoff_A must be a positive number.")
        return None

    if not isinstance(min_residue_index_gap, int) or int(min_residue_index_gap) < 0:
        logger.print("[ERROR] min_residue_index_gap must be a non-negative integer.")
        return None

    try:
        topology, positions_nm, atoms, atom_index, coords_A_prot, bonded_prot = build_openmm_tables(modeller)

    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM tables from modeller: {e}")
        return None

    try:
        prot_centers = collect_protein_ionic_centers(
            topology=topology,
            positions_nm=positions_nm,
            atom_index=atom_index,
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect protein ionic centers: {e}")
        return None

    result_edge_list: List[Dict[str, Any]] = []

    # -----------------------------
    # part 1) intra-protein ionic bonds
    # -----------------------------
    prot_cations = [x for x in prot_centers if x["charge"] == "cation"]
    prot_anions = [x for x in prot_centers if x["charge"] == "anion"]

    if len(prot_cations) > 0 and len(prot_anions) > 0:
        used_cation_centers_per_pair: Set[Tuple[str, str, str, str]] = set()
        used_anion_centers_per_pair: Set[Tuple[str, str, str, str]] = set()

        pp_hits = find_ionic_hits_between_centers(
            center_entries_1=prot_cations,
            center_entries_2=prot_anions,
            distance_cutoff_A=float(ionic_distance_cutoff_A),
        )

        for cation, anion in pp_hits:
            c_res = cation["res_id"]
            a_res = anion["res_id"]

            if c_res == a_res:
                continue

            if isinstance(c_res, int) and isinstance(a_res, int):
                if abs(c_res - a_res) < int(min_residue_index_gap):
                    continue

            pair_a, pair_b = edge_key_pair(c_res, a_res)
            pair_key = (str(pair_a), str(pair_b))

            cation_key = (
                pair_key[0],
                pair_key[1],
                str(c_res),
                str(cation.get("center_id", "")),
            )
            anion_key = (
                pair_key[0],
                pair_key[1],
                str(a_res),
                str(anion.get("center_id", "")),
            )

            if cation_key in used_cation_centers_per_pair:
                continue
            if anion_key in used_anion_centers_per_pair:
                continue

            used_cation_centers_per_pair.add(cation_key)
            used_anion_centers_per_pair.add(anion_key)

            result_edge_list.append(
                {
                    "node1_index": pair_a,
                    "node1_type": "amino_acid",
                    "node2_index": pair_b,
                    "node2_type": "amino_acid",
                    "interaction_type": "IONIC",
                }
            )

    # -----------------------------
    # part 2) protein-ligand ionic bonds
    # -----------------------------
    for ligand_index, lig_mol in enumerate(ligand_mol_list, start=1):
        if lig_mol is None:
            logger.print(f"[WARNING] Skip invalid ligand mol at index {ligand_index}.")
            continue


        if not is_protein_substrate_docked(
            modeller=modeller,
            ligand_mol=lig_mol,
            logger=logger,
            heavy_atom_distance_cutoff_A=float(docked_heavy_atom_distance_cutoff_A),
        ):
            logger.print(
                f"[WARNING] Substrate at index {ligand_index} is not spatially docked to protein. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")

        try:
            lig_coords_A, bonded_lig = build_substrate_tables(lig_mol)

        except Exception as e:
            logger.print(f"[WARNING] Failed to build ligand tables at index {ligand_index}: {e}")
            continue

        try:
            lig_centers = collect_substrate_ionic_centers(
                mol=lig_mol,
                coords_A=lig_coords_A,
            )
        except Exception as e:
            logger.print(f"[WARNING] Failed to collect ligand ionic centers at index {ligand_index}: {e}")
            continue

        if len(lig_centers) == 0 or len(prot_centers) == 0:
            continue

        lig_cations = [x for x in lig_centers if x["charge"] == "cation"]
        lig_anions = [x for x in lig_centers if x["charge"] == "anion"]

        if (len(prot_cations) == 0 and len(prot_anions) == 0) or (len(lig_cations) == 0 and len(lig_anions) == 0):
            continue

        used_protein_residue_set: Set[str] = set()

        # protein cation -> ligand anion
        if len(prot_cations) > 0 and len(lig_anions) > 0:
            hits = find_ionic_hits_between_centers(
                center_entries_1=prot_cations,
                center_entries_2=lig_anions,
                distance_cutoff_A=float(ionic_distance_cutoff_A),
            )
            for prot_center, _ in hits:
                res_id = prot_center["res_id"]
                key = f"{str(res_id)}__{ligand_index}"
                if key in used_protein_residue_set:
                    continue
                used_protein_residue_set.add(key)

                result_edge_list.append(
                    {
                        "node1_index": res_id,
                        "node1_type": "amino_acid",
                        "node2_index": ligand_index,
                        "node2_type": "substrate",
                        "interaction_type": "IONIC",
                    }
                )

        # protein anion -> ligand cation
        if len(prot_anions) > 0 and len(lig_cations) > 0:
            hits = find_ionic_hits_between_centers(
                center_entries_1=prot_anions,
                center_entries_2=lig_cations,
                distance_cutoff_A=float(ionic_distance_cutoff_A),
            )
            for prot_center, _ in hits:
                res_id = prot_center["res_id"]
                key = f"{str(res_id)}__{ligand_index}"
                if key in used_protein_residue_set:
                    continue
                used_protein_residue_set.add(key)

                result_edge_list.append(
                    {
                        "node1_index": res_id,
                        "node1_type": "amino_acid",
                        "node2_index": ligand_index,
                        "node2_type": "substrate",
                        "interaction_type": "IONIC",
                    }
                )

    result_edge_list.sort(key=edge_sort_key)
    return result_edge_list


def calculate_van_der_waals_network(
    modeller: Modeller,
    ligand_mol_list: List[Chem.Mol],
    logger,
    mu: float = 0.01,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
    min_residue_index_gap: int = 3,
) -> List[Dict[str, Any]] | None:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(ligand_mol_list, list):
        logger.print("[ERROR] ligand_mol_list must be a list.")
        return None

    if not isinstance(mu, (int, float)) or float(mu) < 0.0:
        logger.print("[ERROR] mu must be a non-negative number.")
        return None

    if not isinstance(docked_heavy_atom_distance_cutoff_A, (int, float)) or float(docked_heavy_atom_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] docked_heavy_atom_distance_cutoff_A must be a positive number.")
        return None

    if not isinstance(min_residue_index_gap, int) or int(min_residue_index_gap) < 0:
        logger.print("[ERROR] min_residue_index_gap must be a non-negative integer.")
        return None

    try:
        topology, positions_nm, atoms, atom_index, coords_A_prot, bonded_prot = build_openmm_tables(modeller)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM tables from modeller: {e}")
        return None

    try:
        prot_atom_entries = collect_protein_vdw_atoms(
            topology=topology,
            atoms=atoms,
            atom_index=atom_index,
            coords_A=coords_A_prot,
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect protein van der Waals atoms: {e}")
        return None

    result_edge_list: List[Dict[str, Any]] = []

    # -----------------------------
    # part 1) intra-protein vdW
    # -----------------------------
    if len(prot_atom_entries) > 0:
        used_atom_entries_per_pair: Set[Tuple[str, str, str, str]] = set()

        pp_hits = find_vdw_hits_between_atom_entries(
            atom_entries_1=prot_atom_entries,
            atom_entries_2=prot_atom_entries,
            mu=float(mu),
            allow_same_entry_list=True,
        )

        for atom1, atom2 in pp_hits:
            res1 = atom1["res_id"]
            res2 = atom2["res_id"]

            if res1 == res2:
                continue

            if isinstance(res1, int) and isinstance(res2, int):
                if abs(res1 - res2) < int(min_residue_index_gap):
                    continue

            pair_a, pair_b = edge_key_pair(res1, res2)
            pair_key = (str(pair_a), str(pair_b))

            atom1_key = (pair_key[0], pair_key[1], str(res1), str(atom1["atom"]))
            atom2_key = (pair_key[0], pair_key[1], str(res2), str(atom2["atom"]))

            if atom1_key in used_atom_entries_per_pair:
                continue
            if atom2_key in used_atom_entries_per_pair:
                continue

            used_atom_entries_per_pair.add(atom1_key)
            used_atom_entries_per_pair.add(atom2_key)

            result_edge_list.append(
                {
                    "node1_index": pair_a,
                    "node1_type": "amino_acid",
                    "node2_index": pair_b,
                    "node2_type": "amino_acid",
                    "interaction_type": "VDW",
                }
            )

    # -----------------------------
    # part 2) protein-ligand vdW
    # -----------------------------
    for ligand_index, lig_mol in enumerate(ligand_mol_list, start=1):
        if lig_mol is None:
            logger.print(f"[WARNING] Skip invalid ligand mol at index {ligand_index}.")
            continue

        if not is_protein_substrate_docked(
            modeller=modeller,
            ligand_mol=lig_mol,
            logger=logger,
            heavy_atom_distance_cutoff_A=float(docked_heavy_atom_distance_cutoff_A),
        ):
            logger.print(
                f"[WARNING] Substrate at index {ligand_index} is not spatially docked to protein. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")

        try:
            lig_coords_A, bonded_lig = build_substrate_tables(lig_mol)
        except Exception as e:
            logger.print(f"[WARNING] Failed to build ligand tables at index {ligand_index}: {e}")
            continue

        try:
            lig_atom_entries = collect_substrate_vdw_atoms(
                mol=lig_mol,
                coords_A=lig_coords_A,
            )
        except Exception as e:
            logger.print(f"[WARNING] Failed to collect ligand van der Waals atoms at index {ligand_index}: {e}")
            continue

        if len(prot_atom_entries) == 0 or len(lig_atom_entries) == 0:
            continue

        try:
            pl_hits = find_vdw_hits_between_atom_entries(
                atom_entries_1=prot_atom_entries,
                atom_entries_2=lig_atom_entries,
                mu=float(mu),
                allow_same_entry_list=False,
            )
        except Exception as e:
            logger.print(f"[WARNING] Failed to detect van der Waals hits for ligand at index {ligand_index}: {e}")
            continue

        used_protein_residue_set: Set[str] = set()

        for prot_atom, _ in pl_hits:
            res_id = prot_atom["res_id"]
            key = f"{str(res_id)}__{ligand_index}"
            if key in used_protein_residue_set:
                continue

            used_protein_residue_set.add(key)

            result_edge_list.append(
                {
                    "node1_index": res_id,
                    "node1_type": "amino_acid",
                    "node2_index": ligand_index,
                    "node2_type": "substrate",
                    "interaction_type": "VDW",
                }
            )

    result_edge_list.sort(key=edge_sort_key)
    return result_edge_list

def calculate_pipi_stacking_network(
    modeller: Modeller,
    ligand_mol_list: List[Chem.Mol],
    logger,
    ring_center_distance_cutoff_A: float = 6.5,
    min_residue_index_gap: int = 3,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
) -> List[Dict[str, Any]] | None:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(ligand_mol_list, list):
        logger.print("[ERROR] ligand_mol_list must be a list.")
        return None

    if not isinstance(ring_center_distance_cutoff_A, (int, float)) or float(ring_center_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] ring_center_distance_cutoff_A must be a positive number.")
        return None

    if not isinstance(min_residue_index_gap, int) or int(min_residue_index_gap) < 0:
        logger.print("[ERROR] min_residue_index_gap must be a non-negative integer.")
        return None

    if not isinstance(docked_heavy_atom_distance_cutoff_A, (int, float)) or float(docked_heavy_atom_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] docked_heavy_atom_distance_cutoff_A must be a positive number.")
        return None

    try:
        topology, positions_nm, atoms, atom_index, coords_A_prot, bonded_prot = build_openmm_tables(modeller)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM tables from modeller: {e}")
        return None

    try:
        prot_ring_entries = collect_protein_pipi_rings(
            topology=topology,
            positions_nm=positions_nm,
            atom_index=atom_index,
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect protein pi-pi rings: {e}")
        return None

    result_edge_list: List[Dict[str, Any]] = []

    # -----------------------------
    # part 1) intra-protein pi-pi stacking
    # -----------------------------
    if len(prot_ring_entries) > 0:
        used_ring_entries_per_pair: Set[Tuple[str, str, str, str]] = set()

        try:
            pp_hits = find_pipi_hits_between_ring_entries(
                ring_entries_1=prot_ring_entries,
                ring_entries_2=prot_ring_entries,
                ring_center_distance_cutoff_A=float(ring_center_distance_cutoff_A),
                allow_same_entry_list=True,
            )
        except Exception as e:
            logger.print(f"[ERROR] Failed to detect intra-protein pi-pi stacking hits: {e}")
            return None

        for ring1, ring2 in pp_hits:
            res1 = ring1["res_id"]
            res2 = ring2["res_id"]

            if res1 == res2:
                continue

            if isinstance(res1, int) and isinstance(res2, int):
                if abs(res1 - res2) < int(min_residue_index_gap):
                    continue

            center1 = ring1["center"]
            center2 = ring2["center"]
            normal1 = ring1["normal"]
            normal2 = ring2["normal"]

            gamma = min_angle_0_90(angle_deg_between_vectors(normal1, normal2))
            v12 = center2 - center1
            theta = min_angle_0_90(angle_deg_between_vectors(normal1, v12))
            delta = min_angle_0_90(angle_deg_between_vectors(normal2, v12))

            _ = classify_pipi_geometry(theta_deg=theta, delta_deg=delta, gamma_deg=gamma)

            pair_a, pair_b = edge_key_pair(res1, res2)
            pair_key = (str(pair_a), str(pair_b))

            ring1_key = (pair_key[0], pair_key[1], str(res1), str(ring1["ring_id"]))
            ring2_key = (pair_key[0], pair_key[1], str(res2), str(ring2["ring_id"]))

            if ring1_key in used_ring_entries_per_pair:
                continue
            if ring2_key in used_ring_entries_per_pair:
                continue

            used_ring_entries_per_pair.add(ring1_key)
            used_ring_entries_per_pair.add(ring2_key)

            result_edge_list.append(
                {
                    "node1_index": pair_a,
                    "node1_type": "amino_acid",
                    "node2_index": pair_b,
                    "node2_type": "amino_acid",
                    "interaction_type": "PIPISTACK",
                }
            )

    # -----------------------------
    # part 2) protein-ligand pi-pi stacking
    # -----------------------------
    for ligand_index, lig_mol in enumerate(ligand_mol_list, start=1):
        if lig_mol is None:
            logger.print(f"[WARNING] Skip invalid ligand mol at index {ligand_index}.")
            continue

        if not is_protein_substrate_docked(
            modeller=modeller,
            ligand_mol=lig_mol,
            logger=logger,
            heavy_atom_distance_cutoff_A=float(docked_heavy_atom_distance_cutoff_A),
        ):
            logger.print(
                f"[WARNING] Substrate at index {ligand_index} is not spatially docked to protein. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")

        try:
            lig_coords_A, bonded_lig = build_substrate_tables(lig_mol)
        except Exception as e:
            logger.print(f"[WARNING] Failed to build ligand tables at index {ligand_index}: {e}")
            continue

        try:
            lig_ring_entries = collect_substrate_pipi_rings(
                mol=lig_mol,
                coords_A=lig_coords_A,
            )
        except Exception as e:
            logger.print(f"[WARNING] Failed to collect ligand pi-pi rings at index {ligand_index}: {e}")
            continue

        if len(prot_ring_entries) == 0 or len(lig_ring_entries) == 0:
            continue

        try:
            pl_hits = find_pipi_hits_between_ring_entries(
                ring_entries_1=prot_ring_entries,
                ring_entries_2=lig_ring_entries,
                ring_center_distance_cutoff_A=float(ring_center_distance_cutoff_A),
                allow_same_entry_list=False,
            )
        except Exception as e:
            logger.print(f"[WARNING] Failed to detect pi-pi stacking hits for ligand at index {ligand_index}: {e}")
            continue

        used_protein_residue_set: Set[str] = set()

        for prot_ring, lig_ring in pl_hits:
            center1 = prot_ring["center"]
            center2 = lig_ring["center"]
            normal1 = prot_ring["normal"]
            normal2 = lig_ring["normal"]

            gamma = min_angle_0_90(angle_deg_between_vectors(normal1, normal2))
            v12 = center2 - center1
            theta = min_angle_0_90(angle_deg_between_vectors(normal1, v12))
            delta = min_angle_0_90(angle_deg_between_vectors(normal2, v12))

            _ = classify_pipi_geometry(theta_deg=theta, delta_deg=delta, gamma_deg=gamma)

            res_id = prot_ring["res_id"]
            key = f"{str(res_id)}__{ligand_index}"
            if key in used_protein_residue_set:
                continue

            used_protein_residue_set.add(key)

            result_edge_list.append(
                {
                    "node1_index": res_id,
                    "node1_type": "amino_acid",
                    "node2_index": ligand_index,
                    "node2_type": "substrate",
                    "interaction_type": "PIPISTACK",
                }
            )

    result_edge_list.sort(key=edge_sort_key)
    return result_edge_list


def calculate_pication_network(
    modeller: Modeller,
    ligand_mol_list: List[Chem.Mol],
    logger,
    ring_cation_distance_cutoff_A: float = 5.0,
    ring_cation_angle_cutoff_deg: float = 45.0,
    min_residue_index_gap: int = 3,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
) -> List[Dict[str, Any]] | None:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(ligand_mol_list, list):
        logger.print("[ERROR] ligand_mol_list must be a list.")
        return None

    if not isinstance(ring_cation_distance_cutoff_A, (int, float)) or float(ring_cation_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] ring_cation_distance_cutoff_A must be a positive number.")
        return None

    if not isinstance(ring_cation_angle_cutoff_deg, (int, float)) or float(ring_cation_angle_cutoff_deg) < 0.0:
        logger.print("[ERROR] ring_cation_angle_cutoff_deg must be a non-negative number.")
        return None

    if not isinstance(min_residue_index_gap, int) or int(min_residue_index_gap) < 0:
        logger.print("[ERROR] min_residue_index_gap must be a non-negative integer.")
        return None

    if not isinstance(docked_heavy_atom_distance_cutoff_A, (int, float)) or float(docked_heavy_atom_distance_cutoff_A) <= 0.0:
        logger.print("[ERROR] docked_heavy_atom_distance_cutoff_A must be a positive number.")
        return None

    try:
        topology, positions_nm, atoms, atom_index, coords_A_prot, bonded_prot = build_openmm_tables(modeller)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM tables from modeller: {e}")
        return None

    try:
        prot_ring_entries = collect_protein_pipi_rings(
            topology=topology,
            positions_nm=positions_nm,
            atom_index=atom_index,
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect protein aromatic rings for pi-cation: {e}")
        return None

    try:
        prot_cation_entries = collect_protein_pication_centers(
            topology=topology,
            positions_nm=positions_nm,
            atom_index=atom_index,
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect protein cation centers for pi-cation: {e}")
        return None

    result_edge_list: List[Dict[str, Any]] = []

    # -----------------------------
    # part 1) intra-protein pi-cation
    # -----------------------------
    if len(prot_ring_entries) > 0 and len(prot_cation_entries) > 0:
        used_ring_entries_per_pair: Set[Tuple[str, str, str, str]] = set()
        used_cation_entries_per_pair: Set[Tuple[str, str, str, str]] = set()

        try:
            pp_hits = find_pication_hits_between_ring_entries_and_cation_entries(
                ring_entries=prot_ring_entries,
                cation_entries=prot_cation_entries,
                ring_cation_distance_cutoff_A=float(ring_cation_distance_cutoff_A),
                ring_cation_angle_cutoff_deg=float(ring_cation_angle_cutoff_deg),
            )
        except Exception as e:
            logger.print(f"[ERROR] Failed to detect intra-protein pi-cation hits: {e}")
            return None

        for ring_entry, cation_entry in pp_hits:
            res1 = ring_entry["res_id"]
            res2 = cation_entry["res_id"]

            if res1 == res2:
                continue

            if isinstance(res1, int) and isinstance(res2, int):
                if abs(res1 - res2) < int(min_residue_index_gap):
                    continue

            if cation_entry.get("res_name") == "ARG" and "guan_norm" in cation_entry:
                gamma = angle_deg_between_vectors(ring_entry["normal"], cation_entry["guan_norm"])
                _ = classify_pication_arg_geometry(gamma_deg=float(gamma))

            pair_a, pair_b = edge_key_pair(res1, res2)
            pair_key = (str(pair_a), str(pair_b))

            ring_key = (
                pair_key[0],
                pair_key[1],
                str(res1),
                str(ring_entry.get("ring_id", "")),
            )
            cation_key = (
                pair_key[0],
                pair_key[1],
                str(res2),
                str(cation_entry.get("center_id", "")),
            )

            if ring_key in used_ring_entries_per_pair:
                continue
            if cation_key in used_cation_entries_per_pair:
                continue

            used_ring_entries_per_pair.add(ring_key)
            used_cation_entries_per_pair.add(cation_key)

            result_edge_list.append(
                {
                    "node1_index": pair_a,
                    "node1_type": "amino_acid",
                    "node2_index": pair_b,
                    "node2_type": "amino_acid",
                    "interaction_type": "PICATION",
                }
            )

    # -----------------------------
    # part 2) protein-ligand pi-cation
    # -----------------------------
    for ligand_index, lig_mol in enumerate(ligand_mol_list, start=1):
        if lig_mol is None:
            logger.print(f"[WARNING] Skip invalid ligand mol at index {ligand_index}.")
            continue

        if not is_protein_substrate_docked(
            modeller=modeller,
            ligand_mol=lig_mol,
            logger=logger,
            heavy_atom_distance_cutoff_A=float(docked_heavy_atom_distance_cutoff_A),
        ):
            logger.print(
                f"[WARNING] Substrate at index {ligand_index} is not spatially docked to protein. It is recommended to use 'enzywizard substrate' to generate substrate structures and 'enzywizard dock' to generate docked substrate structures.")

        try:
            lig_coords_A, bonded_lig = build_substrate_tables(lig_mol)
        except Exception as e:
            logger.print(f"[WARNING] Failed to build ligand tables at index {ligand_index}: {e}")
            continue

        try:
            lig_ring_entries = collect_substrate_pipi_rings(
                mol=lig_mol,
                coords_A=lig_coords_A,
            )
        except Exception as e:
            logger.print(f"[WARNING] Failed to collect ligand aromatic rings at index {ligand_index}: {e}")
            continue

        try:
            lig_cation_entries = collect_substrate_pication_centers(
                mol=lig_mol,
                coords_A=lig_coords_A,
            )
        except Exception as e:
            logger.print(f"[WARNING] Failed to collect ligand cation centers at index {ligand_index}: {e}")
            continue

        used_protein_residue_set: Set[str] = set()

        # protein aromatic ring vs ligand cation
        if len(prot_ring_entries) > 0 and len(lig_cation_entries) > 0:
            try:
                pr_lc_hits = find_pication_hits_between_ring_entries_and_cation_entries(
                    ring_entries=prot_ring_entries,
                    cation_entries=lig_cation_entries,
                    ring_cation_distance_cutoff_A=float(ring_cation_distance_cutoff_A),
                    ring_cation_angle_cutoff_deg=float(ring_cation_angle_cutoff_deg),
                )
            except Exception as e:
                logger.print(f"[WARNING] Failed to detect protein-ring/ligand-cation pi-cation hits at index {ligand_index}: {e}")
                pr_lc_hits = []

            for prot_ring, lig_cation in pr_lc_hits:
                res_id = prot_ring["res_id"]
                key = f"{str(res_id)}__{ligand_index}"
                if key in used_protein_residue_set:
                    continue

                used_protein_residue_set.add(key)

                result_edge_list.append(
                    {
                        "node1_index": res_id,
                        "node1_type": "amino_acid",
                        "node2_index": ligand_index,
                        "node2_type": "substrate",
                        "interaction_type": "PICATION",
                    }
                )

        # protein cation vs ligand aromatic ring
        if len(prot_cation_entries) > 0 and len(lig_ring_entries) > 0:
            try:
                lr_pc_hits = find_pication_hits_between_ring_entries_and_cation_entries(
                    ring_entries=lig_ring_entries,
                    cation_entries=prot_cation_entries,
                    ring_cation_distance_cutoff_A=float(ring_cation_distance_cutoff_A),
                    ring_cation_angle_cutoff_deg=float(ring_cation_angle_cutoff_deg),
                )
            except Exception as e:
                logger.print(f"[WARNING] Failed to detect ligand-ring/protein-cation pi-cation hits at index {ligand_index}: {e}")
                lr_pc_hits = []

            for lig_ring, prot_cation in lr_pc_hits:
                if prot_cation.get("res_name") == "ARG" and "guan_norm" in prot_cation:
                    gamma = angle_deg_between_vectors(lig_ring["normal"], prot_cation["guan_norm"])
                    _ = classify_pication_arg_geometry(gamma_deg=float(gamma))

                res_id = prot_cation["res_id"]
                key = f"{str(res_id)}__{ligand_index}"
                if key in used_protein_residue_set:
                    continue

                used_protein_residue_set.add(key)

                result_edge_list.append(
                    {
                        "node1_index": res_id,
                        "node1_type": "amino_acid",
                        "node2_index": ligand_index,
                        "node2_type": "substrate",
                        "interaction_type": "PICATION",
                    }
                )

    result_edge_list.sort(key=edge_sort_key)
    return result_edge_list


def calculate_disulfide_bond_network(
    modeller: Modeller,
    logger,
    ss_max_distance_A: float = 2.5,
    min_residue_index_gap: int = 3,
) -> List[Dict[str, Any]] | None:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(ss_max_distance_A, (int, float)) or float(ss_max_distance_A) <= 0.0:
        logger.print("[ERROR] ss_max_distance_A must be a positive number.")
        return None

    if not isinstance(min_residue_index_gap, int) or int(min_residue_index_gap) < 0:
        logger.print("[ERROR] min_residue_index_gap must be a non-negative integer.")
        return None

    try:
        topology, positions_nm, atoms, atom_index, coords_A, bonded = build_openmm_tables(modeller)
    except Exception as e:
        logger.print(f"[ERROR] Failed to build OpenMM tables from modeller: {e}")
        return None

    try:
        cys_sg_atoms = collect_protein_disulfide_sites(
            topology=topology,
            positions_nm=positions_nm,
            atom_index=atom_index,
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to collect protein disulfide sites: {e}")
        return None

    if len(cys_sg_atoms) < 2:
        return []

    try:
        result_edge_list = find_disulfide_bond_hits(
            cys_sg_atoms=cys_sg_atoms,
            ss_max_distance_A=float(ss_max_distance_A),
            min_residue_index_gap=int(min_residue_index_gap),
        )
    except Exception as e:
        logger.print(f"[ERROR] Failed to identify disulfide bonds: {e}")
        return None

    result_edge_list.sort(key=edge_sort_key)
    return result_edge_list

def calculate_all_interaction_network(
    modeller: Modeller,
    ligand_mol_list: List[Chem.Mol],
    substrate_name_list: List[str],
    struct: Structure,
    logger,
    bonded_h_min_distance_A: float = 0.8,
    bonded_h_max_distance_A: float = 1.3,
    da_max_distance_A: float = 3.9,
    ha_max_distance_A: float = 2.5,
    dha_min_angle_deg: float = 90.0,
    ionic_distance_cutoff_A: float = 4.0,
    mu: float = 0.01,
    ring_center_distance_cutoff_A: float = 6.5,
    ring_cation_distance_cutoff_A: float = 5.0,
    ring_cation_angle_cutoff_deg: float = 45.0,
    ss_max_distance_A: float = 2.5,
    docked_heavy_atom_distance_cutoff_A: float = 6.5,
    min_residue_index_gap: int = 3,
) -> List[Dict[str, Any]] | None:
    if not isinstance(modeller, Modeller):
        logger.print("[ERROR] modeller must be an OpenMM Modeller.")
        return None

    if not isinstance(ligand_mol_list, list):
        logger.print("[ERROR] ligand_mol_list must be a list.")
        return None

    if not isinstance(substrate_name_list, list):
        logger.print("[ERROR] substrate_name_list must be a list.")
        return None

    if len(ligand_mol_list) != len(substrate_name_list):
        logger.print("[ERROR] ligand_mol_list and substrate_name_list must have the same length.")
        return None

    if struct is None:
        logger.print("[ERROR] struct is None.")
        return None

    try:
        chain = get_single_chain(struct, logger)
        if chain is None:
            logger.print("[ERROR] Failed to get single chain from structure.")
            return None

        residue_info_list = get_residues_by_chain(chain, logger)
        if residue_info_list is None:
            logger.print("[ERROR] Failed to get residues from structure chain.")
            return None
    except Exception as e:
        logger.print(f"[ERROR] Failed to extract residue information from structure: {e}")
        return None

    aa_index_to_info: Dict[Any, Dict[str, Any]] = {}
    for item in residue_info_list:
        residue_id_tuple, resname, _ = item
        _, resseq, _ = residue_id_tuple

        aa_name = protein_letters_3to1.get(str(resname).capitalize(), "")
        aa_index_to_info[resseq] = {
            "aa_index": resseq,
            "aa_name": normalize_aa_name_to_one_letter(aa_name),
            "node_type": "amino_acid",
        }

    substrate_index_to_info: Dict[int, Dict[str, Any]] = {}
    for i, substrate_name in enumerate(substrate_name_list, start=1):
        substrate_index_to_info[i] = {
            "substrate_index": i,
            "substrate_name": str(substrate_name),
            "node_type": "substrate",
        }

    def build_node(node_index: Any, node_type: str) -> Dict[str, Any] | None:
        if node_type == "amino_acid":
            info = aa_index_to_info.get(node_index)
            if info is None:
                logger.print(f"[ERROR] Amino acid index {node_index} not found in structure residue mapping.")
                return None
            return dict(info)

        if node_type == "substrate":
            info = substrate_index_to_info.get(int(node_index))
            if info is None:
                logger.print(f"[ERROR] Substrate index {node_index} not found in substrate mapping.")
                return None
            return dict(info)

        logger.print(f"[ERROR] Unsupported node_type: {node_type}")
        return None

    all_raw_edges: List[Dict[str, Any]] = []

    hbond_edges = calculate_hydrogen_bond_network(
        modeller=modeller,
        ligand_mol_list=ligand_mol_list,
        logger=logger,
        bonded_h_min_distance_A=bonded_h_min_distance_A,
        bonded_h_max_distance_A=bonded_h_max_distance_A,
        da_max_distance_A=da_max_distance_A,
        ha_max_distance_A=ha_max_distance_A,
        dha_min_angle_deg=dha_min_angle_deg,
        docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
    )
    if hbond_edges is None:
        logger.print("[ERROR] Failed to calculate hydrogen bond network.")
        return None
    all_raw_edges.extend(hbond_edges)

    ionic_edges = calculate_ionic_bond_network(
        modeller=modeller,
        ligand_mol_list=ligand_mol_list,
        logger=logger,
        ionic_distance_cutoff_A=ionic_distance_cutoff_A,
        docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
        min_residue_index_gap=min_residue_index_gap,
    )
    if ionic_edges is None:
        logger.print("[ERROR] Failed to calculate ionic bond network.")
        return None
    all_raw_edges.extend(ionic_edges)

    vdw_edges = calculate_van_der_waals_network(
        modeller=modeller,
        ligand_mol_list=ligand_mol_list,
        logger=logger,
        mu=mu,
        docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
        min_residue_index_gap=min_residue_index_gap,
    )
    if vdw_edges is None:
        logger.print("[ERROR] Failed to calculate van der Waals network.")
        return None
    all_raw_edges.extend(vdw_edges)

    pipi_edges = calculate_pipi_stacking_network(
        modeller=modeller,
        ligand_mol_list=ligand_mol_list,
        logger=logger,
        ring_center_distance_cutoff_A=ring_center_distance_cutoff_A,
        min_residue_index_gap=min_residue_index_gap,
        docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
    )
    if pipi_edges is None:
        logger.print("[ERROR] Failed to calculate pi-pi stacking network.")
        return None
    all_raw_edges.extend(pipi_edges)

    pication_edges = calculate_pication_network(
        modeller=modeller,
        ligand_mol_list=ligand_mol_list,
        logger=logger,
        ring_cation_distance_cutoff_A=ring_cation_distance_cutoff_A,
        ring_cation_angle_cutoff_deg=ring_cation_angle_cutoff_deg,
        min_residue_index_gap=min_residue_index_gap,
        docked_heavy_atom_distance_cutoff_A=docked_heavy_atom_distance_cutoff_A,
    )
    if pication_edges is None:
        logger.print("[ERROR] Failed to calculate pi-cation network.")
        return None
    all_raw_edges.extend(pication_edges)

    ssbond_edges = calculate_disulfide_bond_network(
        modeller=modeller,
        logger=logger,
        ss_max_distance_A=ss_max_distance_A,
        min_residue_index_gap=min_residue_index_gap,
    )
    if ssbond_edges is None:
        logger.print("[ERROR] Failed to calculate disulfide bond network.")
        return None
    all_raw_edges.extend(ssbond_edges)

    merged_result: List[Dict[str, Any]] = []

    for edge in all_raw_edges:
        if not isinstance(edge, dict):
            logger.print("[ERROR] Invalid edge item in raw interaction result.")
            return None

        node1_index = edge.get("node1_index")
        node1_type = edge.get("node1_type")
        node2_index = edge.get("node2_index")
        node2_type = edge.get("node2_type")
        interaction_type = edge.get("interaction_type")

        node1_info = build_node(node1_index, node1_type)
        if node1_info is None:
            return None

        node2_info = build_node(node2_index, node2_type)
        if node2_info is None:
            return None

        merged_result.append(
            {
                "interaction": str(interaction_type),
                "node1": node1_info,
                "node2": node2_info,
            }
        )

    def interaction_sort_key(item: Dict[str, Any]):
        node1 = item.get("node1", {})
        node2 = item.get("node2", {})

        type1 = node1.get("node_type")
        type2 = node2.get("node_type")

        is_intra = (type1 == "amino_acid" and type2 == "amino_acid")

        aa_index = None
        if type1 == "amino_acid":
            aa_index = node1.get("aa_index")
        elif type2 == "amino_acid":
            aa_index = node2.get("aa_index")

        try:
            aa_index = int(aa_index)
        except Exception:
            aa_index = 999999

        return (
            0 if is_intra else 1,
            aa_index,
        )

    merged_result.sort(key=interaction_sort_key)

    return merged_result


def summarize_interaction_counts(interaction_list: List[Dict[str, Any]],logger) -> Dict[str, Any] | None:

    if not isinstance(interaction_list, list):
        logger.print("[ERROR] interaction_list must be a list.")
        return None

    interaction_types = ["HBOND", "IONIC", "VDW", "PIPISTACK", "PICATION", "SSBOND"]
    scope_types = ["overall", "intra_protein", "protein_substrate"]

    result: Dict[str, Dict[str, int]] = {
        scope: {interaction_type: 0 for interaction_type in interaction_types}
        for scope in scope_types
    }

    def get_scope(node1: Dict[str, Any], node2: Dict[str, Any]) -> str | None:
        node1_type = node1.get("node_type")
        node2_type = node2.get("node_type")

        if node1_type == "amino_acid" and node2_type == "amino_acid":
            return "intra_protein"

        if (
            (node1_type == "amino_acid" and node2_type == "substrate")
            or
            (node1_type == "substrate" and node2_type == "amino_acid")
        ):
            return "protein_substrate"

        return None

    for item in interaction_list:
        if not isinstance(item, dict):
            logger.print("[ERROR] Invalid item in interaction_list.")
            return None

        interaction = item.get("interaction")
        node1 = item.get("node1")
        node2 = item.get("node2")

        if interaction not in interaction_types:
            continue

        if not isinstance(node1, dict) or not isinstance(node2, dict):
            logger.print("[ERROR] Invalid node format in interaction item.")
            return None

        scope = get_scope(node1, node2)
        if scope is None:
            logger.print("[ERROR] Unsupported node type combination in interaction item.")
            return None

        result["overall"][interaction] += 1
        result[scope][interaction] += 1

    return result


def generate_interaction_report(
    interaction_list: List[Dict[str, Any]],
    interaction_statistics: Dict[str, Dict[str, int]],
) -> dict:

    raw_report = {
        "output_type": "enzywizard_interaction",
        "interactions": interaction_list,
        "interactions_statistics": interaction_statistics,
    }

    return postprocess_interaction_report_to_schema(raw_report)
