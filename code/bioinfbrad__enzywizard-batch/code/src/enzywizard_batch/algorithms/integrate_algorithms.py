from __future__ import annotations

from typing import Any, Dict, List, Tuple

from ..utils.logging_utils import Logger
from ..utils.integrate_utils import AA_CLASS_COUNT_FIELD_MAPPING,AA_NAME_COUNT_FIELD_MAPPING,AA_SS_COUNT_FIELD_MAPPING,INTERACTION_ORDER,add_count_statistics_to_overall,build_lookup_by_residue,build_protein_derived_statistics,get_clean_new_residue_list,get_disorder_membership_set,get_hydrophobic_cluster_membership_set,get_pocket_membership_set,interaction_one_hot,node_type_one_hot,residue_key,substrate_key
from ..utils.sequence_utils import normalize_aa_name_to_one_letter


def generate_integrate_report(overall_statistics: Dict[str, Any],integrated_graph: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "report_type": "enzywizard_integrate",
        "overall_statistics": overall_statistics,
        "integrated_graph": integrated_graph,
    }


def integrate_reports(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger,) -> Dict[str, Any] | None:
    if not isinstance(report_dict, dict):
        logger.print("[ERROR] report_dict must be a dict.")
        return None

    clean_report = report_dict.get("enzywizard_clean")
    if not isinstance(clean_report, dict):
        logger.print("[ERROR] enzywizard_clean report is required.")
        return None

    if strict:
        required_types = [
            "enzywizard_clean",
            "enzywizard_aaprops",
            "enzywizard_hydrocluster",
            "enzywizard_energy",
            "enzywizard_flexibility",
            "enzywizard_disorder",
            "enzywizard_conservation",
            "enzywizard_embedding",
            "enzywizard_pocket",
            "enzywizard_substrate",
            "enzywizard_dock",
            "enzywizard_interaction",
        ]
        for output_type in required_types:
            if output_type not in report_dict:
                logger.print(f"[ERROR] Missing required report in strict mode: {output_type}")
                return None

    overall_statistics = build_overall_statistics(report_dict, strict, logger)
    if overall_statistics is None:
        return None

    integrated_graph = build_integrated_graph(report_dict, strict, logger)
    if integrated_graph is None:
        return None

    return generate_integrate_report(overall_statistics, integrated_graph)


def build_overall_statistics(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger) -> Dict[str, Any] | None:
    overall_statistics: Dict[str, Any] = {}

    clean_report = report_dict.get("enzywizard_clean")
    aaprops_report = report_dict.get("enzywizard_aaprops")
    hydro_report = report_dict.get("enzywizard_hydrocluster")
    disorder_report = report_dict.get("enzywizard_disorder")
    pocket_report = report_dict.get("enzywizard_pocket")
    energy_report = report_dict.get("enzywizard_energy")
    dock_report = report_dict.get("enzywizard_dock")
    interaction_report = report_dict.get("enzywizard_interaction")

    protein_derived_statistics = build_protein_derived_statistics(clean_report, aaprops_report, logger)
    if protein_derived_statistics is None:
        return None
    overall_statistics.update(protein_derived_statistics)

    if (
        len(overall_statistics) == 0
        and
        aaprops_report is None
        and hydro_report is None
        and disorder_report is None
        and pocket_report is None
        and energy_report is None
        and dock_report is None
        and interaction_report is None
    ):
        if strict:
            logger.print("[ERROR] overall_statistics cannot be empty in strict mode.")
            return None
        return overall_statistics

    if aaprops_report is not None:
        stats = aaprops_report.get("residue_properties_statistics", {})

        if not add_count_statistics_to_overall(
            overall_statistics,
            stats.get("residue_name_statistics"),
            AA_NAME_COUNT_FIELD_MAPPING,
            "residue_name_statistics",
            logger,
        ):
            return None

        if not add_count_statistics_to_overall(
            overall_statistics,
            stats.get("residue_chemical_classification_statistics"),
            AA_CLASS_COUNT_FIELD_MAPPING,
            "residue_chemical_classification_statistics",
            logger,
        ):
            return None

        if not add_count_statistics_to_overall(
            overall_statistics,
            stats.get("residue_secondary_structure_statistics"),
            AA_SS_COUNT_FIELD_MAPPING,
            "residue_secondary_structure_statistics",
            logger,
        ):
            return None
    elif strict:
        logger.print("[ERROR] enzywizard_aaprops is required in strict mode.")
        return None

    if hydro_report is not None:
        hydro_stats = hydro_report.get("hydrophobic_cluster_statistics", {})
        overall_statistics["hydrophobic_cluster_count"] = hydro_stats["hydrophobic_cluster_count"]
        overall_statistics["max_hydrophobic_cluster_area"] = hydro_stats["max_hydrophobic_cluster_area"]
        overall_statistics["total_hydrophobic_cluster_area"] = hydro_stats["total_hydrophobic_cluster_area"]
    elif strict:
        logger.print("[ERROR] enzywizard_hydrocluster is required in strict mode.")
        return None

    if disorder_report is not None:
        disorder_stats = disorder_report.get("disordered_region_statistics", {})
        overall_statistics["disordered_region_count"] = disorder_stats["disordered_region_count"]
        overall_statistics["max_disordered_region_length"] = disorder_stats["max_disordered_region_length"]
        overall_statistics["total_disordered_region_length"] = disorder_stats["total_disordered_region_length"]
    elif strict:
        logger.print("[ERROR] enzywizard_disorder is required in strict mode.")
        return None

    if pocket_report is not None:
        pocket_stats = pocket_report.get("binding_pocket_statistics", {})
        overall_statistics["binding_pocket_count"] = pocket_stats["binding_pocket_count"]
        overall_statistics["max_binding_pocket_volume"] = pocket_stats["max_binding_pocket_volume"]
        overall_statistics["total_binding_pocket_volume"] = pocket_stats["total_binding_pocket_volume"]
    elif strict:
        logger.print("[ERROR] enzywizard_pocket is required in strict mode.")
        return None

    if energy_report is not None:
        energy_terms = energy_report.get("energy_terms", {})
        overall_statistics["total_potential_energy"] = energy_terms["total_potential_energy"]
        overall_statistics["harmonic_bond_potential_energy"] = energy_terms["harmonic_bond_potential_energy"]
        overall_statistics["harmonic_angle_potential_energy"] = energy_terms["harmonic_angle_potential_energy"]
        overall_statistics["custom_bond_potential_energy"] = energy_terms["custom_bond_potential_energy"]
        overall_statistics["custom_torsion_potential_energy"] = energy_terms["custom_torsion_potential_energy"]
        overall_statistics["custom_nonbonded_potential_energy"] = energy_terms["custom_nonbonded_potential_energy"]
        overall_statistics["nonbonded_potential_energy"] = energy_terms["nonbonded_potential_energy"]
        overall_statistics["periodic_torsion_potential_energy"] = energy_terms["periodic_torsion_potential_energy"]
        overall_statistics["cmap_torsion_potential_energy"] = energy_terms["cmap_torsion_potential_energy"]
    elif strict:
        logger.print("[ERROR] enzywizard_energy is required in strict mode.")
        return None

    if dock_report is not None:
        overall_statistics["enzyme_substrate_binding_affinity"] = dock_report["enzyme_substrate_docking_result"]["enzyme_substrate_binding_affinity"]
    elif strict:
        logger.print("[ERROR] enzywizard_dock is required in strict mode.")
        return None

    if interaction_report is not None:
        count_block = interaction_report["molecular_interaction_statistics"]["overall_molecular_interaction_statistics"]

        overall_statistics["hydrogen_bond_count"] = count_block["hydrogen_bond_count"]
        overall_statistics["ionic_bond_count"] = count_block["ionic_bond_count"]
        overall_statistics["van_der_waals_contact_count"] = count_block["van_der_waals_contact_count"]
        overall_statistics["pi_pi_stacking_count"] = count_block["pi_pi_stacking_count"]
        overall_statistics["pi_cation_interaction_count"] = count_block["pi_cation_interaction_count"]
        overall_statistics["disulfide_bond_count"] = count_block["disulfide_bond_count"]
    elif strict:
        logger.print("[ERROR] enzywizard_interaction is required in strict mode.")
        return None

    return overall_statistics


def build_integrated_graph(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger) -> List[Dict[str, Any]] | None:
    interaction_report = report_dict.get("enzywizard_interaction")
    if interaction_report is None:
        if strict:
            logger.print("[ERROR] enzywizard_interaction is required in strict mode.")
            return None
        return []

    interactions = interaction_report.get("molecular_interactions")
    if not isinstance(interactions, list):
        logger.print("[ERROR] molecular_interactions must be a list.")
        return None

    residue_node_list = build_residue_nodes(report_dict, strict, logger)
    if residue_node_list is None:
        return None

    substrate_node_list = build_substrate_nodes(report_dict, strict, logger)
    if substrate_node_list is None:
        return None

    all_node_list = residue_node_list + substrate_node_list

    for i, node in enumerate(all_node_list, start=1):
        node["node_index"] = i

    node_lookup = build_integrated_node_lookup(all_node_list, logger)
    if node_lookup is None:
        return None

    edge_entry_list = build_edge_entries_from_interactions(interactions, node_lookup, strict, logger)
    if edge_entry_list is None:
        return None

    connected_node_index_set = set()
    for edge_entry in edge_entry_list:
        source_node = edge_entry.get("source_node")
        target_node = edge_entry.get("target_node")

        if not isinstance(source_node, dict) or not isinstance(target_node, dict):
            logger.print("[ERROR] Invalid edge entry while collecting connected nodes.")
            return None

        connected_node_index_set.add(source_node["node_index"])
        connected_node_index_set.add(target_node["node_index"])

    isolated_node_entry_list: List[Dict[str, Any]] = []
    for node in all_node_list:
        if node["node_index"] not in connected_node_index_set:
            isolated_node_entry_list.append({
                "isolated_node": node
            })

    return edge_entry_list + isolated_node_entry_list


def build_residue_nodes(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger,) -> List[Dict[str, Any]] | None:
    clean_report = report_dict.get("enzywizard_clean")
    if clean_report is None:
        logger.print("[ERROR] Missing clean report.")
        return None

    aaprops_report = report_dict.get("enzywizard_aaprops")
    hydro_report = report_dict.get("enzywizard_hydrocluster")
    flexibility_report = report_dict.get("enzywizard_flexibility")
    disorder_report = report_dict.get("enzywizard_disorder")
    conservation_report = report_dict.get("enzywizard_conservation")
    embedding_report = report_dict.get("enzywizard_embedding")
    pocket_report = report_dict.get("enzywizard_pocket")

    clean_residue_list = get_clean_new_residue_list(clean_report, logger)
    if clean_residue_list is None:
        return None

    aaprops_lookup = None
    if aaprops_report is not None:
        aaprops_lookup = build_lookup_by_residue(aaprops_report["amino_acid_residue_properties"],"residue_index","residue_name",logger)
        if aaprops_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing aaprops report in strict mode.")
        return None

    flexibility_lookup = None
    if flexibility_report is not None:
        flexibility_lookup = build_lookup_by_residue(flexibility_report["protein_flexibility"],"residue_index","residue_name",logger)
        if flexibility_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing flexibility report in strict mode.")
        return None

    conservation_lookup = None
    if conservation_report is not None:
        conservation_lookup = build_lookup_by_residue(conservation_report["sequence_conservation_scores"],"residue_index","residue_name",logger)
        if conservation_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing conservation report in strict mode.")
        return None

    embedding_lookup = None
    if embedding_report is not None:
        embedding_lookup = build_lookup_by_residue(embedding_report["sequence_embeddings"],"residue_index","residue_name",logger)
        if embedding_lookup is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing embedding report in strict mode.")
        return None

    hydrophobic_membership = None
    if hydro_report is not None:
        hydrophobic_membership = get_hydrophobic_cluster_membership_set(hydro_report, logger)
        if hydrophobic_membership is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing hydrocluster report in strict mode.")
        return None

    disorder_membership = None
    if disorder_report is not None:
        disorder_membership = get_disorder_membership_set(disorder_report, logger)
        if disorder_membership is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing disorder report in strict mode.")
        return None

    pocket_membership = None
    if pocket_report is not None:
        pocket_membership = get_pocket_membership_set(pocket_report, logger)
        if pocket_membership is None:
            return None
    elif strict:
        logger.print("[ERROR] Missing pocket report in strict mode.")
        return None

    node_list: List[Dict[str, Any]] = []

    for residue in clean_residue_list:
        residue_index = residue["residue_index"]
        residue_name = normalize_aa_name_to_one_letter(residue["residue_name"])
        key = residue_key(residue_index, residue_name)

        drop_node = False

        node: Dict[str, Any] = {}
        node["node_index"] = 0
        node["node_type"] = "residue"
        node["node_type_one_hot_encoding"] = node_type_one_hot("residue")
        node["residue_index"] = residue_index
        node["residue_name"] = residue_name

        if aaprops_lookup is not None:
            aa_item = aaprops_lookup.get(key)
            if aa_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing amino_acid_residue_properties entry for residue {residue_index} {residue_name}.")
                    return None
                logger.print(f"[WARNING] Missing amino_acid_residue_properties entry for residue {residue_index} {residue_name}. Node dropped.")
                drop_node = True
            else:
                node["residue_alpha_carbon_coordinate"] = aa_item["residue_alpha_carbon_coordinate"]
                node["residue_chemical_classification"] = aa_item["residue_chemical_classification"]
                node["residue_chemical_classification_multi_hot_encoding"] = aa_item["residue_chemical_classification_multi_hot_encoding"]
                node["residue_secondary_structure"] = aa_item["residue_secondary_structure"]
                node["residue_secondary_structure_one_hot_encoding"] = aa_item["residue_secondary_structure_one_hot_encoding"]
                node["residue_relative_solvent_accessibility"] = aa_item["residue_relative_solvent_accessibility"]
                node["residue_backbone_phi_angle"] = aa_item["residue_backbone_phi_angle"]
                node["residue_backbone_psi_angle"] = aa_item["residue_backbone_psi_angle"]
                node["residue_net_charge"] = aa_item["residue_net_charge"]
                node["residue_pka"] = aa_item["residue_pka"]
                node["residue_volume"] = aa_item["residue_volume"]
                node["residue_hydrophobicity"] = aa_item["residue_hydrophobicity"]
                node["residue_molecular_weight"] = aa_item["residue_molecular_weight"]
                node["residue_isoelectric_point"] = aa_item["residue_isoelectric_point"]
                node["residue_name_one_hot_encoding"] = aa_item["residue_name_one_hot_encoding"]

        if drop_node:
            continue

        if hydrophobic_membership is not None:
            node["is_in_hydrophobic_cluster"] = (key in hydrophobic_membership)

        if flexibility_lookup is not None:
            rmsf_item = flexibility_lookup.get(key)
            if rmsf_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing residue_root_mean_square_fluctuation entry for residue {residue_index} {residue_name}.")
                    return None
                logger.print(f"[WARNING] Missing residue_root_mean_square_fluctuation entry for residue {residue_index} {residue_name}. Node dropped.")
                continue
            node["residue_root_mean_square_fluctuation"] = rmsf_item[ "residue_root_mean_square_fluctuation" ]

        if disorder_membership is not None:
            node["is_in_disordered_region"] = (key in disorder_membership)

        if conservation_lookup is not None:
            cons_item = conservation_lookup.get(key)
            if cons_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing residue_sequence_conservation_score entry for residue {residue_index} {residue_name}.")
                    return None
                logger.print(f"[WARNING] Missing residue_sequence_conservation_score entry for residue {residue_index} {residue_name}. Node dropped.")
                continue
            node["residue_sequence_conservation_score"] = cons_item[ "normalized_shannon_information_content" ]

        if embedding_lookup is not None:
            emb_item = embedding_lookup.get(key)
            if emb_item is None:
                if strict:
                    logger.print(f"[ERROR] Missing residue_embedding entry for residue {residue_index} {residue_name}.")
                    return None
                logger.print(f"[WARNING] Missing residue_embedding entry for residue {residue_index} {residue_name}. Node dropped.")
                continue
            node["residue_embedding"] = emb_item["residue_embedding"]

        if pocket_membership is not None:
            node["is_in_binding_pocket"] = (key in pocket_membership)

        node = reorder_residue_node(node)
        node_list.append(node)

    return node_list


def build_substrate_nodes(report_dict: Dict[str, Dict[str, Any]],strict: bool,logger: Logger) -> List[Dict[str, Any]] | None:
    substrate_report = report_dict.get("enzywizard_substrate")
    dock_report = report_dict.get("enzywizard_dock")

    if substrate_report is None or dock_report is None:
        if strict:
            logger.print("[ERROR] substrate and dock reports are required in strict mode.")
            return None
        return []

    substrate_lookup: Dict[str, Dict[str, Any]] = {}
    for item in substrate_report["substrates"]:
        substrate_name = item["substrate_name"]
        if substrate_name in substrate_lookup:
            logger.print(f"[ERROR] Duplicate substrate_name in substrate report: {substrate_name}")
            return None
        substrate_lookup[substrate_name] = item

    docked_substrate_list = dock_report["enzyme_substrate_docking_result"]["docked_substrates"]
    node_list: List[Dict[str, Any]] = []

    for substrate_index, dock_item in enumerate(docked_substrate_list, start=1):
        substrate_name = dock_item["docked_substrate_name"]

        substrate_item = substrate_lookup.get(substrate_name)
        if substrate_item is None:
            if strict:
                logger.print(f"[ERROR] Missing substrate entry for docked substrate: {substrate_name}")
                return None
            logger.print(f"[WARNING] Missing substrate entry for docked substrate: {substrate_name}. Node dropped.")
            continue

        docked_substrate_structure_name = dock_item["docked_substrate_structure_name"]
        structure_lookup: Dict[str, Dict[str, Any]] = {}
        for structure_item in substrate_item["substrate_possible_structures"]:
            structure_name = structure_item["substrate_structure_name"]
            if structure_name in structure_lookup:
                logger.print(f"[ERROR] Duplicate substrate_structure_name in substrate report: {structure_name}")
                return None
            structure_lookup[structure_name] = structure_item

        structure_item = structure_lookup.get(docked_substrate_structure_name)
        if structure_item is None:
            if strict:
                logger.print(f"[ERROR] Missing substrate structure entry for docked substrate structure: {docked_substrate_structure_name}")
                return None
            logger.print(f"[WARNING] Missing substrate structure entry for docked substrate structure: {docked_substrate_structure_name}. Node dropped.")
            continue

        node: Dict[str, Any] = {}
        node["node_index"] = 0
        node["node_type"] = "substrate"
        node["node_type_one_hot_encoding"] = node_type_one_hot("substrate")
        node["substrate_index"] = substrate_index
        node["substrate_name"] = substrate_name
        node["substrate_smiles"] = substrate_item["substrate_smiles"]
        node["substrate_atom_count"] = substrate_item["substrate_atom_count"]
        node["substrate_molecular_weight"] = substrate_item["substrate_molecular_weight"]
        node["substrate_logp"] = substrate_item["substrate_logp"]
        node["substrate_tpsa"] = substrate_item["substrate_tpsa"]
        node["substrate_heavy_atom_count"] = substrate_item["substrate_heavy_atom_count"]
        node["substrate_hbond_donor_count"] = substrate_item["substrate_hbond_donor_count"]
        node["substrate_hbond_acceptor_count"] = substrate_item["substrate_hbond_acceptor_count"]
        node["substrate_rotatable_bond_count"] = substrate_item["substrate_rotatable_bond_count"]
        node["substrate_molar_refractivity"] = substrate_item["substrate_molar_refractivity"]
        node["substrate_structure_energy"] = structure_item["substrate_structure_energy"]
        node["substrate_structure_max_3d_diameter"] = structure_item["substrate_structure_max_3d_diameter"]
        node["substrate_structure_mean_pairwise_atom_distance"] = structure_item["substrate_structure_mean_pairwise_atom_distance"]
        node["substrate_structure_std_pairwise_atom_distance"] = structure_item["substrate_structure_std_pairwise_atom_distance"]
        node["substrate_structure_asphericity"] = structure_item["substrate_structure_asphericity"]
        node["substrate_structure_spherocity"] = structure_item["substrate_structure_spherocity"]
        node["substrate_structure_principal_moment_ratio"] = structure_item["substrate_structure_principal_moment_ratio"]
        node["substrate_structure_radius_of_gyration"] = structure_item["substrate_structure_radius_of_gyration"]
        node["docked_substrate_center_coordinate"] = dock_item["docked_substrate_center_coordinate"]
        node["substrate_fingerprint_encoding"] = substrate_item["substrate_fingerprint_encoding"]

        node = reorder_substrate_node(node)
        node_list.append(node)

    return node_list


def build_integrated_node_lookup(all_node_list: List[Dict[str, Any]],logger: Logger,) -> Dict[Tuple[Any, ...], Dict[str, Any]] | None:
    lookup: Dict[Tuple[Any, ...], Dict[str, Any]] = {}

    for node in all_node_list:
        node_type = node.get("node_type")
        if node_type == "residue":
            key = ("residue", node["residue_index"], node["residue_name"])
        elif node_type == "substrate":
            key = ("substrate", node["substrate_name"])
        else:
            logger.print("[ERROR] Unknown node_type while building node lookup.")
            return None

        if key in lookup:
            logger.print(f"[ERROR] Duplicate integrated node key: {key}")
            return None

        lookup[key] = node

    return lookup


def build_edge_entries_from_interactions(
    interactions: List[Dict[str, Any]],
    node_lookup: Dict[Tuple[Any, ...], Dict[str, Any]],
    strict: bool,
    logger: Logger
) -> List[Dict[str, Any]] | None:
    merged_edge_count: Dict[Tuple[str, int, int], int] = {}
    merged_edge_nodes: Dict[Tuple[str, int, int], Tuple[Dict[str, Any], Dict[str, Any]]] = {}

    for item in interactions:
        interaction = item["molecular_interaction_type"]
        node1 = item["source_node"]
        node2 = item["target_node"]

        integrated_node_1 = resolve_interaction_node_to_integrated_node(node1, node_lookup, logger)
        if integrated_node_1 is None:
            if strict:
                logger.print("[ERROR] Failed to resolve interaction node1 in strict mode.")
                return None
            logger.print("[WARNING] Interaction node1 cannot be resolved. Edge skipped.")
            continue

        integrated_node_2 = resolve_interaction_node_to_integrated_node(node2, node_lookup, logger)
        if integrated_node_2 is None:
            if strict:
                logger.print("[ERROR] Failed to resolve interaction node2 in strict mode.")
                return None
            logger.print("[WARNING] Interaction node2 cannot be resolved. Edge skipped.")
            continue

        node_index_1 = integrated_node_1["node_index"]
        node_index_2 = integrated_node_2["node_index"]

        if node_index_1 <= node_index_2:
            left_node = integrated_node_1
            right_node = integrated_node_2
            left_index = node_index_1
            right_index = node_index_2
        else:
            left_node = integrated_node_2
            right_node = integrated_node_1
            left_index = node_index_2
            right_index = node_index_1

        merge_key = (interaction, left_index, right_index)

        merged_edge_count[merge_key] = merged_edge_count.get(merge_key, 0) + 1
        merged_edge_nodes[merge_key] = (left_node, right_node)

    edge_entry_list: List[Dict[str, Any]] = []

    for interaction, node_index_1, node_index_2 in sorted(merged_edge_count.keys(),key=lambda x: (x[1], x[2], INTERACTION_ORDER.index(x[0]))):
        source_node, target_node = merged_edge_nodes[(interaction, node_index_1, node_index_2)]

        molecular_interaction = {
            "molecular_interaction_type": interaction,
            "molecular_interaction_one_hot_encoding": interaction_one_hot(interaction),
            "interaction_count": merged_edge_count[(interaction, node_index_1, node_index_2)],
        }

        edge_entry_list.append({
            "molecular_interaction": molecular_interaction,
            "source_node": source_node,
            "target_node": target_node,
        })

    return edge_entry_list


def resolve_interaction_node_to_integrated_node(interaction_node: Dict[str, Any],node_lookup: Dict[Tuple[Any, ...], Dict[str, Any]],logger: Logger,) -> Dict[str, Any] | None:
    node_type = interaction_node.get("node_type")

    if node_type == "residue":
        residue_index = interaction_node.get("residue_index")
        residue_name = interaction_node.get("residue_name")

        if not isinstance(residue_index, int) or not isinstance(residue_name, str):
            logger.print("[ERROR] Invalid residue interaction node.")
            return None

        residue_name_norm = normalize_aa_name_to_one_letter(residue_name)

        return node_lookup.get(("residue", residue_index, residue_name_norm))

    if node_type == "substrate":
        substrate_name = interaction_node.get("substrate_name")
        if not isinstance(substrate_name, str) or substrate_name.strip() == "":
            logger.print("[ERROR] Invalid substrate interaction node.")
            return None

        substrate_name_clean = substrate_name.strip()

        node = node_lookup.get(("substrate", substrate_name_clean))
        if node is not None:
            return node

        if substrate_name_clean.startswith("docked_"):
            substrate_name_clean = substrate_name_clean[len("docked_"):]

        node = node_lookup.get(("substrate", substrate_name_clean))
        if node is not None:
            return node

        logger.print(f"[WARNING] Failed to match substrate node by substrate_name: {substrate_name}")
        return None


def reorder_residue_node(node: Dict[str, Any]) -> Dict[str, Any]:
    ordered: Dict[str, Any] = {}

    ordered["node_index"] = node["node_index"]
    ordered["node_type"] = node["node_type"]
    ordered["node_type_one_hot_encoding"] = node["node_type_one_hot_encoding"]

    ordered["residue_index"] = node["residue_index"]
    ordered["residue_name"] = node["residue_name"]

    optional_fields = [
        "residue_name_one_hot_encoding",
        "residue_alpha_carbon_coordinate",
        "residue_chemical_classification",
        "residue_chemical_classification_multi_hot_encoding",
        "residue_secondary_structure",
        "residue_secondary_structure_one_hot_encoding",
        "residue_relative_solvent_accessibility",
        "residue_backbone_phi_angle",
        "residue_backbone_psi_angle",
        "residue_net_charge",
        "residue_pka",
        "residue_volume",
        "residue_hydrophobicity",
        "residue_molecular_weight",
        "residue_isoelectric_point",
        "residue_root_mean_square_fluctuation",
        "residue_sequence_conservation_score",
        "residue_embedding",
        "is_in_hydrophobic_cluster",
        "is_in_disordered_region",
        "is_in_binding_pocket",
    ]

    for field in optional_fields:
        if field in node:
            ordered[field] = node[field]

    return ordered


def reorder_substrate_node(node: Dict[str, Any]) -> Dict[str, Any]:
    ordered: Dict[str, Any] = {}

    ordered["node_index"] = node["node_index"]
    ordered["node_type"] = node["node_type"]
    ordered["node_type_one_hot_encoding"] = node["node_type_one_hot_encoding"]

    ordered["substrate_index"] = node["substrate_index"]
    ordered["substrate_name"] = node["substrate_name"]

    optional_fields = [
        "substrate_smiles",
        "substrate_atom_count",
        "substrate_molecular_weight",
        "substrate_logp",
        "substrate_tpsa",
        "substrate_heavy_atom_count",
        "substrate_hbond_donor_count",
        "substrate_hbond_acceptor_count",
        "substrate_rotatable_bond_count",
        "substrate_molar_refractivity",
        "substrate_structure_energy",
        "substrate_structure_max_3d_diameter",
        "substrate_structure_mean_pairwise_atom_distance",
        "substrate_structure_std_pairwise_atom_distance",
        "substrate_structure_asphericity",
        "substrate_structure_spherocity",
        "substrate_structure_principal_moment_ratio",
        "substrate_structure_radius_of_gyration",
        "docked_substrate_center_coordinate",
        "substrate_fingerprint_encoding",
    ]

    for field in optional_fields:
        if field in node:
            ordered[field] = node[field]

    return ordered
