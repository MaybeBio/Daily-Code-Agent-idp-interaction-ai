from Bio.PDB.DSSP import DSSP
from typing import Tuple, Any, Dict
from ..resources.aa_resources import DSSP_8STATE
from ..utils.logging_utils import Logger


def get_dssp_fields(dssp: DSSP, dssp_key: Tuple[str, Tuple[str, int, str]], logger: Logger) -> Tuple[str, float, float, float] | None:

    if dssp_key not in dssp:
        logger.print(f"[ERROR] DSSP key {dssp_key} not found.")
        return None

    dssp_res = dssp[dssp_key]

    ss = dssp_res[2]
    rsa = dssp_res[3]
    phi = dssp_res[4]
    psi = dssp_res[5]

    ss = (ss or "-").strip()
    if ss == "" or ss not in DSSP_8STATE:
        ss = "-"

    try:
        return ss, float(rsa), float(phi), float(psi)
    except (TypeError, ValueError) as e:
        logger.print(f"[ERROR] Invalid DSSP numeric values at {dssp_key}: rsa={rsa}, phi={phi}, psi={psi}: {e}")
        return None



def postprocess_aaprops_report_to_schema(raw_report: Dict[str, Any]) -> Dict[str, Any]:
    """
    Map the original EnzyWizard-Aaprops raw report to the new schema-compliant report.
    """

    aa_name_count_field_map: Dict[str, str] = {
        "A": "alanine_count",
        "C": "cysteine_count",
        "D": "aspartic_acid_count",
        "E": "glutamic_acid_count",
        "F": "phenylalanine_count",
        "G": "glycine_count",
        "H": "histidine_count",
        "I": "isoleucine_count",
        "K": "lysine_count",
        "L": "leucine_count",
        "M": "methionine_count",
        "N": "asparagine_count",
        "P": "proline_count",
        "Q": "glutamine_count",
        "R": "arginine_count",
        "S": "serine_count",
        "T": "threonine_count",
        "V": "valine_count",
        "W": "tryptophan_count",
        "Y": "tyrosine_count",
    }

    aa_class_count_field_map: Dict[str, str] = {
        "uncharged_polar": "uncharged_polar_count",
        "positively_charged": "positively_charged_count",
        "negatively_charged": "negatively_charged_count",
        "hydrophobic": "hydrophobic_count",
        "aromatic": "aromatic_count",
        "aliphatic": "aliphatic_count",
        "heterocyclic": "heterocyclic_count",
        "sulfur_containing": "sulfur_containing_count",
    }

    aa_ss_count_field_map: Dict[str, str] = {
        "-": "unassigned_count",
        "H": "alpha_helix_count",
        "B": "beta_bridge_count",
        "E": "extended_strand_count",
        "G": "three_ten_helix_count",
        "I": "pi_helix_count",
        "T": "turn_count",
        "S": "bend_count",
    }

    raw_aa_props = raw_report.get("aa_props", [])
    raw_statistics = raw_report.get("aa_props_statistics", {})

    raw_aa_name_statistics = raw_statistics.get("aa_name_statistics", {})
    raw_aa_class_statistics = raw_statistics.get("aa_class_statistics", {})
    raw_aa_ss_statistics = raw_statistics.get("aa_ss_statistics", {})

    amino_acid_residue_properties = []

    for raw_residue_property in raw_aa_props:
        amino_acid_residue_properties.append(
            {
                "residue_index": raw_residue_property.get("aa_id"),
                "residue_name": raw_residue_property.get("aa_name"),
                "residue_name_one_hot_encoding": raw_residue_property.get("aa_name_one_hot"),
                "residue_chemical_classification": raw_residue_property.get("aa_class"),
                "residue_chemical_classification_multi_hot_encoding": raw_residue_property.get("aa_class_one_hot"),
                "residue_secondary_structure": raw_residue_property.get("aa_ss"),
                "residue_secondary_structure_one_hot_encoding": raw_residue_property.get("aa_ss_one_hot"),
                "residue_relative_solvent_accessibility": raw_residue_property.get("aa_rsa"),
                "residue_backbone_phi_angle": raw_residue_property.get("aa_phi"),
                "residue_backbone_psi_angle": raw_residue_property.get("aa_psi"),
                "residue_net_charge": raw_residue_property.get("aa_net_charge"),
                "residue_pka": raw_residue_property.get("aa_pka"),
                "residue_volume": raw_residue_property.get("aa_volume"),
                "residue_hydrophobicity": raw_residue_property.get("aa_hydrophobicity"),
                "residue_molecular_weight": raw_residue_property.get("aa_molecular_weight"),
                "residue_isoelectric_point": raw_residue_property.get("aa_pi"),
                "residue_alpha_carbon_coordinate": raw_residue_property.get("aa_coord"),
            }
        )

    residue_name_statistics = {
        schema_field_name: raw_aa_name_statistics[raw_field_name]
        for raw_field_name, schema_field_name in aa_name_count_field_map.items()
        if raw_field_name in raw_aa_name_statistics
    }

    residue_chemical_classification_statistics = {
        schema_field_name: raw_aa_class_statistics[raw_field_name]
        for raw_field_name, schema_field_name in aa_class_count_field_map.items()
        if raw_field_name in raw_aa_class_statistics
    }

    residue_secondary_structure_statistics = {
        schema_field_name: raw_aa_ss_statistics[raw_field_name]
        for raw_field_name, schema_field_name in aa_ss_count_field_map.items()
        if raw_field_name in raw_aa_ss_statistics
    }

    return {
        "report_type": raw_report.get("output_type"),
        "amino_acid_residue_properties": amino_acid_residue_properties,
        "residue_properties_statistics": {
            "residue_name_statistics": residue_name_statistics,
            "residue_chemical_classification_statistics": residue_chemical_classification_statistics,
            "residue_secondary_structure_statistics": residue_secondary_structure_statistics,
        },
    }
