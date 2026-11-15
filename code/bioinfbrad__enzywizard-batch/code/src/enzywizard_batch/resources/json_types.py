from __future__ import annotations

from typing import Literal, NotRequired, Required, TypeAlias, TypedDict


CLEAN_TYPE = "enzywizard_clean"
MUTCLEAN_TYPE = "enzywizard_mut_clean"
AAPROPS_TYPE = "enzywizard_aaprops"
HYDROCLUSTER_TYPE = "enzywizard_hydrocluster"
ENERGY_TYPE = "enzywizard_energy"
FLEXIBILITY_TYPE = "enzywizard_flexibility"
DISORDER_TYPE = "enzywizard_disorder"
CONSERVATION_TYPE = "enzywizard_conservation"
EMBEDDING_TYPE = "enzywizard_embedding"
POCKET_TYPE = "enzywizard_pocket"
SUBSTRATE_TYPE = "enzywizard_substrate"
DOCK_TYPE = "enzywizard_dock"
INTERACTION_TYPE = "enzywizard_interaction"
INTEGRATE_TYPE = "enzywizard_integrate"
MUT_INTEGRATE_TYPE = "enzywizard_mut_integrate"


ResidueName: TypeAlias = Literal[
    "A", "C", "D", "E", "F",
    "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R",
    "S", "T", "V", "W", "Y",
]

ResidueSecondaryStructure: TypeAlias = Literal[
    "-",
    "H",
    "B",
    "E",
    "G",
    "I",
    "T",
    "S",
]

MolecularInteractionType: TypeAlias = Literal[
    "HBOND",
    "IONIC",
    "VDW",
    "PIPISTACK",
    "PICATION",
    "SSBOND",
]


'''
common
'''


class ResidueIndexName(TypedDict):
    residue_index: int
    residue_name: str


class CleanResidueInfo(TypedDict):
    residue_index: int
    residue_name: str
    hydrogen_atom_count: int


class ResidueMappingOldToNew(TypedDict):
    old_residue: CleanResidueInfo
    new_residue: CleanResidueInfo


class CleanStatistics(TypedDict):
    removed_heterogen_count: int
    standardized_residue_name_count: int
    repaired_residue_count: int
    added_heavy_atom_count: int
    added_hydrogen_atom_count: int
    retained_residue_count: int


'''
clean_report
'''


class EnzyWizardCleanOutput(TypedDict):
    report_type: Literal["enzywizard_clean"]
    residue_mapping_old_to_new: list[ResidueMappingOldToNew]
    clean_statistics: CleanStatistics


'''
mutclean_report
'''


class EnzyWizardMutCleanOutput(TypedDict):
    report_type: Literal["enzywizard_mut_clean"]

    amino_acid_substitution: str
    cleaned_amino_acid_substitution: str

    wild_type_residue_mapping_old_to_new: list[ResidueMappingOldToNew]
    wild_type_clean_statistics: CleanStatistics

    mutant_residue_mapping_old_to_new: list[ResidueMappingOldToNew]
    mutant_clean_statistics: CleanStatistics


'''
aaprops_report
'''


class AminoAcidResidueProperty(TypedDict):
    residue_index: int
    residue_name: str

    residue_name_one_hot_encoding: list[int]

    residue_chemical_classification: str
    residue_chemical_classification_multi_hot_encoding: list[int]

    residue_secondary_structure: str
    residue_secondary_structure_one_hot_encoding: list[int]

    residue_relative_solvent_accessibility: float
    residue_backbone_phi_angle: float
    residue_backbone_psi_angle: float

    residue_net_charge: float
    residue_pka: float
    residue_volume: float
    residue_hydrophobicity: float
    residue_molecular_weight: float
    residue_isoelectric_point: float

    residue_alpha_carbon_coordinate: list[float]


class ResidueNameStatistics(TypedDict):
    alanine_count: int
    cysteine_count: int
    aspartic_acid_count: int
    glutamic_acid_count: int
    phenylalanine_count: int
    glycine_count: int
    histidine_count: int
    isoleucine_count: int
    lysine_count: int
    leucine_count: int
    methionine_count: int
    asparagine_count: int
    proline_count: int
    glutamine_count: int
    arginine_count: int
    serine_count: int
    threonine_count: int
    valine_count: int
    tryptophan_count: int
    tyrosine_count: int


class ResidueChemicalClassificationStatistics(TypedDict):
    uncharged_polar_count: int
    positively_charged_count: int
    negatively_charged_count: int
    hydrophobic_count: int
    aromatic_count: int
    aliphatic_count: int
    heterocyclic_count: int
    sulfur_containing_count: int


class ResidueSecondaryStructureStatistics(TypedDict):
    unassigned_count: int
    alpha_helix_count: int
    beta_bridge_count: int
    extended_strand_count: int
    three_ten_helix_count: int
    pi_helix_count: int
    turn_count: int
    bend_count: int


class ResiduePropertiesStatistics(TypedDict):
    residue_name_statistics: ResidueNameStatistics
    residue_chemical_classification_statistics: ResidueChemicalClassificationStatistics
    residue_secondary_structure_statistics: ResidueSecondaryStructureStatistics


class EnzyWizardAapropsOutput(TypedDict):
    report_type: Literal["enzywizard_aaprops"]
    amino_acid_residue_properties: list[AminoAcidResidueProperty]
    residue_properties_statistics: ResiduePropertiesStatistics


'''
hydrocluster_report
'''


class HydrophobicClusterEntry(TypedDict):
    hydrophobic_cluster_area: float
    residues: list[ResidueIndexName]


class HydrophobicClusterStatistics(TypedDict):
    hydrophobic_cluster_count: int
    max_hydrophobic_cluster_area: float
    total_hydrophobic_cluster_area: float


class EnzyWizardHydroclusterOutput(TypedDict):
    report_type: Literal["enzywizard_hydrocluster"]
    hydrophobic_cluster_statistics: HydrophobicClusterStatistics
    hydrophobic_clusters: list[HydrophobicClusterEntry]


'''
energy_report
'''


class EnergyTerms(TypedDict):
    total_potential_energy: float
    harmonic_bond_potential_energy: float
    harmonic_angle_potential_energy: float
    custom_bond_potential_energy: float
    custom_torsion_potential_energy: float
    custom_nonbonded_potential_energy: float
    nonbonded_potential_energy: float
    periodic_torsion_potential_energy: float
    cmap_torsion_potential_energy: float


class EnzyWizardEnergyOutput(TypedDict):
    report_type: Literal["enzywizard_energy"]
    energy_terms: EnergyTerms


'''
flexibility_report
'''


class ProteinFlexibilityEntry(TypedDict):
    residue_index: int
    residue_name: str
    residue_root_mean_square_fluctuation: float


class EnzyWizardFlexibilityOutput(TypedDict):
    report_type: Literal["enzywizard_flexibility"]
    protein_flexibility: list[ProteinFlexibilityEntry]


'''
disorder_report
'''


class DisorderedRegionEntry(TypedDict):
    disordered_region_length: int
    residues: list[ResidueIndexName]


class DisorderedRegionStatistics(TypedDict):
    disordered_region_count: int
    max_disordered_region_length: int
    total_disordered_region_length: int


class EnzyWizardDisorderOutput(TypedDict):
    report_type: Literal["enzywizard_disorder"]
    disordered_region_statistics: DisorderedRegionStatistics
    disordered_regions: list[DisorderedRegionEntry]


'''
conservation_report
'''


class SequenceConservationScoreEntry(TypedDict):
    residue_index: int
    residue_name: str
    hmm_profile_raw_score: float
    normalized_emission_probability: float
    normalized_shannon_information_content: float


class EnzyWizardConservationOutput(TypedDict):
    report_type: Literal["enzywizard_conservation"]
    sequence_conservation_scores: list[SequenceConservationScoreEntry]


'''
embedding_report
'''


class SequenceEmbeddingEntry(TypedDict):
    residue_index: int
    residue_name: str
    residue_embedding: list[float]


class EnzyWizardEmbeddingOutput(TypedDict):
    report_type: Literal["enzywizard_embedding"]
    sequence_embeddings: list[SequenceEmbeddingEntry]


'''
pocket_report
'''


class BindingPocketEntry(TypedDict):
    binding_pocket_volume: float
    binding_pocket_sphere_count: int
    residues: list[ResidueIndexName]
    binding_pocket_center_coordinate: list[float]
    binding_pocket_box_size: list[float]


class BindingPocketStatistics(TypedDict):
    binding_pocket_count: int
    max_binding_pocket_volume: float
    total_binding_pocket_volume: float


class EnzyWizardPocketOutput(TypedDict):
    report_type: Literal["enzywizard_pocket"]
    binding_pocket_statistics: BindingPocketStatistics
    binding_pockets: list[BindingPocketEntry]


'''
substrate_report
'''


class SubstratePossibleStructureEntry(TypedDict):
    substrate_structure_name: str
    substrate_structure_energy: float
    substrate_structure_max_3d_diameter: float
    substrate_structure_mean_pairwise_atom_distance: float
    substrate_structure_std_pairwise_atom_distance: float
    substrate_structure_asphericity: float
    substrate_structure_spherocity: float
    substrate_structure_principal_moment_ratio: float
    substrate_structure_radius_of_gyration: float


class SubstrateEntry(TypedDict):
    substrate_name: str
    substrate_smiles: str
    substrate_fingerprint_encoding: list[Literal[0, 1]]

    substrate_atom_count: int
    substrate_molecular_weight: float
    substrate_logp: float
    substrate_tpsa: float
    substrate_heavy_atom_count: int
    substrate_hbond_donor_count: int
    substrate_hbond_acceptor_count: int
    substrate_rotatable_bond_count: int
    substrate_molar_refractivity: float

    substrate_possible_structures: list[SubstratePossibleStructureEntry]


class EnzyWizardSubstrateOutput(TypedDict):
    report_type: Literal["enzywizard_substrate"]
    substrates: list[SubstrateEntry]


'''
dock_report
'''


class DockedSubstrateEntry(TypedDict):
    docked_substrate_name: str
    docked_substrate_structure_name: str
    docked_substrate_center_coordinate: list[float]


class EnzymeSubstrateDockingResult(TypedDict):
    enzyme_substrate_complex_name: str
    enzyme_substrate_binding_affinity: float
    docked_substrate_names: str
    docking_box_center_coordinate: list[float]
    docking_box_size: list[float]
    docked_substrates: list[DockedSubstrateEntry]


class EnzyWizardDockOutput(TypedDict):
    report_type: Literal["enzywizard_dock"]
    enzyme_substrate_docking_result: EnzymeSubstrateDockingResult


'''
interaction_report
'''


class ResidueInteractionNode(TypedDict):
    residue_index: int
    residue_name: str
    node_type: Literal["residue"]


class SubstrateInteractionNode(TypedDict):
    substrate_index: int
    substrate_name: str
    node_type: Literal["substrate"]


InteractionNode: TypeAlias = ResidueInteractionNode | SubstrateInteractionNode


class MolecularInteractionEntry(TypedDict):
    molecular_interaction_type: MolecularInteractionType
    source_node: InteractionNode
    target_node: InteractionNode


class MolecularInteractionTypeCount(TypedDict):
    hydrogen_bond_count: int
    ionic_bond_count: int
    van_der_waals_contact_count: int
    pi_pi_stacking_count: int
    pi_cation_interaction_count: int
    disulfide_bond_count: int


class MolecularInteractionStatistics(TypedDict):
    overall_molecular_interaction_statistics: MolecularInteractionTypeCount
    intra_enzyme_interaction_statistics: MolecularInteractionTypeCount
    enzyme_substrate_interaction_statistics: MolecularInteractionTypeCount


class EnzyWizardInteractionOutput(TypedDict):
    report_type: Literal["enzywizard_interaction"]
    molecular_interactions: list[MolecularInteractionEntry]
    molecular_interaction_statistics: MolecularInteractionStatistics


'''
integrate_report
'''


class IntegratedOverallStatistics(TypedDict, total=False):
    sequence_length: Required[int]
    total_molecular_weight: Required[float]
    total_net_charge: Required[float]
    total_residue_volume: Required[float]
    max_3d_diameter: Required[float]
    radius_of_gyration: Required[float]
    asphericity: float
    spherocity: float
    principal_moment_ratio: float
    bounding_box_volume: Required[float]
    mean_pairwise_ca_distance: Required[float]
    std_pairwise_ca_distance: Required[float]

    residue_name_alanine_count: Required[int]
    residue_name_cysteine_count: Required[int]
    residue_name_aspartic_acid_count: Required[int]
    residue_name_glutamic_acid_count: Required[int]
    residue_name_phenylalanine_count: Required[int]
    residue_name_glycine_count: Required[int]
    residue_name_histidine_count: Required[int]
    residue_name_isoleucine_count: Required[int]
    residue_name_lysine_count: Required[int]
    residue_name_leucine_count: Required[int]
    residue_name_methionine_count: Required[int]
    residue_name_asparagine_count: Required[int]
    residue_name_proline_count: Required[int]
    residue_name_glutamine_count: Required[int]
    residue_name_arginine_count: Required[int]
    residue_name_serine_count: Required[int]
    residue_name_threonine_count: Required[int]
    residue_name_valine_count: Required[int]
    residue_name_tryptophan_count: Required[int]
    residue_name_tyrosine_count: Required[int]

    residue_chemical_classification_uncharged_polar_count: Required[int]
    residue_chemical_classification_positively_charged_count: Required[int]
    residue_chemical_classification_negatively_charged_count: Required[int]
    residue_chemical_classification_hydrophobic_count: Required[int]
    residue_chemical_classification_aromatic_count: Required[int]
    residue_chemical_classification_aliphatic_count: Required[int]
    residue_chemical_classification_heterocyclic_count: Required[int]
    residue_chemical_classification_sulfur_containing_count: Required[int]

    secondary_structure_unassigned_count: Required[int]
    secondary_structure_alpha_helix_count: Required[int]
    secondary_structure_beta_bridge_count: Required[int]
    secondary_structure_extended_strand_count: Required[int]
    secondary_structure_three_ten_helix_count: Required[int]
    secondary_structure_pi_helix_count: Required[int]
    secondary_structure_turn_count: Required[int]
    secondary_structure_bend_count: Required[int]

    hydrophobic_cluster_count: Required[int]
    max_hydrophobic_cluster_area: Required[float]
    total_hydrophobic_cluster_area: Required[float]

    disordered_region_count: Required[int]
    max_disordered_region_length: Required[int]
    total_disordered_region_length: Required[int]

    binding_pocket_count: Required[int]
    max_binding_pocket_volume: Required[float]
    total_binding_pocket_volume: Required[float]

    total_potential_energy: Required[float]
    harmonic_bond_potential_energy: Required[float]
    harmonic_angle_potential_energy: Required[float]
    custom_bond_potential_energy: Required[float]
    custom_torsion_potential_energy: Required[float]
    custom_nonbonded_potential_energy: Required[float]
    nonbonded_potential_energy: Required[float]
    periodic_torsion_potential_energy: Required[float]
    cmap_torsion_potential_energy: Required[float]

    enzyme_substrate_binding_affinity: float

    hydrogen_bond_count: Required[int]
    ionic_bond_count: Required[int]
    van_der_waals_contact_count: Required[int]
    pi_pi_stacking_count: Required[int]
    pi_cation_interaction_count: Required[int]
    disulfide_bond_count: Required[int]


class IntegratedResidueNode(TypedDict):
    node_index: Required[int]
    node_type: Required[Literal["residue"]]
    node_type_one_hot_encoding: Required[list[Literal[0, 1]]]

    residue_index: Required[int]
    residue_name: Required[str]

    residue_name_one_hot_encoding: Required[list[Literal[0, 1]]]
    residue_alpha_carbon_coordinate: Required[list[float]]

    residue_chemical_classification: Required[str]
    residue_chemical_classification_multi_hot_encoding: Required[list[Literal[0, 1]]]

    residue_secondary_structure: Required[str]
    residue_secondary_structure_one_hot_encoding: Required[list[Literal[0, 1]]]

    residue_relative_solvent_accessibility: Required[float]
    residue_backbone_phi_angle: Required[float]
    residue_backbone_psi_angle: Required[float]

    residue_net_charge: Required[float]
    residue_pka: Required[float]
    residue_volume: Required[float]
    residue_hydrophobicity: Required[float]
    residue_molecular_weight: Required[float]
    residue_isoelectric_point: Required[float]

    residue_root_mean_square_fluctuation: Required[float]
    residue_sequence_conservation_score: NotRequired[float]

    residue_embedding: Required[list[float]]

    is_in_hydrophobic_cluster: Required[bool]
    is_in_disordered_region: Required[bool]
    is_in_binding_pocket: Required[bool]


class IntegratedSubstrateNode(TypedDict):
    node_index: Required[int]
    node_type: Required[Literal["substrate"]]
    node_type_one_hot_encoding: Required[list[Literal[0, 1]]]

    substrate_index: Required[int]
    substrate_name: Required[str]

    substrate_smiles: Required[str]
    substrate_atom_count: Required[int]
    substrate_molecular_weight: Required[float]
    substrate_logp: Required[float]
    substrate_tpsa: Required[float]
    substrate_heavy_atom_count: Required[int]
    substrate_hbond_donor_count: Required[int]
    substrate_hbond_acceptor_count: Required[int]
    substrate_rotatable_bond_count: Required[int]
    substrate_molar_refractivity: Required[float]
    substrate_structure_energy: Required[float]
    substrate_structure_max_3d_diameter: Required[float]
    substrate_structure_mean_pairwise_atom_distance: Required[float]
    substrate_structure_std_pairwise_atom_distance: Required[float]
    substrate_structure_asphericity: Required[float]
    substrate_structure_spherocity: Required[float]
    substrate_structure_principal_moment_ratio: Required[float]
    substrate_structure_radius_of_gyration: Required[float]
    docked_substrate_center_coordinate: Required[list[float]]
    substrate_fingerprint_encoding: Required[list[Literal[0, 1]]]


IntegratedNode: TypeAlias = IntegratedResidueNode | IntegratedSubstrateNode


class IntegratedMolecularInteraction(TypedDict):
    molecular_interaction_type: MolecularInteractionType
    molecular_interaction_one_hot_encoding: list[Literal[0, 1]]
    interaction_count: int


class IntegratedInteractionGraphEntry(TypedDict):
    molecular_interaction: IntegratedMolecularInteraction
    source_node: IntegratedNode
    target_node: IntegratedNode


class IntegratedIsolatedNodeGraphEntry(TypedDict):
    isolated_node: IntegratedNode


class IntegratedGraphEdgeNodeReference(TypedDict):
    node_index: int


class IntegratedGraphEdgeEntry(TypedDict):
    molecular_interaction: IntegratedMolecularInteraction
    source_node: IntegratedGraphEdgeNodeReference
    target_node: IntegratedGraphEdgeNodeReference


IntegratedGraphEntry: TypeAlias = (
    IntegratedInteractionGraphEntry | IntegratedIsolatedNodeGraphEntry
)

IntegratedGraphNodeCollection: TypeAlias = list[IntegratedNode]
IntegratedGraphEdgeCollection: TypeAlias = list[IntegratedGraphEdgeEntry]


class EnzyWizardIntegrateOutput(TypedDict):
    report_type: Required[Literal["enzywizard_integrate"]]
    overall_statistics: Required[IntegratedOverallStatistics]
    integrated_graph: Required[list[IntegratedGraphEntry]]


'''
mut_integrate_report
'''


class MutIntegratedOverallStatistics(TypedDict, total=False):
    wild_type_sequence_length: int
    mutant_sequence_length: int
    difference_sequence_length: int

    wild_type_total_molecular_weight: float
    mutant_total_molecular_weight: float
    difference_total_molecular_weight: float

    wild_type_total_net_charge: float
    mutant_total_net_charge: float
    difference_total_net_charge: float

    wild_type_total_residue_volume: float
    mutant_total_residue_volume: float
    difference_total_residue_volume: float

    wild_type_max_3d_diameter: float
    mutant_max_3d_diameter: float
    difference_max_3d_diameter: float

    wild_type_radius_of_gyration: float
    mutant_radius_of_gyration: float
    difference_radius_of_gyration: float

    wild_type_asphericity: float
    mutant_asphericity: float
    difference_asphericity: float

    wild_type_spherocity: float
    mutant_spherocity: float
    difference_spherocity: float

    wild_type_principal_moment_ratio: float
    mutant_principal_moment_ratio: float
    difference_principal_moment_ratio: float

    wild_type_bounding_box_volume: float
    mutant_bounding_box_volume: float
    difference_bounding_box_volume: float

    wild_type_mean_pairwise_ca_distance: float
    mutant_mean_pairwise_ca_distance: float
    difference_mean_pairwise_ca_distance: float

    wild_type_std_pairwise_ca_distance: float
    mutant_std_pairwise_ca_distance: float
    difference_std_pairwise_ca_distance: float

    wild_type_residue_name_count: list[int]
    mutant_residue_name_count: list[int]
    difference_residue_name_count: list[float]

    wild_type_residue_chemical_classification_count: list[int]
    mutant_residue_chemical_classification_count: list[int]
    difference_residue_chemical_classification_count: list[float]

    wild_type_residue_secondary_structure_count: list[int]
    mutant_residue_secondary_structure_count: list[int]
    difference_residue_secondary_structure_count: list[float]

    wild_type_hydrophobic_cluster_count: int
    mutant_hydrophobic_cluster_count: int
    difference_hydrophobic_cluster_count: int

    wild_type_max_hydrophobic_cluster_area: float
    mutant_max_hydrophobic_cluster_area: float
    difference_max_hydrophobic_cluster_area: float

    wild_type_total_hydrophobic_cluster_area: float
    mutant_total_hydrophobic_cluster_area: float
    difference_total_hydrophobic_cluster_area: float

    wild_type_disordered_region_count: int
    mutant_disordered_region_count: int
    difference_disordered_region_count: int

    wild_type_max_disordered_region_length: int
    mutant_max_disordered_region_length: int
    difference_max_disordered_region_length: int

    wild_type_total_disordered_region_length: int
    mutant_total_disordered_region_length: int
    difference_total_disordered_region_length: int

    wild_type_binding_pocket_count: int
    mutant_binding_pocket_count: int
    difference_binding_pocket_count: int

    wild_type_max_binding_pocket_volume: float
    mutant_max_binding_pocket_volume: float
    difference_max_binding_pocket_volume: float

    wild_type_total_binding_pocket_volume: float
    mutant_total_binding_pocket_volume: float
    difference_total_binding_pocket_volume: float

    wild_type_total_potential_energy: float
    mutant_total_potential_energy: float
    difference_total_potential_energy: float

    wild_type_harmonic_bond_potential_energy: float
    mutant_harmonic_bond_potential_energy: float
    difference_harmonic_bond_potential_energy: float

    wild_type_harmonic_angle_potential_energy: float
    mutant_harmonic_angle_potential_energy: float
    difference_harmonic_angle_potential_energy: float

    wild_type_custom_bond_potential_energy: float
    mutant_custom_bond_potential_energy: float
    difference_custom_bond_potential_energy: float

    wild_type_custom_torsion_potential_energy: float
    mutant_custom_torsion_potential_energy: float
    difference_custom_torsion_potential_energy: float

    wild_type_custom_nonbonded_potential_energy: float
    mutant_custom_nonbonded_potential_energy: float
    difference_custom_nonbonded_potential_energy: float

    wild_type_nonbonded_potential_energy: float
    mutant_nonbonded_potential_energy: float
    difference_nonbonded_potential_energy: float

    wild_type_periodic_torsion_potential_energy: float
    mutant_periodic_torsion_potential_energy: float
    difference_periodic_torsion_potential_energy: float

    wild_type_cmap_torsion_potential_energy: float
    mutant_cmap_torsion_potential_energy: float
    difference_cmap_torsion_potential_energy: float

    wild_type_enzyme_substrate_binding_affinity: float
    mutant_enzyme_substrate_binding_affinity: float
    difference_enzyme_substrate_binding_affinity: float

    wild_type_hydrogen_bond_count: int
    mutant_hydrogen_bond_count: int
    difference_hydrogen_bond_count: int

    wild_type_ionic_bond_count: int
    mutant_ionic_bond_count: int
    difference_ionic_bond_count: int

    wild_type_van_der_waals_contact_count: int
    mutant_van_der_waals_contact_count: int
    difference_van_der_waals_contact_count: int

    wild_type_pi_pi_stacking_count: int
    mutant_pi_pi_stacking_count: int
    difference_pi_pi_stacking_count: int

    wild_type_pi_cation_interaction_count: int
    mutant_pi_cation_interaction_count: int
    difference_pi_cation_interaction_count: int

    wild_type_disulfide_bond_count: int
    mutant_disulfide_bond_count: int
    difference_disulfide_bond_count: int


class AminoAcidSubstitutionProperties(TypedDict, total=False):

    wild_type_residue_name_one_hot_encoding: list[float]
    mutant_residue_name_one_hot_encoding: list[float]
    difference_residue_name_one_hot_encoding: list[float]

    wild_type_residue_chemical_classification: str
    mutant_residue_chemical_classification: str

    wild_type_residue_chemical_classification_one_hot_encoding: list[float]
    mutant_residue_chemical_classification_one_hot_encoding: list[float]
    difference_residue_chemical_classification_one_hot_encoding: list[float]


    wild_type_residue_secondary_structure_one_hot_encoding: list[float]
    mutant_residue_secondary_structure_one_hot_encoding: list[float]
    difference_residue_secondary_structure_one_hot_encoding: list[float]

    wild_type_residue_relative_solvent_accessibility: float
    mutant_residue_relative_solvent_accessibility: float
    difference_residue_relative_solvent_accessibility: float

    wild_type_residue_backbone_phi_angle: float
    mutant_residue_backbone_phi_angle: float
    difference_residue_backbone_phi_angle: float

    wild_type_residue_backbone_psi_angle: float
    mutant_residue_backbone_psi_angle: float
    difference_residue_backbone_psi_angle: float

    wild_type_residue_net_charge: float
    mutant_residue_net_charge: float
    difference_residue_net_charge: float

    wild_type_residue_pka: float
    mutant_residue_pka: float
    difference_residue_pka: float

    wild_type_residue_volume: float
    mutant_residue_volume: float
    difference_residue_volume: float

    wild_type_residue_hydrophobicity: float
    mutant_residue_hydrophobicity: float
    difference_residue_hydrophobicity: float

    wild_type_residue_molecular_weight: float
    mutant_residue_molecular_weight: float
    difference_residue_molecular_weight: float

    wild_type_residue_isoelectric_point: float
    mutant_residue_isoelectric_point: float
    difference_residue_isoelectric_point: float

    wild_type_residue_root_mean_square_fluctuation: float
    mutant_residue_root_mean_square_fluctuation: float
    difference_residue_root_mean_square_fluctuation: float

    wild_type_residue_sequence_conservation_score: float
    mutant_residue_sequence_conservation_score: float
    difference_residue_sequence_conservation_score: float

    wild_type_mutation_site_distance_to_centroid: float
    mutant_mutation_site_distance_to_centroid: float
    difference_mutation_site_distance_to_centroid: float

    wild_type_mutation_site_distance_to_nearest_binding_pocket: float
    mutant_mutation_site_distance_to_nearest_binding_pocket: float
    difference_mutation_site_distance_to_nearest_binding_pocket: float

    wild_type_mutation_site_distance_to_nearest_hydrophobic_cluster: float
    mutant_mutation_site_distance_to_nearest_hydrophobic_cluster: float
    difference_mutation_site_distance_to_nearest_hydrophobic_cluster: float

    wild_type_mutation_site_distance_to_nearest_disordered_region: float
    mutant_mutation_site_distance_to_nearest_disordered_region: float
    difference_mutation_site_distance_to_nearest_disordered_region: float

    wild_type_mutation_site_distance_to_nearest_substrate: float
    mutant_mutation_site_distance_to_nearest_substrate: float
    difference_mutation_site_distance_to_nearest_substrate: float



class EnzyWizardMutIntegrateOutput(TypedDict):
    report_type: Required[Literal["enzywizard_mut_integrate"]]
    cleaned_amino_acid_substitution: Required[str]

    overall_statistics: Required[MutIntegratedOverallStatistics]
    amino_acid_substitution_properties: Required[AminoAcidSubstitutionProperties]

    wild_type_integrated_graph: Required[list[IntegratedGraphEntry]]
    mutant_integrated_graph: Required[list[IntegratedGraphEntry]]
