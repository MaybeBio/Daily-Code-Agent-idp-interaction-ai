"""
Core constants for IDPInterp.
"""

# ESM-2 max sequence length
ESM2_MAX_SEQ_LEN = 1022

# Standard MobiDB disorder triplets we care about (original generic-disorder set)
TARGET_TRIPLETS = [
    "curated-disorder-merge",
    "homology-disorder-merge",
    "prediction-disorder-th_50",
    "curated-disorder-priority",
    "curated-lip-merge",
    "derived-lip-merge",
]

# Sub-state axes on top of generic disorder: MoRF, fuzzy, molten globule, random coil,
# and compositional classes. All are residue-level region triplets served by the MobiDB
# JSON download API, so they flow through process.py unchanged. Each protein is parsed
# for every triplet key in its MobiDB record, so extending this list only requires
# re-running 01_fetch_mobidb and process.py.
#
# Continuous annotations (e.g. prediction-plddt-alphafold) and DisProt IDPO structural
# states are not region triplets and are handled separately.

# Binding-mode sub-states
BINDING_MODE_TRIPLETS = [
    "derived-binding_mode_disorder_to_order-mobi",       # MoRFs (disorder→order on binding)
    "derived-binding_mode_disorder_to_disorder-mobi",     # fuzzy complexes (stay disordered)
]

# Compositional sub-types (validate the 9-class biochemical taxonomy; also
# molten-globule (compact) vs. random-coil (extended) proxies)
COMPOSITIONAL_TRIPLETS = [
    "prediction-compact-mobidb_lite_sub",                 # molten-globule proxy
    "prediction-extended-mobidb_lite_sub",                # random-coil proxy
    "prediction-polyampholyte-mobidb_lite_sub",
    "prediction-low_complexity-mobidb_lite_sub",
    "prediction-proline_rich-mobidb_lite_sub",
    "prediction-glycine_rich-mobidb_lite_sub",
]

# Structural evidence + negative control
STRUCTURAL_EVIDENCE_TRIPLETS = [
    "derived-missing_residues-th_90",                     # X-ray missing residues (disorder proxy)
    "derived-fixed-th_90",                                # ordered crystal residues (NEGATIVE CONTROL)
]

# AlphaFold-derived disorder
ALPHAFOLD_TRIPLETS = [
    "prediction-disorder-alphafold",
]

# Amyloid / aggregation axis
AMYLOID_TRIPLETS = [
    "curated-lip-mfib",                                   # MFIB amyloid fibrils
]

# DisProt experimental structural states. Not MobiDB keys: these columns are written by
# scripts/08_disprot_substates.py (+disprot_mode=build), which adds mask_disprot_* to
# gold_curated.parquet and appends region rows. Without that build step the evaluator
# skips them.
DISPROT_TRIPLETS = [
    "disprot-molten_globule",        # experimental molten globule (vs random coil)
    "disprot-pre_molten_globule",
    "disprot-disorder_to_order",     # experimental MoRF
    "disprot-order_to_disorder",
]

# Additional triplets for phase separation, binding, etc. (original extras)
EXTRA_TRIPLETS = [
    "curated-phase_separation-merge",
    "curated-binding_mode_disorder_to_disorder-fuzdb",
]

# Full evaluation set used by process.py and 04_evaluate_atlas.py
ALL_TARGET_TRIPLETS = (
    TARGET_TRIPLETS
    + EXTRA_TRIPLETS
    + BINDING_MODE_TRIPLETS
    + COMPOSITIONAL_TRIPLETS
    + STRUCTURAL_EVIDENCE_TRIPLETS
    + ALPHAFOLD_TRIPLETS
    + AMYLOID_TRIPLETS
    + DISPROT_TRIPLETS
)
