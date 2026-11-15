####################################
################ BioPython functions
####################################
### Import dependencies
import os
import math
import numpy as np
from pathlib import Path
from typing import Optional, List, Dict
from collections import defaultdict
from scipy.spatial import cKDTree
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import is_aa

# analyze sequence composition of design
def validate_design_sequence(sequence, num_clashes, advanced_settings):
    note_array = []

    # Check if protein contains clashes after relaxation
    if num_clashes > 0:
        note_array.append('Relaxed structure contains clashes.')

    # Check if the sequence contains disallowed amino acids
    if advanced_settings["omit_AAs"]:
        restricted_AAs = advanced_settings["omit_AAs"].split(',')
        for restricted_AA in restricted_AAs:
            if restricted_AA in sequence:
                note_array.append('Contains: '+restricted_AA+'!')

    # Analyze the protein
    analysis = ProteinAnalysis(sequence)

    # Calculate the reduced extinction coefficient per 1% solution
    extinction_coefficient_reduced = analysis.molar_extinction_coefficient()[0]
    molecular_weight = round(analysis.molecular_weight() / 1000, 2)
    extinction_coefficient_reduced_1 = round(extinction_coefficient_reduced / molecular_weight * 0.01, 2)

    # Check if the absorption is high enough
    if extinction_coefficient_reduced_1 <= 2:
        note_array.append(f'Absorption value is {extinction_coefficient_reduced_1}, consider adding tryptophane to design.')

    # Join the notes into a single string
    notes = ' '.join(note_array)

    return notes

# temporary function, calculate RMSD of input PDB and trajectory target
def target_pdb_rmsd(trajectory_pdb, starting_pdb, chain_ids_string):
    # Parse the PDB files
    parser = PDBParser(QUIET=True)
    structure_trajectory = parser.get_structure('trajectory', trajectory_pdb)
    structure_starting = parser.get_structure('starting', starting_pdb)
    
    # Extract chain A from trajectory_pdb
    chain_trajectory = structure_trajectory[0]['A']
    
    # Extract the specified chains from starting_pdb
    chain_ids = chain_ids_string.split(',')
    residues_starting = []
    for chain_id in chain_ids:
        chain_id = chain_id.strip()
        chain = structure_starting[0][chain_id]
        for residue in chain:
            if is_aa(residue, standard=True):
                residues_starting.append(residue)
    
    # Extract residues from chain A in trajectory_pdb
    residues_trajectory = [residue for residue in chain_trajectory if is_aa(residue, standard=True)]
    
    # Ensure that both structures have the same number of residues
    min_length = min(len(residues_starting), len(residues_trajectory))
    residues_starting = residues_starting[:min_length]
    residues_trajectory = residues_trajectory[:min_length]
    
    # Collect CA atoms from the two sets of residues
    atoms_starting = [residue['CA'] for residue in residues_starting if 'CA' in residue]
    atoms_trajectory = [residue['CA'] for residue in residues_trajectory if 'CA' in residue]
    
    # Calculate RMSD using structural alignment
    sup = Superimposer()
    sup.set_atoms(atoms_starting, atoms_trajectory)
    rmsd = sup.rms
    
    return round(rmsd, 2)

# detect C alpha clashes for deformed trajectories
def calculate_clash_score(pdb_file, threshold=2.4, only_ca=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    atoms = []
    atom_info = []  # Detailed atom info for debugging and processing

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'H':  # Skip hydrogen atoms
                        continue
                    if only_ca and atom.get_name() != 'CA':
                        continue
                    atoms.append(atom.coord)
                    atom_info.append((chain.id, residue.id[1], atom.get_name(), atom.coord))

    tree = cKDTree(atoms)
    pairs = tree.query_pairs(threshold)

    valid_pairs = set()
    for (i, j) in pairs:
        chain_i, res_i, name_i, coord_i = atom_info[i]
        chain_j, res_j, name_j, coord_j = atom_info[j]

        # Exclude clashes within the same residue
        if chain_i == chain_j and res_i == res_j:
            continue

        # Exclude directly sequential residues in the same chain for all atoms
        if chain_i == chain_j and abs(res_i - res_j) == 1:
            continue

        # If calculating sidechain clashes, only consider clashes between different chains
        if not only_ca and chain_i == chain_j:
            continue

        valid_pairs.add((i, j))

    return len(valid_pairs)

three_to_one_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# identify interacting residues at the binder interface
def hotspot_residues(trajectory_pdb, binder_chain="B", atom_distance_cutoff=4.0):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", trajectory_pdb)

    # Get the specified chain
    binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A')
    binder_coords = np.array([atom.coord for atom in binder_atoms])

    # Get atoms and coords for the target chain
    target_atoms = Selection.unfold_entities(structure[0]['A'], 'A')
    target_coords = np.array([atom.coord for atom in target_atoms])

    # Build KD trees for both chains
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    # Prepare to collect interacting residues
    interacting_residues = {}

    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    # Process each binder atom's interactions
    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        # Convert three-letter code to single-letter code using the manual dictionary
        if binder_resname in three_to_one_map:
            aa_single_letter = three_to_one_map[binder_resname]
            for close_idx in close_indices:
                target_residue = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter

    return interacting_residues


# Return 0-indexed residue positions on target (chain A) and binder (binder_chain) that
# have any atom within `cutoff` Å of an atom from the other chain (project convention: 4 Å all-atom).
def interface_residue_indices(pdb_file, binder_chain="B", cutoff=4.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    model = structure[0]

    target_residues = list(model["A"].get_residues())
    binder_residues = list(model[binder_chain].get_residues())

    target_atoms = Selection.unfold_entities(model["A"], "A")
    binder_atoms = Selection.unfold_entities(model[binder_chain], "A")
    target_coords = np.array([a.coord for a in target_atoms])
    binder_coords = np.array([a.coord for a in binder_atoms])
    target_resid_per_atom = [a.get_parent().id[1] for a in target_atoms]
    binder_resid_per_atom = [a.get_parent().id[1] for a in binder_atoms]

    target_tree = cKDTree(target_coords)
    binder_tree = cKDTree(binder_coords)
    pairs = binder_tree.query_ball_tree(target_tree, cutoff)

    t_iface, b_iface = set(), set()
    for b_idx, close in enumerate(pairs):
        if close:
            b_iface.add(binder_resid_per_atom[b_idx])
            for t_idx in close:
                t_iface.add(target_resid_per_atom[t_idx])

    target_id_to_pos = {res.id[1]: i for i, res in enumerate(target_residues)}
    binder_id_to_pos = {res.id[1]: i for i, res in enumerate(binder_residues)}
    target_idx = sorted(target_id_to_pos[r] for r in t_iface if r in target_id_to_pos)
    binder_idx = sorted(binder_id_to_pos[r] for r in b_iface if r in binder_id_to_pos)
    return target_idx, binder_idx, len(target_residues), len(binder_residues)


# calculate secondary structure percentage of design
def calc_ss_percentage(pdb_file, advanced_settings, chain_id="B", atom_distance_cutoff=4.0):
    # Parse the structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]  # Consider only the first model in the structure

    # Calculate DSSP for the model
    dssp = DSSP(model, pdb_file, dssp=advanced_settings["dssp_path"])

    # Prepare to count residues
    ss_counts = defaultdict(int)
    ss_interface_counts = defaultdict(int)
    plddts_interface = []
    plddts_ss = []

    # Get chain and interacting residues once
    chain = model[chain_id]
    interacting_residues = set(hotspot_residues(pdb_file, chain_id, atom_distance_cutoff).keys())

    for residue in chain:
        residue_id = residue.id[1]
        if (chain_id, residue_id) in dssp:
            ss = dssp[(chain_id, residue_id)][2]  # Get the secondary structure
            ss_type = 'loop'
            if ss in ['H', 'G', 'I']:
                ss_type = 'helix'
            elif ss == 'E':
                ss_type = 'sheet'

            ss_counts[ss_type] += 1

            if ss_type != 'loop':
                # calculate secondary structure normalised pLDDT
                avg_plddt_ss = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_ss.append(avg_plddt_ss)

            if residue_id in interacting_residues:
                ss_interface_counts[ss_type] += 1

                # calculate interface pLDDT
                avg_plddt_residue = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_interface.append(avg_plddt_residue)

    # Calculate percentages
    total_residues = sum(ss_counts.values())
    total_interface_residues = sum(ss_interface_counts.values())

    percentages = calculate_percentages(total_residues, ss_counts['helix'], ss_counts['sheet'])
    interface_percentages = calculate_percentages(total_interface_residues, ss_interface_counts['helix'], ss_interface_counts['sheet'])

    i_plddt = round(sum(plddts_interface) / len(plddts_interface) / 100, 2) if plddts_interface else 0
    ss_plddt = round(sum(plddts_ss) / len(plddts_ss) / 100, 2) if plddts_ss else 0

    return (*percentages, *interface_percentages, i_plddt, ss_plddt)

def calculate_percentages(total, helix, sheet):
    helix_percentage = round((helix / total) * 100,2) if total > 0 else 0
    sheet_percentage = round((sheet / total) * 100,2) if total > 0 else 0
    loop_percentage = round(((total - helix - sheet) / total) * 100,2) if total > 0 else 0

    return helix_percentage, sheet_percentage, loop_percentage


####################################
###### Full-target validation (ported from Bindome/full-target fork)
###### Superimpose crop onto full-length target, flag confident binder clashes
####################################

def _residues(chain):
    return [r for r in Selection.unfold_entities(chain, "R") if is_aa(r, standard=True)]

def _ca_for_resnums(chain, resnums_set):
    ca = []
    for r in _residues(chain):
        if r.id[1] in resnums_set and "CA" in r:
            ca.append(r["CA"])
    return ca

def _chain_ca_atoms(structure, chain_id):
    out = []
    for a in Selection.unfold_entities(structure, "A"):
        if a.get_parent().get_parent().id == chain_id and a.id == "CA":
            out.append(a)
    return out

def save_superposition(fixed_struct, moved_struct, out_path, save_full_structure=False):
    builder = StructureBuilder()
    builder.init_structure("superposition")

    builder.init_model(0)
    for ch in moved_struct[0]:
        builder.structure[0].add(ch.copy())

    if save_full_structure:
        builder.init_model(1)
        for ch in fixed_struct[0]:
            builder.structure[1].add(ch.copy())

    io = PDBIO()
    io.set_structure(builder.get_structure())
    io.save(str(out_path))

def calculate_clashes_after_superpostion_on_full_target(
    model_pdb_path: str,
    target_full_structure: Structure,
    target_full_pae: np.ndarray,
    target_crop_indices: List[int],     # PDB residue numbers (1-based)
    superimposed_filepath: Optional[str] = None,
    clash_radius: float = 3.0,             # Clash distance cutoff (Å)
    target_chain_id: str = "A",                         # Target chain ID
    binder_chain_id: str = "B",
    clash_plddt_threshold: float = 70.0,
    clash_pae_threshold: float = 15.0,
) -> Dict[str, float]:
    """
    Superimpose a cropped target structure onto the full target,
    then detect CA-based clashes between the binder and the
    full target outside the crop region. A clash only counts (fails)
    if the clashing full-target residue is confidently placed
    (high pLDDT AND low PAE relative to the crop).
    """

    crop_resnums = set(target_crop_indices)

    parser = PDBParser(QUIET=True)
    target_crop_binder = parser.get_structure("protein", model_pdb_path)

    target_full_chain = target_full_structure[0][target_chain_id]
    target_crop_chain = target_crop_binder[0][target_chain_id]

    # 1) Superimpose crop target onto full target (CA atoms on crop resnums)
    fixed_ca = _ca_for_resnums(target_full_chain, crop_resnums)
    moving_ca = [r["CA"] for r in _residues(target_crop_chain) if "CA" in r]

    if len(fixed_ca) != len(moving_ca):
        raise ValueError(f"Mismatch in number of CA atoms for alignment: {model_pdb_path}")

    sup = Superimposer()
    sup.set_atoms(fixed=fixed_ca, moving=moving_ca)
    sup.apply(target_crop_binder.get_atoms())

    # Save combined PDB (optional but useful)
    if superimposed_filepath:
        Path(superimposed_filepath).parent.mkdir(parents=True, exist_ok=True)
        save_superposition(target_full_structure, target_crop_binder, superimposed_filepath)

    # 2) Build CA lists for full target excluding crop, and binder CA
    full_res = _residues(target_full_chain)

    crop_zero_idx = np.array([i for i, r in enumerate(full_res) if r.id[1] in crop_resnums], dtype=int)

    full_ex_crop_zero_idx = []
    full_ex_crop_ca_atoms = []
    for i, r in enumerate(full_res):
        if r.id[1] in crop_resnums:
            continue
        if "CA" in r:
            full_ex_crop_zero_idx.append(i)     # maps KD-tree index -> residue index in full_res
            full_ex_crop_ca_atoms.append(r["CA"])

    binder_ca_atoms = _chain_ca_atoms(target_crop_binder, binder_chain_id)

    coords_target = np.asarray([a.get_coord() for a in full_ex_crop_ca_atoms], float)
    coords_binder = np.asarray([a.get_coord() for a in binder_ca_atoms], float)

    if coords_target.size == 0:
        # The crop covers the entire target
        return {"clash_pass": True}

    # 3) KD-tree query: any target CA within threshold of any binder CA
    tree = cKDTree(coords_target)
    neighbor_lists = tree.query_ball_point(coords_binder, r=clash_radius, workers=4)

    if not any(neighbor_lists):
        return {"clash_pass": True}

    hit_idx = np.unique(np.concatenate([np.asarray(x, dtype=int) for x in neighbor_lists if x]))
    clash_res_indices = np.unique(np.asarray(full_ex_crop_zero_idx, dtype=int)[hit_idx])  # indices into full_res, zero-indexed

    # 4) Clash confidence: a clash only fails if the clashing residue is confidently
    #    placed — high pLDDT (max over clashers) AND low PAE (min over clashers) vs the crop.
    plddts = [full_res[i]["CA"].get_bfactor() for i in clash_res_indices if "CA" in full_res[i]]
    max_plddt = float(np.max(plddts)) if plddts else float("nan")

    sub1 = target_full_pae[np.ix_(clash_res_indices, crop_zero_idx)]
    sub2 = target_full_pae[np.ix_(crop_zero_idx, clash_res_indices)].T
    min_pae = float(((sub1 + sub2) * 0.5).min())

    # Pass if clashing residues are "not confident" i.e. low pLDDT OR high PAE
    clash_pass = bool(max_plddt <= clash_plddt_threshold or min_pae >= clash_pae_threshold)

    return {"clash_pass": clash_pass}

def superimpose_on_full_target_and_save_pdb(best_model_pdb, target_full_structure, target_crop_indices, target_settings, mpnn_design_name, best_model_number, design_paths):
    crop_resnums = set(target_crop_indices)

    parser = PDBParser(QUIET=True)
    target_crop_binder = parser.get_structure("protein", best_model_pdb)

    target_full_chain = target_full_structure[0][target_settings["chains"]]
    target_crop_chain = target_crop_binder[0][target_settings["chains"]]

    fixed_ca = _ca_for_resnums(target_full_chain, crop_resnums)
    moving_ca = [r["CA"] for r in _residues(target_crop_chain) if "CA" in r]

    if len(fixed_ca) != len(moving_ca):
        raise ValueError(f"Mismatch in number of CA atoms for alignment: {best_model_pdb}")

    sup = Superimposer()
    sup.set_atoms(fixed=fixed_ca, moving=moving_ca)
    sup.apply(target_crop_binder.get_atoms())

    superimposed_filepath = os.path.join(design_paths["Accepted/Superimposed"], f"{mpnn_design_name}_model{best_model_number}.pdb")

    Path(superimposed_filepath).parent.mkdir(parents=True, exist_ok=True)
    save_superposition(target_full_structure, target_crop_binder, superimposed_filepath)


def binder_crop_terminus_engagement(
    model_pdb_path,
    target_crop_indices=None,           # original PDB resnums of the crop (from starting_pdb)
    full_target_structure=None,         # full-length target; enables cut-vs-real-end detection
    terminus_len: int = 5,
    target_chain_id: str = "A",
    binder_chain_id: str = "B",
    cutoff: float = 4.0,
    only_cut_ends: bool = True,
):
    """How much does the binder's interface engage the CROP's terminal residues?

    A gap crop is a fragment carved from a full protein, with hydrophilic caps
    snapped onto both termini. Under full-length superimposition, a *cut* terminus
    is exactly where the real chain continues, so a binder gripping it is an
    artifact of cropping (it grips the arbitrary padding, and in the real protein
    that space is occupied by the flanking chain). This flags such designs.

    - `terminus_len` K: the first K and last K crop residues (by chain order) are
      "termini". The check is positional, so it is independent of how chain A is
      numbered (original vs renumbered 1..L).
    - `only_cut_ends`: when `full_target_structure` + `target_crop_indices` are
      given, a terminus counts ONLY if the full chain extends past it (the crop's
      min/max resnum is not the protein's first/last residue). A genuine protein
      terminus is exempt — binding there is not a cropping artifact. Without the
      full structure both ends are treated as termini (conservative).

    Returns a dict of counts + booleans; does not itself reject anything.
    """
    result = {
        "n_terminus_iface": 0,
        "terminus_iface_fraction": 0.0,
    }
    try:
        t_idx, _b_idx, target_len, _blen = interface_residue_indices(
            model_pdb_path, binder_chain=binder_chain_id, cutoff=cutoff)
    except Exception:
        return result
    if target_len <= 0:
        return result

    # decide which ends are "penalised" termini
    n_end_active, c_end_active = True, True
    if only_cut_ends and full_target_structure is not None and target_crop_indices:
        try:
            crop_sorted = sorted(int(x) for x in target_crop_indices)
            full_nums = [r.id[1] for r in _residues(full_target_structure[0][target_chain_id])]
            if full_nums:
                n_end_active = crop_sorted[0] > min(full_nums)    # chain continues before the crop
                c_end_active = crop_sorted[-1] < max(full_nums)   # chain continues after the crop
        except Exception as _e:
            # fall back to treating both ends as termini -- but SAY so, since silent fallback
            # was the pipeline bug on 2026-10-05: real-terminus contacts were being counted
            # as cropping artifacts, which with a tight gate rejects legitimate designs.
            print(f"[binder_crop_terminus_engagement] WARNING: only_cut_ends setup failed "
                  f"({type(_e).__name__}: {_e}); both ends treated as cut termini. "
                  f"Pass full_target_structure + target_crop_indices with original crop "
                  f"numbering to enable the real-terminus exemption.")
    elif only_cut_ends:
        # requested but missing inputs -- same fallback, same warning (prints once per call)
        print(f"[binder_crop_terminus_engagement] NOTE: only_cut_ends=True but "
              f"full_target_structure/target_crop_indices not provided; both ends treated as "
              f"cut termini. This counts the REAL protein terminus as a cropping artifact.")

    K = min(terminus_len, target_len)
    n_positions = set(range(0, K)) if n_end_active else set()
    c_positions = set(range(target_len - K, target_len)) if c_end_active else set()
    term_positions = n_positions | c_positions

    n_iface = len(t_idx)
    term_iface = [i for i in t_idx if i in term_positions]
    result["n_terminus_iface"] = len(term_iface)
    result["terminus_iface_fraction"] = round(len(term_iface) / n_iface, 3) if n_iface else 0.0
    return result

def binder_terminus_clearance(
    model_pdb_path,
    target_chain_id: str = "A",
    binder_chain_id: str = "B",
    iface_cutoff: float = 5.0,
):
    """How much room does each binder terminus have from its own interface?

    A purification tag or fusion partner on a binder terminus needs that terminus
    to sit away from the binding surface. Long binders against short, extended IDP
    crops tend to *wrap* the target, which leaves both termini buried in the
    interface and the design untaggable even when every confidence metric is good.

    For the binder's first and last residue this measures the minimum heavy-atom
    distance to the nearest OTHER binder residue that contacts the target
    (a binder residue within `iface_cutoff` of any target heavy atom).

    Returns the two clearances and their max -- the usable tag end. A terminus
    that is itself an interface residue still reports its distance to the rest of
    the interface, so `*_is_iface` is given separately for an unambiguous reject.
    """
    result = {
        "binder_nterm_clearance": 0.0,
        "binder_cterm_clearance": 0.0,
        "binder_term_clearance": 0.0,
        "binder_nterm_is_iface": False,
        "binder_cterm_is_iface": False,
    }
    try:
        model = PDBParser(QUIET=True).get_structure("complex", model_pdb_path)[0]
        if target_chain_id not in model or binder_chain_id not in model:
            return result
        bres = [r for r in model[binder_chain_id] if r.id[0] == " "]
        if len(bres) < 3:
            return result

        tcoord = np.array([a.coord for a in model[target_chain_id].get_atoms()
                           if a.element != "H"])
        bxyz, brid = [], []
        for r in bres:
            for a in r:
                if a.element != "H":
                    bxyz.append(a.coord)
                    brid.append(r.id[1])
        if not len(tcoord) or not bxyz:
            return result
        bxyz = np.array(bxyz)
        brid = np.array(brid)

        d_to_target, _ = cKDTree(tcoord).query(bxyz, workers=1)
        iface = set(brid[d_to_target < iface_cutoff].tolist())
        if not iface:
            return result

        def clearance(rid):
            others = [i for i in iface if i != rid]
            if not others:
                return 0.0
            sel = bxyz[brid == rid]
            pool = bxyz[np.isin(brid, list(others))]
            return float(cKDTree(pool).query(sel, workers=1)[0].min())

        n_id, c_id = bres[0].id[1], bres[-1].id[1]
        n_cl, c_cl = clearance(n_id), clearance(c_id)
        result["binder_nterm_clearance"] = round(n_cl, 2)
        result["binder_cterm_clearance"] = round(c_cl, 2)
        result["binder_term_clearance"] = round(max(n_cl, c_cl), 2)
        result["binder_nterm_is_iface"] = n_id in iface
        result["binder_cterm_is_iface"] = c_id in iface
    except Exception:
        return result
    return result
