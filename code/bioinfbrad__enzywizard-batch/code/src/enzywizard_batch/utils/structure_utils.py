from __future__ import annotations

from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

from Bio.PDB import Atom
from Bio.PDB.Structure import Structure
from ..utils.logging_utils import Logger
from Bio.Data.IUPACData import protein_letters_3to1
from openmm.app import Modeller
from openmm.app.element import hydrogen as ELEMENT_H

from typing import List, Tuple, Dict, Any
import numpy as np

def get_first_model(struct: Structure, logger: Logger) -> Model | None:
    for m in struct:
        return m
    logger.print(f"[ERROR] No model found in structure")
    return None

def get_single_chain(struct: Structure, logger:Logger) -> Chain | None:
    m = get_first_model(struct, logger)
    if not m:
        return None
    for chain in m.get_chains():
        return chain
    logger.print(f"[ERROR] No chain found in structure")
    return None

def get_chain_length(chain: Chain, logger:Logger) -> int | None:
    if chain is None:
        logger.print(f"[ERROR] Bad chain input")
        return None

    length = 0

    for residue in chain:
        if residue.id[0] != " ":
            continue
        length += 1
    if not length:
        logger.print(f"[ERROR] Invalid sequence length")
        return None
    return length



def get_residues_by_chain(chain: Chain, logger: Logger) -> List[Tuple[Tuple[str,int,str], str, Tuple[float, float, float]]] | None:
    result: List[Tuple[Tuple[str, int, str], str, Tuple[float, float, float]]] = []

    for res in chain.get_residues():
        hetflag, resseq, icode = res.id

        if hetflag != " ":
            continue

        resname = res.get_resname().strip()

        # 必须有 CA 原子
        if "CA" not in res:
            logger.print(f"[ERROR] Residue {resseq} {resname} missing CA atom, skipped.")
            return None

        ca: Atom = res["CA"]
        coord = tuple(ca.get_coord())  # (x, y, z)

        result.append(((hetflag,resseq,icode), resname, coord))

    return result

def get_residue_ca_coord_by_aa_id(struct: Structure, aa_id: int, logger: Logger) -> List[float] | None:
    if struct is None:
        logger.print("[ERROR] Structure input is None.")
        return None

    if not isinstance(aa_id, int) or aa_id <= 0:
        logger.print("[ERROR] catalytic_residue must be a positive integer aa_id.")
        return None

    chain = get_single_chain(struct, logger)
    if chain is None:
        return None

    residue_list = get_residues_by_chain(chain, logger)
    if residue_list is None or len(residue_list) == 0:
        logger.print("[ERROR] Failed to get cleaned protein residues.")
        return None

    valid_aa_id_list = [int(residue_key[1]) for residue_key, _, _ in residue_list]
    for residue_key, _, coord in residue_list:
        if int(residue_key[1]) == aa_id:
            return [float(coord[0]), float(coord[1]), float(coord[2])]

    logger.print(
        f"[ERROR] catalytic_residue aa_id {aa_id} is out of range. "
        f"Valid aa_id range: {min(valid_aa_id_list)}-{max(valid_aa_id_list)}."
    )
    return None

def get_sequence(residues: List[Tuple[Tuple[str, int, str], str, Tuple[float, float, float]]],logger: Logger) -> str | None:

    if residues is None:
        logger.print("[ERROR] Residues input is None")
        return None

    if len(residues) == 0:
        logger.print("[ERROR] Empty residues list")
        return None

    seq_chars: List[str] = []

    for residue_key, resname, _ in residues:
        aa = protein_letters_3to1.get(resname.capitalize())
        if aa is None:
            logger.print(f"[ERROR] Unsupported residue name: {resname}")
            return None
        seq_chars.append(aa)

    return "".join(seq_chars)

def get_structure_box(struct: Structure, logger: Logger) -> Dict[str, Any] | None:
    try:
        coord_list = []

        for atom in struct.get_atoms():
            coord = atom.get_coord()
            if coord is None or len(coord) != 3:
                continue

            try:
                x = float(coord[0])
                y = float(coord[1])
                z = float(coord[2])
            except Exception:
                continue

            coord_list.append([x, y, z])

        if len(coord_list) == 0:
            logger.print("[ERROR] No valid atom coordinates found in Structure.")
            return None

        coords = np.asarray(coord_list, dtype=float)

        mn = np.min(coords, axis=0)
        mx = np.max(coords, axis=0)

        center_coord = ((mn + mx) / 2.0).tolist()
        box_boundaries = (mx - mn).tolist()

        return {
            "center_coord": [float(x) for x in center_coord],
            "box_boundaries": [float(x) for x in box_boundaries],
        }

    except Exception as e:
        logger.print(f"[ERROR] Failed to compute box space for Structure. Reason: {e}")
        return None


def structure_has_hydrogen(struct: Structure, logger) -> bool:
    if not isinstance(struct, Structure):
        logger.print("[ERROR] struct must be a Bio.PDB Structure.")
        return False

    try:
        for atom in struct.get_atoms():
            element = getattr(atom, "element", None)

            if element == "H" or (element is None and atom.get_name().startswith("H")):
                return True

        return False

    except Exception as e:
        logger.print(f"[ERROR] Failed to inspect hydrogens in structure. Reason: {e}")
        return False

def structure_has_too_few_hydrogens(struct: Structure,logger,min_hydrogen_count: int = 10,min_hydrogen_ratio: float = 0.05) -> bool:
    if not isinstance(struct, Structure):
        logger.print("[ERROR] struct must be a Bio.PDB Structure.")
        return True

    if not isinstance(min_hydrogen_count, int) or min_hydrogen_count < 0:
        logger.print("[ERROR] min_hydrogen_count must be a non-negative integer.")
        return True

    if not isinstance(min_hydrogen_ratio, (int, float)) or float(min_hydrogen_ratio) < 0.0:
        logger.print("[ERROR] min_hydrogen_ratio must be a non-negative number.")
        return True

    try:
        atom_count = 0
        h_count = 0

        for atom in struct.get_atoms():
            atom_count += 1

            element = getattr(atom, "element", None)
            if element == "H" or (element is None and atom.get_name().startswith("H")):
                h_count += 1

        if atom_count <= 0:
            logger.print("[ERROR] Structure contains no atoms.")
            return True

        h_ratio = float(h_count) / float(atom_count)

        if h_count < min_hydrogen_count:
            return True

        if h_ratio < float(min_hydrogen_ratio):
            return True

        return False

    except Exception as e:
        logger.print(f"[ERROR] Failed to count hydrogens in structure. Reason: {e}")
        return True


def get_fasta_dict_from_structure(
    struct: Structure,
    logger: Logger,
    header: str = "query",
) -> Dict[str, str] | None:
    if not isinstance(struct, Structure):
        logger.print("[ERROR] struct must be a Bio.PDB Structure.")
        return None

    if not isinstance(header, str) or not header.strip():
        logger.print("[ERROR] header must be a non-empty string.")
        return None

    chain = get_single_chain(struct, logger)
    if chain is None:
        logger.print("[ERROR] Failed to get single chain from structure.")
        return None

    residues = get_residues_by_chain(chain, logger)
    if residues is None:
        logger.print("[ERROR] Failed to get residues from structure chain.")
        return None

    sequence = get_sequence(residues, logger)
    if sequence is None:
        logger.print("[ERROR] Failed to get sequence from structure residues.")
        return None

    return {
        "header": header.strip(),
        "sequence": sequence,
    }
