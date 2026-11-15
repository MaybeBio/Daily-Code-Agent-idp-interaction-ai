from __future__ import annotations

from Bio.PDB import MMCIFParser, PDBParser, MMCIFIO, PDBIO
from Bio.PDB.Structure import Structure
from pathlib import Path
from Bio.PDB.DSSP import DSSP
from ..utils.logging_utils import Logger
import gzip
import json
import tempfile
from ..utils.common_utils import convert_to_json_serializable, InlineJSONEncoder, wrap_leaf_lists_as_rawjson, get_clean_filename, get_optimized_filename
from openmm.app import PDBFile,PDBxFile, Modeller
from ..utils.structure_utils import get_single_chain,get_residues_by_chain,get_sequence
from ..utils.conservation_utils import load_msa_sto,load_msa_aligned_fasta,load_msa_a3m, write_sto,write_aligned_fasta,write_a3m
from typing import List, Dict,Any, Tuple
import subprocess
from rdkit import Chem
from ..utils.substrate_utils import is_valid_mol_3d, build_docked_mol_from_atom_info
from Bio.PDB import StructureBuilder
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
import copy
import numpy as np
from pdbfixer import PDBFixer
import re



def file_exists(path: str | Path) -> bool:
    p = Path(path)
    return p.exists() and p.is_file()

def get_stem(input_path: str | Path) -> str:
    p = Path(input_path)
    if is_fasta_gz_path(p):
        return p.name[:-len(".fasta.gz")]
    return p.stem

def is_fasta_gz_path(path: str | Path) -> bool:
    return Path(path).name.lower().endswith(".fasta.gz")

def open_fasta_text(path: str | Path):
    p = Path(path)
    if is_fasta_gz_path(p):
        return gzip.open(p, "rt", encoding="utf-8")
    return p.open("r", encoding="utf-8")

MAXFILENAME=150

def check_filename_length(name: str, logger: Logger) -> bool:
    if len(name) > MAXFILENAME:
        logger.print(f"[ERROR] Filename too long (>{MAXFILENAME}): {name}")
        return False
    return True

def load_protein_structure(path: str | Path, protein_name:str, logger: Logger) -> Structure | None:
    p = Path(path)

    try:
        if p.suffix.lower() in {".cif", ".mmcif"}:
            parser = MMCIFParser(QUIET=True)
        elif p.suffix.lower() == ".pdb":
            parser = PDBParser(QUIET=True)
        else:
            logger.print(f"[ERROR] Unsupported format: {p}")
            return None

        return parser.get_structure(protein_name, str(p))

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading structure for {str(p)}: {e}")
        return None

def load_openmm_structure(path: str | Path, logger: Logger) -> PDBFile | PDBxFile | None:
    p = Path(path)

    try:
        if p.suffix.lower() in {".cif", ".mmcif"}:
            return PDBxFile(str(p))
        elif p.suffix.lower() == ".pdb":
            return PDBFile(str(p))
        else:
            logger.print(f"[ERROR] Unsupported format: {p}")
            return None

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading OpenMM structure for {str(p)}: {e}")
        return None

def structure_to_pdbfile(struct: Structure, logger: Logger, protein_name: str = "structure") -> PDBFile | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            io = PDBIO()
            io.set_structure(struct)
            io.save(str(pdb_path))

            pdb = PDBFile(str(pdb_path))
            return pdb

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert Structure to PDBFile: {e}")
        return None

def pdbfile_to_structure(pdb: PDBFile, logger: Logger, protein_name: str = "structure") -> Structure | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            with open(pdb_path, "w") as f:
                PDBFile.writeFile(pdb.topology, pdb.positions, f)

            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(protein_name, str(pdb_path))

            return struct

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert PDBFile to Structure: {e}")
        return None

def modeller_to_structure(modeller: Modeller, logger: Logger, protein_name: str = "structure") -> Structure | None:
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            pdb_path = Path(tmp.name)

            with open(pdb_path, "w") as f:
                PDBFile.writeFile(modeller.topology, modeller.positions, f)

            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(protein_name, str(pdb_path))

            return struct

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert Modeller to Structure: {e}")
        return None

def load_dssp(struct: Structure, logger: Logger) -> DSSP | None:
    try:
        model = next(struct.get_models())

        with tempfile.TemporaryDirectory() as td:
            cif_path = Path(td) / "tmp.cif"

            io = MMCIFIO()
            io.set_structure(model)
            io.save(str(cif_path))

            dssp = DSSP(model,str(cif_path),dssp="mkdssp")
            return dssp

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading dssp: {e}")
        return None

def write_cif(struct: Structure, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = MMCIFIO()
    io.set_structure(struct)
    io.save(str(output_path))

def write_pdb(struct: Structure, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(struct)
    io.save(str(output_path))

def write_json_from_dict(dict_data: dict, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dict_data=convert_to_json_serializable(dict_data)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(dict_data, f, indent=2, ensure_ascii=False)

def write_json_from_dict_inline_leaf_lists(dict_data: dict, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dict_data = convert_to_json_serializable(dict_data)
    dict_data = wrap_leaf_lists_as_rawjson(dict_data)

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(
            dict_data,
            f,
            cls=InlineJSONEncoder,
            indent=2,
            ensure_ascii=False
        )

def write_fasta(struct: Structure, output_path: str | Path, logger: Logger) -> bool:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    chain = get_single_chain(struct, logger)
    if chain is None:
        return False

    residues = get_residues_by_chain(chain, logger)
    if residues is None:
        return False

    seq = get_sequence(residues, logger)
    if seq is None:
        return False

    header = str(struct.id).strip() if getattr(struct, "id", None) else None
    if header is None:
        logger.print("[ERROR] Structure header is None")
        return False

    try:
        with output_path.open("w", encoding="utf-8") as f:
            f.write(f">{header}\n")
            f.write(f"{seq}\n")
        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to write FASTA to {output_path}: {e}")
        return False

def load_msa(path: str | Path, logger: Logger) -> List[Dict[str, str]] | None:
    p = Path(path)

    try:
        suffix = p.suffix.lower()

        if suffix in {".sto", ".stockholm"}:
            return load_msa_sto(p, logger)

        elif suffix in {".fa", ".fasta", ".afa"} or is_fasta_gz_path(p):
            return load_msa_aligned_fasta(p, logger)

        elif suffix == ".a3m":
            return load_msa_a3m(p, logger)

        else:
            logger.print(f"[ERROR] Unsupported MSA format: {str(p)}")
            return None

    except Exception as e:
        logger.print(f"[ERROR] Exception in load_msa from {str(p)}: {e}")
        return None

def load_fasta(path: str | Path, logger: Logger) -> Dict[str, str] | None:
    p = Path(path)

    try:
        header: str | None = None
        seq_parts: list[str] = []

        with open_fasta_text(p) as f:
            for raw_line in f:
                line = raw_line.rstrip("\n").strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if header is not None:
                        logger.print(f"[ERROR] Multiple sequences found in FASTA file: {str(p)}")
                        return None

                    header = line[1:].strip()
                else:
                    if header is None:
                        logger.print(f"[ERROR] Invalid FASTA format in {str(p)}: sequence line appears before header.")
                        return None
                    seq_parts.append(line)

        if header is None:
            logger.print(f"[ERROR] No header found in FASTA file: {str(p)}")
            return None

        sequence = "".join(seq_parts)

        if sequence.strip() == "":
            logger.print(f"[ERROR] Empty sequence found in FASTA file: {str(p)}")
            return None

        return {"header": header, "sequence": sequence}

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading FASTA from {str(p)}: {e}")
        return None

def write_msa(msa_list: List[Dict[str, str]], output_path: str | Path, logger: Logger) -> bool:
    p = Path(output_path)

    try:
        suffix = p.suffix.lower()

        if suffix in {".sto", ".stockholm"}:
            return write_sto(msa_list, p, logger)

        elif suffix in {".fa", ".fasta", ".afa"}:
            return write_aligned_fasta(msa_list, p, logger)

        elif suffix == ".a3m":
            return write_a3m(msa_list, p, logger)

        else:
            logger.print(f"[ERROR] Unsupported MSA format: {str(p)}")
            return False

    except Exception as e:
        logger.print(f"[ERROR] Exception in write_msa: {e}")
        return False

def write_hmm(sto_path: str | Path, output_path: str | Path, logger: Logger) -> bool:
    sto_file = Path(sto_path)
    hmm_file = Path(output_path)

    try:
        if not sto_file.exists():
            logger.print(f"[ERROR] Input Stockholm file not found: {str(sto_file)}")
            return False

        hmm_file.parent.mkdir(parents=True, exist_ok=True)

        p = subprocess.run(
            [
                "hmmbuild",
                "--hand",
                str(hmm_file),
                str(sto_file),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if p.returncode != 0:
            logger.print(f"[ERROR] hmmbuild failed for {str(sto_file)}: {p.stderr.strip()}")
            return False

        if (not hmm_file.exists()) or hmm_file.stat().st_size == 0:
            logger.print(f"[ERROR] Failed to generate HMM file: {str(hmm_file)}")
            return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in write_hmm: {e}")
        return False


def write_sdf(mol_3d: Chem.Mol, sdf_path: str | Path, logger: Logger,) -> bool:
    if not is_valid_mol_3d(mol_3d, logger):
        return False

    try:
        sdf_path = Path(sdf_path)
        sdf_path.parent.mkdir(parents=True, exist_ok=True)

        writer = Chem.SDWriter(str(sdf_path))
        conf_id = mol_3d.GetConformer().GetId()
        writer.write(mol_3d, confId=conf_id)
        writer.close()

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to save SDF file.")
            return False

        return True
    except Exception as e:
        logger.print(f"[ERROR] Failed to save Mol(3D) to SDF file: {sdf_path}. Reason: {e}")
        return False


def save_substrate_structures(substrate_feature_list: List[Dict[str, Any]],output_dir: str | Path,logger: Logger) -> bool:

    if not isinstance(substrate_feature_list, list):
        logger.print("[ERROR] substrate_feature_list must be a list.")
        return False

    try:
        output_dir = Path(output_dir)

        tasks: List[Tuple[Chem.Mol, Path]] = []

        for item in substrate_feature_list:
            structures = item.get("structures", [])

            for s in structures:
                mol = s.get("structure_mol")
                name = s.get("structure_name")

                if not mol or not name:
                    logger.print("[ERROR] Invalid structure entry.")
                    return False

                clean_name = get_clean_filename(name)
                path = output_dir / get_optimized_filename(f"{clean_name}.sdf")

                tasks.append((mol, path))

        for mol, path in tasks:
            if not write_sdf(mol, path, logger):
                return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Failed to save substrate structures. Reason: {e}")
        return False

def write_protein_pdbqt(struct: Structure,pdbqt_path: str | Path,logger: Logger) -> bool:
    def parse_meeko_bad_residues_from_stderr(stderr: str) -> Tuple[List[str], List[str]]:
        direct_residue_set = set()
        range_residue_set = set()

        if not isinstance(stderr, str) or not stderr.strip():
            return [], []

        # Case 1:
        # matched with excess inter-residue bond(s): A:39
        for m in re.finditer(
                r"matched with excess inter-residue bond\(s\):\s*([A-Za-z0-9_]*:\d+)",
                stderr
        ):
            direct_residue_set.add(m.group(1))

        # Case 2:
        # No template matched for residue_key='A:836'
        for m in re.finditer(
                r"No template matched for residue_key='?([A-Za-z0-9_]*:\d+)'?",
                stderr
        ):
            direct_residue_set.add(m.group(1))

        # Case 3:
        # Expected 2 paddings for (A:37, A:39)
        for m in re.finditer(
                r"Expected\s+\d+\s+paddings\s+for\s+\(([A-Za-z0-9_]*:\d+),\s*([A-Za-z0-9_]*:\d+)\)",
                stderr
        ):
            res1 = m.group(1)
            res2 = m.group(2)

            chain1, num1_text = res1.split(":")
            chain2, num2_text = res2.split(":")

            try:
                num1 = int(num1_text)
                num2 = int(num2_text)
            except Exception:
                range_residue_set.add(res1)
                range_residue_set.add(res2)
                continue

            if chain1 == chain2:
                start = min(num1, num2)
                end = max(num1, num2)

                # 保守限制，避免一次删太大一段
                if end - start <= 5:
                    for n in range(start, end + 1):
                        range_residue_set.add(f"{chain1}:{n}")
                else:
                    range_residue_set.add(res1)
                    range_residue_set.add(res2)
            else:
                range_residue_set.add(res1)
                range_residue_set.add(res2)

        direct_residue_list = sorted(direct_residue_set)
        range_residue_list = sorted(range_residue_set)

        return direct_residue_list, range_residue_list

    def build_meeko_delete_residues_text(residue_list: List[str]) -> str:
        chain_to_resseq_list: Dict[str, List[int]] = {}

        for item in residue_list:
            if ":" not in item:
                continue

            chain_id, resseq_text = item.split(":", 1)

            try:
                resseq = int(resseq_text)
            except Exception:
                continue

            if chain_id not in chain_to_resseq_list:
                chain_to_resseq_list[chain_id] = []

            chain_to_resseq_list[chain_id].append(resseq)

        part_list: List[str] = []

        for chain_id in sorted(chain_to_resseq_list.keys()):
            unique_resseq_list = sorted(set(chain_to_resseq_list[chain_id]))
            if len(unique_resseq_list) == 0:
                continue

            resseq_text = ",".join(str(x) for x in unique_resseq_list)
            part_list.append(f"{chain_id}:{resseq_text}")

        return ",".join(part_list)

    try:
        pdbqt_path = Path(pdbqt_path)
        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_pdb = Path(tmp_dir) / "protein.pdb"

            # 写 PDB
            write_pdb(struct, tmp_pdb)

            if not tmp_pdb.exists() or tmp_pdb.stat().st_size <= 0:
                logger.print("[ERROR] Failed to write temporary PDB file.")
                return False

            # 调用 Meeko CLI
            p = subprocess.run(
                [
                    "mk_prepare_receptor.py",
                    "--read_pdb", str(tmp_pdb),
                    "--write_pdbqt", str(pdbqt_path),
                    "--allow_bad_res"
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )

            if p.returncode != 0:
                direct_residue_list, range_residue_list = parse_meeko_bad_residues_from_stderr(p.stderr)
                if len(direct_residue_list) > 0:
                    delete_text = build_meeko_delete_residues_text(direct_residue_list)
                    p2 = subprocess.run(
                        [
                            "mk_prepare_receptor.py",
                            "--read_pdb", str(tmp_pdb),
                            "--write_pdbqt", str(pdbqt_path),
                            "--allow_bad_res",
                            "--delete_residues", str(delete_text)
                        ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=False,
                    )
                    if p2.returncode == 0:
                        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
                            logger.print("[ERROR] Failed to generate PDBQT file.")
                            return False
                        return True

                if len(range_residue_list) > 0:
                    delete_text = build_meeko_delete_residues_text(range_residue_list)
                    p3 = subprocess.run(
                        [
                            "mk_prepare_receptor.py",
                            "--read_pdb", str(tmp_pdb),
                            "--write_pdbqt", str(pdbqt_path),
                            "--allow_bad_res",
                            "--delete_residues", str(delete_text)
                        ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=False,
                    )
                    if p3.returncode == 0:
                        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
                            logger.print("[ERROR] Failed to generate PDBQT file.")
                            return False

                        return True

                receptor_output = "\n".join(((p.stderr or p.stdout or "").strip()).splitlines()[-20:])
                if receptor_output:
                    logger.print(f"[ERROR] mk_prepare_receptor.py failed with return code {p.returncode}. Output tail: {receptor_output}")
                else:
                    logger.print(f"[ERROR] mk_prepare_receptor.py failed with return code {p.returncode}.")
                return False

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to generate PDBQT file.")
            return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert Structure to PDBQT: {pdbqt_path}. Reason: {e}")
        return False

def write_substrate_pdbqt_from_sdf(sdf_path: str | Path,pdbqt_path: str | Path,logger: Logger) -> bool:
    try:
        sdf_path = Path(sdf_path)
        pdbqt_path = Path(pdbqt_path)
        pdbqt_path.parent.mkdir(parents=True, exist_ok=True)

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Invalid input SDF file.")
            return False

        # 调用 Meeko CLI
        p = subprocess.run(
            [
                "mk_prepare_ligand.py",
                "-i", str(sdf_path),
                "-o", str(pdbqt_path),
                "--add_index_map",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        if p.returncode != 0:
            ligand_output = "\n".join(((p.stderr or p.stdout or "").strip()).splitlines()[-20:])
            if ligand_output:
                logger.print(f"[ERROR] mk_prepare_ligand.py failed with return code {p.returncode}. Output tail: {ligand_output}")
            else:
                logger.print(f"[ERROR] mk_prepare_ligand.py failed with return code {p.returncode}.")
            return False

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to generate PDBQT file.")
            return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert SDF to PDBQT: {sdf_path}. Reason: {e}")
        return False

def load_sdf_mol_3d(sdf_path: str | Path, logger: Logger) -> Chem.Mol | None:
    try:
        sdf_path = Path(sdf_path)

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Invalid input SDF file.")
            return None

        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        if supplier is None or len(supplier) == 0:
            logger.print("[ERROR] Failed to load SDF file.")
            return None

        mol = supplier[0]
        if mol is None:
            logger.print("[ERROR] Failed to parse Mol from SDF file.")
            return None

        if mol.GetNumConformers() <= 0:
            logger.print("[ERROR] Input SDF does not contain 3D coordinates.")
            return None

        return mol

    except Exception as e:
        logger.print(f"[ERROR] Failed to read Mol(3D) from SDF file: {sdf_path}. Reason: {e}")
        return None

def write_docked_sdf_from_atom_info(original_mol_3d: Chem.Mol,docked_atom_info_list: List[Dict[str, Any]],sdf_path: str | Path,logger: Logger) -> Chem.Mol | None:

    if original_mol_3d is None or original_mol_3d.GetNumConformers() <= 0:
        logger.print("[ERROR] Invalid original Mol(3D).")
        return None

    if not isinstance(docked_atom_info_list, list) or len(docked_atom_info_list) == 0:
        logger.print("[ERROR] Invalid docked_atom_info_list.")
        return None

    try:
        sdf_path = Path(sdf_path)
        sdf_path.parent.mkdir(parents=True, exist_ok=True)

        mol = build_docked_mol_from_atom_info(
            original_mol_3d,
            docked_atom_info_list,
            logger
        )
        if mol is None:
            return None

        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        if not sdf_path.exists() or sdf_path.stat().st_size <= 0:
            logger.print("[ERROR] Failed to save docked SDF file.")
            return None

        return mol

    except Exception as e:
        logger.print(f"[ERROR] Failed to write docked atom information to SDF file: {sdf_path}. Reason: {e}")
        return None

def write_docked_complex_from_mol_list(
    struct: Structure,
    docked_mol_list: List[Chem.Mol],
    protein_name: str,
    substrate_names: str,
    output_dir: str | Path,
    logger: Logger
) -> str | None:
    if struct is None:
        logger.print("[ERROR] struct is None.")
        return None

    if not isinstance(docked_mol_list, list) or len(docked_mol_list) == 0:
        logger.print("[ERROR] Invalid docked_mol_list.")
        return None

    if not isinstance(protein_name, str) or not protein_name.strip():
        logger.print("[ERROR] Invalid protein_name.")
        return None

    if not isinstance(substrate_names, str) or not substrate_names.strip():
        logger.print("[ERROR] Invalid substrate_names.")
        return None

    if not isinstance(output_dir, (str, Path)):
        logger.print("[ERROR] output_dir must be a str or Path.")
        return None

    try:
        protein_chain = get_single_chain(struct, logger)
        if protein_chain is None:
            return None

        residue_info_list = get_residues_by_chain(protein_chain, logger)
        if residue_info_list is None or len(residue_info_list) == 0:
            logger.print("[ERROR] Failed to get valid protein residues.")
            return None

        max_protein_resseq = max(res_id[1] for res_id, _, _ in residue_info_list)
        ligand_resseq_start = max_protein_resseq + 1

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        complex_name = f"docked_{protein_name}_{substrate_names}"
        complex_name = get_optimized_filename(complex_name)
        cif_path = output_dir / f"{complex_name}.cif"
        pdb_path = output_dir / f"{complex_name}.pdb"

        builder = StructureBuilder.StructureBuilder()
        builder.init_structure("complex")
        builder.init_model(0)
        builder.init_chain("A")

        new_struct = builder.get_structure()
        new_model: Model = new_struct[0]
        new_protein_chain: Chain = new_model["A"]

        for residue in protein_chain.get_residues():
            new_protein_chain.add(copy.deepcopy(residue))

        # === 计算 protein 当前最大 atom serial，ligand 从后面继续编号 ===
        max_serial = 0
        for atom in new_protein_chain.get_atoms():
            try:
                serial = int(atom.serial_number)
                if serial > max_serial:
                    max_serial = serial
            except Exception:
                continue

        current_serial = max_serial + 1

        ligand_chain = Chain("L")
        new_model.add(ligand_chain)

        for ligand_i, mol in enumerate(docked_mol_list, start=1):
            if mol is None or mol.GetNumConformers() <= 0:
                logger.print(f"[ERROR] Invalid docked Mol for ligand index {ligand_i}.")
                return None

            conf = mol.GetConformer()
            ligand_res_id = ("H", ligand_resseq_start + ligand_i - 1, " ")
            ligand_residue = Residue(ligand_res_id, "LIG", " ")

            element_count_dict: Dict[str, int] = {}

            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol().upper().strip()
                if not symbol:
                    symbol = "X"

                element_count_dict[symbol] = element_count_dict.get(symbol, 0) + 1
                atom_name = f"{symbol}{element_count_dict[symbol]}"
                atom_name = atom_name[:4]
                fullname = atom_name.rjust(4)

                pos = conf.GetAtomPosition(atom.GetIdx())
                serial_number = current_serial
                current_serial += 1

                pdb_atom = Atom(
                    name=atom_name,
                    coord=np.array([float(pos.x), float(pos.y), float(pos.z)], dtype=float),
                    bfactor=1.0,
                    occupancy=1.0,
                    altloc=" ",
                    fullname=fullname,
                    serial_number=serial_number,
                    element=symbol,
                )
                ligand_residue.add(pdb_atom)

            ligand_chain.add(ligand_residue)

        write_cif(new_struct, cif_path)
        write_pdb(new_struct, pdb_path)

        if not cif_path.exists() or cif_path.stat().st_size == 0:
            logger.print("[ERROR] Failed to write complex CIF.")
            return None

        if not pdb_path.exists() or pdb_path.stat().st_size == 0:
            logger.print("[ERROR] Failed to write complex PDB.")
            return None

        return str(cif_path)

    except Exception as e:
        logger.print(f"[ERROR] Failed to build docked complex CIF/PDB from Mol list. Reason: {e}")
        return None

def load_openmm_modeller(path: str | Path, logger) -> Modeller | None:
    if not isinstance(path, (str, Path)):
        logger.print("[ERROR] path must be a str or Path.")
        return None

    try:
        p = Path(path)
    except Exception as e:
        logger.print(f"[ERROR] Failed to parse path: {e}")
        return None

    if not p.exists() or p.stat().st_size <= 0:
        logger.print(f"[ERROR] Invalid input structure file: {p}")
        return None

    try:
        suffix = p.suffix.lower()

        if suffix in {".cif", ".mmcif"}:
            obj = PDBxFile(str(p))
        elif suffix == ".pdb":
            obj = PDBFile(str(p))
        else:
            logger.print(f"[ERROR] Unsupported structure format: {p}")
            return None

        return Modeller(obj.topology, obj.positions)

    except Exception as e:
        logger.print(f"[ERROR] Failed to load OpenMM Modeller from {str(p)}: {e}")
        return None

def load_substrate_name_and_mol_3d_list(substrate_names: str,substrate_dir: str | Path,logger: Logger) -> Tuple[List[str], List[Chem.Mol]] | None:
    if not isinstance(substrate_names, str) or not substrate_names.strip():
        logger.print("[ERROR] substrate_names is empty.")
        return None

    if not isinstance(substrate_dir, (str, Path)):
        logger.print("[ERROR] substrate_dir must be a str or Path.")
        return None

    try:
        substrate_dir = Path(substrate_dir)
    except Exception as e:
        logger.print(f"[ERROR] Failed to parse substrate_dir: {e}")
        return None

    if not substrate_dir.exists() or not substrate_dir.is_dir():
        logger.print(f"[ERROR] Invalid substrate_dir: {substrate_dir}")
        return None

    substrate_name_list = [x.strip() for x in str(substrate_names).split(";")]
    if len(substrate_name_list) == 0:
        logger.print("[ERROR] substrate_names is empty.")
        return None

    if any(not x for x in substrate_name_list):
        logger.print("[ERROR] substrate_names contains empty substrate name.")
        return None

    if len(set(substrate_name_list)) != len(substrate_name_list):
        logger.print("[ERROR] Duplicate substrate names are not allowed.")
        return None

    mol_3d_list: List[Chem.Mol] = []

    for substrate_name in substrate_name_list:
        substrate_file_stem = get_optimized_filename(substrate_name)
        if not substrate_file_stem:
            logger.print(f"[ERROR] Invalid substrate filename generated from substrate: {substrate_name}")
            return None

        sdf_path = substrate_dir / f"{substrate_file_stem}.sdf"

        if not sdf_path.exists() or not sdf_path.is_file():
            logger.print(f"[ERROR] Substrate SDF not found: {sdf_path}")
            return None

        mol_3d = load_sdf_mol_3d(sdf_path, logger)
        if mol_3d is None:
            logger.print(f"[ERROR] Failed to load Mol(3D) from SDF: {sdf_path}")
            return None

        mol_3d_list.append(mol_3d)

    return substrate_name_list, mol_3d_list

def load_pdbfixer(path: str | Path,logger: Logger) -> PDBFixer | None:

    p = Path(path)

    try:
        if not p.exists():
            logger.print(f"[ERROR] File not found: {str(p)}")
            return None

        if p.stat().st_size == 0:
            logger.print(f"[ERROR] File is empty: {str(p)}")
            return None

        suffix = p.suffix.lower()

        if suffix not in {".pdb", ".cif", ".mmcif"}:
            logger.print(f"[ERROR] Unsupported format: {str(p)}")
            return None

        fixer = PDBFixer(filename=str(p))

        if fixer.topology is None:
            logger.print(f"[ERROR] Failed to load topology: {str(p)}")
            return None

        return fixer

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading PDBFixer for {str(p)}: {e}")
        return None


def write_cif_from_pdbfixer(fixer: PDBFixer, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        PDBxFile.writeFile(
            fixer.topology,
            fixer.positions,
            f,
            keepIds=True
        )

def write_pdb_from_pdbfixer(fixer: PDBFixer, output_path: str | Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            f,
            keepIds=True
        )

def pdbfixer_to_modeller(fixer: PDBFixer,logger: Logger) -> Modeller | None:

    if fixer is None:
        logger.print("[ERROR] fixer is None.")
        return None

    if not isinstance(fixer, PDBFixer):
        logger.print("[ERROR] fixer must be a PDBFixer instance.")
        return None

    try:
        topology = fixer.topology
        positions = fixer.positions

        if topology is None:
            logger.print("[ERROR] fixer.topology is None.")
            return None

        if positions is None:
            logger.print("[ERROR] fixer.positions is None.")
            return None
        modeller = Modeller(topology, positions)

        return modeller

    except Exception as e:
        logger.print(f"[ERROR] Failed to convert PDBFixer to Modeller: {e}")
        return None
