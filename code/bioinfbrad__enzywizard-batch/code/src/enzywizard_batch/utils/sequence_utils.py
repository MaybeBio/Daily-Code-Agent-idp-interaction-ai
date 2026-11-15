from __future__ import annotations
from pathlib import Path
from typing import List, Dict, Any

from Bio.Data.IUPACData import protein_letters_3to1

from ..utils.logging_utils import Logger
from ..utils.conservation_utils import check_msa_sto,check_msa_aligned_fasta,check_msa_a3m, clean_sto,clean_aligned_fasta,clean_a3m, remove_a3m_insertions, is_all_gap, postprocess_conservation_report_to_schema

def is_fasta_gz_path(path: str | Path) -> bool:
    return Path(path).name.lower().endswith(".fasta.gz")

def check_msa(path: str | Path, sequence_dict: Dict[str, str], msa_list: List[Dict[str, str]], logger: Logger) -> bool:
    p = Path(path)
    query_sequence = sequence_dict.get("sequence",None)

    try:
        suffix = p.suffix.lower()

        if suffix in {".sto", ".stockholm"}:
            return check_msa_sto(query_sequence, msa_list, logger)

        elif suffix in {".fa", ".fasta", ".afa"} or is_fasta_gz_path(p):
            return check_msa_aligned_fasta(query_sequence, msa_list, logger)

        elif suffix == ".a3m":
            return check_msa_a3m(query_sequence, msa_list, logger)

        else:
            logger.print(f"[ERROR] Unsupported MSA format: {str(p)}")
            return False

    except Exception as e:
        logger.print(f"[ERROR] Exception in check_msa for {str(p)}: {e}")
        return False

def clean_msa(path: str | Path, msa_list: List[Dict[str, str]], logger: Logger) -> List[Dict[str, str]] | None:
    p = Path(path)

    try:
        suffix = p.suffix.lower()

        if suffix in {".sto", ".stockholm"}:
            return clean_sto(msa_list, logger)

        elif suffix in {".fa", ".fasta", ".afa"} or is_fasta_gz_path(p):
            return clean_aligned_fasta(msa_list, logger)

        elif suffix == ".a3m":
            return clean_a3m(msa_list, logger)

        else:
            logger.print(f"[ERROR] Unsupported MSA format: {suffix}")
            return None

    except Exception as e:
        logger.print(f"[ERROR] Exception in clean_msa: {e}")
        return None

def clean_msa_to_sto(msa_list: List[Dict[str, str]], logger: Logger) -> List[Dict[str, str]] | None:
    try:
        if not isinstance(msa_list, list) or len(msa_list) == 0:
            logger.print("[ERROR] msa_list is empty or invalid.")
            return None

        valid_chars_sto = set("ACDEFGHIKLMNPQRSTVWY-.")
        valid_chars_a3m = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.")
        cleaned: List[Dict[str, str]] = []

        first_header = msa_list[0].get("header")
        first_seq = msa_list[0].get("sequence")

        if first_header is None or first_seq is None:
            logger.print("[ERROR] First record is invalid.")
            return None

        first_header = first_header.strip()
        first_seq = first_seq.strip()

        if first_header == "" or first_seq == "":
            logger.print("[ERROR] First record is invalid.")
            return None

        first_header = first_header.split()[0]
        first_seq = "".join(ch for ch in first_seq if ch in valid_chars_a3m)
        first_seq = remove_a3m_insertions(first_seq).upper()

        if first_seq == "":
            logger.print("[ERROR] First sequence is invalid after cleaning.")
            return None

        expected_len = len(first_seq)
        if expected_len == 0:
            logger.print("[ERROR] First sequence is invalid after cleaning.")
            return None

        if is_all_gap(first_seq):
            logger.print("[ERROR] First sequence becomes all-gap after cleaning.")
            return None

        header_set = set()

        for record in msa_list:
            if not isinstance(record, dict):
                continue

            header = record.get("header")
            seq = record.get("sequence")

            if header is None or seq is None:
                continue

            header = header.strip()
            seq = seq.strip()

            if header == "" or seq == "":
                continue

            header = header.split()[0]

            if header in header_set:
                continue
            header_set.add(header)

            seq = "".join(ch for ch in seq if ch in valid_chars_a3m)
            seq = remove_a3m_insertions(seq).upper()

            if seq == "":
                continue

            for ch in seq:
                if ch not in valid_chars_sto:
                    seq = ""
                    break

            if seq == "":
                continue

            if len(seq) != expected_len:
                continue

            if is_all_gap(seq):
                continue

            cleaned.append({"header": header, "sequence": seq})

        if len(cleaned) == 0:
            logger.print("[ERROR] No valid sequences remain after cleaning.")
            return None

        if first_header not in [r["header"] for r in cleaned]:
            logger.print("[ERROR] Query sequence removed during cleaning.")
            return None

        return cleaned

    except Exception as e:
        logger.print(f"[ERROR] Exception in clean_msa_to_sto: {e}")
        return None


def normalize_aa_name_to_one_letter(aa_name: Any) -> str:

    aa_name_clean = aa_name.strip()

    if len(aa_name_clean) == 1:
        return aa_name_clean.upper()

    if len(aa_name_clean) == 3:
        aa_name_3 = aa_name_clean.upper().capitalize()
        aa_name_1 = protein_letters_3to1.get(aa_name_3)
        if isinstance(aa_name_1, str) and aa_name_1 != "":
            return aa_name_1.upper()
    return "X"