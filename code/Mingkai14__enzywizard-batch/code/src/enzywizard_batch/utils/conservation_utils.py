from __future__ import annotations

import gzip
from pathlib import Path
from typing import List, Dict, Any
import math

from ..utils.logging_utils import Logger

def is_fasta_gz_path(path: str | Path) -> bool:
    return Path(path).name.lower().endswith(".fasta.gz")

def open_text_or_fasta_gz(path: str | Path):
    p = Path(path)
    if is_fasta_gz_path(p):
        return gzip.open(p, "rt", encoding="utf-8")
    return p.open("r", encoding="utf-8")

def remove_gaps(seq: str) -> str:
    return seq.replace("-", "").replace(".", "")

def remove_a3m_insertions(seq: str) -> str:
    return "".join(ch for ch in seq if not ch.islower())

def remove_a3m_insertions_and_gaps(seq: str) -> str:
    return "".join(ch for ch in seq if (not ch.islower()) and ch not in {"-", "."})

def is_all_gap(seq: str) -> bool:
    return all(ch in {"-", "."} for ch in seq)


def load_msa_sto(path: str | Path, logger: Logger) -> List[Dict[str, str]] | None:
    p = Path(path)

    try:
        msa_list: List[Dict[str, str]] = []
        seq_dict: Dict[str, str] = {}
        header_order: List[str] = []

        with p.open("r", encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.rstrip("\n").strip()

                if not line:
                    continue

                if line == "//":
                    break

                if line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) < 2:
                    logger.print(f"[ERROR] Invalid Stockholm line in {str(p)}: {line}")
                    return None

                header = parts[0]
                seq_part = parts[1]

                if header not in seq_dict:
                    seq_dict[header] = seq_part
                    header_order.append(header)
                else:
                    seq_dict[header] += seq_part

        if len(header_order) == 0:
            logger.print(f"[ERROR] No sequence found in Stockholm file: {str(p)}")
            return None

        for header in header_order:
            msa_list.append({"header": header, "sequence": seq_dict[header]})

        return msa_list

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading Stockholm MSA from {str(p)}: {e}")
        return None

def load_msa_aligned_fasta(path: str | Path, logger: Logger) -> List[Dict[str, str]] | None:
    p = Path(path)

    try:
        msa_list: List[Dict[str, str]] = []
        current_header: str | None = None
        current_seq_parts: List[str] = []

        with open_text_or_fasta_gz(p) as f:
            for raw_line in f:
                line = raw_line.rstrip("\n").strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if current_header is not None:
                        msa_list.append({
                            "header": current_header,
                            "sequence": "".join(current_seq_parts)
                        })

                    current_header = line[1:].strip()
                    current_seq_parts = []
                else:
                    if current_header is None:
                        logger.print(f"[ERROR] Invalid FASTA format in {str(p)}: sequence line appears before header.")
                        return None
                    current_seq_parts.append(line)

        if current_header is not None:
            msa_list.append({
                "header": current_header,
                "sequence": "".join(current_seq_parts)
            })

        if len(msa_list) == 0:
            logger.print(f"[ERROR] No sequence found in aligned FASTA file: {str(p)}")
            return None

        return msa_list
    except Exception as e:
        logger.print(f"[ERROR] Exception in loading aligned FASTA MSA from {str(p)}: {e}")
        return None

def load_msa_a3m(path: str | Path, logger: Logger) -> List[Dict[str, str]] | None:
    p = Path(path)

    try:
        msa_list: List[Dict[str, str]] = []
        current_header: str | None = None
        current_seq_parts: List[str] = []

        with p.open("r", encoding="utf-8") as f:
            for raw_line in f:
                line = raw_line.rstrip("\n").strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if current_header is not None:
                        msa_list.append({
                            "header": current_header,
                            "sequence": "".join(current_seq_parts)
                        })

                    current_header = line[1:].strip()
                    current_seq_parts = []
                else:
                    if current_header is None:
                        logger.print(f"[ERROR] Invalid A3M format in {str(p)}: sequence line appears before header.")
                        return None
                    current_seq_parts.append(line)

        if current_header is not None:
            msa_list.append({
                "header": current_header,
                "sequence": "".join(current_seq_parts)
            })

        if len(msa_list) == 0:
            logger.print(f"[ERROR] No sequence found in A3M file: {str(p)}")
            return None

        return msa_list

    except Exception as e:
        logger.print(f"[ERROR] Exception in loading A3M MSA from {str(p)}: {e}")
        return None

def check_msa_sto(query_sequence: str, msa_list: List[Dict[str, str]], logger: Logger) -> bool:
    if not isinstance(msa_list, list):
        logger.print("[ERROR] msa_list is not a list.")
        return False

    if len(msa_list) == 0:
        logger.print("[ERROR] msa_list is empty.")
        return False

    if query_sequence is None:
        logger.print("[ERROR] query_sequence is empty.")
        return False

    query_sequence = query_sequence.strip()
    if query_sequence == "":
        logger.print("[ERROR] query_sequence is empty.")
        return False
    try:

        valid_chars = set("ACDEFGHIKLMNPQRSTVWYBXZUO-.")

        first_seq = msa_list[0]["sequence"]
        first_seq_ungapped = remove_gaps(first_seq)

        if first_seq_ungapped != query_sequence:
            logger.print("[ERROR] The first Stockholm MSA sequence does not match query_sequence after gap removal.")
            return False

        expected_len = len(first_seq)
        for record in msa_list:
            header = record.get("header")
            seq = record.get("sequence")

            if header is None or header.strip() == "":
                logger.print("[ERROR] Empty header found in MSA.")
                return False

            if seq is None or seq.strip() == "":
                logger.print(f"[ERROR] Empty sequence found for header {header}.")
                return False
            for ch in seq:
                if ch not in valid_chars:
                    logger.print(f"[ERROR] Invalid character '{ch}' in sequence for header {header}")
                    return False

            if len(seq) != expected_len:
                logger.print(f"[ERROR] Inconsistent sequence length in Stockholm MSA for header {header}.")
                return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in checking Stockholm MSA: {e}")
        return False


def check_msa_aligned_fasta(query_sequence: str, msa_list: List[Dict[str, str]], logger: Logger) -> bool:
    if not isinstance(msa_list, list):
        logger.print("[ERROR] msa_list is not a list.")
        return False

    if len(msa_list) == 0:
        logger.print("[ERROR] msa_list is empty.")
        return False

    if query_sequence is None:
        logger.print("[ERROR] query_sequence is empty.")
        return False

    query_sequence = query_sequence.strip()
    if query_sequence == "":
        logger.print("[ERROR] query_sequence is empty.")
        return False
    try:

        valid_chars = set("ACDEFGHIKLMNPQRSTVWYBXZUO-.")

        first_seq = msa_list[0]["sequence"]
        first_seq_ungapped = remove_gaps(first_seq)

        if first_seq_ungapped != query_sequence:
            logger.print("[ERROR] The first aligned FASTA MSA sequence does not match query_sequence after gap removal.")
            return False

        expected_len = len(first_seq)
        for record in msa_list:
            header = record.get("header")
            seq = record.get("sequence")

            if header is None or header.strip() == "":
                logger.print("[ERROR] Empty header found in MSA.")
                return False

            if seq is None or seq.strip() == "":
                logger.print(f"[ERROR] Empty sequence found for header {header}.")
                return False

            for ch in seq:
                if ch not in valid_chars:
                    logger.print(f"[ERROR] Invalid character '{ch}' in sequence for header {header}")
                    return False

            if len(seq) != expected_len:
                logger.print(f"[ERROR] Inconsistent sequence length in aligned FASTA MSA for header {header}.")
                return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in checking aligned FASTA MSA: {e}")
        return False


def check_msa_a3m(query_sequence: str, msa_list: List[Dict[str, str]], logger: Logger) -> bool:
    if not isinstance(msa_list, list):
        logger.print("[ERROR] msa_list is not a list.")
        return False

    if len(msa_list) == 0:
        logger.print("[ERROR] msa_list is empty.")
        return False

    if query_sequence is None:
        logger.print("[ERROR] query_sequence is empty.")
        return False

    query_sequence = query_sequence.strip()
    if query_sequence == "":
        logger.print("[ERROR] query_sequence is empty.")
        return False
    try:

        valid_chars = set("ACDEFGHIKLMNPQRSTVWYBXZUO-.abcdefghijklmnopqrstuvwxyz")

        first_seq = msa_list[0]["sequence"]
        first_seq_query_like = remove_a3m_insertions_and_gaps(first_seq)

        if first_seq_query_like != query_sequence:
            logger.print("[ERROR] The first A3M MSA sequence does not match query_sequence after removing lowercase insertions and gaps.")
            return False

        expected_aligned_len = len(remove_a3m_insertions(first_seq))
        for record in msa_list:
            header = record.get("header")
            seq = record.get("sequence")

            if header is None or header.strip() == "":
                logger.print("[ERROR] Empty header found in MSA.")
                return False

            if seq is None or seq.strip() == "":
                logger.print(f"[ERROR] Empty sequence found for header {header}.")
                return False
            if is_all_gap(seq):
                logger.print(f"[ERROR] Sequence contains only gaps for header {header}.")
                return False
            for ch in seq:
                if ch not in valid_chars:
                    logger.print(f"[ERROR] Invalid character '{ch}' in sequence for header {header}")
                    return False

            aligned_seq = remove_a3m_insertions(seq)

            if len(aligned_seq) != expected_aligned_len:
                logger.print(f"[ERROR] Inconsistent aligned length in A3M MSA for header {header} after removing lowercase insertions.")
                return False

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in checking A3M MSA: {e}")
        return False



# clean


def write_sto(msa_list: List[Dict[str, str]], output_path: str | Path, logger: Logger) -> bool:
    p = Path(output_path)

    try:
        if not isinstance(msa_list, list):
            logger.print("[ERROR] msa_list is not a list.")
            return False

        if len(msa_list) == 0:
            logger.print("[ERROR] msa_list is empty.")
            return False

        valid_chars = set("ACDEFGHIKLMNPQRSTVWY-.")
        first_record = msa_list[0]

        first_header = first_record.get("header")
        first_seq = first_record.get("sequence")

        if first_header is None or not isinstance(first_header, str) or first_header.strip() == "":
            logger.print("[ERROR] First header is empty or invalid.")
            return False

        if any(ch.isspace() for ch in first_header):
            logger.print(f"[ERROR] Header contains whitespace: {first_header}")
            return False

        if first_seq is None or not isinstance(first_seq, str) or first_seq.strip() == "":
            logger.print("[ERROR] First sequence is empty or invalid.")
            return False

        expected_len = len(first_seq)
        if expected_len == 0:
            logger.print("[ERROR] First sequence is empty.")
            return False

        for ch in first_seq:
            if ch not in valid_chars:
                logger.print(f"[ERROR] Invalid character '{ch}' in first sequence for header {first_header}")
                return False

        header_set = set()
        rf_chars: List[str] = []
        for ch in first_seq:
            if ch in {"-", "."}:
                rf_chars.append(".")
            else:
                rf_chars.append("x")
        rf = "".join(rf_chars)

        with p.open("w", encoding="utf-8") as f:
            f.write("# STOCKHOLM 1.0\n")
            f.write("\n")

            for record in msa_list:
                if not isinstance(record, dict):
                    logger.print("[ERROR] Invalid record found in msa_list.")
                    return False

                header = record.get("header")
                seq = record.get("sequence")

                if header is None or not isinstance(header, str) or header.strip() == "":
                    logger.print("[ERROR] Empty header found in msa_list.")
                    return False

                if any(ch.isspace() for ch in header):
                    logger.print(f"[ERROR] Header contains whitespace: {header}")
                    return False

                if header in header_set:
                    logger.print(f"[ERROR] Duplicate header found in msa_list: {header}")
                    return False
                header_set.add(header)

                if seq is None or not isinstance(seq, str) or seq.strip() == "":
                    logger.print(f"[ERROR] Empty sequence found for header {header}.")
                    return False

                for ch in seq:
                    if ch not in valid_chars:
                        logger.print(f"[ERROR] Invalid character '{ch}' in sequence for header {header}")
                        return False

                if len(seq) != expected_len:
                    logger.print(f"[ERROR] Inconsistent sequence length for header {header}.")
                    return False

                f.write(f"{header} {seq}\n")

            f.write(f"#=GC RF {rf}\n")
            f.write("//\n")

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in write_sto: {e}")
        return False


def write_aligned_fasta(msa_list: List[Dict[str, str]], output_path: str | Path, logger: Logger) -> bool:
    p = Path(output_path)

    try:
        if not isinstance(msa_list, list):
            logger.print("[ERROR] msa_list is not a list.")
            return False

        if len(msa_list) == 0:
            logger.print("[ERROR] msa_list is empty.")
            return False

        valid_chars = set("ACDEFGHIKLMNPQRSTVWY-.")
        first_record = msa_list[0]

        first_header = first_record.get("header")
        first_seq = first_record.get("sequence")

        if first_header is None or not isinstance(first_header, str) or first_header.strip() == "":
            logger.print("[ERROR] First header is empty or invalid.")
            return False

        if first_seq is None or not isinstance(first_seq, str) or first_seq.strip() == "":
            logger.print("[ERROR] First sequence is empty or invalid.")
            return False

        expected_len = len(first_seq)
        if expected_len == 0:
            logger.print("[ERROR] First sequence is empty.")
            return False

        for ch in first_seq:
            if ch not in valid_chars:
                logger.print(f"[ERROR] Invalid character '{ch}' in first sequence for header {first_header}")
                return False

        header_set = set()

        with p.open("w", encoding="utf-8") as f:
            for record in msa_list:
                if not isinstance(record, dict):
                    logger.print("[ERROR] Invalid record found in msa_list.")
                    return False

                header = record.get("header")
                seq = record.get("sequence")

                if header is None or not isinstance(header, str) or header.strip() == "":
                    logger.print("[ERROR] Empty header found in msa_list.")
                    return False

                if header in header_set:
                    logger.print(f"[ERROR] Duplicate header found in msa_list: {header}")
                    return False
                header_set.add(header)

                if seq is None or not isinstance(seq, str) or seq.strip() == "":
                    logger.print(f"[ERROR] Empty sequence found for header {header}.")
                    return False

                for ch in seq:
                    if ch not in valid_chars:
                        logger.print(f"[ERROR] Invalid character '{ch}' in sequence for header {header}")
                        return False

                if len(seq) != expected_len:
                    logger.print(f"[ERROR] Inconsistent sequence length for header {header}.")
                    return False

                f.write(f">{header}\n")
                f.write(f"{seq}\n")

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in write_aligned_fasta: {e}")
        return False


def write_a3m(msa_list: List[Dict[str, str]], output_path: str | Path, logger: Logger) -> bool:
    p = Path(output_path)

    try:
        if not isinstance(msa_list, list):
            logger.print("[ERROR] msa_list is not a list.")
            return False

        if len(msa_list) == 0:
            logger.print("[ERROR] msa_list is empty.")
            return False

        valid_chars = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.")
        first_record = msa_list[0]

        first_header = first_record.get("header")
        first_seq = first_record.get("sequence")

        if first_header is None or not isinstance(first_header, str) or first_header.strip() == "":
            logger.print("[ERROR] First header is empty or invalid.")
            return False

        if first_seq is None or not isinstance(first_seq, str) or first_seq.strip() == "":
            logger.print("[ERROR] First sequence is empty or invalid.")
            return False

        first_aligned_seq = "".join(ch for ch in first_seq if not ch.islower())
        expected_aligned_len = len(first_aligned_seq)
        if expected_aligned_len == 0:
            logger.print("[ERROR] First aligned sequence is empty after removing lowercase insertions.")
            return False

        for ch in first_seq:
            if ch not in valid_chars:
                logger.print(f"[ERROR] Invalid character '{ch}' in first sequence for header {first_header}")
                return False

        header_set = set()

        with p.open("w", encoding="utf-8") as f:
            for record in msa_list:
                if not isinstance(record, dict):
                    logger.print("[ERROR] Invalid record found in msa_list.")
                    return False

                header = record.get("header")
                seq = record.get("sequence")

                if header is None or not isinstance(header, str) or header.strip() == "":
                    logger.print("[ERROR] Empty header found in msa_list.")
                    return False

                if header in header_set:
                    logger.print(f"[ERROR] Duplicate header found in msa_list: {header}")
                    return False
                header_set.add(header)

                if seq is None or not isinstance(seq, str) or seq.strip() == "":
                    logger.print(f"[ERROR] Empty sequence found for header {header}.")
                    return False

                for ch in seq:
                    if ch not in valid_chars:
                        logger.print(f"[ERROR] Invalid character '{ch}' in sequence for header {header}")
                        return False

                aligned_seq = "".join(ch for ch in seq if not ch.islower())
                if len(aligned_seq) != expected_aligned_len:
                    logger.print(f"[ERROR] Inconsistent aligned length in A3M for header {header}.")
                    return False

                f.write(f">{header}\n")
                f.write(f"{seq}\n")

        return True

    except Exception as e:
        logger.print(f"[ERROR] Exception in write_a3m: {e}")
        return False


def clean_sto(msa_list: List[Dict[str, str]], logger: Logger) -> List[Dict[str, str]] | None:
    try:
        if not isinstance(msa_list, list) or len(msa_list) == 0:
            logger.print("[ERROR] msa_list is empty or invalid.")
            return None

        valid_chars = set("ACDEFGHIKLMNPQRSTVWY-.")
        cleaned: List[Dict[str, str]] = []

        first_seq = msa_list[0].get("sequence")
        if first_seq is None:
            logger.print("[ERROR] First sequence invalid.")
            return None

        first_seq = first_seq.strip().upper()
        if first_seq == "":
            logger.print("[ERROR] First sequence invalid.")
            return None

        first_seq = "".join(ch for ch in first_seq if ch in valid_chars)
        if first_seq == "":
            logger.print("[ERROR] First sequence invalid after cleaning.")
            return None

        expected_len = len(first_seq)

        header_set = set()

        for record in msa_list:
            header = record.get("header")
            seq = record.get("sequence")

            if header is None or seq is None:
                continue

            header = header.strip()
            seq = seq.strip().upper()

            if header == "" or seq == "":
                continue

            header = header.split()[0]

            if header in header_set:
                continue
            header_set.add(header)

            seq = "".join(ch for ch in seq if ch in valid_chars)

            if len(seq) != expected_len:
                continue

            if all(ch in {"-", "."} for ch in seq):
                continue

            cleaned.append({"header": header, "sequence": seq})

        if len(cleaned) == 0:
            logger.print("[ERROR] No valid sequences remain after cleaning.")
            return None

        if cleaned[0]["sequence"].strip() == "":
            logger.print("[ERROR] Query sequence lost after cleaning.")
            return None

        query_header = msa_list[0].get("header")
        if query_header is None:
            logger.print("[ERROR] First header invalid.")
            return None

        query_header = query_header.strip()
        if query_header == "":
            logger.print("[ERROR] First header invalid.")
            return None

        query_header = query_header.split()[0]

        if query_header not in [r["header"] for r in cleaned]:
            logger.print("[ERROR] Query sequence removed during cleaning.")
            return None

        return cleaned

    except Exception as e:
        logger.print(f"[ERROR] Exception in clean_sto: {e}")
        return None


def clean_aligned_fasta(msa_list: List[Dict[str, str]], logger: Logger) -> List[Dict[str, str]] | None:
    try:
        if not isinstance(msa_list, list) or len(msa_list) == 0:
            logger.print("[ERROR] msa_list is empty or invalid.")
            return None

        valid_chars = set("ACDEFGHIKLMNPQRSTVWY-.")
        cleaned: List[Dict[str, str]] = []

        first_seq = msa_list[0].get("sequence")
        if first_seq is None:
            logger.print("[ERROR] First sequence invalid.")
            return None

        first_seq = first_seq.strip().upper()
        if first_seq == "":
            logger.print("[ERROR] First sequence invalid.")
            return None

        first_seq = "".join(ch for ch in first_seq if ch in valid_chars)
        if first_seq == "":
            logger.print("[ERROR] First sequence invalid after cleaning.")
            return None

        expected_len = len(first_seq)

        header_set = set()

        for record in msa_list:
            header = record.get("header")
            seq = record.get("sequence")

            if header is None or seq is None:
                continue

            header = header.strip()
            seq = seq.strip().upper()

            if header == "" or seq == "":
                continue

            header = header.split()[0]

            if header in header_set:
                continue
            header_set.add(header)

            seq = "".join(ch for ch in seq if ch in valid_chars)

            if len(seq) != expected_len:
                continue

            if all(ch in {"-", "."} for ch in seq):
                continue

            cleaned.append({"header": header, "sequence": seq})

        if len(cleaned) == 0:
            logger.print("[ERROR] No valid sequences remain after cleaning.")
            return None

        if cleaned[0]["sequence"].strip() == "":
            logger.print("[ERROR] Query sequence lost after cleaning.")
            return None

        query_header = msa_list[0].get("header")
        if query_header is None:
            logger.print("[ERROR] First header invalid.")
            return None

        query_header = query_header.strip()
        if query_header == "":
            logger.print("[ERROR] First header invalid.")
            return None

        query_header = query_header.split()[0]

        if query_header not in [r["header"] for r in cleaned]:
            logger.print("[ERROR] Query sequence removed during cleaning.")
            return None

        return cleaned

    except Exception as e:
        logger.print(f"[ERROR] Exception in clean_aligned_fasta: {e}")
        return None


def clean_a3m(msa_list: List[Dict[str, str]], logger: Logger) -> List[Dict[str, str]] | None:
    try:
        if not isinstance(msa_list, list) or len(msa_list) == 0:
            logger.print("[ERROR] msa_list is empty or invalid.")
            return None

        valid_chars = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.")
        cleaned: List[Dict[str, str]] = []

        first_seq = msa_list[0].get("sequence")
        if first_seq is None:
            logger.print("[ERROR] First sequence invalid.")
            return None

        first_seq = first_seq.strip()
        if first_seq == "":
            logger.print("[ERROR] First sequence invalid.")
            return None

        first_seq = "".join(ch for ch in first_seq if ch in valid_chars)
        if first_seq == "":
            logger.print("[ERROR] First sequence invalid after cleaning.")
            return None

        expected_len = len(remove_a3m_insertions(first_seq))

        header_set = set()

        for record in msa_list:
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

            seq = "".join(ch for ch in seq if ch in valid_chars)

            aligned_seq = remove_a3m_insertions(seq)

            if len(aligned_seq) != expected_len:
                continue

            if all(ch in {"-", "."} for ch in aligned_seq):
                continue

            cleaned.append({"header": header, "sequence": seq})

        if len(cleaned) == 0:
            logger.print("[ERROR] No valid sequences remain after cleaning.")
            return None

        if cleaned[0]["sequence"].strip() == "":
            logger.print("[ERROR] Query sequence lost after cleaning.")
            return None

        query_header = msa_list[0].get("header")
        if query_header is None:
            logger.print("[ERROR] First header invalid.")
            return None

        query_header = query_header.strip()
        if query_header == "":
            logger.print("[ERROR] First header invalid.")
            return None

        query_header = query_header.split()[0]

        if query_header not in [r["header"] for r in cleaned]:
            logger.print("[ERROR] Query sequence removed during cleaning.")
            return None

        return cleaned

    except Exception as e:
        logger.print(f"[ERROR] Exception in clean_a3m: {e}")
        return None


def postprocess_conservation_report_to_schema(
    raw_report: Dict[str, Any],
) -> Dict[str, Any] | None:
    """
    Postprocess the raw EnzyWizard-Conservation report into the official JSON Schema field names.
    """

    if not isinstance(raw_report, dict):
        return None

    raw_conservation_scores = raw_report.get("conservation_scores")
    if not isinstance(raw_conservation_scores, list):
        return None

    sequence_conservation_scores: List[Dict[str, Any]] = []

    for raw_item in raw_conservation_scores:
        if not isinstance(raw_item, dict):
            return None

        sequence_conservation_scores.append(
            {
                "residue_index": raw_item.get("aa_id"),
                "residue_name": raw_item.get("aa_name"),
                "hmm_profile_raw_score": raw_item.get("hmm_emission_log_score"),
                "normalized_emission_probability": raw_item.get("emission_probability"),
                "normalized_shannon_information_content": raw_item.get("conservation_score"),
            }
        )

    return {
        "report_type": raw_report.get("output_type"),
        "sequence_conservation_scores": sequence_conservation_scores,
    }