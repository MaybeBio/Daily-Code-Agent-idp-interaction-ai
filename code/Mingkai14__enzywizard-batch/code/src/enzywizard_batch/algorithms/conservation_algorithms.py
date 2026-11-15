from __future__ import annotations
from ..utils.logging_utils import Logger
from typing import List, Dict, Any
from pathlib import Path
import math
from Bio.Data.IUPACData import protein_letters_1to3
from ..utils.sequence_utils import normalize_aa_name_to_one_letter, postprocess_conservation_report_to_schema


def get_emission_probabilities_from_hmm(hmm_file: str | Path, logger: Logger) -> List[Dict[str, Any]] | None:
    p = Path(hmm_file)

    try:
        if not p.exists():
            logger.print(f"[ERROR] HMM file not found: {str(p)}")
            return None

        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()

        alphabet: List[str] | None = None
        emission_list: List[Dict[str, Any]] = []
        in_hmm_block = False

        for line in lines:
            stripped = line.strip()

            valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
            if stripped.startswith("HMM "):
                toks = stripped.split()
                if toks[0] != "HMM":
                    continue

                alphabet = [x.upper() for x in toks[1:]]
                if len(alphabet) == 0:
                    logger.print(f"[ERROR] Empty HMM alphabet line in {str(p)}")
                    return None

                if any(aa not in valid_aa for aa in alphabet):
                    logger.print(f"[ERROR] Invalid HMM alphabet line in {str(p)}")
                    return None

                in_hmm_block = True
                continue

            if not in_hmm_block:
                continue

            if alphabet is None:
                logger.print(f"[ERROR] Missing HMM alphabet in {str(p)}")
                return None

            toks = stripped.split()
            if len(toks) == 0:
                continue

            if not toks[0].isdigit():
                continue

            pos = int(toks[0])

            if len(toks) < 1 + len(alphabet):
                logger.print(f"[ERROR] Invalid HMM match emission row at position {pos} in {str(p)}")
                return None

            raw_score_dict: Dict[str, float | None] = {}
            raw_prob_list: List[float] = []

            for aa, raw in zip(alphabet, toks[1:1 + len(alphabet)]):
                aa = str(aa).upper()

                if raw == "*":
                    raw_score_dict[aa] = 0.0
                    raw_prob_list.append(0.0)
                else:
                    try:
                        raw_value = float(raw)
                    except Exception as e:
                        logger.print(f"[ERROR] Invalid HMM emission value '{raw}' at position {pos} in {str(p)}: {e}")
                        return None

                    raw_score_dict[aa] = raw_value
                    raw_prob_list.append(math.exp(-raw_value))

            total = sum(raw_prob_list)
            if total <= 0.0:
                logger.print(f"[ERROR] Invalid HMM emission probability sum at position {pos} in {str(p)}")
                return None

            norm_prob_dict: Dict[str, float] = {}
            for aa, prob in zip(alphabet, raw_prob_list):
                norm_prob_dict[str(aa).upper()] = float(prob / total)

            emission_list.append({
                "position": pos,
                "raw_scores": raw_score_dict,
                "probabilities": norm_prob_dict
            })

        if len(emission_list) == 0:
            logger.print(f"[ERROR] No match emission rows found in HMM file: {str(p)}")
            return None

        return emission_list

    except Exception as e:
        logger.print(f"[ERROR] Exception in get_emission_probabilities_from_hmm: {e}")
        return None

def compute_shannon_entropy_from_emission_prob(prob_list: List[float], logger: Logger) -> float | None:
    try:
        if not isinstance(prob_list, list):
            logger.print("[ERROR] prob_list is not a list.")
            return None

        if len(prob_list) == 0:
            logger.print("[ERROR] prob_list is empty.")
            return None

        total = 0.0
        cleaned_prob_list: List[float] = []

        for p in prob_list:
            if not isinstance(p, (int, float)):
                logger.print("[ERROR] Invalid probability value found.")
                return None

            p = float(p)
            if math.isnan(p) or p < 0.0:
                logger.print("[ERROR] Invalid probability value found.")
                return None

            cleaned_prob_list.append(p)
            total += p

        if total <= 0.0:
            logger.print("[ERROR] Sum of probabilities must be positive.")
            return None

        entropy = 0.0
        for p in cleaned_prob_list:
            p_norm = p / total
            if p_norm > 0.0:
                entropy -= p_norm * math.log(p_norm)

        return float(entropy)

    except Exception as e:
        logger.print(f"[ERROR] Exception in compute_shannon_entropy_from_emission_prob: {e}")
        return None

def compute_conservation_scores(hmm_file: str | Path, sequence_dict: Dict[str,str], logger: Logger) -> List[Dict[str, Any]] | None:
    p = Path(hmm_file)
    query_sequence = sequence_dict.get("sequence",None)

    try:
        if query_sequence is None:
            logger.print("[ERROR] query_sequence is empty.")
            return None

        query_sequence = query_sequence.strip().upper()
        if query_sequence == "":
            logger.print("[ERROR] query_sequence is empty.")
            return None

        valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
        for ch in query_sequence:
            if ch not in valid_chars:
                logger.print(f"[ERROR] Invalid character '{ch}' in query_sequence.")
                return None

        emission_list = get_emission_probabilities_from_hmm(p, logger)
        if emission_list is None:
            return None

        if len(emission_list) != len(query_sequence):
            logger.print(
                f"[ERROR] HMM length {len(emission_list)} does not match query_sequence length {len(query_sequence)}."
            )
            return None

        result_list: List[Dict[str, Any]] = []

        for i, query_aa in enumerate(query_sequence, start=1):
            emission_prob_dict = emission_list[i - 1]["probabilities"]
            raw_score_dict = emission_list[i - 1]["raw_scores"]

            if query_aa not in emission_prob_dict:
                logger.print(f"[ERROR] Amino acid '{query_aa}' not found in HMM alphabet at position {i}.")
                return None

            prob_list = list(emission_prob_dict.values())
            entropy = compute_shannon_entropy_from_emission_prob(prob_list, logger)
            if entropy is None:
                logger.print(f"[ERROR] Failed to compute Shannon entropy at position {i}.")
                return None

            aa_name = protein_letters_1to3.get(query_aa)
            if aa_name is None:
                logger.print(f"[ERROR] Failed to map amino acid '{query_aa}' to 3-letter name.")
                return None

            result_list.append({
                "aa_id": i,
                "aa_name": normalize_aa_name_to_one_letter(aa_name.upper()),
                "hmm_emission_log_score": raw_score_dict.get(query_aa),
                "emission_probability": float(emission_prob_dict[query_aa]),
                "conservation_score": float((math.log(20.0) - entropy) / math.log(20.0)),
            })

        return result_list

    except Exception as e:
        logger.print(f"[ERROR] Exception in compute_conservation_scores: {e}")
        return None

def generate_conservation_report(conservation_scores: List[Dict[str, Any]]) -> dict | None:
    raw_report = {
        "output_type": "enzywizard_conservation",
        "conservation_scores": conservation_scores,
    }

    return postprocess_conservation_report_to_schema(raw_report)
