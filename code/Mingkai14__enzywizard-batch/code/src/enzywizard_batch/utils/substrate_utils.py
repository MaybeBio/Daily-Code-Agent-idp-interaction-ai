from __future__ import annotations

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from ..utils.logging_utils import Logger
from ..utils.common_utils import get_clean_filename

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


import random
import time
import math
import numpy as np
from typing import Any, Dict, List, Optional, Set, Tuple
from urllib.parse import quote
from ..resources.substrate_resources import *

import requests

# Validation
def is_valid_smiles(smiles: str) -> bool:
    try:
        if not isinstance(smiles, str):
            return False
        smiles = smiles.strip()
        if not smiles:
            return False
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False

def is_valid_mol_2d(mol: Chem.Mol, logger: Logger) -> bool:
    if mol is None:
        logger.print("[ERROR] Input Mol(2D) is None.")
        return False

    try:
        if not isinstance(mol, Chem.Mol):
            logger.print("[ERROR] Input object is not an RDKit Mol(2D).")
            return False

        if mol.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input Mol(2D) contains no atoms.")
            return False

        Chem.SanitizeMol(mol)
        return True
    except Exception as e:
        logger.print(f"[ERROR] Input Mol(2D) is invalid or failed sanitization. Reason: {e}")
        return False


def is_valid_mol_h(mol_h: Chem.Mol, logger: Logger) -> bool:
    if mol_h is None:
        logger.print("[ERROR] Input Mol(H) is None.")
        return False

    try:
        if not isinstance(mol_h, Chem.Mol):
            logger.print("[ERROR] Input object is not an RDKit Mol(H).")
            return False

        if mol_h.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input Mol(H) contains no atoms.")
            return False

        Chem.SanitizeMol(mol_h)

        has_h = any(atom.GetAtomicNum() == 1 for atom in mol_h.GetAtoms())
        if not has_h:
            logger.print("[WARNING] Input Mol(H) does not contain explicit hydrogen atoms.")


        return True
    except Exception as e:
        logger.print(f"[ERROR] Input Mol(H) is invalid or failed sanitization. Reason: {e}")
        return False


def is_valid_conf_3d(conf: Chem.Conformer, logger: Logger) -> bool:
    if conf is None:
        logger.print("[ERROR] Input conformer is None.")
        return False

    try:
        if not isinstance(conf, Chem.Conformer):
            logger.print("[ERROR] Input object is not an RDKit Conformer.")
            return False

        if not conf.Is3D():
            logger.print("[ERROR] Input conformer is not 3D.")
            return False

        if conf.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input conformer contains no atoms.")
            return False

        for atom_idx in range(conf.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            if any(math.isnan(v) or math.isinf(v) for v in [pos.x, pos.y, pos.z]):
                logger.print("[ERROR] Input conformer contains invalid 3D coordinates.")
                return False

        return True
    except Exception as e:
        logger.print(f"[ERROR] Input conformer(3D) is invalid. Reason: {e}")
        return False


def is_valid_mol_3d(mol_3d: Chem.Mol, logger: Logger) -> bool:
    if mol_3d is None:
        logger.print("[ERROR] Input Mol(3D) is None.")
        return False

    try:
        if not isinstance(mol_3d, Chem.Mol):
            logger.print("[ERROR] Input object is not an RDKit Mol(3D).")
            return False

        if mol_3d.GetNumAtoms() <= 0:
            logger.print("[ERROR] Input Mol(3D) contains no atoms.")
            return False

        Chem.SanitizeMol(mol_3d)

        if mol_3d.GetNumConformers() <= 0:
            logger.print("[ERROR] Input Mol(3D) contains no conformer.")
            return False

        conf = mol_3d.GetConformer()
        if not is_valid_conf_3d(conf, logger):
            return False

        if conf.GetNumAtoms() != mol_3d.GetNumAtoms():
            logger.print("[ERROR] Atom count mismatch between Mol(3D) and conformer.")
            return False

        return True
    except Exception as e:
        logger.print(f"[ERROR] Input Mol(3D) is invalid. Reason: {e}")
        return False

# fetch substrate
def clean_compound_name(raw_name: str | None) -> str:
    try:
        if raw_name is None:
            return ""

        name = str(raw_name).strip()

        if len(name) >= 2 and name[0] == name[-1] and name[0] in {"'", '"'}:
            name = name[1:-1].strip()

        name = re.sub(r"\s+", " ", name)
        name = TRAILING_PUNCT_RE.sub("", name).strip()
        name = name.replace("_", " ")
        name = re.sub(r"\s+", " ", name).strip()
        return name
    except Exception:
        return ""


def casefold_compound_name(name: str) -> str:
    try:
        return clean_compound_name(name).casefold()
    except Exception:
        return ""


def expand_synonyms_with_brackets(name: str) -> List[str]:
    try:
        cleaned_name = clean_compound_name(name)
        if not cleaned_name:
            return []

        result_set: Set[str] = {cleaned_name}

        name_without_brackets = re.sub(r"[\(\[\{][^)\]\}]+[\)\]\}]", "", cleaned_name)
        name_without_brackets = clean_compound_name(name_without_brackets)
        if name_without_brackets:
            result_set.add(name_without_brackets)

        for match in BRACKET_RE.finditer(cleaned_name):
            inner = clean_compound_name(match.group(1))
            if inner and len(inner) > 1:
                result_set.add(inner)

        result_list = [x for x in result_set if len(x) > 1]
        result_list.sort(key=lambda x: (-len(x), x))
        return result_list
    except Exception:
        return []


def light_normalize_for_match(name: str) -> str:
    try:
        normalized_name = clean_compound_name(name).casefold()
        if not normalized_name:
            return ""

        normalized_name = SEP_RE.sub(" ", normalized_name)
        normalized_name = MID_PUNCT_RE.sub(" ", normalized_name)
        normalized_name = re.sub(r"\s+", " ", normalized_name).strip()
        return normalized_name
    except Exception:
        return ""


def strip_nonstructural_tail(name: str) -> str:
    try:
        current_name = name.strip()
        if not current_name:
            return ""

        for pattern in PREFIX_PATTERNS:
            stripped_name = re.sub(pattern, "", current_name).strip()
            if stripped_name != current_name:
                current_name = stripped_name

        changed = True
        while changed:
            changed = False
            for phrase in SUFFIX_PHRASES:
                pattern = re.compile(rf"\b{phrase}\b\s*$", re.IGNORECASE)
                stripped_name = pattern.sub("", current_name).strip()
                if stripped_name != current_name and stripped_name:
                    current_name = stripped_name
                    changed = True
                    break

        current_name = re.sub(r"\s+", " ", current_name).strip()
        return current_name
    except Exception:
        return ""


def build_normalized_key_set(name: str) -> Set[str]:
    try:
        key_set: Set[str] = set()

        for form in expand_synonyms_with_brackets(name):
            normalized_form = light_normalize_for_match(form)
            if normalized_form:
                key_set.add(normalized_form)

            stripped_form = strip_nonstructural_tail(normalized_form)
            stripped_form = light_normalize_for_match(stripped_form)
            if stripped_form:
                key_set.add(stripped_form)

        return key_set
    except Exception:
        return set()


def build_retry_session(session: Optional[requests.Session] = None) -> requests.Session:
    try:
        used_session = session if session is not None else requests.Session()

        retry = Retry(
            total=5,
            connect=5,
            read=5,
            status=5,
            backoff_factor=0.5,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET", "POST"]),
            raise_on_status=False,
            respect_retry_after_header=True,
        )

        adapter = HTTPAdapter(max_retries=retry, pool_connections=20, pool_maxsize=20)
        used_session.mount("http://", adapter)
        used_session.mount("https://", adapter)
        return used_session
    except Exception:
        return requests.Session()


def request_json(
    session: requests.Session,
    method: str,
    url: str,
    logger: Logger,
    timeout: int = 15,
    json_body: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    log_http_error: bool = True,
) -> Optional[Dict[str, Any]]:
    try:
        time.sleep(0.2)
        response = session.request(
            method=method,
            url=url,
            timeout=timeout,
            json=json_body,
            headers=headers,
        )

        if response.status_code == 429:
            retry_after = response.headers.get("Retry-After")
            try:
                wait_time = float(retry_after) if retry_after else 1.0
            except Exception:
                wait_time = 1.0

            wait_time = min(30.0, max(0.5, wait_time)) + random.uniform(0.0, 0.5)
            time.sleep(wait_time)

            response = session.request(
                method=method,
                url=url,
                timeout=timeout,
                json=json_body,
                headers=headers,
            )

        if response.status_code < 200 or response.status_code >= 300:
            if log_http_error:
                logger.print(f"[WARNING] Request failed: {method} {url}, status_code={response.status_code}.")
            return None

        try:
            data = response.json()
        except Exception as e:
            logger.print(f"[WARNING] Failed to parse JSON response: {method} {url}. Reason: {e}")
            return None

        if isinstance(data, dict):
            return data

        logger.print(f"[WARNING] Unexpected JSON response type: {method} {url}, type={type(data).__name__}.")
        return None

    except Exception as e:
        logger.print(f"[WARNING] Request failed: {method} {url}. Reason: {e}")
        return None


def pubchem_name_to_cid(compound_name: str,session: requests.Session,logger: Logger,timeout: int = 15) -> Optional[int]:
    try:
        quoted_name = quote(compound_name, safe="")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quoted_name}/cids/JSON"
        data = request_json(session=session, method="GET", url=url, logger=logger, timeout=timeout)
        if data is None:
            return None

        cids = data.get("IdentifierList", {}).get("CID")
        if isinstance(cids, list) and len(cids) > 0:
            return int(cids[0])

        return None
    except Exception:
        return None


def pubchem_cid_to_smiles(cid: int,session: requests.Session,logger: Logger,timeout: int = 15) -> Optional[str]:
    try:
        props = "CanonicalSMILES,IsomericSMILES,ConnectivitySMILES,SMILES"
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{props}/JSON"
        data = request_json(session=session, method="GET", url=url, logger=logger, timeout=timeout)
        if data is None:
            return None

        rows = data.get("PropertyTable", {}).get("Properties")
        if not isinstance(rows, list) or len(rows) == 0:
            return None

        for row in rows:
            if not isinstance(row, dict):
                continue
            for key in ("IsomericSMILES", "CanonicalSMILES", "SMILES", "ConnectivitySMILES"):
                value = row.get(key)
                if isinstance(value, str) and value.strip():
                    value = value.strip()
                    if is_valid_smiles(value):
                        return value

        return None
    except Exception:
        return None


def pubchem_cid_to_synonyms(cid: int,session: requests.Session,logger: Logger,timeout: int = 15) -> List[str]:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        data = request_json(session=session, method="GET", url=url, logger=logger, timeout=timeout)
        if data is None:
            return []

        info_list = data.get("InformationList", {}).get("Information")
        if not isinstance(info_list, list) or len(info_list) == 0:
            return []

        synonym_list = info_list[0].get("Synonym")
        if not isinstance(synonym_list, list):
            return []

        output_list: List[str] = []
        seen_keys: Set[str] = set()

        for synonym in synonym_list:
            cleaned_synonym = clean_compound_name(synonym)
            if not cleaned_synonym:
                continue

            synonym_key = cleaned_synonym.casefold()
            if synonym_key in seen_keys:
                continue

            seen_keys.add(synonym_key)
            output_list.append(cleaned_synonym)

        return output_list
    except Exception:
        return []


def extract_smiles_from_chebi_payload(payload: Any) -> Optional[str]:
    try:
        if payload is None:
            return None

        if isinstance(payload, dict):
            default_structure = payload.get("default_structure")
            if isinstance(default_structure, dict):
                smiles = default_structure.get("smiles")
                if isinstance(smiles, str) and smiles.strip():
                    return smiles.strip()

            for key in (
                "smiles",
                "SMILES",
                "smilesString",
                "smiles_string",
                "canonicalSmiles",
                "canonical_smiles",
                "canonicalSMILES",
                "CanonicalSMILES",
            ):
                value = payload.get(key)
                if isinstance(value, str) and value.strip():
                    return value.strip()

            for key in (
                "structure",
                "chemicalStructure",
                "chemical_structure",
                "mol",
                "molecule",
                "defaultStructure",
                "default_structure",
            ):
                nested_value = payload.get(key)
                nested_smiles = extract_smiles_from_chebi_payload(nested_value)
                if nested_smiles:
                    return nested_smiles

            for key in ("structures", "chemicalStructures", "chemical_structures"):
                nested_list = payload.get(key)
                if isinstance(nested_list, list):
                    for item in nested_list:
                        nested_smiles = extract_smiles_from_chebi_payload(item)
                        if nested_smiles:
                            return nested_smiles

            for value in payload.values():
                nested_smiles = extract_smiles_from_chebi_payload(value)
                if nested_smiles:
                    return nested_smiles

        elif isinstance(payload, list):
            for item in payload:
                nested_smiles = extract_smiles_from_chebi_payload(item)
                if nested_smiles:
                    return nested_smiles

        return None
    except Exception:
        return None


def extract_id_name_from_chebi_hit(hit: Dict[str, Any]) -> Tuple[Optional[str], Optional[str]]:
    try:
        compound_id: Optional[str] = None
        compound_name: Optional[str] = None

        for key in ("chebiId", "chebi_id", "id", "compoundId", "compound_id"):
            value = hit.get(key)
            if isinstance(value, (str, int)):
                compound_id = str(value)
                break

        for key in ("chebiAsciiName", "name", "compoundName", "label", "title"):
            value = hit.get(key)
            if isinstance(value, str) and value.strip():
                compound_name = value.strip()
                break

        if compound_id and compound_id.isdigit():
            compound_id = f"CHEBI:{compound_id}"

        return compound_id, compound_name
    except Exception:
        return None, None


def chebi_fetch_compound(chebi_id: str,session: requests.Session,logger: Logger,timeout: int = 15) -> Optional[Dict[str, Any]]:
    try:
        chebi_number = chebi_id.replace("CHEBI:", "")
        url = f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{chebi_number}/"
        return request_json(session=session, method="GET", url=url, logger=logger, timeout=timeout)
    except Exception:
        return None


def chebi_search_exact(compound_name: str,session: requests.Session,logger: Logger,timeout: int = 15) -> List[Dict[str, Any]]:
    try:
        candidates: List[Dict[str, Any]] = []

        endpoints = [
            ("POST", "https://www.ebi.ac.uk/chebi/backend/api/public/search/"),
            ("POST", "https://www.ebi.ac.uk/chebi/backend/api/public/advanced_search/"),
        ]

        payloads = [
            {"term": compound_name},
            {"text": compound_name},
            {"query": compound_name},
            {"search": compound_name},
        ]

        headers = {"Content-Type": "application/json", "Accept": "application/json"}

        for method, url in endpoints:
            for body in payloads:
                data = request_json(
                    session=session,
                    method=method,
                    url=url,
                    logger=logger,
                    timeout=timeout,
                    json_body=body,
                    headers=headers,
                    log_http_error=False,
                )
                if data is None:
                    continue

                for key in ("results", "result", "hits", "compounds", "entities", "data"):
                    value = data.get(key)
                    if isinstance(value, list):
                        candidates.extend([x for x in value if isinstance(x, dict)])

                if "page" in data and isinstance(data["page"], dict):
                    content = data["page"].get("content")
                    if isinstance(content, list):
                        candidates.extend([x for x in content if isinstance(x, dict)])

                if len(candidates) > 0:
                    return candidates

        ols_url = f"https://www.ebi.ac.uk/ols4/api/search?q={quote(compound_name, safe='')}&ontology=chebi&exact=true"
        data = request_json(session=session, method="GET", url=ols_url, logger=logger, timeout=timeout)
        if data is None:
            return []

        docs = data.get("response", {}).get("docs")
        if not isinstance(docs, list):
            return []

        for doc in docs:
            if not isinstance(doc, dict):
                continue

            iri = doc.get("iri")
            label = doc.get("label")
            chebi_id: Optional[str] = None

            if isinstance(iri, str):
                match = re.search(r"CHEBI[_:](\d+)", iri)
                if match:
                    chebi_id = f"CHEBI:{match.group(1)}"

            if chebi_id:
                candidates.append({"chebiId": chebi_id, "name": label})

        return candidates
    except Exception:
        return []


def choose_best_chebi_smiles(query_name: str,hit_list: List[Dict[str, Any]],session: requests.Session,logger: Logger,timeout: int = 15) -> Optional[str]:
    try:
        query_casefold = casefold_compound_name(query_name)

        level1_matches: List[Tuple[str, str]] = []
        level2_candidates: List[Tuple[str, str]] = []

        for hit in hit_list:
            compound_id, compound_name = extract_id_name_from_chebi_hit(hit)
            if not compound_id or not compound_name:
                continue

            if casefold_compound_name(compound_name) == query_casefold:
                level1_matches.append((compound_id, compound_name))
            else:
                level2_candidates.append((compound_id, compound_name))

        matched_list = level1_matches

        if len(matched_list) == 0:
            query_key_set = build_normalized_key_set(query_name)
            if len(query_key_set) > 0:
                for compound_id, compound_name in level2_candidates:
                    name_key_set = build_normalized_key_set(compound_name)
                    if len(name_key_set) > 0 and len(query_key_set & name_key_set) > 0:
                        matched_list.append((compound_id, compound_name))

        if len(matched_list) == 0:
            return None

        for compound_id, _ in matched_list:
            detail = chebi_fetch_compound(
                chebi_id=compound_id,
                session=session,
                logger=logger,
                timeout=timeout,
            )
            if detail is None:
                continue

            smiles = extract_smiles_from_chebi_payload(detail)
            if isinstance(smiles, str) and smiles.strip():
                smiles = smiles.strip()
                if is_valid_smiles(smiles):
                    return smiles

        return None
    except Exception:
        return None

# process substrate
def get_mol_from_smiles(smiles: str, logger: Logger) -> Chem.Mol | None:
    if not isinstance(smiles, str):
        logger.print("[ERROR] Input SMILES must be a string.")
        return None

    smiles = smiles.strip()
    if not smiles:
        logger.print("[ERROR] Input SMILES is empty.")
        return None

    if not is_valid_smiles(smiles):
        logger.print(f"[ERROR] Invalid SMILES: {smiles}")
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.print(f"[ERROR] Failed to convert SMILES to Mol(2D): {smiles}")
            return None

        if not is_valid_mol_2d(mol, logger):
            return None

        return mol
    except Exception as e:
        logger.print(f"[ERROR] Failed to convert SMILES to Mol(2D): {smiles}. Reason: {e}")
        return None

def get_fingerprint_from_mol_2d(mol_2d: Chem.Mol,logger: Logger,radius: int = 2,n_bits: int = 512) -> List[int] | None:
    if not is_valid_mol_2d(mol_2d, logger):
        return None

    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol_2d,
            radius=int(radius),
            nBits=int(n_bits),
        )
        arr = np.zeros((int(n_bits),), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return [int(x) for x in arr]
    except Exception as e:
        logger.print(f"[ERROR] Failed to generate fingerprint from Mol(2D). Reason: {e}")
        return None

def get_2d_descriptor_dict_from_mol_2d(mol_2d: Chem.Mol,logger: Logger) -> Dict[str, Any] | None:
    if not is_valid_mol_2d(mol_2d, logger):
        return None

    try:
        descriptor_dict: Dict[str, Any] = {}
        descriptor_dict["substrate_num_atoms"] = int(mol_2d.GetNumAtoms())
        descriptor_dict["substrate_num_heavy_atoms"] = int(mol_2d.GetNumHeavyAtoms())
        descriptor_dict["substrate_num_hetero_atoms"] = int(
            sum(1 for atom in mol_2d.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
        )
        descriptor_dict["substrate_num_rings"] = int(rdMolDescriptors.CalcNumRings(mol_2d))
        descriptor_dict["substrate_num_aromatic_rings"] = int(
            rdMolDescriptors.CalcNumAromaticRings(mol_2d)
        )
        descriptor_dict["substrate_molecular_weight"] = float(Descriptors.MolWt(mol_2d))
        descriptor_dict["substrate_mol_logp"] = float(Crippen.MolLogP(mol_2d))
        descriptor_dict["substrate_molar_refractivity"] = float(Crippen.MolMR(mol_2d))
        descriptor_dict["substrate_tpsa"] = float(rdMolDescriptors.CalcTPSA(mol_2d))
        descriptor_dict["substrate_hbond_donor_count"] = int(
            rdMolDescriptors.CalcNumHBD(mol_2d)
        )
        descriptor_dict["substrate_hbond_acceptor_count"] = int(
            rdMolDescriptors.CalcNumHBA(mol_2d)
        )
        descriptor_dict["substrate_rotatable_bond_count"] = int(
            rdMolDescriptors.CalcNumRotatableBonds(mol_2d)
        )
        descriptor_dict["substrate_formal_charge"] = int(Chem.GetFormalCharge(mol_2d))
        return descriptor_dict
    except Exception as e:
        logger.print(f"[ERROR] Failed to generate 2D descriptors from Mol(2D). Reason: {e}")
        return None

def get_3d_descriptor_dict_from_mol_3d(mol_3d: Chem.Mol,logger: Logger) -> Dict[str, Any] | None:
    if not is_valid_mol_3d(mol_3d, logger):
        return None

    try:
        conf = mol_3d.GetConformer()
        coords = np.array(
            [
                [
                    conf.GetAtomPosition(atom_idx).x,
                    conf.GetAtomPosition(atom_idx).y,
                    conf.GetAtomPosition(atom_idx).z,
                ]
                for atom_idx in range(conf.GetNumAtoms())
            ],
            dtype=float,
        )

        pairwise_distances: List[float] = []
        for i in range(len(coords)):
            diff = coords[i + 1:] - coords[i]
            if len(diff) > 0:
                pairwise_distances.extend(np.linalg.norm(diff, axis=1).tolist())

        if pairwise_distances:
            distance_array = np.array(pairwise_distances, dtype=float)
            max_3d_diameter = float(np.max(distance_array))
            mean_pairwise_atom_distance = float(np.mean(distance_array))
            std_pairwise_atom_distance = float(np.std(distance_array))
        else:
            max_3d_diameter = 0.0
            mean_pairwise_atom_distance = 0.0
            std_pairwise_atom_distance = 0.0

        pmi1 = float(rdMolDescriptors.CalcPMI1(mol_3d))
        pmi3 = float(rdMolDescriptors.CalcPMI3(mol_3d))
        if pmi1 > 0:
            principal_moment_ratio = float(pmi3 / pmi1)
        else:
            principal_moment_ratio = 0.0

        return {
            "structure_max_3d_diameter": max_3d_diameter,
            "structure_mean_pairwise_atom_distance": mean_pairwise_atom_distance,
            "structure_std_pairwise_atom_distance": std_pairwise_atom_distance,
            "structure_asphericity": float(rdMolDescriptors.CalcAsphericity(mol_3d)),
            "structure_spherocity": float(rdMolDescriptors.CalcSpherocityIndex(mol_3d)),
            "structure_principal_moment_ratio": principal_moment_ratio,
            "structure_radius_of_gyration": float(
                rdMolDescriptors.CalcRadiusOfGyration(mol_3d)
            ),
        }
    except Exception as e:
        logger.print(f"[ERROR] Failed to generate 3D descriptors from Mol(3D). Reason: {e}")
        return None

def get_mol_h_from_mol_2d(mol_2d: Chem.Mol, logger: Logger) -> Chem.Mol | None:
    if not is_valid_mol_2d(mol_2d, logger):
        return None

    try:
        mol_h = Chem.AddHs(mol_2d, addCoords=False)
        if mol_h is None:
            logger.print("[ERROR] Failed to add hydrogens to Mol(2D).")
            return None

        if not is_valid_mol_h(mol_h, logger):
            return None

        return mol_h
    except Exception as e:
        logger.print(f"[ERROR] Failed to add hydrogens to Mol(2D). Reason: {e}")
        return None

def get_mol_3d_list_from_embedded_mol(mol_with_confs: Chem.Mol,logger: Logger,) -> List[Chem.Mol] | None:
    if mol_with_confs is None:
        logger.print("[ERROR] Embedded Mol is None.")
        return None

    try:
        if mol_with_confs.GetNumConformers() <= 0:
            logger.print("[ERROR] Embedded Mol contains no conformers.")
            return None

        mol_3d_list: List[Chem.Mol] = []

        for conf in mol_with_confs.GetConformers():
            if not is_valid_conf_3d(conf, logger):
                continue

            new_mol = Chem.Mol(mol_with_confs)
            conf_id_list = [c.GetId() for c in new_mol.GetConformers()]

            for conf_id in conf_id_list:
                if conf_id != conf.GetId():
                    new_mol.RemoveConformer(conf_id)

            if is_valid_mol_3d(new_mol, logger):
                mol_3d_list.append(new_mol)

        if len(mol_3d_list) == 0:
            logger.print("[WARNING] No valid 3D conformers were obtained from embedded Mol.")
            return []

        return mol_3d_list
    except Exception as e:
        logger.print(f"[ERROR] Failed to split embedded Mol into Mol(3D) list. Reason: {e}")
        return None


def get_mol_3d_list_from_mol_h(
    mol_h: Chem.Mol,
    logger: Logger,
    num_confs: int = 5,
    prune_rms: float = 0.5,
    random_seed: int = 202602,
) -> List[Chem.Mol] | None:
    if not is_valid_mol_h(mol_h, logger):
        return None

    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = int(random_seed)
        params.pruneRmsThresh = float(prune_rms)

        params.useRandomCoords = True

        mol_h = Chem.Mol(mol_h)
        conf_ids = AllChem.EmbedMultipleConfs(
            mol_h,
            numConfs=int(num_confs),
            params=params,
        )

        if conf_ids is None or len(conf_ids) == 0:
            logger.print("[WARNING] No conformers were generated by EmbedMultipleConfs.")
            return []

        logger.print(f"[INFO] Successfully generated {len(conf_ids)} conformers.")

        mol_3d_list = get_mol_3d_list_from_embedded_mol(mol_h, logger)
        if mol_3d_list is None:
            return None

        return mol_3d_list
    except Exception as e:
        logger.print(f"[ERROR] Failed to generate Mol(3D) list from Mol(H). Reason: {e}")
        return None

def get_mmff_energy_from_mol_3d(mol_3d: Chem.Mol, logger: Logger) -> float | None:
    if not is_valid_mol_3d(mol_3d, logger):
        return None

    try:
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol_3d, mmffVariant="MMFF94s")
        if mmff_props is None:
            logger.print("[ERROR] MMFF properties are unavailable for this Mol(3D).")
            return None

        conf_id = mol_3d.GetConformer().GetId()
        ff = AllChem.MMFFGetMoleculeForceField(mol_3d, mmff_props, confId=conf_id)

        if ff is None:
            logger.print(f"[ERROR] MMFF force field is unavailable for this Mol(3D): conf_id={conf_id}.")
            return None

        energy = ff.CalcEnergy()
        if energy is None or math.isnan(float(energy)) or math.isinf(float(energy)):
            logger.print(f"[ERROR] MMFF energy is invalid for this Mol(3D): conf_id={conf_id}.")
            return None

        return float(energy)
    except Exception as e:
        try:
            conf_id = mol_3d.GetConformer().GetId()
        except Exception:
            conf_id = "unknown"
        logger.print(f"[ERROR] Failed to calculate MMFF energy from Mol(3D): conf_id={conf_id}. Reason: {e}")
        return None


def get_uff_energy_from_mol_3d(mol_3d: Chem.Mol, logger: Logger) -> float | None:
    if not is_valid_mol_3d(mol_3d, logger):
        return None

    try:
        conf_id = mol_3d.GetConformer().GetId()
        ff = AllChem.UFFGetMoleculeForceField(mol_3d, confId=conf_id)

        if ff is None:
            logger.print(f"[ERROR] UFF force field is unavailable for this Mol(3D): conf_id={conf_id}.")
            return None

        energy = ff.CalcEnergy()
        if energy is None or math.isnan(float(energy)) or math.isinf(float(energy)):
            logger.print(f"[ERROR] UFF energy is invalid for this Mol(3D): conf_id={conf_id}.")
            return None

        logger.print(f"[INFO] UFF energy calculated, conf_id={conf_id}, energy={energy:.4f}")

        return float(energy)
    except Exception as e:
        try:
            conf_id = mol_3d.GetConformer().GetId()
        except Exception:
            conf_id = "unknown"
        logger.print(f"[ERROR] Failed to calculate UFF energy from Mol(3D): conf_id={conf_id}. Reason: {e}")
        return None


def get_energy_from_mol_3d(mol_3d: Chem.Mol, logger: Logger) -> float | None:
    if not is_valid_mol_3d(mol_3d, logger):
        return None

    energy = get_uff_energy_from_mol_3d(mol_3d, logger)
    if energy is not None:
        return energy

    logger.print("[ERROR] Failed to get UFF energy from Mol(3D).")
    return None

def get_minimized_mol_3d_list_from_mol_3d_list(mol_3d_list: List[Chem.Mol],logger: Logger) -> List[Chem.Mol] | None:

    if mol_3d_list is None or not isinstance(mol_3d_list, list):
        logger.print("[ERROR] Input Mol(3D) list is invalid.")
        return None

    if not mol_3d_list:
        logger.print("[ERROR] Input Mol(3D) list is empty.")
        return None

    minimized_list: List[Chem.Mol] = []

    try:
        for mol_3d in mol_3d_list:
            if not is_valid_mol_3d(mol_3d, logger):
                continue

            try:
                mol_copy = Chem.Mol(mol_3d)

                # UFF 最小化
                conf_id = mol_copy.GetConformer().GetId()
                result = AllChem.UFFOptimizeMolecule(mol_copy, confId=conf_id)

                # result=0 收敛，>0 未收敛（但通常也可以用）
                if result != 0:
                    logger.print(f"[WARNING] UFF optimization did not converge, conf_id={conf_id}")
                else:
                    logger.print(f"[INFO] UFF optimization success, conf_id={conf_id}")

                if not is_valid_mol_3d(mol_copy, logger):
                    continue

                minimized_list.append(mol_copy)

            except Exception as e:
                try:
                    conf_id = mol_3d.GetConformer().GetId()
                except Exception:
                    conf_id = "unknown"
                logger.print(f"[WARNING] Failed to minimize a Mol(3D): conf_id={conf_id}. Reason: {e}")
                continue

        if len(minimized_list) == 0:
            logger.print("[WARNING] No valid minimized Mol(3D) was obtained.")
            return []

        return minimized_list

    except Exception as e:
        logger.print(f"[ERROR] Failed to minimize Mol(3D) list. Reason: {e}")
        return None

def get_substrate_report_suffix_from_feature_list(substrate_feature_list: List[Dict[str, Any]],logger: Logger,) -> str | None:
    if not isinstance(substrate_feature_list, list):
        logger.print("[ERROR] substrate_feature_list must be a list.")
        return None

    try:
        substrate_name_list: List[str] = []

        for item in substrate_feature_list:
            if not isinstance(item, dict):
                logger.print("[ERROR] Invalid substrate entry in substrate_feature_list.")
                return None

            substrate_name = item.get("substrate_name")
            if not isinstance(substrate_name, str) or not substrate_name.strip():
                logger.print("[ERROR] Invalid substrate_name in substrate_feature_list.")
                return None

            substrate_name_list.append(substrate_name.strip())

        if len(substrate_name_list) == 0:
            logger.print("[ERROR] No substrate_name was found in substrate_feature_list.")
            return None

        suffix = "_".join(substrate_name_list)
        suffix = get_clean_filename(suffix)

        return suffix

    except Exception as e:
        logger.print(f"[ERROR] Failed to generate substrate report suffix. Reason: {e}")
        return None


def build_docked_mol_from_atom_info(
    original_mol_3d: Chem.Mol,
    docked_atom_info_list: List[Dict[str, Any]],
    logger: Logger,
) -> Chem.Mol | None:

    if original_mol_3d is None or original_mol_3d.GetNumConformers() <= 0:
        logger.print("[ERROR] Invalid original Mol(3D).")
        return None

    if not isinstance(docked_atom_info_list, list) or len(docked_atom_info_list) == 0:
        logger.print("[ERROR] Invalid docked_atom_info_list.")
        return None

    try:
        atom_num = original_mol_3d.GetNumAtoms()

        used_original_atom_index_set = set()
        kept_original_atom_index_list: List[int] = []

        for item in docked_atom_info_list:
            original_atom_index = int(item.get("original_atom_index", 0))

            if original_atom_index <= 0 or original_atom_index > atom_num:
                logger.print("[ERROR] Invalid original atom index in docked_atom_info_list.")
                return None

            if original_atom_index in used_original_atom_index_set:
                logger.print("[ERROR] Duplicate original atom index.")
                return None

            used_original_atom_index_set.add(original_atom_index)
            kept_original_atom_index_list.append(original_atom_index)

        kept_original_atom_index_list.sort()

        old_to_new_index_dict = {
            old_idx: new_idx
            for new_idx, old_idx in enumerate(kept_original_atom_index_list)
        }

        rw_mol = Chem.RWMol()
        new_conf = Chem.Conformer(len(kept_original_atom_index_list))

        # ===== build atoms =====
        for old_index in kept_original_atom_index_list:
            old_atom = original_mol_3d.GetAtomWithIdx(old_index - 1)

            new_atom = Chem.Atom(old_atom.GetAtomicNum())
            new_atom.SetFormalCharge(old_atom.GetFormalCharge())
            new_atom.SetIsAromatic(old_atom.GetIsAromatic())
            new_atom.SetChiralTag(old_atom.GetChiralTag())
            new_atom.SetNoImplicit(old_atom.GetNoImplicit())
            new_atom.SetNumExplicitHs(old_atom.GetNumExplicitHs())
            new_atom.SetNumRadicalElectrons(old_atom.GetNumRadicalElectrons())

            rw_mol.AddAtom(new_atom)

        kept_old_index_set = set(kept_original_atom_index_list)

        # ===== build bonds =====
        for bond in original_mol_3d.GetBonds():
            b = bond.GetBeginAtomIdx() + 1
            e = bond.GetEndAtomIdx() + 1

            if b in kept_old_index_set and e in kept_old_index_set:
                rw_mol.AddBond(
                    old_to_new_index_dict[b],
                    old_to_new_index_dict[e],
                    bond.GetBondType(),
                )

        # ===== set coordinates =====
        for item in docked_atom_info_list:
            old_idx = int(item["original_atom_index"])
            x, y, z = float(item["x"]), float(item["y"]), float(item["z"])

            new_idx = old_to_new_index_dict[old_idx]
            new_conf.SetAtomPosition(new_idx, (x, y, z))

        mol = rw_mol.GetMol()
        mol.RemoveAllConformers()
        mol.AddConformer(new_conf, assignId=True)

        Chem.SanitizeMol(mol)

        return mol

    except Exception as e:
        logger.print(f"[ERROR] Failed to build docked Mol. Reason: {e}")
        return None

def postprocess_substrate_report_to_schema(
    raw_report: Dict[str, Any],
    logger: Logger,
) -> Dict[str, Any] | None:
    """
    Postprocess the raw EnzyWizard-Substrate report into the official JSON Schema field names.
    """

    if not isinstance(raw_report, dict):
        logger.print("[ERROR] raw_report must be a dictionary.")
        return None

    def is_integer(value: Any) -> bool:
        return isinstance(value, int) and not isinstance(value, bool)

    def is_number(value: Any) -> bool:
        return (
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            and math.isfinite(float(value))
        )

    def is_fingerprint(value: Any) -> bool:
        return (
            isinstance(value, list)
            and all(is_integer(x) and x in (0, 1) for x in value)
        )

    try:
        substrates = raw_report.get("substrates")
        if not isinstance(substrates, list):
            logger.print("[ERROR] raw_report['substrates'] must be a list.")
            return None

        schema_substrate_list: List[Dict[str, Any]] = []

        for item in substrates:
            if not isinstance(item, dict):
                logger.print("[ERROR] Invalid substrate entry in raw report.")
                return None

            substrate_name = item.get("substrate_name")
            smiles = item.get("smiles")
            if not isinstance(substrate_name, str) or not substrate_name.strip():
                logger.print("[ERROR] Invalid substrate_name in raw report.")
                return None

            if not isinstance(smiles, str) or not smiles.strip():
                logger.print(f"[ERROR] Invalid SMILES in raw report for substrate: {substrate_name}")
                return None

            fingerprint = item.get("fingerprint")
            if not is_fingerprint(fingerprint):
                logger.print(f"[ERROR] Invalid fingerprint in raw report for substrate: {substrate_name}")
                return None

            integer_field_list = [
                "num_atoms",
                "heavy_atom_count",
                "hbond_donor_count",
                "hbond_acceptor_count",
                "rotatable_bond_count",
            ]
            for field_name in integer_field_list:
                if not is_integer(item.get(field_name)):
                    logger.print(f"[ERROR] Invalid substrate integer field '{field_name}' for substrate: {substrate_name}")
                    return None

            number_field_list = [
                "mol_weight",
                "logp",
                "tpsa",
                "molar_refractivity",
            ]
            for field_name in number_field_list:
                if not is_number(item.get(field_name)):
                    logger.print(f"[ERROR] Invalid substrate numeric field '{field_name}' for substrate: {substrate_name}")
                    return None

            raw_structure_list = item.get("structures", [])
            if not isinstance(raw_structure_list, list):
                logger.print("[ERROR] Invalid structures field in raw substrate report.")
                return None

            schema_structure_list: List[Dict[str, Any]] = []

            for structure_item in raw_structure_list:
                if not isinstance(structure_item, dict):
                    logger.print("[ERROR] Invalid structure entry in raw substrate report.")
                    return None

                structure_name = structure_item.get("structure_name")
                if not isinstance(structure_name, str) or not structure_name.strip():
                    logger.print(f"[ERROR] Invalid structure_name in raw report for substrate: {substrate_name}")
                    return None

                structure_number_field_list = [
                    "structure_energy",
                    "structure_max_3d_diameter",
                    "structure_mean_pairwise_atom_distance",
                    "structure_std_pairwise_atom_distance",
                    "structure_asphericity",
                    "structure_spherocity",
                    "structure_principal_moment_ratio",
                    "structure_radius_of_gyration",
                ]
                for field_name in structure_number_field_list:
                    if not is_number(structure_item.get(field_name)):
                        logger.print(f"[ERROR] Invalid substrate structure numeric field '{field_name}' for structure: {structure_name}")
                        return None

                schema_structure_list.append(
                    {
                        "substrate_structure_name": structure_name,
                        "substrate_structure_energy": structure_item["structure_energy"],
                        "substrate_structure_max_3d_diameter": structure_item["structure_max_3d_diameter"],
                        "substrate_structure_mean_pairwise_atom_distance": structure_item["structure_mean_pairwise_atom_distance"],
                        "substrate_structure_std_pairwise_atom_distance": structure_item["structure_std_pairwise_atom_distance"],
                        "substrate_structure_asphericity": structure_item["structure_asphericity"],
                        "substrate_structure_spherocity": structure_item["structure_spherocity"],
                        "substrate_structure_principal_moment_ratio": structure_item["structure_principal_moment_ratio"],
                        "substrate_structure_radius_of_gyration": structure_item["structure_radius_of_gyration"],
                    }
                )

            schema_substrate_list.append(
                {
                    "substrate_name": substrate_name,
                    "substrate_smiles": smiles,
                    "substrate_fingerprint_encoding": fingerprint,
                    "substrate_atom_count": item["num_atoms"],
                    "substrate_molecular_weight": item["mol_weight"],
                    "substrate_logp": item["logp"],
                    "substrate_tpsa": item["tpsa"],
                    "substrate_heavy_atom_count": item["heavy_atom_count"],
                    "substrate_hbond_donor_count": item["hbond_donor_count"],
                    "substrate_hbond_acceptor_count": item["hbond_acceptor_count"],
                    "substrate_rotatable_bond_count": item["rotatable_bond_count"],
                    "substrate_molar_refractivity": item["molar_refractivity"],
                    "substrate_possible_structures": schema_structure_list,
                }
            )

        return {
            "report_type": raw_report.get("output_type", "enzywizard_substrate"),
            "substrates": schema_substrate_list,
        }

    except Exception as e:
        logger.print(f"[ERROR] Failed to postprocess substrate report to schema. Reason: {e}")
        return None
