"""
DisProt client for experimental disorder sub-state annotations.

MobiDB gives us generic disorder + binding modes + compositional proxies, but
the experimentally-curated *structural states* (molten globule, pre-molten
globule) and disorder→order transitions (MoRFs) live in DisProt's IDPO
annotations. These axes separate a molten globule from a random coil, which a
scalar disorder score does not.

DisProt is small (~3k entries) and these states are sparse, so the workflow is to
probe the counts first and build masks only if they are sufficient.

Verified entry schema (https://disprot.org/api/<DPxxxxx>):
    {"acc": "P03265", "sequence": "...",
     "regions": [{"start": 294, "end": 334, "term_id": "IDPO:0000002",
                  "term_name": "disorder", "term_namespace": "Structural state",
                  "region_id": "..."}]}
Regions are 1-indexed inclusive.
"""
import logging
import time

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

DISPROT_SEARCH = "https://disprot.org/api/search"

# Target structural states / transitions, matched case-insensitively on
# region["term_name"] (robust to IDPO id-format changes).
DEFAULT_TARGET_STATES = {
    "molten_globule": ["molten globule"],
    "pre_molten_globule": ["pre-molten globule", "premolten globule"],
    "disorder_to_order": ["disorder to order"],   # MoRF
    "order_to_disorder": ["order to disorder"],
    "disorder": ["disorder"],                      # reference / sanity
}


def _session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": "IDPInterp/0.1 (research)"})
    s.mount("https://", HTTPAdapter(
        max_retries=Retry(total=5, backoff_factor=1,
                          status_forcelist=[429, 500, 502, 503, 504])))
    return s


def fetch_released(release: str = "current") -> list[dict]:
    """Fetch all released DisProt entries (list of dicts with acc/sequence/regions)."""
    s = _session()
    params = {"release": release, "show_ambiguous": "false",
              "show_obsolete": "false", "format": "json"}
    for attempt in range(3):
        try:
            r = s.get(DISPROT_SEARCH, params=params, timeout=120)
            r.raise_for_status()
            payload = r.json()
            data = payload.get("data", payload if isinstance(payload, list) else [])
            logger.info(f"DisProt: fetched {len(data)} released entries.")
            return data
        except Exception as e:
            logger.warning(f"DisProt fetch attempt {attempt+1} failed: {e}")
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("Could not fetch DisProt released set.")


def region_table(entries: list[dict]) -> pd.DataFrame:
    """Flatten entries → one row per region with acc/start/end/term_name/namespace."""
    rows = []
    for e in entries:
        acc = e.get("acc")
        if not acc:
            continue
        for reg in e.get("regions", []):
            rows.append({
                "acc": acc,
                "start": int(reg["start"]),
                "end": int(reg["end"]),
                "term_name": (reg.get("term_name") or "").strip().lower(),
                "term_namespace": reg.get("term_namespace", ""),
            })
    return pd.DataFrame(rows)


def _matches(term_name: str, needles: list[str]) -> bool:
    return any(n in term_name for n in needles)


def probe_counts(
    regions: pd.DataFrame,
    eval_accs: set[str] | None = None,
    target_states: dict[str, list[str]] | None = None,
) -> dict:
    """
    Report, per target state, how many proteins/regions exist DB-wide and (if
    eval_accs given) within our evaluation corpus. This is the viability gate:
    a state with too few regions in-corpus (<~50) is not worth an F1 axis.
    """
    target_states = target_states or DEFAULT_TARGET_STATES
    out = {}
    for state, needles in target_states.items():
        m = regions[regions["term_name"].apply(lambda t: _matches(t, needles))]
        rec = {
            "db_proteins": int(m["acc"].nunique()),
            "db_regions": int(len(m)),
            "db_residues": int((m["end"] - m["start"] + 1).clip(lower=0).sum()),
        }
        if eval_accs is not None:
            mi = m[m["acc"].isin(eval_accs)]
            rec.update({
                "corpus_proteins": int(mi["acc"].nunique()),
                "corpus_regions": int(len(mi)),
                "viable": bool(len(mi) >= 50),
            })
        out[state] = rec
    return out
