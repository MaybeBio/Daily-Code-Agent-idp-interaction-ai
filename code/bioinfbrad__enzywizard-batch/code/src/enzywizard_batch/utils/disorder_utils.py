from __future__ import annotations
from ..utils.logging_utils import Logger
from typing import Any, Dict, List

def moving_average(values: List[float], window_size: int, logger: Logger) -> List[float] | None:
    if window_size <= 0:
        logger.print("[ERROR] window_size must be positive")
        return None

    n = len(values)
    if n == 0:
        logger.print("[ERROR] Empty input for moving_average")
        return None

    half = window_size // 2
    result: List[float] = []

    for i in range(n):
        left = max(0, i - half)
        right = min(n, i + half + 1)
        sub = values[left:right]
        result.append(sum(sub) / len(sub))

    return result



def postprocess_disorder_report_to_schema(
    raw_report: Dict[str, Any],
    logger: Logger,
) -> Dict[str, Any] | None:
    """
    Postprocess the raw EnzyWizard-Disorder report into the official JSON Schema field names.
    """

    try:
        raw_statistics = raw_report.get("disorder_region_statistics", {})
        raw_regions = raw_report.get("disorder_regions", [])

        disordered_regions: List[Dict[str, Any]] = []

        for raw_region in raw_regions:
            residues: List[Dict[str, Any]] = []

            for raw_residue in raw_region.get("residues", []):
                residues.append({
                    "residue_index": raw_residue.get("aa_id"),
                    "residue_name": raw_residue.get("aa_name"),
                })

            disordered_regions.append({
                "disordered_region_length": raw_region.get("length"),
                "residues": residues,
            })

        schema_report: Dict[str, Any] = {
            "report_type": raw_report.get("output_type"),
            "disordered_region_statistics": {
                "disordered_region_count": raw_statistics.get("region_num"),
                "max_disordered_region_length": raw_statistics.get("max_region_length"),
                "total_disordered_region_length": raw_statistics.get("total_region_length"),
            },
            "disordered_regions": disordered_regions,
        }

        return schema_report

    except Exception as e:
        logger.print(f"[ERROR] Failed to postprocess disorder report to schema: {e}")
        return None