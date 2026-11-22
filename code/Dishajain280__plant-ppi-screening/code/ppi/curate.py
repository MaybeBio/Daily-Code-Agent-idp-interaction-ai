"""Validate a candidate list before any predictions are scored.

The failure mode this guards against is a label set that looks complete but is
quietly wrong: duplicate pairs counted twice, "sources" that say TODO, or so many
easy decoys that the metrics measure nothing. Every check here corresponds to a
mistake that produces a confident, meaningless number.
"""

from __future__ import annotations

from dataclasses import dataclass

ALLOWED_LABELS = {"0", "1"}
ALLOWED_VERIFIED = {"yes", "no"}
ALLOWED_TYPES = {"positive_curated", "negative_curated", "decoy_localization", "calibration"}
CALIBRATION = "calibration"
PLACEHOLDER_SOURCES = {"", "todo", "tbd", "?", "unknown", "none", "n/a", "-"}


@dataclass
class Finding:
    level: str  # "error" | "warning"
    where: str
    message: str

    def __str__(self) -> str:
        mark = "ERROR  " if self.level == "error" else "warning"
        return f"{mark}  {self.where}: {self.message}"


def _key(row: dict) -> tuple[str, str]:
    """Unordered protein pair, so A-B and B-A are recognised as the same test."""
    a = (row.get("acc_a") or row.get("gene_a") or "").strip().upper()
    b = (row.get("acc_b") or row.get("gene_b") or "").strip().upper()
    return tuple(sorted((a, b)))  # type: ignore[return-value]


def validate(rows: list[dict], min_pairs: int = 15, max_decoy_share: float = 0.4) -> list[Finding]:
    findings: list[Finding] = []
    seen_ids: dict[str, int] = {}
    seen_pairs: dict[tuple[str, str], str] = {}

    for index, row in enumerate(rows, start=2):  # row 1 is the header
        pair_id = (row.get("pair_id") or "").strip()
        where = pair_id or f"row {index}"

        if not pair_id:
            findings.append(Finding("error", f"row {index}", "pair_id is empty"))
            continue

        if pair_id in seen_ids:
            findings.append(Finding("error", where,
                                    f"duplicate pair_id (also on row {seen_ids[pair_id]})"))
        seen_ids[pair_id] = index

        label = (row.get("label") or "").strip()
        if label not in ALLOWED_LABELS:
            findings.append(Finding("error", where, f"label must be 0 or 1, got {label!r}"))

        verified = (row.get("verified") or "").strip().lower()
        if verified not in ALLOWED_VERIFIED:
            findings.append(Finding("error", where, f"verified must be yes or no, got {verified!r}"))

        label_type = (row.get("label_type") or "").strip()
        if label_type not in ALLOWED_TYPES:
            findings.append(Finding("error", where,
                                    f"label_type must be one of {sorted(ALLOWED_TYPES)}"))
        elif label == "0" and label_type == "positive_curated":
            findings.append(Finding("error", where,
                                    "a curated positive is marked negative"))
        elif label == "1" and label_type.startswith("decoy"):
            findings.append(Finding("error", where, "a decoy is marked positive"))

        source = (row.get("source") or "").strip()
        if source.lower() in PLACEHOLDER_SOURCES:
            findings.append(Finding("warning", where, "no source recorded; this label is unusable"))
        elif len(source) < 20:
            findings.append(Finding("warning", where,
                                    "source is very short; record author, year, method or PMID"))

        key = _key(row)
        is_calibration = label_type == CALIBRATION
        if key[0] and key[1] and not is_calibration:
            if key in seen_pairs:
                findings.append(Finding("error", where,
                                        f"same protein pair already listed as {seen_pairs[key]}; "
                                        f"counting it twice biases AUC"))
            else:
                seen_pairs[key] = pair_id
        if (row.get("gene_a") or "").strip().upper() == (row.get("gene_b") or "").strip().upper() != "":
            findings.append(Finding("warning", where,
                                    "same gene on both sides; homodimers are legitimate, but "
                                    "confirm that is what you meant"))

    findings.extend(_composition_findings(rows, min_pairs, max_decoy_share))
    return findings


def _composition_findings(rows: list[dict], min_pairs: int,
                          max_decoy_share: float) -> list[Finding]:
    findings: list[Finding] = []
    # Calibration pairs exist to check the pipeline, so they never enter the metrics.
    rows = [r for r in rows if (r.get("label_type") or "").strip() != CALIBRATION]
    total = len(rows)
    if total == 0:
        return [Finding("error", "file", "no benchmark rows (only calibration pairs)")]
    if total < min_pairs:
        findings.append(Finding(
            "warning", "file",
            f"only {total} pairs; with fewer than {min_pairs} the confidence interval will be "
            f"too wide to distinguish a real signal from luck"))

    positives = sum(1 for r in rows if (r.get("label") or "").strip() == "1")
    decoys = sum(1 for r in rows if (r.get("label_type") or "").strip() == "decoy_localization")
    unverified = sum(1 for r in rows if (r.get("verified") or "").strip().lower() != "yes")

    if positives < 8:
        findings.append(Finding("warning", "file",
                                f"only {positives} positives; AUC on very few positives is unstable"))
    if not 0 < positives < total:
        findings.append(Finding("error", "file", "the file needs both positives and negatives"))

    decoy_share = decoys / total
    if decoy_share > max_decoy_share:
        findings.append(Finding(
            "warning", "file",
            f"{decoy_share:.0%} of pairs are decoys; easy negatives inflate every metric. "
            f"Report decoy-only scores separately from curated negatives"))
    if unverified:
        findings.append(Finding(
            "warning", "file",
            f"{unverified} of {total} pairs are unverified; scoring is blocked until they are "
            f"checked or excluded"))
    return findings


def summary(rows: list[dict]) -> dict[str, int]:
    return {
        "pairs": len(rows),
        "calibration": sum(1 for r in rows
                           if (r.get("label_type") or "").strip() == CALIBRATION),
        "positives": sum(1 for r in rows if (r.get("label") or "").strip() == "1"),
        "negatives": sum(1 for r in rows if (r.get("label") or "").strip() == "0"),
        "curated_negatives": sum(1 for r in rows
                                 if (r.get("label_type") or "").strip() == "negative_curated"),
        "decoys": sum(1 for r in rows
                      if (r.get("label_type") or "").strip() == "decoy_localization"),
        "verified": sum(1 for r in rows if (r.get("verified") or "").strip().lower() == "yes"),
    }
