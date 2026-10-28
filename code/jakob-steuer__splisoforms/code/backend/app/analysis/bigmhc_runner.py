"""
bigmhc_runner.py
=================
Subprocess wrapper for BigMHC EL (presentation) and IM (immunogenicity) models.

BigMHC — Karchin Lab, Johns Hopkins.
  Paper: Albert et al., Nature Machine Intelligence 2023
  Repo:  https://github.com/KarchinLab/bigmhc
  License: BigMHC Academic License

BigMHC is a script-based tool (not a PyPI library).  We invoke it via
subprocess using its predict.py CLI.  The repo is expected at:
  {PROJECT_ROOT}/tools/bigmhc/

  git clone --depth=1 https://github.com/KarchinLab/bigmhc.git tools/bigmhc

Two models (selected via -m flag):
  el  — BigMHC EL: MHC-I eluted ligand presentation score (0–1).
         Outperforms NetMHCpan-4.1 on EL benchmarks (AUROC 0.9733 vs 0.9496).
  im  — BigMHC IM: immunogenicity from transfer learning on T-cell assay data.
         PPVn 0.4375 vs HLAthena 0.2638.  This is the key new signal not in MHCflurry.

Usage:
    from app.analysis.bigmhc_runner import run_bigmhc, is_available

    if is_available():
        results = run_bigmhc(peptides, alleles)
        # results[(peptide, allele)] = BigMHCResult(el_score, im_score)
"""

import csv
import logging
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ─── Paths ────────────────────────────────────────────────────────────────────

def _bigmhc_src() -> str:
    """Return path to the BigMHC src/ directory."""
    env_dir = os.environ.get("BIGMHC_DIR")
    if env_dir:
        return os.path.join(env_dir, "src")
    # Fallback: walk up from this file to project root, then tools/bigmhc/src
    here       = os.path.dirname(os.path.abspath(__file__))
    proj_root  = os.path.dirname(os.path.dirname(os.path.dirname(here)))
    return os.path.join(proj_root, "tools", "bigmhc", "src")


def _bigmhc_predict_script() -> str:
    return os.path.join(_bigmhc_src(), "predict.py")


def _detect_device() -> str:
    """Return device flag for BigMHC CLI: GPU index (e.g. '0') or 'cpu'."""
    try:
        import torch
        if torch.cuda.is_available():
            return str(torch.cuda.current_device())
    except ImportError:
        pass
    return "cpu"


# ─── Availability check ────────────────────────────────────────────────────────

_bigmhc_available: Optional[bool] = None


def is_available() -> bool:
    """Return True if the BigMHC predict.py script is present."""
    global _bigmhc_available
    if _bigmhc_available is None:
        _bigmhc_available = os.path.isfile(_bigmhc_predict_script())
        if not _bigmhc_available:
            logger.debug(
                f"BigMHC predict.py not found at {_bigmhc_predict_script()}. "
                "Clone repo: git clone --depth=1 "
                "https://github.com/KarchinLab/bigmhc.git tools/bigmhc"
            )
    return _bigmhc_available


# ─── Result dataclass ──────────────────────────────────────────────────────────

@dataclass
class BigMHCResult:
    el_score: Optional[float]   # EL presentation score (0–1); None = not scored by the tool
    im_score: Optional[float]   # IM immunogenicity score (0–1); None = not scored by the tool


# ─── Internal runner ───────────────────────────────────────────────────────────

def _run_model(
    input_csv : str,
    model     : str,   # "el" or "im"
    output_csv: str,
) -> None:
    """
    Run BigMHC predict.py for one model.

    Input CSV columns (zero-indexed):
      col 0 = allele  (-a=0)
      col 1 = peptide (-p=1)
      1 header row    (-c=1)

    Output: same CSV with prediction column appended:
      BigMHC_EL  (when model="el")
      BigMHC_IM  (when model="im")
    """
    cmd = [
        sys.executable,              # same Python interpreter / venv
        _bigmhc_predict_script(),
        f"-i={input_csv}",
        f"-m={model}",
        f"-a=0",                     # allele in column 0
        f"-p=1",                     # peptide in column 1
        f"-c=1",                     # 1 header row
        f"-d={_detect_device()}",
        f"-o={output_csv}",
        f"-v=0",                     # suppress verbose output
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=_bigmhc_src(),           # predict.py imports sibling modules via cwd
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"BigMHC predict.py ({model}) failed (exit {result.returncode}):\n"
            f"  stdout: {result.stdout[-500:]}\n"
            f"  stderr: {result.stderr[-500:]}"
        )


# ─── Main runner ───────────────────────────────────────────────────────────────

def run_bigmhc(
    peptides: List[str],
    alleles : List[str],
) -> Dict[Tuple[str, str], BigMHCResult]:
    """
    Run BigMHC EL + IM for all (peptide × allele) combinations.

    Writes a temporary CSV, calls predict.py twice (EL and IM), reads
    both output CSVs, and returns a results dict.

    Args:
        peptides: Unique peptide sequences (8–11 aa, standard amino acids).
        alleles:  HLA alleles, e.g. ["HLA-A*02:01", "HLA-B*07:02"].

    Returns:
        Dict mapping (peptide, allele) → BigMHCResult.
        Keys use the original allele strings as passed in.
    """
    if not peptides or not alleles:
        return {}

    # Build Cartesian product: allele col 0, peptide col 1 (matches CLI defaults)
    rows = [(a, p) for p in peptides for a in alleles]

    with tempfile.TemporaryDirectory(prefix="bigmhc_") as tmpdir:
        input_csv  = os.path.join(tmpdir, "input.csv")
        el_out_csv = os.path.join(tmpdir, "el_out.csv")
        im_out_csv = os.path.join(tmpdir, "im_out.csv")

        # Write input
        with open(input_csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["mhc", "pep"])
            writer.writerows(rows)

        logger.info(
            f"  BigMHC: {len(rows)} peptide×allele pairs "
            f"({len(peptides)} peptides × {len(alleles)} alleles)"
        )

        # Run EL model
        try:
            _run_model(input_csv, "el", el_out_csv)
        except RuntimeError as exc:
            logger.error(f"BigMHC EL failed: {exc}")
            el_out_csv = None

        # Run IM model
        try:
            _run_model(input_csv, "im", im_out_csv)
        except RuntimeError as exc:
            logger.error(f"BigMHC IM failed: {exc}")
            im_out_csv = None

        if el_out_csv is None and im_out_csv is None:
            return {}

        # Parse EL scores: output column is "BigMHC_EL"
        el_scores: Dict[Tuple[str, str], float] = {}
        if el_out_csv and os.path.isfile(el_out_csv):
            with open(el_out_csv, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    # Column names: mhc, pep, BigMHC_EL
                    allele  = row.get("mhc", "").strip()
                    peptide = row.get("pep", "").strip()
                    score   = row.get("BigMHC_EL") or row.get("bigmhc_el") or "0"
                    try:
                        el_scores[(peptide, allele)] = float(score)
                    except (ValueError, TypeError):
                        pass

        # Parse IM scores: output column is "BigMHC_IM"
        im_scores: Dict[Tuple[str, str], float] = {}
        if im_out_csv and os.path.isfile(im_out_csv):
            with open(im_out_csv, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    allele  = row.get("mhc", "").strip()
                    peptide = row.get("pep", "").strip()
                    score   = row.get("BigMHC_IM") or row.get("bigmhc_im") or "0"
                    try:
                        im_scores[(peptide, allele)] = float(score)
                    except (ValueError, TypeError):
                        pass

        # Merge into results dict; allele key uses original input strings
        all_results: Dict[Tuple[str, str], BigMHCResult] = {}
        # A pair the tool did not score must stay None (abstain from the vote),
        # not 0.0 (a confident negative that also demotes P1→P2 via is_immunogenic
        # and is indistinguishable from a real tiny score after rounding) — review NEO-7.
        for orig_allele, peptide in rows:
            key = (peptide, orig_allele)
            all_results[key] = BigMHCResult(
                el_score=el_scores.get(key),
                im_score=im_scores.get(key),
            )

    logger.info(f"  BigMHC: {len(all_results)} results returned")
    return all_results


# ─── Percentile rank conversion ────────────────────────────────────────────────

def el_score_to_percentile_rank(
    score    : float,
    allele   : str,
    reference: Optional[Dict[str, List[float]]] = None,
) -> float:
    """
    Convert a BigMHC EL score (0–1) to an approximate %rank.

    If no reference is provided, a linear approximation is used:
        %rank ≈ (1 - score) * 100
    This is used as a proxy in consensus_rank computation.
    """
    if reference and allele in reference:
        bg   = reference[allele]
        rank = sum(s >= score for s in bg) / len(bg) * 100
        return round(rank, 3)
    return round(max(0.0, min(100.0, (1.0 - score) * 100.0)), 3)
