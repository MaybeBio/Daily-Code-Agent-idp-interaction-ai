"""
netmhcpan_runner.py
====================
Thin subprocess wrapper for NetMHCpan-4.2.

NetMHCpan is a binary-only tool distributed by DTU Health Tech:
  https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/

Installation:
  1. Register and download the Linux/macOS tarball.
  2. Extract and set NETMHCPAN_BIN env-var to the `netMHCpan` binary path,
     OR place it on $PATH.

Usage:
    from app.analysis.netmhcpan_runner import run_netmhcpan, is_available

    if is_available():
        results = run_netmhcpan(peptides, allele)
        # results[(peptide, allele)] = NetMHCpanResult(affinity_nm, rank_el, rank_ba)

NetMHCpan allele format: "HLA-A02:01" (no asterisk), e.g. from "HLA-A*02:01".
This runner normalises automatically.
"""

import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ─── Configuration ─────────────────────────────────────────────────────────────
# Override by setting the NETMHCPAN_BIN environment variable.
_DEFAULT_BINARY_NAMES = ["netMHCpan", "netmhcpan", "NetMHCpan"]


@dataclass
class NetMHCpanResult:
    affinity_nm : float         # BA IC50 in nM
    rank_el     : float         # EL %rank (lower = more likely presented)
    rank_ba     : float         # BA %rank (binding affinity percentile)


# ─── Availability check ────────────────────────────────────────────────────────

def _find_binary() -> Optional[str]:
    """Return path to the netMHCpan binary, or None if not found."""
    env_path = os.environ.get("NETMHCPAN_BIN")
    if env_path and os.path.isfile(env_path) and os.access(env_path, os.X_OK):
        return env_path
    for name in _DEFAULT_BINARY_NAMES:
        found = shutil.which(name)
        if found:
            return found
    return None


def _resolve_native_binary(wrapper_path: str) -> Optional[Tuple[str, str]]:
    """Find the actual compiled binary inside a NetMHCpan installation.

    The tcsh wrapper calls $NMHOME/Linux_x86_64/bin/netMHCpan-4.2 (or similar).
    When NMHOME has a long path, the binary crashes with buffer overflow.

    To work around both issues this function:
    1. Creates a short symlink /tmp/nmhc -> real NMHOME (avoids buffer overflow)
    2. Creates NMHOME/bin -> NMHOME/Linux_x86_64/bin if missing (binary expects
       $NETMHCpan/bin/ for helpers like nnalign)

    Returns (native_binary_via_short_path, short_nmhome) or None if not found.
    """
    real_nmhome = os.path.dirname(os.path.abspath(wrapper_path))
    arch_bin_dir = None
    native_name = None

    for arch_dir in ("Linux_x86_64", "Linux_ia64"):
        bin_dir = os.path.join(real_nmhome, arch_dir, "bin")
        if not os.path.isdir(bin_dir):
            continue
        for name in os.listdir(bin_dir):
            if name.lower().startswith("netmhcpan"):
                candidate = os.path.join(bin_dir, name)
                if os.access(candidate, os.X_OK):
                    arch_bin_dir = os.path.join(arch_dir, "bin")
                    native_name = os.path.join(arch_bin_dir, name)
                    break
        if native_name:
            break

    if not native_name:
        return None

    # NMHOME/bin must resolve to the arch bin dir: the binary execs its helper
    # (nnalign_…) from $NMHOME/bin. A DANGLING link (left behind when the install
    # directory was moved — the 2026-07 run lost NetMHCpan on 18/40 shards this
    # way) passes os.path.exists() checks on the link itself but not on its target,
    # so repair it rather than leave it. Use a RELATIVE link so a future move keeps
    # working.
    root_bin = os.path.join(real_nmhome, "bin")
    try:
        if os.path.islink(root_bin) and not os.path.isdir(root_bin):
            logger.warning(f"Repairing dangling symlink {root_bin} -> {os.readlink(root_bin)}")
            os.unlink(root_bin)
        if not os.path.exists(root_bin):
            os.symlink(arch_bin_dir, root_bin)
            logger.info(f"Created symlink {root_bin} -> {arch_bin_dir}")
    except OSError as e:
        logger.warning(f"Could not create bin symlink: {e}")

    short_link = os.path.join(tempfile.gettempdir(), "nmhc_root")
    try:
        if os.path.islink(short_link):
            if os.readlink(short_link) != real_nmhome:
                os.unlink(short_link)
                os.symlink(real_nmhome, short_link)
        elif os.path.isdir(short_link):
            import shutil as _sh
            _sh.rmtree(short_link)
            os.symlink(real_nmhome, short_link)
        elif not os.path.exists(short_link):
            os.symlink(real_nmhome, short_link)
    except OSError as e:
        logger.warning(f"Could not create short symlink: {e}; using real path")
        short_link = real_nmhome

    binary = os.path.join(short_link, native_name)
    return binary, short_link


_netmhcpan_available: Optional[bool] = None  # None = not yet probed


def is_available() -> bool:
    """
    Return True if a usable NetMHCpan binary is found AND actually runs.
    Result is cached after the first call.
    """
    global _netmhcpan_available
    if _netmhcpan_available is not None:
        return _netmhcpan_available

    binary = _find_binary()
    if binary is None:
        _netmhcpan_available = False
        return False

    # Probe: run with -version to verify the binary works on this OS/arch.
    # (-h causes netMHCpan to wait for stdin; -version exits immediately.)
    try:
        result = subprocess.run(
            [binary, "-version"],
            capture_output=True, text=True, timeout=10,
        )
        output = result.stdout + result.stderr
        if "no binaries found" in output or "not found" in output.lower():
            logger.warning(
                f"NetMHCpan binary found at {binary} but does not run on this "
                f"platform ({output.strip().splitlines()[0]}). "
                "NetMHCpan is Linux-only; it will be skipped on macOS."
            )
            _netmhcpan_available = False
        else:
            # `-version` succeeds even when the helper binaries are unreachable
            # (the per-allele run then fails with "nnalign_…: not found"). Check
            # that $NMHOME/bin actually resolves to a directory holding a helper.
            resolved = _resolve_native_binary(binary)
            helper_ok = True
            if resolved:
                root_bin = os.path.join(os.path.dirname(os.path.abspath(binary)), "bin")
                helper_ok = os.path.isdir(root_bin) and any(
                    n.startswith("nnalign") for n in os.listdir(root_bin))
            if not helper_ok:
                logger.warning(
                    f"NetMHCpan wrapper at {binary} runs, but $NMHOME/bin does not "
                    f"contain the nnalign helper — predictions would fail per allele; "
                    "marking unavailable.")
            _netmhcpan_available = helper_ok
    except Exception as exc:
        logger.warning(f"NetMHCpan probe failed ({exc}); marking unavailable.")
        _netmhcpan_available = False

    return _netmhcpan_available


# ─── Allele normalisation ──────────────────────────────────────────────────────

def _normalise_allele(allele: str) -> str:
    """
    Convert MHCflurry-style allele names to NetMHCpan format.
    "HLA-A*02:01"  →  "HLA-A02:01"
    "HLA-B*07:02"  →  "HLA-B07:02"
    """
    return allele.replace("*", "")


# ─── Output parser ─────────────────────────────────────────────────────────────

def _parse_netmhcpan_output(stdout: str, allele: str) -> Dict[str, NetMHCpanResult]:
    """
    Parse the text output from `netMHCpan -BA` (combined BA+EL run).

    NetMHCpan-4.2 combined output columns (space-separated, fixed-width):
      Pos  MHC  Peptide  Core  Of  Gp  Gl  Ip  Il  Icore  Identity
      Score_EL  %Rank_EL  Score_BA  %Rank_BA  Affinity(nM)  BindLevel

    Note on column order: NetMHCpan-4.2's combined output places %Rank_BA
    BEFORE Affinity(nM) — opposite to what an earlier version of this parser
    assumed.  Migration g8e1d4c7b2a9 fixed the swapped historical values.

    Returns dict: peptide → NetMHCpanResult
    """
    results: Dict[str, NetMHCpanResult] = {}

    for line in stdout.splitlines():
        line = line.strip()
        # Skip blank, header, separator, and summary lines
        if not line or line.startswith("#") or line.startswith("-") or line.startswith("Pos"):
            continue
        # Data lines start with a numeric position
        parts = line.split()
        if len(parts) < 16:
            continue
        try:
            int(parts[0])  # first column is position (integer)
        except ValueError:
            continue

        try:
            peptide     = parts[2]
            score_el    = float(parts[11])   # Score_EL
            rank_el     = float(parts[12])   # %Rank_EL
            score_ba    = float(parts[13])   # Score_BA
            rank_ba     = float(parts[14])   # %Rank_BA   (was incorrectly read as parts[15])
            affinity_nm = float(parts[15])   # Affinity(nM) (was incorrectly read as parts[14])
            results[peptide] = NetMHCpanResult(
                affinity_nm=affinity_nm,
                rank_el=rank_el,
                rank_ba=rank_ba,
            )
        except (IndexError, ValueError) as exc:
            logger.debug(f"Could not parse NetMHCpan line: {line!r} — {exc}")

    return results


# ─── Main runner ───────────────────────────────────────────────────────────────

def run_netmhcpan(
    peptides : List[str],
    allele   : str,
    *,
    batch_size: int = 500,
) -> Dict[Tuple[str, str], NetMHCpanResult]:
    """
    Run NetMHCpan-4.2 in combined BA+EL mode for one allele against a list of peptides.

    Args:
        peptides:   Unique peptide sequences (8–11 aa).
        allele:     HLA allele in MHCflurry format (e.g. "HLA-A*02:01").
        batch_size: Max peptides per subprocess call (avoids arg-length limits).

    Returns:
        Dict mapping (peptide, allele) → NetMHCpanResult.
        Peptides that fail prediction are omitted.
    """
    wrapper = _find_binary()
    if wrapper is None:
        raise RuntimeError(
            "NetMHCpan binary not found. Set NETMHCPAN_BIN or add netMHCpan to $PATH."
        )

    norm_allele = _normalise_allele(allele)
    all_results: Dict[Tuple[str, str], NetMHCpanResult] = {}

    resolved = _resolve_native_binary(wrapper)
    if resolved:
        native_bin, nmhome = resolved
        short_tmpdir = os.path.join(tempfile.gettempdir(), "nmhc_tmp")
        os.makedirs(short_tmpdir, exist_ok=True)
        run_env = os.environ.copy()
        run_env["NMHOME"] = nmhome
        run_env["NETMHCpan"] = nmhome
        run_env["TMPDIR"] = short_tmpdir
        binary = native_bin
        logger.debug(f"Using native binary {native_bin} with TMPDIR={short_tmpdir}")
    else:
        run_env = None
        binary = wrapper

    # Process in batches to keep individual calls manageable
    for batch_start in range(0, len(peptides), batch_size):
        batch = peptides[batch_start : batch_start + batch_size]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pep", delete=False) as fh:
            for pep in batch:
                fh.write(pep + "\n")
            tmp_path = fh.name

        try:
            cmd = [
                binary,
                "-a", norm_allele,
                "-f", tmp_path,
                "-BA",           # combined BA+EL output
                "-s",            # sort by %Rank_EL
                "-inptype", "1", # input is peptide list (not FASTA)
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
                env=run_env,
            )
            if result.returncode != 0:
                logger.warning(
                    f"NetMHCpan returned exit {result.returncode} for allele {allele}: "
                    f"{result.stderr[:200]}"
                )
                continue

            parsed = _parse_netmhcpan_output(result.stdout, norm_allele)
            if not parsed:
                if result.stdout.strip():
                    stdout_lines = result.stdout.strip().splitlines()
                    logger.warning(
                        f"NetMHCpan produced output ({len(stdout_lines)} lines) but "
                        f"0 parsed for {allele}. First 5 lines:\n"
                        + "\n".join(stdout_lines[:5])
                    )
                else:
                    logger.warning(
                        f"NetMHCpan returned empty stdout for {allele} "
                        f"(exit={result.returncode}, stderr={result.stderr[:500]})"
                    )
            for peptide, res in parsed.items():
                all_results[(peptide, allele)] = res

        except subprocess.TimeoutExpired:
            logger.error(f"NetMHCpan timed out for allele {allele}, batch starting at {batch_start}")
        except Exception as exc:
            logger.error(f"NetMHCpan error for allele {allele}: {exc}")
        finally:
            os.unlink(tmp_path)

    return all_results


def _netmhcpan_threads() -> int:
    """Alleles are independent subprocess calls, so run them concurrently.
    NETMHCPAN_THREADS > SLURM_CPUS_PER_TASK > 4 (capped at 8: the binary is
    single-threaded and the box's other work needs air)."""
    for var in ("NETMHCPAN_THREADS", "SLURM_CPUS_PER_TASK"):
        v = os.environ.get(var)
        if v and v.isdigit() and int(v) > 0:
            return max(1, min(int(v), 8))
    return 4


def run_netmhcpan_multi_allele(
    peptides : List[str],
    alleles  : List[str],
    *,
    batch_size: int = 500,
    threads: Optional[int] = None,
) -> Dict[Tuple[str, str], NetMHCpanResult]:
    """
    Run NetMHCpan for multiple alleles concurrently and merge results.

    Sequential allele calls were the dominant cost of the HPC pipeline (~50 s per
    isoform for 35 alleles on an 8-CPU allocation that mostly idled; shard 1 of the
    2026-10 restage hit the 36 h wall limit). Each allele is one single-threaded
    subprocess, so a small thread pool scales almost linearly. Logs progress every
    5 alleles.
    """
    all_results: Dict[Tuple[str, str], NetMHCpanResult] = {}
    if not peptides or not alleles:
        return all_results
    n_threads = threads or _netmhcpan_threads()

    # Resolve/repair symlinks once, in the calling thread, before fanning out.
    wrapper = _find_binary()
    if wrapper:
        _resolve_native_binary(wrapper)

    if n_threads <= 1 or len(alleles) == 1:
        for i, allele in enumerate(alleles, 1):
            try:
                all_results.update(run_netmhcpan(peptides, allele, batch_size=batch_size))
                if i % 5 == 0:
                    logger.info(f"  NetMHCpan: {i}/{len(alleles)} alleles done, {len(all_results)} results")
            except Exception as exc:
                logger.warning(f"  NetMHCpan skipped allele {allele}: {exc}")
        return all_results

    from concurrent.futures import ThreadPoolExecutor, as_completed
    done = 0
    with ThreadPoolExecutor(max_workers=min(n_threads, len(alleles))) as pool:
        futs = {pool.submit(run_netmhcpan, peptides, a, batch_size=batch_size): a for a in alleles}
        for fut in as_completed(futs):
            allele = futs[fut]
            try:
                all_results.update(fut.result())
            except Exception as exc:
                logger.warning(f"  NetMHCpan skipped allele {allele}: {exc}")
            done += 1
            if done % 5 == 0:
                logger.info(f"  NetMHCpan: {done}/{len(alleles)} alleles done, "
                            f"{len(all_results)} results ({n_threads} threads)")
    return all_results
