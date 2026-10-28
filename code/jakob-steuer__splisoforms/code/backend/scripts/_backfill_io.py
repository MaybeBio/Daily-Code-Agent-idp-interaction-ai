"""
_backfill_io.py
===============
Shared terminal / logfile output utilities for backfill scripts.

Output contract every backfill script should follow
---------------------------------------------------
  [<script_name>] <human title>
    scope: <N items> | <key params>
  <progress bar showing X / N + running counts>
  ─────────────────────────────────────────
    <TITLE> BACKFILL COMPLETE
    <key>           <value>  (<percent>)
    ...
  ─────────────────────────────────────────
    Detailed log: <logfile path>

  • stdout      header, single tqdm line, summary block, logfile path
  • stderr      progress bar only (tqdm default; off when redirected)
  • logfile     per-isoform detail at INFO+, rotated at 10 MB × 2 backups,
                stored under backend/logs/ (which is .gitignored)

Four entry points
-----------------
  setup_logging(name)         configures root logger + returns logfile Path
  setup_worker_logging(path)  ProcessPoolExecutor `initializer=` callback that
                              attaches the same logfile handler in each worker
  Progress(total, desc)       context-managed progress wrapper (tqdm w/ fallback)
  print_summary(...)          bordered final summary block
"""

from __future__ import annotations

import logging
import logging.handlers
import sys
from pathlib import Path
from typing import Iterable, Optional


# Default log directory: backend/logs/ at the project root (ignored in .gitignore)
_DEFAULT_LOG_DIR = Path(__file__).resolve().parent.parent / "logs"


def setup_logging(name: str, log_dir: Optional[Path] = None,
                  stdout_level: int = logging.WARNING,
                  file_level: int = logging.INFO,
                  max_bytes: int = 10 * 1024 * 1024,
                  backup_count: int = 2) -> Path:
    """
    Configure root logger with a quiet stdout handler and a rotating file handler.

    stdout      WARNING+ only (errors, startup failures) — does not compete with the
                tqdm progress bar.
    logfile     INFO+ per-item details, rotated at max_bytes × backup_count copies.
                Path: <log_dir>/<name>.log  (default log_dir is backend/logs/).

    Returns the absolute path to the logfile so the caller can print it in the
    final summary.
    """
    log_dir = log_dir or _DEFAULT_LOG_DIR
    log_dir.mkdir(parents=True, exist_ok=True)
    logfile = log_dir / f"{name}.log"

    root = logging.getLogger()
    # Reset any pre-existing handlers so re-invocations in the same process
    # don't double-attach (matters for tests, smoke runs, etc).
    for h in list(root.handlers):
        root.removeHandler(h)
    root.setLevel(logging.DEBUG)

    # Stdout: quiet
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(stdout_level)
    sh.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    root.addHandler(sh)

    # File: rotating, verbose
    fh = logging.handlers.RotatingFileHandler(
        logfile, maxBytes=max_bytes, backupCount=backup_count,
    )
    fh.setLevel(file_level)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(processName)s: %(message)s",
        datefmt="%H:%M:%S",
    ))
    root.addHandler(fh)

    return logfile


def setup_worker_logging(logfile_path: str, level: int = logging.INFO,
                          max_bytes: int = 10 * 1024 * 1024,
                          backup_count: int = 2) -> None:
    """
    ProcessPoolExecutor `initializer=` callback.

    On macOS / Windows the workers spawn fresh (no shared state with the
    parent), so the parent's setup_logging() handlers are not inherited.
    This callback runs once per worker process and attaches the same
    rotating file handler so worker log records reach the same logfile.

    Errors and warnings still go to stderr at WARNING+ (the default).
    """
    root = logging.getLogger()
    # Avoid duplicate handlers if this worker process is reused
    for h in list(root.handlers):
        if isinstance(h, logging.handlers.RotatingFileHandler) and \
           getattr(h, "baseFilename", None) == str(logfile_path):
            return
    root.setLevel(level)
    fh = logging.handlers.RotatingFileHandler(
        logfile_path, maxBytes=max_bytes, backupCount=backup_count,
    )
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(processName)s: %(message)s",
        datefmt="%H:%M:%S",
    ))
    root.addHandler(fh)


class Progress:
    """
    Progress wrapper that prefers tqdm but gracefully falls back to periodic
    stderr ticks (one line per ~1% of work, capped at ~100 updates total).

    Usage:
        with Progress(total=10000, desc="Pfam") as p:
            for ... :
                p.update(1, status=f"ok={s} skip={k}")
    """

    def __init__(self, total: int, desc: str, unit: str = "iso"):
        self._total = max(total, 1)
        self._desc = desc
        self._n = 0
        self._last_tick = 0
        try:
            from tqdm import tqdm
            self._bar = tqdm(
                total=total, desc=desc, unit=unit, ncols=90,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}] {postfix}",
            )
        except ImportError:
            self._bar = None
            sys.stderr.write(f"(install tqdm for a progress bar: pip install tqdm)\n")

    def update(self, n: int = 1, status: str = "") -> None:
        if self._bar is not None:
            if status:
                self._bar.set_postfix_str(status, refresh=False)
            self._bar.update(n)
            return

        self._n += n
        # In fallback mode, tick every ~1% (min step 100, max step 10000).
        step = min(10000, max(100, self._total // 100))
        if self._n - self._last_tick >= step or self._n >= self._total:
            pct = self._n / self._total * 100
            sys.stderr.write(f"\r{self._desc}: {self._n:,}/{self._total:,} ({pct:5.1f}%) {status}")
            sys.stderr.flush()
            self._last_tick = self._n

    def close(self) -> None:
        if self._bar is not None:
            self._bar.close()
        else:
            sys.stderr.write("\n")

    def __enter__(self) -> "Progress":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def print_summary(title: str,
                  counts: list[tuple[str, int]],
                  total: Optional[int] = None,
                  logfile: Optional[Path] = None,
                  extras: Optional[list[tuple[str, list[tuple[str, int]]]]] = None,
                  notes: Optional[list[str]] = None,
                  errors: Optional[list[tuple]] = None,
                  errors_preview: int = 10) -> None:
    """
    Print a bordered, two-column summary block.

    counts     ordered list of (label, count) rows; labels can start with
               ✓ / – / ✗ for visual status.  If `total` is provided, percentages
               are appended.
    extras     optional named sub-tables — each is (heading, list[(label, count)]).
               Sub-rows get the same count+percentage rendering as the main counts.
    notes      optional list of free-form lines (no count, no percentage); rendered
               under a "Notes:" heading.  Use this for output paths, parameters, etc.
    errors     optional list of (id, msg) tuples; first `errors_preview` are
               shown inline, rest are summarised as "… N more (see logfile)".
    logfile    optional logfile path appended as the last line.
    """
    if total is None:
        total = sum(n for _, n in counts) or 1

    label_w = max((len(lbl) for lbl, _ in counts), default=8) + 2
    border = "─" * 60

    print()
    print(border)
    print(f"  {title}")
    print(border)
    for label, n in counts:
        pct = n / total * 100
        print(f"  {label:<{label_w}} {n:>8,}  ({pct:5.1f}%)")

    if extras:
        for heading, rows in extras:
            if not rows:
                continue
            print()
            print(f"  {heading}")
            sub_w = max(len(r[0]) for r in rows) + 2
            for label, n in rows:
                pct = n / total * 100
                print(f"    {label:<{sub_w}} {n:>8,}  ({pct:5.1f}%)")

    if notes:
        print()
        print("  Notes:")
        for line in notes:
            print(f"    {line}")

    if errors:
        print()
        print(f"  Errors ({len(errors)}):")
        for entry in errors[:errors_preview]:
            print(f"    {entry[0]}: {entry[1] if len(entry) > 1 else ''}")
        if len(errors) > errors_preview:
            print(f"    … and {len(errors) - errors_preview} more (see logfile)")

    print(border)
    if logfile:
        print(f"  Detailed log: {logfile}")


def exit_status(total: int, errors: int, max_error_fraction: float = 0.5) -> int:
    """
    Process exit code for a backfill: 0 when the run did useful work, 1 when it
    did not. "PASS" in run_local_backfills.sh is keyed on the exit code, and the
    2026-10 review found several steps that errored on every item yet exited 0
    (review ENG-2). A run is a failure when there was work to do and every item
    errored, or more than `max_error_fraction` of items did.
    """
    if total <= 0:
        return 0
    if errors >= total:
        return 1
    return 1 if errors > max_error_fraction * total else 0
