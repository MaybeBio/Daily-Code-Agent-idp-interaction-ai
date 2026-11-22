"""Build AlphaFold Server batch files from a candidate list.

The server accepts a JSON array of jobs, at most 100 per file and 5,000 tokens
per job (one token per residue). This module enforces both limits and keeps a
local record of what was submitted, because the free quota is 30 jobs a day.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

MAX_JOBS_PER_FILE = 100
MAX_TOKENS_PER_JOB = 5000
DAILY_JOB_QUOTA = 30


@dataclass
class Job:
    pair_id: str
    chain_a: str
    chain_b: str
    tokens: int

    def as_server_job(self, seed: int | None = None) -> dict:
        # The server docs recommend an empty seed list, which lets it pick one seed.
        return {
            "name": self.pair_id,
            "modelSeeds": [] if seed is None else [seed],
            "sequences": [
                {"proteinChain": {"sequence": self.chain_a, "count": 1}},
                {"proteinChain": {"sequence": self.chain_b, "count": 1}},
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }


@dataclass
class BatchPlan:
    jobs: list[Job] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)

    @property
    def total_tokens(self) -> int:
        return sum(job.tokens for job in self.jobs)

    def files(self, max_jobs: int = MAX_JOBS_PER_FILE) -> list[list[Job]]:
        return [self.jobs[i:i + max_jobs] for i in range(0, len(self.jobs), max_jobs)]


def plan(candidates: list[dict], sequences: dict[str, str]) -> BatchPlan:
    """Turn candidate rows into jobs, dropping anything the server would reject.

    ``candidates`` rows need pair_id, acc_a, acc_b, and optional range_a/range_b
    strings like "1-120". ``sequences`` maps accession to full sequence.
    """
    from .uniprot import slice_sequence

    result = BatchPlan()
    for row in candidates:
        pair_id = row["pair_id"]
        try:
            seq_a = slice_sequence(sequences[row["acc_a"]],
                                   *_parse_range(row.get("range_a", "")))
            seq_b = slice_sequence(sequences[row["acc_b"]],
                                   *_parse_range(row.get("range_b", "")))
        except KeyError as exc:
            result.skipped.append((pair_id, f"missing sequence for {exc}"))
            continue
        except ValueError as exc:
            result.skipped.append((pair_id, str(exc)))
            continue

        tokens = len(seq_a) + len(seq_b)
        if tokens > MAX_TOKENS_PER_JOB:
            result.skipped.append(
                (pair_id, f"{tokens} tokens exceeds the {MAX_TOKENS_PER_JOB} cap; "
                          f"truncate to domains in candidates.csv"))
            continue
        result.jobs.append(Job(pair_id=pair_id, chain_a=seq_a, chain_b=seq_b, tokens=tokens))
    return result


def _parse_range(text: str | None) -> tuple[int | None, int | None]:
    if not text or not str(text).strip():
        return None, None
    start, _, end = str(text).partition("-")
    return int(start), (int(end) if end.strip() else None)


def write_batches(plan_obj: BatchPlan, out_dir: Path, seed: int | None = None,
                  max_jobs: int = MAX_JOBS_PER_FILE) -> list[Path]:
    """Write one JSON file per group of jobs and return the written paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for index, group in enumerate(plan_obj.files(max_jobs), start=1):
        path = out_dir / f"batch_{index:03d}.json"
        path.write_text(json.dumps([job.as_server_job(seed) for job in group], indent=2),
                        encoding="utf-8")
        written.append(path)
    return written


class QuotaTracker:
    """Counts jobs submitted per calendar day against the server's free quota."""

    def __init__(self, path: Path, daily_limit: int = DAILY_JOB_QUOTA):
        self.path = Path(path)
        self.daily_limit = daily_limit

    def _load(self) -> dict:
        if self.path.exists():
            return json.loads(self.path.read_text(encoding="utf-8"))
        return {}

    def remaining(self) -> int:
        used = self._load().get(date.today().isoformat(), 0)
        return max(0, self.daily_limit - used)

    def record(self, jobs: int) -> None:
        data = self._load()
        key = date.today().isoformat()
        data[key] = data.get(key, 0) + jobs
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(data, indent=2), encoding="utf-8")
