"""Named stages with a run log: `bsc fetch`, `bsc score`, `bsc coordrg`, `bsc analyse`, `bsc select_eval`,
`bsc predict`, `bsc envelopes`, `bsc coverage`, `bsc report`, `bsc all`."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

from bsc import config as C

STAGES = ["fetch", "score", "coordrg", "analyse", "select_eval", "predict", "envelopes", "coverage", "report"]


def _git() -> dict:
    try:
        rev = subprocess.run(["git", "-C", str(C.REPO), "rev-parse", "HEAD"], capture_output=True, text=True)
        dirty = subprocess.run(["git", "-C", str(C.REPO), "status", "--porcelain", "--", ".", ":!results"],
                               capture_output=True, text=True)
        return {"commit": rev.stdout.strip(), "dirty": bool(dirty.stdout.strip())}
    except OSError:
        return {}


def _inputs_sha() -> dict:
    rec = C.DATA / "checksums.json"
    return json.load(open(rec)) if rec.exists() else {}


def run(stage: str, argv: list[str]) -> None:
    t0 = time.time()
    if stage == "fetch":
        from bsc import fetch
        ap = argparse.ArgumentParser(prog="bsc fetch")
        ap.add_argument("--archive-dir", type=Path, default=None)
        ap.add_argument("--models", nargs="*", default=None)
        a = ap.parse_args(argv)
        fetch.main(a.models, a.archive_dir)
    elif stage == "score":
        from bsc import score
        ap = argparse.ArgumentParser(prog="bsc score")
        ap.add_argument("--models", nargs="*", default=None)
        a = ap.parse_args(argv)
        score.main(a.models)
    elif stage == "coordrg":
        from bsc import coordrg
        ap = argparse.ArgumentParser(prog="bsc coordrg")
        ap.add_argument("--archive", type=Path, default=None,
                        help="PeptoneBench Predictions.tar.gz; missing ensembles are extracted in one pass")
        a = ap.parse_args(argv)
        coordrg.main(a.archive)
    elif stage == "analyse":
        from bsc import analyse
        analyse.main()
    elif stage == "select_eval":
        from bsc import select_eval
        select_eval.main()
    elif stage == "predict":
        from bsc import predict
        predict.main()
    elif stage == "envelopes":
        from bsc import envelopes
        envelopes.main()
    elif stage == "coverage":
        from bsc import coverage
        coverage.main()
    elif stage == "report":
        from bsc import report
        report.main()
    else:
        sys.exit(f"unknown stage {stage!r}; stages: {', '.join(STAGES)}")
    C.RESULTS.mkdir(parents=True, exist_ok=True)
    log = C.RESULTS / "run_log.json"
    hist = json.load(open(log)) if log.exists() else []
    hist.append({"stage": stage, "args": argv, "seconds": round(time.time() - t0, 1),
                 "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t0)),
                 **_git(), "inputs_sha256": _inputs_sha(),
                 "config_sha256": hashlib.sha256(Path(C.__file__).read_bytes()).hexdigest()[:16]})
    json.dump(hist, open(log, "w"), indent=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=[*STAGES, "all", "list"])
    ap.add_argument("rest", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    rest = [x for x in a.rest if x != "--"]
    if a.stage == "list":
        for s in STAGES:
            print(f"  {s}")
        return
    for s in (STAGES if a.stage == "all" else [a.stage]):
        print(f"▶ {s}", flush=True)
        run(s, rest if a.stage != "all" else [])
        print("✓ done", flush=True)


if __name__ == "__main__":
    main()
