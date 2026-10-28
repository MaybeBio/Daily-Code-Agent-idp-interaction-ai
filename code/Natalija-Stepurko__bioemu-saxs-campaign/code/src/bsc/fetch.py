"""Download the public inputs from Zenodo with checksum verification, and unpack what the
analysis needs. Nothing here is redistributed by this repository."""
from __future__ import annotations

import hashlib
import json
import tarfile
from pathlib import Path

import requests

from bsc import config as C


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def zenodo_files() -> dict[str, dict]:
    r = requests.get(C.ZENODO_API, timeout=60)
    r.raise_for_status()
    return {f["key"]: f for f in r.json()["files"]}


def download(name: str, dest: Path, remote: dict) -> Path:
    """Stream one Zenodo file to dest.part, check its md5 against the record, then rename."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    expected = remote["checksum"].split(":")[-1]
    if dest.exists() and _md5(dest) == expected:
        return dest
    part = dest.with_suffix(dest.suffix + ".part")
    with requests.get(remote["links"]["self"], stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(part, "wb") as f:
            for chunk in r.iter_content(1 << 22):
                f.write(chunk)
    if _md5(part) != expected:
        part.unlink()
        raise RuntimeError(f"{name}: md5 mismatch after download; delete and rerun")
    part.rename(dest)
    return dest


def record_checksums(paths: dict[str, Path]) -> None:
    rec = C.DATA / "checksums.json"
    C.DATA.mkdir(parents=True, exist_ok=True)
    old = json.load(open(rec)) if rec.exists() else {}
    new = {k: _sha256(p) for k, p in paths.items()}
    for k, v in new.items():
        if k in old and old[k] != v:
            raise RuntimeError(f"{k}: sha256 changed since first download ({old[k][:12]} -> {v[:12]})")
    old.update(new)
    json.dump(old, open(rec, "w"), indent=2)


def extract_saxs(archive: Path) -> None:
    """Unpack the SAXS table and curves only."""
    with tarfile.open(archive) as t:
        members = [m for m in t.getmembers() if m.name.lstrip("./").startswith("PeptoneDB-SAXS/")]
        t.extractall(C.DATA, members=members, filter="data")


def extract_predictions(archive: Path, models: list[str]) -> None:
    """Unpack the back-calculated SAXS curves for the named models only (the conformers are not needed)."""
    want = tuple(f"Predictions/PeptoneDB-SAXS-expt/{m}/" for m in models)
    with tarfile.open(archive) as t:
        members = [m for m in t.getmembers() if m.name.startswith(want)]
        t.extractall(C.DATA, members=members, filter="data")


def fetch_sasbdb_summaries(labels: list[str]) -> None:
    """Entry summaries (measured molecular weight, Guinier Rg, project title) from the SASBDB REST API."""
    import time
    out = json.load(open(C.SASBDB_SUMMARY)) if C.SASBDB_SUMMARY.exists() else {}
    keys = ("experimental_mw", "guinier_i0_mw", "porod_mw", "guinier_rg", "pddf_dmax", "symmetry")
    for lab in labels:
        if lab in out:
            continue
        r = requests.get(C.SASBDB_API.format(label=lab), timeout=60)
        r.raise_for_status()
        d = r.json()
        out[lab] = {k: d.get(k) for k in keys}
        out[lab]["title"] = (d.get("project") or {}).get("title")
        time.sleep(0.1)
    json.dump(out, open(C.SASBDB_SUMMARY, "w"), indent=1)


def main(models: list[str] | None = None, archive_dir: Path | None = None) -> None:
    """Fetch and unpack. `archive_dir` can point at already-downloaded archives."""
    models = models or [C.MODEL, *C.COMPARATORS]
    remote = zenodo_files()
    paths = {}
    for name in C.ZENODO_FILES:
        dest = (archive_dir or C.DATA / "archives") / name
        paths[name] = download(name, dest, remote[name])
        print(f"  {name}: ok ({dest.stat().st_size / 1e9:.2f} GB)", flush=True)
    record_checksums(paths)
    extract_saxs(paths["PeptoneDBs.tar.gz"])
    extract_predictions(paths["Predictions.tar.gz"], models)
    n = len(list(C.SAXS_DIR.glob("*.dat")))
    print(f"  {n} experimental curves; predictions for {', '.join(models)}", flush=True)
    import pandas as pd
    fetch_sasbdb_summaries(pd.read_csv(C.SAXS_TABLE)["label"].tolist())
    print(f"  SASBDB summaries for {len(json.load(open(C.SASBDB_SUMMARY)))} entries", flush=True)
