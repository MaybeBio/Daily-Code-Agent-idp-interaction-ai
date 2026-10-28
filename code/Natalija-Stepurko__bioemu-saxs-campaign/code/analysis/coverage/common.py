"""Shared paths and helpers for the SAXS coverage analysis."""
import datetime
import io
import json
import os
import pathlib
import struct
import time
import zipfile

import requests

# Locations. COVERAGE_WORK holds the downloads and the intermediate tables (about 0.9 GB, not tracked);
# COVERAGE_DOWNLOADS can point the download cache elsewhere. The SAXS table is read from the main
# pipeline's data directory (BSC_ROOT, as for `bsc fetch`); the clean-entry flags from the tracked results.
REPO = pathlib.Path(__file__).resolve().parents[2]
ROOT = pathlib.Path(os.environ.get("COVERAGE_WORK", REPO / "analysis" / "coverage" / "work"))
DL = pathlib.Path(os.environ.get("COVERAGE_DOWNLOADS", ROOT / "downloads"))
TAB = ROOT / "tables"
FIG = ROOT / "figures"
BSC_ROOT = pathlib.Path(os.environ.get("BSC_ROOT", REPO))
SAXS_TABLE = BSC_ROOT / "data" / "PeptoneDB-SAXS" / "PeptoneDB-SAXS.csv"
RESOLV = REPO / "results" / "resolvability.csv"
ENTRIES = REPO / "results" / "entries.csv"
EXPORT = REPO / "results" / "coverage"   # the small outputs the page and the `coverage` stage read
VERSIONS = DL / "versions.json"
# BioEmu-1 training data (scripts 11-17): cluster definitions, Zenodo metadata and the topology and metadata
# members read out of the remote archives. MMseqs2 is a static binary (MMSEQS, default under COVERAGE_TOOLS).
BT = DL / "bioemu_training"
TOOLS = pathlib.Path(os.environ.get("COVERAGE_TOOLS", ROOT / "tools"))
MMSEQS = os.environ.get("MMSEQS", str(TOOLS / "mmseqs" / "bin" / "mmseqs"))
THREADS = int(os.environ.get("COVERAGE_THREADS", "6"))
for _d in [TAB, FIG] + [DL / x for x in ("cath", "sasbdb", "uniprot", "rcsb", "interpro", "pdbe", "ted")] + [BT]:
    _d.mkdir(parents=True, exist_ok=True)

S = requests.Session()
S.headers["User-Agent"] = "saxs-coverage-analysis/1.0 (python-requests)"


def get(url, params=None, tries=5, timeout=60, ok404=False, **kw):
    for i in range(tries):
        try:
            r = S.get(url, params=params, timeout=timeout, **kw)
            if r.status_code == 404 and ok404:
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2 * (i + 1)); continue
            r.raise_for_status()
            return r
        except requests.RequestException:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"failed {url}")


def cached_json(path, url, params=None, ok404=False, **kw):
    path = pathlib.Path(path)
    if path.exists():
        return json.loads(path.read_text())
    r = get(url, params=params, ok404=ok404, **kw)
    if r.status_code == 404 or r.status_code == 204 or not r.content:
        data = {"__status__": r.status_code}
    else:
        data = r.json()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return data


def record_version(key, info):
    """Record the release of a source the first time it is fetched. An existing record is kept: the
    versions describe the cached downloads, and deleting the download directory removes both."""
    v = json.loads(VERSIONS.read_text()) if VERSIONS.exists() else {}
    if key in v:
        return
    info = dict(info)
    info.setdefault("accessed", datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"))
    v[key] = info
    VERSIONS.write_text(json.dumps(v, indent=2))


ORG_CATS = ["Human", "Mouse", "E. coli", "Yeast", "Other bacteria", "Viruses", "Plants", "Other"]


def org_category(taxid, sci_name, lineage_names):
    """Map an organism to the eight comparison groups. lineage_names: iterable of names."""
    name = sci_name.strip() if isinstance(sci_name, str) else ""
    lin = {x.strip() for x in lineage_names if isinstance(x, str) and x}
    try:
        taxid = int(taxid)
    except (TypeError, ValueError):
        taxid = None
    if taxid == 9606 or name.startswith("Homo sapiens"):
        return "Human"
    if taxid == 10090 or name.startswith("Mus musculus"):
        return "Mouse"
    if name.startswith("Escherichia coli") or "Escherichia coli" in lin:
        return "E. coli"
    if (name.startswith("Saccharomyces cerevisiae") or name.startswith("Schizosaccharomyces pombe")
            or "Saccharomyces cerevisiae" in lin or "Schizosaccharomyces pombe" in lin):
        return "Yeast"
    if "Bacteria" in lin:
        return "Other bacteria"
    if "Viruses" in lin:
        return "Viruses"
    if "Viridiplantae" in lin:
        return "Plants"
    return "Other"


class HTTPRangeFile(io.RawIOBase):
    """Seekable read-only view of a remote file through HTTP range requests, enough for zipfile.ZipFile to
    list an archive and extract single members without downloading it. The central directory is fetched in
    one request; other reads fetch the requested bytes plus a small read-ahead."""

    def __init__(self, url, readahead=1 << 16):
        self.url, self.pos, self.ra, self.regions, self.nreq = url, 0, readahead, [], 0
        self.size = None
        self._get(0, 0)

    def _get(self, a, b):
        for t in range(8):
            try:
                r = S.get(self.url, headers={"Range": f"bytes={a}-{b}"}, timeout=(20, 180))
                self.nreq += 1
                if r.status_code == 206:
                    self.size = int(r.headers["Content-Range"].split("/")[1])
                    return r.content
                wait = int(r.headers.get("Retry-After", 0) or 0) or 5 * (t + 1)
                print(f"  HTTP {r.status_code} for bytes {a}-{b}; waiting {wait} s", flush=True)
                time.sleep(wait)
            except requests.RequestException as e:
                print(f"  {type(e).__name__}; retrying", flush=True)
                time.sleep(5 * (t + 1))
        raise RuntimeError(f"range {a}-{b} of {self.url} failed")

    def prefetch_central_directory(self):
        a = max(0, self.size - 65558)
        tail = self._get(a, self.size - 1)
        i = tail.rfind(b"PK\x05\x06")
        cd_size, cd_off = struct.unpack("<II", tail[i + 12:i + 20])
        if cd_off == 0xFFFFFFFF or cd_size == 0xFFFFFFFF:          # zip64 end-of-central-directory record
            j = tail.rfind(b"PK\x06\x06")
            cd_size, cd_off = struct.unpack("<QQ", tail[j + 40:j + 56])
        self.regions = [(a, tail), (cd_off, self._get(cd_off, cd_off + cd_size - 1))]

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = max(0, min(n, self.size - self.pos))
        for a, buf in self.regions:
            if a <= self.pos and self.pos + n <= a + len(buf):
                out = buf[self.pos - a:self.pos - a + n]
                self.pos += n
                return out
        buf = self._get(self.pos, min(self.size, self.pos + max(n, self.ra)) - 1)
        self.regions = self.regions[:2] + [(self.pos, buf)]      # the central directory and one read-ahead
        out = buf[:n]
        self.pos += len(out)
        return out

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def open_remote_zip(url):
    """zipfile.ZipFile over HTTP range requests; returns (ZipFile, HTTPRangeFile)."""
    f = HTTPRangeFile(url)
    f.prefetch_central_directory()
    return zipfile.ZipFile(f), f


def read_remote_member(url, info, tries=8):
    """One member of a remote zip in a single range request, from its ZipInfo (stored or deflated)."""
    import zlib
    a = info.header_offset
    b = a + 30 + len(info.filename.encode()) + 4096 + info.compress_size
    for t in range(tries):
        try:
            r = S.get(url, headers={"Range": f"bytes={a}-{b}"}, timeout=(20, 180))
            if r.status_code == 206:
                break
            time.sleep(int(r.headers.get("Retry-After", 0) or 0) or 5 * (t + 1))
        except requests.RequestException:
            time.sleep(5 * (t + 1))
    else:
        raise RuntimeError(f"{info.filename}: range request failed")
    buf = r.content
    assert buf[:4] == b"PK\x03\x04", info.filename
    n, e = struct.unpack("<HH", buf[26:30])
    data = buf[30 + n + e:30 + n + e + info.compress_size]
    if info.compress_type == zipfile.ZIP_DEFLATED:
        data = zlib.decompress(data, -15)
    elif info.compress_type != zipfile.ZIP_STORED:
        raise ValueError(f"{info.filename}: compression {info.compress_type}")
    assert len(data) == info.file_size and zlib.crc32(data) == info.CRC, info.filename
    return data
