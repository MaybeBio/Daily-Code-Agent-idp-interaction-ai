#!/usr/bin/env python3
"""
fetch_gencode_reference.py
==========================
Phase 0 of the GENCODE v47 reference build (see audit-trail / commit history
around 2026-05-20 for the design rationale).

Downloads the GENCODE 47 (= Ensembl 114) annotation + protein FASTA + MANE
Select summary and parses them into a normalised TSV with one row per
protein-coding transcript.  This TSV is the authoritative reference that
all subsequent phases (cleanup, rename, predict-queue) align against.

Outputs
-------
    data/gencode/v47/
        gencode.v47.annotation.gtf.gz                  # source GTF (kept)
        gencode.v47.pc_translations.fa.gz              # source protein FASTA
        MANE.GRCh38.v1.4.summary.txt.gz                # MANE Select table
        CHECKSUMS                                      # source server's md5s
        source_urls.txt                                # provenance
        gencode_v47.tsv                                # parsed reference (one row/transcript)

The parsed TSV columns:
    ensembl_transcript_id     e.g. ENST00000280772     (no version suffix)
    ensembl_protein_id        e.g. ENSP00000284984     (no version suffix)
    ensembl_gene_id           e.g. ENSG00000151150     (no version suffix)
    gene_symbol               e.g. ANK3
    transcript_name           e.g. ANK3-219            (Ensembl's stable display name)
    transcript_type           e.g. protein_coding
    protein_length            integer (aa count)
    protein_md5               32-char md5 of the protein sequence
    chromosome                e.g. chr10               (kept verbatim from GTF)
    strand                    "+" or "-"
    exon_count                integer
    cds_start_min, cds_end_max    genomic min/max across CDS exons (1-based, inclusive)
    is_mane_select            "1" | "0"
    is_mane_plus_clinical     "1" | "0"
    appris_principal          "" | "principal1" .. "principal5"  (APPRIS principal-level tag verbatim)
    has_non_canonical_aa      "1" | "0"     (X / B / Z / U present in sequence)
    protein_sequence          raw amino-acid string

Run
---
    cd backend
    python scripts/fetch_gencode_reference.py
    python scripts/fetch_gencode_reference.py --force   # re-download even if cached
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import logging
import os
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

# ── Path setup ────────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent))

from _backfill_io import Progress, setup_logging  # noqa: E402

log = logging.getLogger("fetch_gencode_reference")


# ── Constants ─────────────────────────────────────────────────────────────────

# Versioned subdirectory; bump alongside the GENCODE release.
GENCODE_VERSION = "v47"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_OUT_DIR = REPO_ROOT / "data" / "gencode" / GENCODE_VERSION

GENCODE_BASE = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47"
MANE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/"
    "MANE.GRCh38.v1.4.summary.txt.gz"
)

SOURCES = [
    (f"{GENCODE_BASE}/gencode.v47.annotation.gtf.gz",       "gencode.v47.annotation.gtf.gz"),
    (f"{GENCODE_BASE}/gencode.v47.pc_translations.fa.gz",   "gencode.v47.pc_translations.fa.gz"),
    (f"{GENCODE_BASE}/MD5SUMS",                             "MD5SUMS_gencode"),
    (MANE_URL,                                              "MANE.GRCh38.v1.4.summary.txt.gz"),
]

NON_CANONICAL_AA = set("XBZU*")   # AF3 cannot reliably handle these


# ── HTTP helpers ──────────────────────────────────────────────────────────────

def _download(url: str, dest: Path, force: bool = False) -> None:
    """Download `url` → `dest` with simple progress; skip if already present."""
    if dest.exists() and not force:
        size_mb = dest.stat().st_size / 1024 / 1024
        log.info("  cached: %s (%.1f MB)", dest.name, size_mb)
        return
    log.info("  downloading: %s", url)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "splisoforms/1.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            total = int(resp.headers.get("Content-Length", "0"))
            with open(dest, "wb") as f:
                downloaded = 0
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total > 0 and downloaded % (32 * 1024 * 1024) < 1024 * 1024:
                        pct = downloaded / total * 100
                        log.info("    %.0f / %.0f MB (%.0f%%)",
                                 downloaded / 1024 / 1024, total / 1024 / 1024, pct)
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code} downloading {url}") from exc


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


# ── GTF parsing ───────────────────────────────────────────────────────────────

@dataclass
class TranscriptRecord:
    """One protein-coding transcript, populated incrementally from the GTF."""
    ensembl_transcript_id: str          # without version
    ensembl_protein_id: str = ""        # filled from protein FASTA pairing
    ensembl_gene_id: str = ""
    gene_symbol: str = ""
    transcript_name: str = ""
    transcript_type: str = ""
    chromosome: str = ""
    strand: str = ""
    cds_starts: List[int] = field(default_factory=list)
    cds_ends:   List[int] = field(default_factory=list)
    exon_count: int = 0
    is_mane_select: bool = False
    is_mane_plus_clinical: bool = False
    # Empty string if not principal; else the APPRIS tag verbatim
    # ("principal1" .. "principal5").  Phase 3 parses the integer suffix as a
    # tie-breaker rank — principal1 outranks principal5.
    appris_principal: str = ""
    # Filled later from the protein FASTA
    protein_sequence: str = ""

    @property
    def cds_start_min(self) -> Optional[int]:
        return min(self.cds_starts) if self.cds_starts else None

    @property
    def cds_end_max(self) -> Optional[int]:
        return max(self.cds_ends) if self.cds_ends else None

    @property
    def has_non_canonical_aa(self) -> bool:
        return any(c in NON_CANONICAL_AA for c in self.protein_sequence)


def _parse_attrs(attr_str: str) -> Dict[str, str]:
    """GTF attribute parser: `key "value"; key "value"; ...` → dict."""
    out: Dict[str, str] = {}
    tags: List[str] = []
    for part in attr_str.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, val = part.split(" ", 1)
        val = val.strip().strip('"')
        if key == "tag":
            tags.append(val)
        else:
            out[key] = val
    out["_tags"] = ",".join(tags)
    return out


def _strip_version(ensembl_id: str) -> str:
    """ENST00000280772.5 → ENST00000280772.  Idempotent on already-stripped IDs."""
    return ensembl_id.split(".", 1)[0] if ensembl_id else ensembl_id


def parse_gtf(gtf_path: Path) -> Dict[str, TranscriptRecord]:
    """
    Stream the GENCODE GTF and emit a TranscriptRecord per protein-coding transcript.
    Only `transcript`, `exon`, and `CDS` lines are inspected (skipping `gene`, `UTR`, etc.).
    """
    records: Dict[str, TranscriptRecord] = {}

    log.info("parsing GTF: %s", gtf_path.name)
    n_lines = 0
    n_transcripts_seen = 0

    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            n_lines += 1
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in {"transcript", "exon", "CDS"}:
                continue

            attrs = _parse_attrs(parts[8])
            biotype = attrs.get("transcript_type") or attrs.get("transcript_biotype")
            if biotype != "protein_coding":
                continue

            tid = _strip_version(attrs.get("transcript_id", ""))
            if not tid:
                continue

            rec = records.get(tid)
            if rec is None:
                rec = TranscriptRecord(ensembl_transcript_id=tid)
                records[tid] = rec
                n_transcripts_seen += 1

            if feature == "transcript":
                rec.ensembl_gene_id   = _strip_version(attrs.get("gene_id", ""))
                rec.gene_symbol       = attrs.get("gene_name", "")
                rec.transcript_name   = attrs.get("transcript_name", "")
                rec.transcript_type   = biotype
                rec.chromosome        = parts[0]
                rec.strand            = parts[6]
                tags = attrs.get("_tags", "")
                rec.is_mane_select        = "MANE_Select" in tags
                rec.is_mane_plus_clinical = "MANE_Plus_Clinical" in tags
                # Capture the principal level itself ("principal1" .. "5"),
                # not just the boolean — Phase 3 uses the number as a rank.
                rec.appris_principal      = next(
                    (t for t in tags.split(",") if t.startswith("appris_principal")),
                    "",
                )

            elif feature == "exon":
                rec.exon_count += 1

            elif feature == "CDS":
                try:
                    rec.cds_starts.append(int(parts[3]))
                    rec.cds_ends.append(int(parts[4]))
                except ValueError:
                    pass

    log.info("GTF scan: %d lines processed, %d protein-coding transcripts captured",
             n_lines, n_transcripts_seen)
    return records


# ── Protein FASTA parsing ─────────────────────────────────────────────────────

def parse_protein_fasta(
    fa_path: Path,
    records: Dict[str, TranscriptRecord],
) -> None:
    """
    Read the GENCODE protein-translations FASTA and attach the sequence + ENSP
    to each TranscriptRecord.  Header format:

        >ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.7|...|OR4F5-201|OR4F5|305
    """
    log.info("parsing protein FASTA: %s", fa_path.name)
    cur_tid: Optional[str] = None
    cur_pid: str = ""
    cur_seq: List[str] = []
    n_matched = 0
    n_skipped_no_record = 0

    def _flush() -> None:
        nonlocal n_matched, n_skipped_no_record
        if cur_tid is None:
            return
        rec = records.get(cur_tid)
        if rec is None:
            n_skipped_no_record += 1
            return
        rec.ensembl_protein_id = cur_pid
        rec.protein_sequence   = "".join(cur_seq)
        n_matched += 1

    with gzip.open(fa_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                cur_seq = []
                # Header tokens are pipe-delimited.  Layout:
                #   ENSP|ENST|ENSG|OTT_gene|OTT_transcript|transcript_name|gene_symbol|length
                tokens = line[1:].split("|")
                cur_pid = _strip_version(tokens[0]) if len(tokens) >= 1 else ""
                cur_tid = _strip_version(tokens[1]) if len(tokens) >= 2 else None
            else:
                cur_seq.append(line.strip())
    _flush()

    log.info("protein FASTA: matched %d transcripts; %d sequences had no GTF record",
             n_matched, n_skipped_no_record)


# ── MANE Select cross-check ───────────────────────────────────────────────────

def load_mane_set(mane_path: Path) -> set[str]:
    """Return the set of MANE-Select ENST IDs (without version).

    MANE summary file is tab-delimited; ENST is in the `Ensembl_nuc` column.
    The header line is the first line and starts with '#NCBI_GeneID'.
    """
    mane_enst: set[str] = set()
    with gzip.open(mane_path, "rt") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        try:
            ensembl_nuc_col = header.index("Ensembl_nuc")
            mane_status_col = header.index("MANE_status")
        except ValueError:
            log.warning("MANE summary header missing expected columns; got %s", header)
            return mane_enst
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(ensembl_nuc_col, mane_status_col):
                continue
            if row[mane_status_col] != "MANE Select":
                continue
            mane_enst.add(_strip_version(row[ensembl_nuc_col]))
    log.info("MANE Select: %d transcripts in MANE v1.4", len(mane_enst))
    return mane_enst


# ── Output ────────────────────────────────────────────────────────────────────

TSV_COLUMNS = [
    "ensembl_transcript_id", "ensembl_protein_id", "ensembl_gene_id",
    "gene_symbol", "transcript_name", "transcript_type",
    "protein_length", "protein_md5",
    "chromosome", "strand", "exon_count",
    "cds_start_min", "cds_end_max",
    "is_mane_select", "is_mane_plus_clinical", "appris_principal",
    "has_non_canonical_aa", "protein_sequence",
]


def write_tsv(out_path: Path, records: Iterable[TranscriptRecord]) -> int:
    n = 0
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        w.writerow(TSV_COLUMNS)
        for rec in records:
            if not rec.protein_sequence:
                # No protein sequence in the FASTA → skip the row.
                continue
            md5 = hashlib.md5(rec.protein_sequence.encode()).hexdigest()
            w.writerow([
                rec.ensembl_transcript_id,
                rec.ensembl_protein_id,
                rec.ensembl_gene_id,
                rec.gene_symbol,
                rec.transcript_name,
                rec.transcript_type,
                len(rec.protein_sequence),
                md5,
                rec.chromosome,
                rec.strand,
                rec.exon_count,
                rec.cds_start_min if rec.cds_start_min is not None else "",
                rec.cds_end_max   if rec.cds_end_max   is not None else "",
                "1" if rec.is_mane_select        else "0",
                "1" if rec.is_mane_plus_clinical else "0",
                rec.appris_principal,
                "1" if rec.has_non_canonical_aa  else "0",
                rec.protein_sequence,
            ])
            n += 1
    return n


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch GENCODE v47 reference and emit normalised TSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help=f"Output directory (default: {DEFAULT_OUT_DIR})")
    parser.add_argument("--force", action="store_true",
                        help="Re-download source files even if already cached")
    args = parser.parse_args()

    logfile = setup_logging("fetch_gencode_reference")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    log.info("output directory: %s", args.out_dir)

    # ── Download ─────────────────────────────────────────────────────────────
    log.info("─" * 60)
    log.info("downloading GENCODE %s sources", GENCODE_VERSION)
    log.info("─" * 60)
    for url, name in SOURCES:
        _download(url, args.out_dir / name, force=args.force)

    # ── Provenance ───────────────────────────────────────────────────────────
    src_urls = args.out_dir / "source_urls.txt"
    with open(src_urls, "w") as fh:
        from datetime import datetime, timezone
        fh.write(f"# Fetched at: {datetime.now(timezone.utc).isoformat()}Z\n")
        fh.write(f"# GENCODE release: {GENCODE_VERSION}\n")
        for url, name in SOURCES:
            fh.write(f"{name}\t{url}\n")
    log.info("provenance written: %s", src_urls)

    # ── Checksums (local) ────────────────────────────────────────────────────
    checks = args.out_dir / "CHECKSUMS"
    with open(checks, "w") as fh:
        for _, name in SOURCES:
            p = args.out_dir / name
            if p.exists():
                fh.write(f"{_md5_file(p)}  {name}\n")
    log.info("local checksums: %s", checks)

    # ── Parse ────────────────────────────────────────────────────────────────
    log.info("─" * 60)
    gtf_path = args.out_dir / "gencode.v47.annotation.gtf.gz"
    fa_path  = args.out_dir / "gencode.v47.pc_translations.fa.gz"
    mane_path = args.out_dir / "MANE.GRCh38.v1.4.summary.txt.gz"

    records = parse_gtf(gtf_path)
    parse_protein_fasta(fa_path, records)

    # Cross-check MANE Select tag against the official MANE summary; the GTF
    # tag should already be accurate but the cross-check guards against
    # version drift.
    mane_set = load_mane_set(mane_path)
    upgrades = 0
    for rec in records.values():
        if rec.ensembl_transcript_id in mane_set and not rec.is_mane_select:
            rec.is_mane_select = True
            upgrades += 1
    if upgrades:
        log.warning("MANE cross-check: upgraded %d transcripts whose GTF tag was missing", upgrades)

    # ── Write parsed TSV ────────────────────────────────────────────────────
    tsv_path = args.out_dir / "gencode_v47.tsv"
    n_rows = write_tsv(tsv_path, records.values())

    # ── Summary ─────────────────────────────────────────────────────────────
    log.info("─" * 60)
    n_mane    = sum(1 for r in records.values() if r.is_mane_select and r.protein_sequence)
    n_appris  = sum(1 for r in records.values() if r.appris_principal and r.protein_sequence)
    n_noncanon = sum(1 for r in records.values()
                     if r.protein_sequence and r.has_non_canonical_aa)
    # Distribution of protein-md5 (unique vs total)
    md5_set: set[str] = set()
    for r in records.values():
        if r.protein_sequence:
            md5_set.add(hashlib.md5(r.protein_sequence.encode()).hexdigest())

    log.info("Output: %s", tsv_path)
    log.info("  rows                            %s", f"{n_rows:>7,}")
    log.info("  unique protein_md5 values        %s", f"{len(md5_set):>7,}")
    log.info("  MANE Select                      %s", f"{n_mane:>7,}")
    log.info("  APPRIS principal                 %s", f"{n_appris:>7,}")
    log.info("  with non-canonical AAs (X/B/Z/U) %s", f"{n_noncanon:>7,}")
    log.info("─" * 60)
    log.info("detailed log: %s", logfile)


if __name__ == "__main__":
    main()
