#!/usr/bin/env python3
"""
build_source_metadata.py — generate unified per-source metadata.tsv files.

Produces one `metadata.tsv` per data source under `data/sources/<source>/`, all
sharing the same 31-column schema (see docs/unified_metadata_schema.md). The
universal key is `protein_md5`; cross-source dedup happens downstream at ingest.

Grain: one row per *native transcript ID* (design principle #3 of the schema) —
  - gencode      : one row per Ensembl transcript (ENST); GENCODE transcript tags
                   (MANE_Select, readthrough_transcript, NMD_exception, cds_end_NF …)
                   are carried in `gencode_tags`
  - ccrcc        : one row per (PacBio isoform, sample). PB.* IDs are numbered
                   independently per SQANTI sample run, so the SAME id can name
                   different transcripts in different samples (review EXP-1). Each
                   ORF's origin sample is known from the author pep file it came
                   from (`isoform_metadata.tsv:source_file`), and tissue flags,
                   SQANTI class and exon coords are taken from THAT sample only.
                   `sample_source` records it. A pbid that occurs in more than one
                   sample with a different protein gets `<pbid>@<sample>` as its
                   transcript_id (orig_transcript_id stays the bare pbid).
  - breast_cancer: one row per PacBio isoform (pb_id); the same protein md5 is
                   split into tissue-specific pb_ids upstream, so tissue is taken
                   from each row directly (`sample_source` = tissue)

Usage:
    python scripts/build_source_metadata.py --source gencode
    python scripts/build_source_metadata.py --source ccrcc
    python scripts/build_source_metadata.py --source bc
    python scripts/build_source_metadata.py --source all

Outputs (TSV): data/sources/<source>/metadata.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

# csv field-size limit — sequence/exons columns can be long
csv.field_size_limit(10_000_000)

REPO_ROOT = Path(__file__).resolve().parent.parent          # …/splisoforms
SOURCES_DIR = REPO_ROOT / "data" / "sources"

# ── Unified schema (31 columns; 30–31 are optional, nullable extensions) ─────
COLUMNS = [
    "protein_md5",          # 1  dedup key
    "source",               # 2  gencode_v47 | ccrcc | breast_cancer
    "transcript_id",        # 3  native ID for this source
    "orig_transcript_id",   # 4  upstream ID (ENST / ESPRESSO / PB.*)
    "gene_symbol",          # 5
    "gene_id",              # 6
    "seq_len",              # 7  protein length (aa)
    "sequence",             # 8  protein AA sequence
    "chrom",                # 9
    "strand",               # 10
    "genomic_start",        # 11 transcript span min
    "genomic_end",          # 12 transcript span max
    "exon_count",           # 13
    "biotype",              # 14
    "structural_category",  # 15 FSM/ISM/NIC/NNC or 'reference'
    "subcategory",          # 16 SQANTI3 subcategory (long-read); "" otherwise
    "is_canonical",         # 17
    "is_mane_select",       # 18 GENCODE MANE Select flag
    "appris_tag",           # 19 GENCODE/APPRIS principal tag
    "fate",                 # 20 reference | af3 | afdb | already_in_db
    "structure_provenance", # 21 af3_local | afdb | none
    "has_structure",        # 22 computed
    "exons",                # 23 JSON [[start,end],…] transcript order
    "in_tumor",             # 24
    "in_normal",            # 25
    "in_metastasis",        # 26
    "in_fibroblast",        # 27 BC cultured-fibroblast control
    "best_hit_acc",         # 28 closest UniProt accession
    "best_hit_identity",    # 29 % identity to best UniProt hit
    "sample_source",        # 30 origin sample / tissue of this transcript ("" if n/a)
    "gencode_tags",         # 31 comma-joined GENCODE transcript tags ("" if n/a)
]

# GENCODE transcript tags worth carrying (ingest priority, NMD, QC)
GENCODE_TAGS_KEPT = {
    "MANE_Select", "MANE_Plus_Clinical", "Ensembl_canonical", "basic", "CCDS",
    "readthrough_transcript", "NMD_exception", "NMD_likely_if_extended",
    "cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF",
}


# ── Shared helpers ─────────────────────────────────────────────────────────────
def _open(path: Path):
    """Open a plain or gzipped text file transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def structure_exists(source: str, md5: str) -> bool:
    """True iff data/sources/<source>/structures/<md5>/<md5>_model.cif exists."""
    if not md5:
        return False
    return (SOURCES_DIR / source / "structures" / md5 / f"{md5}_model.cif").exists()


def _gtf_attrs(attr_field: str) -> Dict[str, str]:
    """Parse a GTF/GFF column-9 attribute string into a dict."""
    attrs: Dict[str, str] = {}
    for chunk in attr_field.rstrip(";").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        # GENCODE/Ensembl: key "value";   PacBio GFF: key "value";
        if " " in chunk:
            key, _, val = chunk.partition(" ")
            attrs[key] = val.strip().strip('"')
    return attrs


def parse_exons_by_transcript(
    gtf_path: Path,
    key_attr: str = "transcript_id",
    strip_version: bool = False,
    wanted: Optional[set] = None,
) -> Dict[str, Dict]:
    """
    Parse exon features from a GTF/GFF into {tx_id: {chrom, strand, exons}}.

    `exons` is a list of [start, end] in transcript order (ascending for +,
    descending for −). genomic_start/end and exon_count are derived by callers.
    If `wanted` is given, only those tx_ids are retained (saves memory on the
    large PacBio GFFs). `strip_version` drops the ".N" suffix from the key.
    """
    raw: Dict[str, Dict] = {}
    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attrs = _gtf_attrs(f[8])
            txid = attrs.get(key_attr)
            if txid is None:
                continue
            if strip_version:
                txid = txid.split(".")[0]
            if wanted is not None and txid not in wanted:
                continue
            rec = raw.get(txid)
            if rec is None:
                rec = {"chrom": f[0], "strand": f[6], "exons": []}
                raw[txid] = rec
            rec["exons"].append([int(f[3]), int(f[4])])

    # finalise: sort exons in transcript order
    for rec in raw.values():
        rec["exons"].sort(key=lambda e: e[0], reverse=(rec["strand"] == "-"))
    return raw


def parse_transcript_tags(gtf_path: Path, strip_version: bool = True,
                          kept: Optional[set] = None) -> Dict[str, List[str]]:
    """{tx_id: sorted list of GENCODE `tag` values} from the `transcript` features
    of a GENCODE GTF. Only tags in `kept` (default GENCODE_TAGS_KEPT) are retained."""
    kept = GENCODE_TAGS_KEPT if kept is None else kept
    tag_re = re.compile(r'tag "([^"]+)"')
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    out: Dict[str, List[str]] = {}
    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "transcript":
                continue
            m = tid_re.search(f[8])
            if not m:
                continue
            txid = m.group(1).split(".")[0] if strip_version else m.group(1)
            tags = sorted({t for t in tag_re.findall(f[8]) if t in kept})
            if tags:
                out[txid] = tags
    return out


def _coords_from_exons(rec: Optional[Dict]) -> Dict:
    """Derive genomic_start/end, exon_count, exons-json from a parsed exon rec."""
    if not rec or not rec["exons"]:
        return {"chrom": "", "strand": "", "genomic_start": "",
                "genomic_end": "", "exon_count": "", "exons": ""}
    starts = [e[0] for e in rec["exons"]]
    ends = [e[1] for e in rec["exons"]]
    return {
        "chrom": rec["chrom"],
        "strand": rec["strand"],
        "genomic_start": min(starts),
        "genomic_end": max(ends),
        "exon_count": len(rec["exons"]),
        "exons": json.dumps(rec["exons"], separators=(",", ":")),
    }


def parse_fasta(path: Path, key: str = "first_token") -> Dict[str, str]:
    """Parse a FASTA into {id: sequence}. key='first_token' uses the header up to
    the first whitespace or pipe."""
    seqs: Dict[str, str] = {}
    cur_id: Optional[str] = None
    buf: List[str] = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf)
                header = line[1:].strip()
                # split on whitespace or ' | '
                cur_id = header.replace(" |", " ").split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id is not None:
            seqs[cur_id] = "".join(buf)
    return seqs


def _blank_row(source: str) -> Dict[str, str]:
    row = {c: "" for c in COLUMNS}
    row["source"] = source
    return row


def write_tsv(rows: List[Dict], out_path: Path) -> None:
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def summarise(rows: List[Dict], source: str) -> None:
    n = len(rows)
    print(f"\n── {source}: {n} rows ──")
    if not n:
        return
    # null rate for key columns
    for col in ("protein_md5", "sequence", "chrom", "exons",
                "structural_category", "fate"):
        filled = sum(1 for r in rows if r.get(col) not in ("", None))
        print(f"  {col:22s} filled {filled:>7}/{n}  ({100*filled//n}%)")
    with_struct = sum(1 for r in rows if r.get("has_structure") == "True")
    print(f"  has_structure=True     {with_struct:>7}/{n}")
    uniq_md5 = len({r['protein_md5'] for r in rows if r['protein_md5']})
    print(f"  unique protein_md5     {uniq_md5:>7}")


# ── GENCODE builder ─────────────────────────────────────────────────────────────
def build_gencode() -> List[Dict]:
    src = "gencode_v47"
    base = SOURCES_DIR / src
    tsv = base / "gencode_v47.tsv"
    gtf = base / "gencode.v47.annotation.gtf.gz"

    print(f"[gencode] parsing exons from {gtf.name} …", file=sys.stderr)
    exons = parse_exons_by_transcript(gtf, key_attr="transcript_id",
                                      strip_version=True)
    print(f"[gencode] {len(exons)} transcripts with exon coords", file=sys.stderr)
    print(f"[gencode] parsing transcript tags from {gtf.name} …", file=sys.stderr)
    tags = parse_transcript_tags(gtf, strip_version=True)
    print(f"[gencode] {len(tags)} transcripts with kept tags", file=sys.stderr)

    rows: List[Dict] = []
    with tsv.open() as fh:
        for m in csv.DictReader(fh, delimiter="\t"):
            enst = m["ensembl_transcript_id"]
            md5 = m["protein_md5"]
            enst_bare = enst.split(".")[0]
            tx_tags = tags.get(enst_bare, [])
            # MANE from the TSV flag OR the GTF tag (belt and braces)
            is_mane = (m.get("is_mane_select") == "1") or ("MANE_Select" in tx_tags)
            is_canon = is_mane or bool(m.get("appris_principal"))
            row = _blank_row(src)
            row.update({
                "protein_md5": md5,
                "transcript_id": enst,
                "orig_transcript_id": enst,
                "gene_symbol": m.get("gene_symbol", ""),
                "gene_id": m.get("ensembl_gene_id", ""),
                "seq_len": m.get("protein_length", ""),
                "sequence": m.get("protein_sequence", ""),
                "biotype": m.get("transcript_type", ""),
                "structural_category": "reference",
                "subcategory": "",
                "is_canonical": str(is_canon),
                "is_mane_select": str(is_mane),
                "appris_tag": m.get("appris_principal", "") or "",
                "fate": "reference",
                "structure_provenance": "af3_local",
                "has_structure": str(structure_exists(src, md5)),
                "sample_source": "",
                "gencode_tags": ",".join(tx_tags),
            })
            row.update(_coords_from_exons(exons.get(enst_bare)))
            rows.append(row)
    return rows


# ── ccRCC builder ───────────────────────────────────────────────────────────────
SQANTI_SAMPLE_TISSUE = {
    "Normal": "in_normal",
    "ccRCC2": "in_tumor", "ccRCC3": "in_tumor",
    "ccRCC4": "in_tumor", "ccRCC5": "in_tumor",
}
SQANTI_FILENAME = "scisoseq_transcripts.sorted.filtered_lite_classification.txt"

# The author-curated ORF files (author_pep/) are one per sample. Verified by
# matching every header's (pbid, genomic span) against the per-sample GFFs: file
# _1_ → Normal (14,969/14,969), _2_ → ccRCC2, _3_ → ccRCC3, _4_ → ccRCC4,
# _5_ → ccRCC5, with zero ambiguous matches. `isoform_metadata.tsv:source_file`
# names that pep file, which is how each row's origin sample is recovered.
CCRCC_PEP_FILE_SAMPLE = {"1": "Normal", "2": "ccRCC2", "3": "ccRCC3",
                         "4": "ccRCC4", "5": "ccRCC5"}
_CCRCC_PEP_IDX = re.compile(r"^p28443_o31598_([1-5])_")


def ccrcc_sample_from_source_file(source_file: str) -> str:
    """Origin sample of a ccRCC ORF from its author pep file name ("" if unknown)."""
    m = _CCRCC_PEP_IDX.match(source_file or "")
    return CCRCC_PEP_FILE_SAMPLE.get(m.group(1), "") if m else ""


def _load_ccrcc_sqanti(sqanti_root: Path) -> Dict[tuple, Dict]:
    """{(sample, pbid): {structural_category, subcategory}} — per sample, NOT merged:
    PB.* ids are sample-local, so the same id in two samples is two transcripts."""
    rec: Dict[tuple, Dict] = {}
    for sample_dir in sorted(p for p in sqanti_root.iterdir() if p.is_dir()):
        if sample_dir.name not in SQANTI_SAMPLE_TISSUE:
            print(f"[ccrcc] unmapped SQANTI dir {sample_dir.name} — skipping",
                  file=sys.stderr)
            continue
        path = sample_dir / SQANTI_FILENAME
        if not path.exists():
            print(f"[ccrcc] missing {path} — skipping", file=sys.stderr)
            continue
        with path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rec[(sample_dir.name, r["isoform"])] = {
                    "structural_category": r.get("structural_category") or "",
                    "subcategory": r.get("subcategory") or "",
                }
    return rec


def build_ccrcc() -> List[Dict]:
    src = "ccrcc"
    base = SOURCES_DIR / src
    meta = base / "isoform_metadata.tsv"

    # 1. sequences — the AUTHORITATIVE source is af3_targets.fa, keyed by
    #    protein_md5. It is the frozen, provenance-tracked set of author-curated
    #    "best ORFsequences forAF" proteins that were actually folded by AF3, so
    #    every structure (dir = <md5>) matches it exactly. Sourcing the sequence
    #    here guarantees md5(sequence) == protein_md5 == structure.
    #
    #    Do NOT use transdecoder/*.pep for the sequence: those are a *later*
    #    TransDecoder re-run that produced different ORF calls (e.g. PB.133394.56
    #    .p1 = 293 aa in the re-run vs the 132 aa protein that was folded). Reading
    #    them mismatched sequence↔structure for ~75% of ccRCC isoforms. The pep is
    #    still loaded as a fallback for the rare row whose md5 is absent from
    #    af3_targets (non-af3-folded; skipped at ingest when it has no structure).
    print("[ccrcc] loading af3_targets.fa (folded proteins, by md5) …", file=sys.stderr)
    af3_by_md5: Dict[str, str] = parse_fasta(base / "af3_targets.fa")  # key = first token = md5
    print(f"[ccrcc] {len(af3_by_md5)} folded proteins", file=sys.stderr)

    print("[ccrcc] loading transdecoder pep (fallback only) …", file=sys.stderr)
    seqs: Dict[str, str] = {}
    for pep in sorted((base / "transdecoder").glob("*.pep")):
        seqs.update(parse_fasta(pep))
    print(f"[ccrcc] {len(seqs)} pep sequences", file=sys.stderr)

    # 2. tissue + structural_category from SQANTI
    print("[ccrcc] loading SQANTI per-sample classification …", file=sys.stderr)
    sqanti = _load_ccrcc_sqanti(base / "sqanti3")
    print(f"[ccrcc] {len(sqanti)} pbids in SQANTI", file=sys.stderr)

    # 3. exon coords: PacBio GFF keyed by transcript_id == pbid, PER SAMPLE
    #    (<sample>.filtered_lite.gff). Collect wanted pbids first to bound memory.
    wanted = set()
    with meta.open() as fh:
        for m in csv.DictReader(fh, delimiter="\t"):
            wanted.add(m["pbid"])
    print(f"[ccrcc] parsing exon GFFs for {len(wanted)} pbids …", file=sys.stderr)
    exons_by_sample: Dict[str, Dict[str, Dict]] = {}
    for gff in sorted((base / "exon_gff").glob("*.gff")):
        sample = gff.name.split(".")[0]
        exons_by_sample[sample] = parse_exons_by_transcript(
            gff, key_attr="transcript_id", wanted=wanted)
        print(f"[ccrcc]   {sample}: {len(exons_by_sample[sample])} pbids with exons",
              file=sys.stderr)

    def _lookup_any_sample(pbid: str):
        """Fallback when the origin sample is unknown: first sample carrying the id."""
        for sample in sorted(exons_by_sample):
            if pbid in exons_by_sample[sample]:
                return sample
        return ""

    rows: List[Dict] = []
    seen = set()            # (pbid, sample)
    tid_taken = set()       # bare pbids already used as transcript_id
    n_suffixed = n_no_sample = 0
    with meta.open() as fh:
        for m in csv.DictReader(fh, delimiter="\t"):
            pbid = m["pbid"]
            sample = ccrcc_sample_from_source_file(m.get("source_file", ""))
            if not sample:
                sample = _lookup_any_sample(pbid)
                n_no_sample += 1
            if (pbid, sample) in seen:
                continue
            seen.add((pbid, sample))
            # transcript_id must be unique DB-wide; a pbid reused by another sample
            # (for a different protein) is suffixed with its sample.
            if pbid in tid_taken:
                tid = f"{pbid}@{sample}"
                n_suffixed += 1
            else:
                tid = pbid
                tid_taken.add(pbid)
            md5 = m["md5"]
            fate = m.get("fate", "")
            prov = {"af3": "af3_local", "afdb": "afdb"}.get(fate, "none")
            sq = sqanti.get((sample, pbid), {})
            row = _blank_row(src)
            row.update({
                "protein_md5": md5,
                "transcript_id": tid,
                "orig_transcript_id": pbid,
                "gene_symbol": m.get("gene_name", ""),
                "gene_id": m.get("gene_pbid", ""),
                "seq_len": m.get("seq_len", ""),
                # Folded protein (matches structure) first; pep only if md5 unknown.
                "sequence": af3_by_md5.get(md5) or seqs.get(m.get("protein_id", ""), ""),
                "biotype": "protein_coding",
                "structural_category": sq.get("structural_category", ""),
                "subcategory": sq.get("subcategory", ""),
                "is_canonical": "False",
                "is_mane_select": "False",
                "appris_tag": "",
                "fate": fate,
                "structure_provenance": prov,
                "has_structure": str(structure_exists(src, md5)),
                "in_tumor": str(sample.startswith("ccRCC")),
                "in_normal": str(sample == "Normal"),
                "in_metastasis": "False",
                "in_fibroblast": "False",
                "best_hit_acc": m.get("best_hit_acc", ""),
                "best_hit_identity": m.get("best_hit_identity", ""),
                "sample_source": sample,
                "gencode_tags": "",
            })
            row.update(_coords_from_exons(exons_by_sample.get(sample, {}).get(pbid)))
            rows.append(row)
    print(f"[ccrcc] {len(rows)} rows; {n_suffixed} pbids suffixed with @sample; "
          f"{n_no_sample} rows without a pep-file sample (fallback lookup)",
          file=sys.stderr)
    return rows


# ── Breast-cancer builder ────────────────────────────────────────────────────────
BC_TISSUE_TO_FLAG = {
    "breast": "in_normal",        # healthy breast control
    "tumor": "in_tumor",
    "fibroblast": "in_fibroblast",
}


def build_bc() -> List[Dict]:
    src = "breast_cancer"
    base = SOURCES_DIR / src
    meta = base / "isoform_metadata.tsv"

    # 1. sequences: af3_targets.fa keyed by md5 (af3-fated only)
    print("[bc] loading af3_targets.fa sequences …", file=sys.stderr)
    seqs = parse_fasta(base / "af3_targets.fa")
    print(f"[bc] {len(seqs)} af3 sequences (by md5)", file=sys.stderr)

    # 2. exon coords: 3 raw GTFs keyed by transcript_id == orig_transcript_id
    print("[bc] parsing raw GTF exons …", file=sys.stderr)
    exons: Dict[str, Dict] = {}
    for gtf in sorted((base / "raw_gtfs").glob("*.gtf")):
        exons.update(parse_exons_by_transcript(gtf, key_attr="transcript_id"))
    print(f"[bc] {len(exons)} transcripts with exon coords", file=sys.stderr)

    rows: List[Dict] = []
    with meta.open() as fh:
        for m in csv.DictReader(fh, delimiter="\t"):
            md5 = m["md5"]
            fate = m.get("fate", "")
            prov = "af3_local" if fate == "af3" else (
                "afdb" if fate == "afdb" else "none")
            tissue = m.get("tissue", "")
            flag = BC_TISSUE_TO_FLAG.get(tissue)
            structural = m.get("structural_class", "")
            # exon coords: orig_transcript_id may be comma-separated; first match
            orig = m.get("orig_transcript_id", "")
            ex_rec = None
            for cand in orig.split(","):
                ex_rec = exons.get(cand.strip())
                if ex_rec:
                    break
            row = _blank_row(src)
            row.update({
                "protein_md5": md5,
                "transcript_id": m.get("pb_id", ""),
                "orig_transcript_id": orig,
                "gene_symbol": "",
                "gene_id": m.get("gene_pb_id", ""),
                "seq_len": m.get("seq_len", ""),
                "sequence": seqs.get(md5, ""),
                "biotype": "protein_coding",
                "structural_category": structural,
                "subcategory": "",
                "is_canonical": str(structural == "FSM"),
                "is_mane_select": "False",
                "appris_tag": "",
                "fate": fate,
                "structure_provenance": prov,
                "has_structure": str(structure_exists(src, md5)),
                "in_tumor": str(flag == "in_tumor"),
                "in_normal": str(flag == "in_normal"),
                "in_metastasis": "False",
                "in_fibroblast": str(flag == "in_fibroblast"),
                "best_hit_acc": m.get("best_hit_acc", ""),
                "best_hit_identity": m.get("best_hit_identity", ""),
                "sample_source": tissue,
                "gencode_tags": "",
            })
            row.update(_coords_from_exons(ex_rec))
            rows.append(row)
    return rows


BUILDERS = {"gencode": build_gencode, "ccrcc": build_ccrcc, "bc": build_bc}
OUT_DIR = {"gencode": "gencode_v47", "ccrcc": "ccrcc", "bc": "breast_cancer"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True,
                    choices=["gencode", "ccrcc", "bc", "all"])
    ap.add_argument("--out", default=None,
                    help="override output path (single-source only)")
    args = ap.parse_args()

    targets = ["gencode", "ccrcc", "bc"] if args.source == "all" else [args.source]
    for src in targets:
        rows = BUILDERS[src]()
        out = (Path(args.out) if args.out
               else SOURCES_DIR / OUT_DIR[src] / "metadata.tsv")
        write_tsv(rows, out)
        summarise(rows, src)
        print(f"  → wrote {out}")


if __name__ == "__main__":
    main()
