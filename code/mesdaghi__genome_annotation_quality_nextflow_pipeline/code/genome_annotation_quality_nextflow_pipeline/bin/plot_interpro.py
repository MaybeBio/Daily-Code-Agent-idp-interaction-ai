#!/usr/bin/env python3
"""
Plot a query proteome's InterProScan coverage alongside precomputed
model organisms loaded from ipr_coverage.json.

This script is self-contained — it does not require ipr_coverage_reporter.py
to be importable. It includes its own InterProScan XML parser and FASTA
length reader, identical in behaviour to the reporter script.

Usage
-----
    python3 ipr_coverage_query.py QUERY.fasta QUERY.xml \
        [--json ipr_coverage.json] \
        [--name "My species"] \
        [--outdir .]

Produces three files in --outdir:
    query_ipr_coverage.png                       - per-library bar chart
    query_ipr_merged_coverage_distribution.png   - merged-coverage violin/box
    query_ipr_summary.tsv                        - tab-separated summary table

The query is appended at the end of every plot and rendered with a thicker
border and bold italic label so it stands out from the precomputed organisms.
"""
from __future__ import annotations

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ─── Shared constants ────────────────────────────────────────────────────────
LIBRARY_COLOURS = {
    "PFAM":             "#4E79A7",
    "GENE3D":           "#F28E2B",
    "SUPERFAMILY":      "#E15759",
    "PANTHER":          "#76B7B2",
    "PROSITE_PROFILES": "#59A14F",
    "PROSITE_PATTERNS": "#EDC948",
    "CDD":              "#B07AA1",
    "SMART":            "#FF9DA7",
    "PRINTS":           "#9C755F",
    "HAMAP":            "#BAB0AC",
    "PIRSF":            "#D4A6C8",
    "PIRSR":            "#86BCB6",
    "SFLD":             "#F1CE63",
    "NCBIFAM":          "#A0CBE8",
    "FUNFAM":           "#FFBE7D",
    "MOBIDB_LITE":      "#8CD17D",
    "COILS":            "#B6992D",
    "ANTIFAM":          "#D37295",
}
DEFAULT_COLOUR = "#AAAAAA"
TOTAL_COLOUR = "#2B2B2B"

# Libraries that carry an InterPro (IPR) accession — used for "with IPR" count
IPR_LIBRARIES = {
    "PFAM", "GENE3D", "SUPERFAMILY", "PANTHER", "PROSITE_PROFILES",
    "PROSITE_PATTERNS", "CDD", "SMART", "PRINTS", "HAMAP", "PIRSF",
    "PIRSR", "SFLD", "NCBIFAM", "FUNFAM",
}

SPECIES_PALETTE = ["#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                   "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"]


# ─── FASTA / XML parsing (mirrors ipr_coverage_reporter.py) ──────────────────
def parse_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """
    Return {protein_id: sequence_length} for every record in the FASTA.
    Uses the first whitespace-delimited token of the header as the ID.
    """
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_len = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_id is not None:
            lengths[current_id] = current_len
    return lengths


def parse_xml(xml_path: Path):
    """
    Returns
    -------
    xml_proteins   : int   — unique protein IDs appearing in the XML
    with_ipr       : int   — proteins with ≥1 IPR-bearing library hit
    lib_counts     : dict  — {library: number of proteins with ≥1 hit from it}
    protein_ipr    : dict  — {protein_id: {"length": int,
                                           "intervals": [(start, end), ...]}}
                             Only intervals from matches that carry an
                             <entry ac="IPR..."> element are included.
    """
    def _strip_ns(tag: str) -> str:
        return tag.split("}")[-1] if "}" in tag else tag

    protein_libs: dict[str, set] = defaultdict(set)
    protein_ipr: dict[str, dict] = {}

    current_ids: list[str] = []
    current_libs: set[str] = set()
    current_seq_len: int = 0
    current_intervals: list[tuple[int, int]] = []

    match_has_ipr: bool = False
    match_intervals: list[tuple[int, int]] = []
    inside_match: bool = False
    inside_protein: bool = False

    for event, elem in ET.iterparse(xml_path, events=("start", "end")):
        tag = _strip_ns(elem.tag)

        if event == "start":
            if tag == "protein":
                inside_protein = True
                current_ids = []
                current_libs = set()
                current_seq_len = 0
                current_intervals = []
            elif inside_protein and tag == "xref":
                pid = elem.get("id") or elem.get("name")
                if pid:
                    current_ids.append(pid)
            elif inside_protein and tag == "signature-library-release":
                lib = elem.get("library", "UNKNOWN")
                current_libs.add(lib)
            elif inside_protein and tag.endswith("-match"):
                inside_match = True
                match_has_ipr = False
                match_intervals = []
            elif inside_match and tag == "entry":
                ac = elem.get("ac", "")
                if ac.startswith("IPR"):
                    match_has_ipr = True
            elif inside_match and tag.endswith("-location") \
                    and not tag.endswith("-location-fragment"):
                s = elem.get("start")
                e = elem.get("end")
                if s is not None and e is not None:
                    try:
                        s_i, e_i = int(s), int(e)
                        if e_i >= s_i:
                            match_intervals.append((s_i, e_i))
                    except ValueError:
                        pass

        elif event == "end":
            if tag == "sequence" and inside_protein:
                seq = (elem.text or "").strip()
                current_seq_len = len(seq)
                elem.clear()
            elif inside_match and tag.endswith("-match"):
                if match_has_ipr:
                    current_intervals.extend(match_intervals)
                inside_match = False
                match_has_ipr = False
                match_intervals = []
            elif tag == "protein":
                for pid in current_ids:
                    protein_libs[pid] |= current_libs
                    protein_ipr[pid] = {
                        "length": current_seq_len,
                        "intervals": list(current_intervals),
                    }
                inside_protein = False
                elem.clear()

    xml_proteins = len(protein_libs)
    with_ipr = sum(1 for libs in protein_libs.values() if libs & IPR_LIBRARIES)

    lib_counts: dict[str, int] = defaultdict(int)
    for libs in protein_libs.values():
        for lib in libs:
            lib_counts[lib] += 1

    return xml_proteins, with_ipr, dict(lib_counts), protein_ipr


# ─── Coverage calculations ──────────────────────────────────────────────────
def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sort and merge overlapping/adjacent closed intervals."""
    if not intervals:
        return []
    s = sorted(intervals)
    merged = [s[0]]
    for start, end in s[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            if end > last_end:
                merged[-1] = (last_start, end)
        else:
            merged.append((start, end))
    return merged


def compute_coverages(
    protein_ipr: dict[str, dict],
    fasta_lengths: dict[str, int],
) -> tuple[list[float], list[float]]:
    """
    Return two parallel coverage vectors over every protein in the FASTA:

    merged_cov  — fraction of residues covered by the union of all IPR hits
                  (overlaps counted once).
    longest_cov — length of the single longest IPR hit divided by protein length.

    Proteins with zero IPR hits contribute 0.0 to both.
    """
    merged_cov: list[float] = []
    longest_cov: list[float] = []

    for pid, length in fasta_lengths.items():
        if length <= 0:
            continue
        info = protein_ipr.get(pid)
        intervals = info["intervals"] if info else []

        if not intervals:
            merged_cov.append(0.0)
            longest_cov.append(0.0)
            continue

        merged = merge_intervals(intervals)
        covered = sum(e - s + 1 for s, e in merged)
        merged_cov.append(min(covered / length, 1.0))

        longest = max(e - s + 1 for s, e in intervals)
        longest_cov.append(min(longest / length, 1.0))

    return merged_cov, longest_cov

QUERY_COLOUR = "#222222"          # dark, distinct from the species palette
QUERY_HIGHLIGHT = "#C8102E"       # accent for borders / label colour


# ── Loading ──────────────────────────────────────────────────────────────────
def load_models(json_path: Path) -> dict:
    """Load the precomputed IPR_COVERAGE dict produced by ipr_coverage_reporter."""
    with open(json_path) as fh:
        data = json.load(fh)
    if not data:
        raise ValueError(f"{json_path} is empty.")
    # Normalise legacy 'Genus species' keys to 'Genus_species' so labels and
    # cross-tool comparisons are consistent with the rest of the pipeline.
    return {k.replace(" ", "_"): v for k, v in data.items()}


def compute_query(fasta_path: Path, xml_path: Path) -> dict:
    """Parse the query and return a record matching the JSON schema."""
    fasta_lengths = parse_fasta_lengths(fasta_path)
    total = len(fasta_lengths)
    if total == 0:
        raise ValueError(f"No proteins found in {fasta_path}")

    xml_proteins, with_ipr, lib_counts, protein_ipr = parse_xml(xml_path)

    merged_cov, longest_cov = compute_coverages(protein_ipr, fasta_lengths)
    merged_arr = np.asarray(merged_cov, dtype=float)
    longest_arr = np.asarray(longest_cov, dtype=float)

    lib_pcts = {
        lib: round(n / total * 100, 3) for lib, n in lib_counts.items()
    }

    def _stats(a: np.ndarray) -> dict:
        if a.size == 0:
            return {"n": 0, "mean": 0.0, "median": 0.0, "std": 0.0}
        return {
            "n": int(a.size),
            "mean": round(float(a.mean()), 4),
            "median": round(float(np.median(a)), 4),
            "std": round(float(a.std(ddof=1)), 4) if a.size > 1 else 0.0,
        }

    return {
        "total":        total,
        "xml_proteins": xml_proteins,
        "with_ipr":     with_ipr,
        "no_hit":       max(total - xml_proteins, 0),
        "ipr_pct":      round(with_ipr / total * 100, 3) if total else 0.0,
        "any_pct":      round(xml_proteins / total * 100, 3) if total else 0.0,
        "lib_counts":   lib_counts,
        "lib_pcts":     lib_pcts,
        "merged_coverage_stats":  _stats(merged_arr),
        "longest_coverage_stats": _stats(longest_arr),
        "merged_coverage":  merged_cov,
        "longest_coverage": longest_cov,
    }


# ── Bar chart ────────────────────────────────────────────────────────────────
def make_bar_chart(records: list[tuple[str, dict, bool]], output_path: Path):
    """records: list of (species, record_dict, is_query)."""
    # Library order across all species (IPR-bearing libraries only)
    all_libs: dict[str, int] = {}
    for _, rec, _ in records:
        for lib, n in rec.get("lib_counts", {}).items():
            if lib in IPR_LIBRARIES:
                all_libs[lib] = all_libs.get(lib, 0) + n
    lib_order = sorted(all_libs, key=all_libs.get, reverse=True)

    species_labels = [s for s, _, _ in records]
    n_species = len(records)
    n_bars = len(lib_order) + 1  # +1 for Total

    cluster_width = 0.86
    bar_width = cluster_width / n_bars
    offsets = (np.arange(n_bars) - (n_bars - 1) / 2) * bar_width

    fig_w = max(10, n_species * n_bars * 0.28 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, 9))
    x = np.arange(n_species)
    handles = []

    # Per-library bars
    for j, lib in enumerate(lib_order):
        pcts = np.array(
            [rec.get("lib_pcts", {}).get(lib, 0.0) for _, rec, _ in records],
            dtype=float,
        )
        colour = LIBRARY_COLOURS.get(lib, DEFAULT_COLOUR)
        # Highlight query bars with a thicker red border
        for i, (_, _, is_q) in enumerate(records):
            edge = QUERY_HIGHLIGHT if is_q else "white"
            lw = 1.2 if is_q else 0.3
            ax.bar(x[i] + offsets[j], pcts[i], width=bar_width,
                   color=colour, edgecolor=edge, linewidth=lw)
        handles.append(mpatches.Patch(color=colour, label=lib))

    # Total bar (any IPR hit)
    totals_pct = np.array(
        [rec.get("ipr_pct", 0.0) for _, rec, _ in records],
        dtype=float,
    )
    total_xs = x + offsets[-1]
    for i, (_, _, is_q) in enumerate(records):
        edge = QUERY_HIGHLIGHT if is_q else "white"
        lw = 1.2 if is_q else 0.3
        ax.bar(total_xs[i], totals_pct[i], width=bar_width,
               color=TOTAL_COLOUR, edgecolor=edge, linewidth=lw)
    handles.append(mpatches.Patch(color=TOTAL_COLOUR, label="Total (any IPR hit)"))
    handles.append(
        mpatches.Patch(facecolor="white", edgecolor=QUERY_HIGHLIGHT, linewidth=1.2,
                       label="Query (highlighted)")
    )

    # Total-bar labels
    for xt, pct in zip(total_xs, totals_pct):
        ax.text(xt, pct + 1, f"{pct:.1f}%", ha="center", va="bottom",
                fontsize=8, fontweight="bold")

    # Per-species summary annotation above each cluster
    ymax = 118
    for i, (_, rec, _) in enumerate(records):
        label = (f"IPR {rec.get('with_ipr', 0):,} / {rec.get('total', 0):,}   "
                 f"any-hit {rec.get('xml_proteins', 0):,}")
        ax.text(i, ymax - 4, label, ha="center", va="top", fontsize=8,
                color="#555555")

    # Visual separator between models and query
    if any(is_q for _, _, is_q in records):
        first_q = next(i for i, (_, _, is_q) in enumerate(records) if is_q)
        if first_q > 0:
            ax.axvline(first_q - 0.5, color="#999999", linestyle="--",
                       linewidth=0.8, alpha=0.6)

    ax.set_xticks(x)
    xtick_labels = ax.set_xticklabels(species_labels, fontsize=10, style="italic",
                                       rotation=20, ha="right")
    for tick, (_, _, is_q) in zip(xtick_labels, records):
        if is_q:
            tick.set_fontweight("bold")
            tick.set_color(QUERY_HIGHLIGHT)

    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_ylabel("Percentage of input proteins (%)", fontsize=10)
    ax.set_title("InterProScan per-library coverage "
                 "(denominator = FASTA input)", fontsize=12, fontweight="bold")
    ax.tick_params(axis="x", labelsize=10)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    n_cols = min(len(handles), 7)
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(0.5, -0.10), ncol=n_cols,
              fontsize=8, frameon=False, title="Signature library",
              title_fontsize=9)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved → {output_path}")


# ── Distribution plot ────────────────────────────────────────────────────────
def make_distribution_plot(records: list[tuple[str, dict, bool]],
                           series_key: str,
                           output_path: Path,
                           title: str,
                           ylabel: str):
    """
    series_key: 'merged_coverage' or 'longest_coverage'
    Each record must contain the per-protein vector under series_key.
    """
    species = [s for s, _, _ in records]
    data = [np.asarray(rec.get(series_key, []), dtype=float)
            for _, rec, _ in records]
    is_query = [is_q for _, _, is_q in records]

    # Sanity check: missing vectors → skip with a clear error
    for s, d in zip(species, data):
        if d.size == 0:
            print(f"  WARNING: no '{series_key}' data for {s}; "
                  f"that species will appear as an empty slot.")

    n_species = len(species)
    positions = np.arange(1, n_species + 1)

    fig, ax = plt.subplots(figsize=(max(8, n_species * 1.4 + 2), 7))

    # Violins
    nonempty_idx = [i for i, d in enumerate(data) if d.size]
    nonempty_data = [data[i] for i in nonempty_idx]
    nonempty_pos = positions[nonempty_idx]

    if nonempty_data:
        parts = ax.violinplot(nonempty_data, positions=nonempty_pos,
                              widths=0.85, showmeans=False,
                              showmedians=False, showextrema=False)
        for body, idx in zip(parts["bodies"], nonempty_idx):
            colour = (QUERY_COLOUR if is_query[idx]
                      else SPECIES_PALETTE[idx % len(SPECIES_PALETTE)])
            body.set_facecolor(colour)
            body.set_edgecolor(QUERY_HIGHLIGHT if is_query[idx] else colour)
            body.set_alpha(0.45 if is_query[idx] else 0.35)
            body.set_linewidth(1.6 if is_query[idx] else 0.8)

        # Inner box plots
        bp = ax.boxplot(nonempty_data, positions=nonempty_pos,
                        widths=0.18, patch_artist=True, showfliers=False,
                        medianprops={"color": "white", "linewidth": 1.5},
                        whiskerprops={"color": "#222222", "linewidth": 1.0},
                        capprops={"color": "#222222", "linewidth": 1.0},
                        boxprops={"edgecolor": "#222222", "linewidth": 1.0})
        for box, idx in zip(bp["boxes"], nonempty_idx):
            colour = (QUERY_COLOUR if is_query[idx]
                      else SPECIES_PALETTE[idx % len(SPECIES_PALETTE)])
            box.set_facecolor(colour)
            box.set_alpha(0.95)
            if is_query[idx]:
                box.set_edgecolor(QUERY_HIGHLIGHT)
                box.set_linewidth(1.6)

        means = [float(d.mean()) for d in nonempty_data]
        ax.scatter(nonempty_pos, means, marker="D", s=28, color="white",
                   edgecolor="#222222", linewidth=0.8, zorder=5)

    # Visual separator between models and query
    if any(is_query):
        first_q = is_query.index(True)
        if first_q > 0:
            ax.axvline(first_q + 0.5, color="#999999", linestyle="--",
                       linewidth=0.8, alpha=0.6)

    # Stats annotation
    stats_lines = []
    for s, d, is_q in zip(species, data, is_query):
        n = d.size
        mean = float(d.mean()) if n else 0.0
        median = float(np.median(d)) if n else 0.0
        n_zero = int((d == 0).sum())
        zero_pct = (n_zero / n * 100) if n else 0.0
        prefix = "▶ " if is_q else "  "
        stats_lines.append(
            f"{prefix}{s}: n={n:,}  mean={mean:.3f}  median={median:.3f}  "
            f"zeros={n_zero:,} ({zero_pct:.1f}%)"
        )

    ax.set_xticks(positions)
    xtick_labels = ax.set_xticklabels(species, fontsize=10, style="italic",
                                       rotation=20, ha="right")
    for tick, is_q in zip(xtick_labels, is_query):
        if is_q:
            tick.set_fontweight("bold")
            tick.set_color(QUERY_HIGHLIGHT)

    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    legend_handles = [
        mpatches.Patch(facecolor="#888888", edgecolor="#888888", alpha=0.35,
                       label="Density (violin)"),
        mpatches.Patch(facecolor="#888888", edgecolor="#222222",
                       label="IQR (box) + median (white line)"),
        plt.Line2D([0], [0], marker="D", color="w", markerfacecolor="white",
                   markeredgecolor="#222222", markersize=7, label="Mean"),
        mpatches.Patch(facecolor=QUERY_COLOUR, edgecolor=QUERY_HIGHLIGHT,
                       linewidth=1.6, label="Query (highlighted)"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False,
              fontsize=8)

    ax.text(0.0, -0.18, "\n".join(stats_lines), transform=ax.transAxes,
            fontsize=8, family="monospace", va="top", ha="left",
            color="#333333")

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved → {output_path}")


# ── Summary table ────────────────────────────────────────────────────────────
def write_summary_table(records: list[tuple[str, dict, bool]], output_path: Path):
    """Tab-separated summary of headline numbers for every species."""
    cols = [
        "species", "is_query", "total", "xml_proteins", "with_ipr",
        "no_hit", "ipr_pct", "any_pct",
        "merged_mean", "merged_median",
        "longest_mean", "longest_median",
    ]
    rows = []
    for species, rec, is_q in records:
        ms = rec.get("merged_coverage_stats", {})
        ls = rec.get("longest_coverage_stats", {})
        rows.append([
            species,
            "yes" if is_q else "no",
            rec.get("total", 0),
            rec.get("xml_proteins", 0),
            rec.get("with_ipr", 0),
            rec.get("no_hit", 0),
            rec.get("ipr_pct", 0.0),
            rec.get("any_pct", 0.0),
            ms.get("mean", 0.0), ms.get("median", 0.0),
            ls.get("mean", 0.0), ls.get("median", 0.0),
        ])

    with open(output_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(c) for c in row) + "\n")

    # Pretty-print to stdout too
    print("\nSummary:")
    widths = [max(len(str(r[i])) for r in [cols] + rows) for i in range(len(cols))]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*cols))
    print("  ".join("-" * w for w in widths))
    for row in rows:
        print(fmt.format(*[str(c) for c in row]))
    print(f"\nTable saved → {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description="Plot a query proteome alongside precomputed model organisms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fasta", type=Path, help="Query FASTA file")
    p.add_argument("xml", type=Path, help="Query InterProScan XML file")
    p.add_argument("--json", type=Path, default=Path("ipr_coverage.json"),
                   help="Precomputed model-organism JSON")
    p.add_argument("--name", type=str, default=None,
                   help="Display name for the query "
                        "(default: derived from the XML filename)")
    p.add_argument("--outdir", type=Path, default=Path("."),
                   help="Output directory for plots and summary table")
    args = p.parse_args()

    for f in (args.fasta, args.xml, args.json):
        if not f.exists():
            sys.exit(f"ERROR: file not found: {f}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    query_name = args.name or args.xml.stem.replace("_", " ")

    print(f"Loading model organisms from {args.json} ...")
    models = load_models(args.json)
    print(f"  → {len(models)} species: {', '.join(models.keys())}")

    print(f"\nProcessing query: {query_name}")
    print(f"  FASTA: {args.fasta}")
    print(f"  XML:   {args.xml}")
    query_record = compute_query(args.fasta, args.xml)
    print(f"  → total={query_record['total']:,}  "
          f"with_ipr={query_record['with_ipr']:,}  "
          f"ipr_pct={query_record['ipr_pct']}%")

    # Compose record list: models first (alphabetical), query last
    records: list[tuple[str, dict, bool]] = [
        (sp, rec, False) for sp, rec in sorted(models.items())
    ]
    records.append((query_name, query_record, True))

    print()
    make_bar_chart(records, args.outdir / "query_ipr_coverage.png")
    make_distribution_plot(
        records, "merged_coverage",
        args.outdir / "query_ipr_merged_coverage_distribution.png",
        title="Distribution of IPR domain coverage per protein "
              "(merged, overlaps counted once)",
        ylabel="Fraction of protein residues covered by IPR domains",
    )
    write_summary_table(records, args.outdir / "query_ipr_summary.tsv")


if __name__ == "__main__":
    main()
