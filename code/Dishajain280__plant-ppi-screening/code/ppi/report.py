"""Render screening results as a markdown table and an SVG figure.

No plotting library: the figure is written as plain SVG so it renders in a
GitHub README and in any browser without installing anything.
"""

from __future__ import annotations

from pathlib import Path

WIDTH, HEIGHT = 720, 320
PAD_L, PAD_R, PAD_T, PAD_B = 56, 20, 24, 44


def markdown_table(rows: list[dict], decay: int = 3) -> str:
    """Rows sorted by score, positive/negative marked, ready for a README."""
    ordered = sorted(rows, key=lambda r: -r["score"])
    lines = [
        "| Rank | Pair | Label | ipTM | min chain-pair PAE (A) |",
        "| ---: | --- | --- | ---: | ---: |",
    ]
    for rank, row in enumerate(ordered, start=1):
        label = "positive" if row["label"] == 1 else "negative"
        pae = row.get("min_chain_pair_pae")
        pae_text = f"{pae:.2f}" if isinstance(pae, (int, float)) and pae == pae else "n/a"
        lines.append(
            f"| {rank} | {row['pair_id']} | {label} | {row['score']:.{decay}f} | {pae_text} |"
        )
    return "\n".join(lines)


def svg_chart(rows: list[dict], threshold: float | None = None,
              title: str = "Ranked interaction scores") -> str:
    """Ranked-score dot plot: filled dots positive, hollow dots negative."""
    ordered = sorted(rows, key=lambda r: -r["score"])
    if not ordered:
        return '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"/>'

    scores = [r["score"] for r in ordered]
    lo, hi = min(scores), max(scores)
    if hi - lo < 1e-9:
        hi = lo + 1.0
    plot_w = WIDTH - PAD_L - PAD_R
    plot_h = HEIGHT - PAD_T - PAD_B

    def x(i: int) -> float:
        return PAD_L + (plot_w * (i + 0.5) / len(ordered))

    def y(value: float) -> float:
        return PAD_T + plot_h * (1 - (value - lo) / (hi - lo))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="#ffffff"/>',
        f'<text x="{PAD_L}" y="16" font-size="13" font-weight="bold" fill="#1f2328">{title}</text>',
        f'<line x1="{PAD_L}" y1="{PAD_T + plot_h}" x2="{WIDTH - PAD_R}" y2="{PAD_T + plot_h}" '
        f'stroke="#d0d7de"/>',
        f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T + plot_h}" stroke="#d0d7de"/>',
    ]

    for tick in (0.0, 0.5, 1.0):
        value = lo + (hi - lo) * tick
        yy = y(value)
        parts.append(
            f'<line x1="{PAD_L}" y1="{yy:.1f}" x2="{WIDTH - PAD_R}" y2="{yy:.1f}" '
            f'stroke="#f0f3f6"/>'
        )
        parts.append(
            f'<text x="{PAD_L - 8}" y="{yy + 4:.1f}" font-size="10" fill="#59636e" '
            f'text-anchor="end">{value:.2f}</text>'
        )

    if threshold is not None and lo <= threshold <= hi:
        ty = y(threshold)
        parts.append(
            f'<line x1="{PAD_L}" y1="{ty:.1f}" x2="{WIDTH - PAD_R}" y2="{ty:.1f}" '
            f'stroke="#b54708" stroke-dasharray="5 4"/>'
        )
        parts.append(
            f'<text x="{WIDTH - PAD_R}" y="{ty - 5:.1f}" font-size="10" fill="#b54708" '
            f'text-anchor="end">threshold {threshold:.2f}</text>'
        )

    for i, row in enumerate(ordered):
        cx, cy = x(i), y(row["score"])
        if row["label"] == 1:
            parts.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.5" fill="#1a7f37"/>')
        else:
            parts.append(
                f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4" fill="none" stroke="#8250df" '
                f'stroke-width="1.6"/>'
            )

    parts.append(
        f'<text x="{PAD_L + plot_w / 2}" y="{HEIGHT - 12}" font-size="10" fill="#59636e" '
        f'text-anchor="middle">pairs, ranked by ipTM (green filled = known interactor, '
        f'purple hollow = negative)</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def write_report(rows: list[dict], metrics: dict, out_dir: Path,
                 threshold: float | None = None, note: str | None = None) -> tuple[Path, Path]:
    """Write report.md and figure.svg; returns both paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    md = ["# Screening results", ""]
    if note:
        md += [f"> {note}", ""]
    md += ["## Metrics", ""]
    for name, metric in metrics.items():
        if name in {"top_k", "selected_at_threshold"}:
            continue
        md.append(f"- **{name}**: {metric}")
    md += ["", "## Ranked pairs", "", markdown_table(rows), ""]

    report_path = out_dir / "report.md"
    report_path.write_text("\n".join(md), encoding="utf-8")

    figure_path = out_dir / "figure.svg"
    figure_path.write_text(svg_chart(rows, threshold), encoding="utf-8")
    return report_path, figure_path
