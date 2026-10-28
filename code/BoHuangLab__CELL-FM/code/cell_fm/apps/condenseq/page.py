"""The Condensate Titration tab's panels, captions and figures, in the Virtual OpenCell tab's look.

The two tabs read as one site: the same flat OpenCell styling (viewer.scoped_css), the same layout
-- a navbar, the input and the plot on the left, a black image viewer with a sample strip and
contrast settings on the right -- and the same channel colours, the nucleus blue and the protein
grey. Generated images go through viewer.render, so both viewers draw
alike.

Nothing here touches the models: the app hands in a pipeline.Result and the reference nucleus.
"""

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from ..virtual_opencell import viewer
from . import metrics

SCOPE = "#condensate-titration"

# the plot palette: categorical slots 1 and 2 (validated, CVD dE 24.7 / normal dE 33.6)
CURVE = "#2a78d6"     # the titration curve and the area under it
GAP = "#eb6834"       # the area above the curve
INK = "#0b0b0b"
MUTED = "#52514e"
BACKGROUND = "#fff"   # plots sit on white, like the Virtual OpenCell map

NUM_THUMBNAILS = 64   # the sample strip: evenly spaced along the ladder; the slider reaches the rest

# OpenCell's channel switch. viewer.render calls the protein channel the "target"; here it is the IDP
CHANNELS = [("Nucleus", "Nucleus"), ("IDP", "Target"), ("Both", "Both")]

def _item(label: str, value: str) -> str:
    return (f'<div class="vo-metadata-item"><div class="vo-metadata-label">{label}</div>'
            f'<div class="vo-metadata-value">{value}</div></div>')


# ---------------------------------------------------------------------------
# A run's numbers and curve
# ---------------------------------------------------------------------------

def metrics_html(result=None) -> str:
    """AUC and AAC as OpenCell metadata; dashes before any run."""
    values = ["—", "—"] if result is None else [f"{result.auc:.4f}", f"{result.aac:.4f}"]
    before = "Run the analysis to measure them. " if result is None else ""
    return f"""
<div class="vo-target cs-metrics">
  <div class="vo-section-header">Titration</div>
  <div class="vo-section-caption">{before}AUC is how much of the scanned range the IDP spends
  condensed; AAC how much of that it loses to re-dissolving at high expression. Both are normalised
  to [0, 1].</div>
  <div class="vo-metadata">{''.join(_item(l, v) for l, v in zip(["AUC", "AAC"], values))}</div>
</div>
"""


CURVE_HTML = """
<div class="vo-target">
  <div class="vo-section-header">Titration curve</div>
  <div class="vo-section-caption">One grey dot per generated image, called condensed (1) or diffuse (0)
  by the ViT; the line is their moving average, the condensate probability. Blue is the area under
  it (AUC), orange the condensation lost to re-dissolving past the peak (AAC).</div>
</div>
"""


def _style(ax):
    ax.set_facecolor(BACKGROUND)
    ax.tick_params(labelsize=10, colors=MUTED)
    ax.grid(True, alpha=0.18, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
        ax.spines[side].set_linewidth(0.8)


def curve_figure(result, ladder) -> Figure:
    """The titration curve, with the two integrated areas drawn as the areas they are.

    ladder is the (lowest, highest) expression level; before any run the axes are drawn empty
    over it, so the panel has its final shape from the start.
    """
    fig = Figure(figsize=(7.4, 4.4), dpi=100, facecolor=BACKGROUND)
    ax = fig.add_subplot()
    _style(ax)

    if result is None:
        ax.text(0.5, 0.5, "Run the analysis to draw the curve", transform=ax.transAxes,
                ha="center", va="center", fontsize=11, color=MUTED)
    else:
        x, y, curve = result.intensities, result.predictions, result.curve
        flat = metrics.no_reentrant_curve(curve)

        # the raw per-image calls, kept visible but recessive: the curve is a summary of these
        ax.scatter(x, y, s=7, alpha=0.16, color=MUTED, linewidths=0, zorder=1)
        # AUC: what the sequence actually does
        ax.fill_between(x, 0, curve, color=CURVE, alpha=0.18, linewidth=0, zorder=2,
                        label=f"AUC {result.auc:.3f}")
        # AAC: the condensation lost to reentrant dissolution at high concentration
        if result.aac > 1e-9:
            ax.fill_between(x, curve, flat, where=flat > curve, color=GAP, alpha=0.30,
                            linewidth=0, zorder=3, label=f"AAC {result.aac:.3f}")
            ax.plot(x, flat, color=GAP, linewidth=1.2, linestyle="--", alpha=0.8, zorder=4)
        ax.plot(x, curve, color=CURVE, linewidth=2.0, zorder=5)

        if result.c_sat != float("inf"):
            ax.axvline(result.c_sat, color=MUTED, linewidth=1.0, linestyle=":", zorder=4)
            ax.annotate(f"$c_{{sat}}$ {result.c_sat:.0f}", xy=(result.c_sat, 1.0), xytext=(6, -13),
                        textcoords="offset points", fontsize=9, color=MUTED, ha="left", va="top")
        ax.legend(loc="upper left", frameon=False, fontsize=10, labelcolor=INK)

    ax.set_xscale("log")
    ax.set_xlim(*ladder)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Protein expression level (a.u., log scale)", fontsize=11, color=INK)
    ax.set_ylabel("Condensate probability", fontsize=11, color=INK)
    fig.tight_layout()
    return fig


def results_table(result) -> pd.DataFrame:
    """The numbers behind the plot, so the chart is never the only way to read them."""
    return pd.DataFrame({
        "protein_intensity_level": np.round(result.intensities, 4),
        "predicted_class": result.predictions,
        "condensate_probability_ma": np.round(result.curve, 6),
    })


# ---------------------------------------------------------------------------
# The image viewer and sample strip
# ---------------------------------------------------------------------------

def _sample(result, nucleus: np.ndarray, index: int) -> np.ndarray:
    """One image as viewer.render takes it: (2, H, W), nucleus first, on the uint16 scale its
    auto-contrast expects. Before any run the IDP channel is empty."""
    protein = np.zeros_like(nucleus) if result is None else result.images[index].astype(np.float32)
    return np.stack([nucleus, protein]) * 65535.0


def view(result, nucleus, index: int, channel: str = "Both",
         nucleus_settings=(0, 100, 1.0), target_settings=(0, 100, 1.0)):
    """Image `index` of the ladder as an RGB array; the reference nucleus alone before a run."""
    if nucleus is None:
        return None
    if result is not None:
        index = int(np.clip(index, 0, len(result.images) - 1))
    return viewer.render(_sample(result, nucleus, index), channel, nucleus_settings, target_settings)


def caption_html(result=None, index: int = 0) -> str:
    if result is None:
        return ('<div class="cs-viewer-caption">The reference nucleus every image is generated in. '
                'Run the analysis to add the IDP.</div>')
    index = int(np.clip(index, 0, len(result.images) - 1))
    call = "condensed" if result.predictions[index] else "diffuse"
    return (f'<div class="cs-viewer-caption">Image {index + 1} of {len(result.images)} &middot; '
            f'expression level <b>{result.image_intensities[index]:.0f}</b> a.u. &middot; '
            f'called <b>{call}</b> &middot; condensate probability {result.curve[index]:.2f}</div>')


def thumbnail_indices(n: int) -> np.ndarray:
    """Which images of an n-image ladder the strip shows."""
    return np.unique(np.linspace(0, n - 1, min(n, NUM_THUMBNAILS)).round().astype(int))


def thumbnails(result, nucleus) -> list:
    """(image, caption) per strip entry, both channels at default settings, captioned with the level."""
    return [(view(result, nucleus, i), f"{result.image_intensities[i]:.0f}")
            for i in thumbnail_indices(len(result.images))]


# ---------------------------------------------------------------------------
# The substitution scan
# ---------------------------------------------------------------------------

def scan_header_html(scan) -> str:
    return f"""
<div class="vo-target">
  <div class="vo-section-header">Amino acid substitution scan</div>
  <div class="vo-section-caption">Every residue type in {len(scan.idr_index)} held-out IDPs replaced by
  each other type, each mutant titrated as above, and its AUC and AAC compared with the wild type's.
  A cell is the mean change over the IDPs that contain the FROM residue: left the paper's Fig. 5d
  (AUC), right its Supplementary Fig. 5c (AAC). Click any cell for the distribution behind it, or a
  residue on either axis to compare whole rows (Mutate FROM) or columns (Mutate TO) below.</div>
</div>
"""


def substitution_html(scan, i: int, j: int) -> str:
    """The selected substitution, named large like a target, with its summary statistics."""
    source, target = scan.amino_acids[i], scan.amino_acids[j]
    items, n = [], 0
    for metric in ("auc", "aac"):
        values = scan.delta[metric][i, j]
        values = values[np.isfinite(values)]
        n = len(values)
        label = metric.upper()
        items += [_item(f"Mean Δ {label}", f"{values.mean():+.4f}"),
                  _item(f"Median Δ {label}", f"{np.median(values):+.4f}")]
    return f"""
<div class="vo-target">
  <div class="cs-substitution">{source} → {target}</div>
  <div class="vo-status">Every {source} in the IDP replaced by {target} &middot; {n} of
  {len(scan.idr_index)} IDPs contain {source}</div>
  <div class="vo-metadata cs-metrics">{''.join(items)}</div>
</div>
"""


COMPARE_HTML = """
<div class="vo-target">
  <div class="vo-section-header">Compare rows or columns</div>
  <div class="vo-section-caption">A row is everything one residue becomes (Mutate FROM), a column
  everything that becomes one residue (Mutate TO). Click residues on a matrix's axis to add or remove
  them, up to eight; clicking on the other axis starts a comparison there. Each point is a matrix
  cell, with its 95% confidence interval over the IDPs; hover for the numbers.</div>
</div>
"""


def per_idp_html(repo: str, filename: str) -> str:
    return f"""
<div class="vo-target">
  <div class="vo-section-header">Per-IDP values</div>
  <div class="vo-section-caption">One row per IDP behind the selected cell; click a header to sort.
  Precomputed: the whole scan is <a href="https://huggingface.co/{repo}/blob/main/{filename}"
  target="_blank">{filename}</a> in the model repo, checked against the matrices the paper's figures
  were drawn from (<code>ana_all_mutation_log_scale.py</code>, AAC being
  <code>compute_area_if_no_reentrant_minus_auc</code>).</div>
</div>
"""


# ---------------------------------------------------------------------------
# Static text
# ---------------------------------------------------------------------------

NAVBAR_HTML = """
<div class="vo-navbar-brand">
  <span class="vo-navbar-title">Condensate Titration</span>
  <span class="vo-navbar-menu">CELL-FM virtual CondenSeq experiments on 66-aa IDPs</span>
</div>
"""

HOW_HTML = """
<div class="vo-target">
  <div class="vo-section-header">How a titration is computed</div>
  <div class="vo-section-caption">CELL-FM generates one image per expression level along a log-spaced
  ladder, a ViT classifier calls each image condensed or diffuse, a moving average over the calls is
  the titration curve, and the curve is integrated against log expression into AUC and AAC.</div>
</div>
"""

FOOTER_HTML = """
<div class="vo-footer">
  Each image is CELL-FM's prediction from the IDP sequence and a protein expression level alone,
  drawn in one shared cell: the nucleus channel is the same real CondenSeq nucleus for every image,
  so differences come from the sequence and the level. Sampling is fixed at 512 levels from 20 to
  2048 a.u., 100 ODE steps and seed 6, as in the offline runs. The IDP library and the condensate
  phenotypes the models learned from are CondenSeq's
  (Kappel et al., <i>Nat. Methods</i> 2025).
</div>
"""

# Additions to the OpenCell look for what the Virtual OpenCell tab does not have
CSS = """
#condensate-titration .cs-sequence textarea {
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace !important;
  font-size: 15px !important; line-height: 1.55; background: #fff;
}
#condensate-titration .cs-run {
  background: var(--vo-blue) !important; color: #fff !important; border: none !important;
  border-radius: 5px !important; font-family: "Lato", sans-serif; font-size: 17px; letter-spacing: 0.5px;
}
#condensate-titration .cs-run:hover { background: #0189bd !important; }
#condensate-titration .cs-metrics .vo-metadata-value { font-family: "Lato", sans-serif; font-size: 26px; font-weight: 300; }
#condensate-titration .cs-viewer-caption { font-family: "Lato", sans-serif; font-size: 14px; color: #555; margin: 2px 0; }
#condensate-titration .cs-viewer-caption b { color: #182026; font-weight: 600; }
#condensate-titration .cs-substitution { font-family: "Lato", sans-serif; font-size: 44px; font-weight: 400; line-height: 1.1; }
#condensate-titration .cs-diagram img { max-width: 1100px; margin: 0 auto; }
"""
