"""The amino-acid substitution scan: mean change in AUC (Fig. 5d) and AAC (Supplementary
Fig. 5c), and the per-IDP changes behind every cell of the two matrices.

The scan replaced every residue type by every other type in 100 held-out CondenSeq IDPs and
ran each mutant through the same generate -> classify -> integrate path as pipeline.py: 512
images per sequence on a log-spaced expression ladder, integrated as in metrics.py. That is
millions of generated images, far past a demo's GPU budget, so the results are precomputed by
notebooks/tools/build_substitution_scan_assets.py and published in the model repo as
condenseq/substitution_scan.npz; this module only reads and draws them.

A matrix cell is the mean change over the IDPs that carry the FROM residue -- an IDP without
it is left out rather than counted as a zero change, as the analysis script did. The heatmaps
are drawn the way that script drew the paper's (CELL-Diff2-Dev
cell_diff/tasks/analysis/ana_all_mutation_log_scale.py): bwr, centred on zero over each
matrix's own range, rows FROM and columns TO, in its residue order.
"""

import os
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

MODEL_REPO = os.environ.get("CELLFM_MODEL_REPO", "BoHuangLab/CELL-FM")
LOCAL_WEIGHTS = os.environ.get("CELLFM_LOCAL_WEIGHTS", "")
SCAN_FILE = "condenseq/substitution_scan.npz"

METRICS = ("auc", "aac")
LABELS = {"auc": "AUC", "aac": "AAC"}
PAPER_FIGURES = {"auc": "Fig. 5d", "aac": "Supplementary Fig. 5c"}

# The app's palette: AUC in the titration curve's blue, AAC in the orange of the gap above it
COLORS = {"auc": "#2a78d6", "aac": "#eb6834"}
INK = "#0b0b0b"
MUTED = "#52514e"
SURFACE = "#fff"  # white, like the tab's other plots and the Virtual OpenCell map

# Shown before any click: F -> S, the substitution the paper's NUP98 validation makes one F at a time
DEFAULT_SELECTION = ("F", "S")
# and compared: the rows of the residues the paper finds drive condensation, aromatics and arginine
DEFAULT_COMPARISON = ("row", ("F", "Y", "W", "R"))

# Series colours for compared rows or columns: the categorical palette's eight slots in fixed order
# (validated on white with the dataviz checks; three are under 3:1 against it, hence the direct
# labels and the table view that go with the plot)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def mean_matrix(delta: np.ndarray) -> np.ndarray:
    """(from, to, idp) changes -> (from, to) mean over the scored IDPs; 0 where none are scored,
    which is the diagonal."""
    scored = np.isfinite(delta)
    total = np.where(scored, delta, 0.0).sum(axis=2)
    count = scored.sum(axis=2)
    return np.divide(total, count, out=np.zeros_like(total), where=count > 0)


@dataclass(frozen=True)
class Scan:
    amino_acids: tuple       # matrix order, R H K D E S T N Q G A V I L M C F Y W P
    idr_index: np.ndarray    # (100,) CondenSeq index of each wild-type IDP
    wt_sequence: np.ndarray  # (100,) its 66-aa sequence
    wt: dict                 # metric -> (100,) wild-type value
    delta: dict              # metric -> (20, 20, 100) mutant minus wild type, NaN where not scored
    mean: dict               # metric -> (20, 20) the matrix the paper plots

    def index(self, residue: str) -> int:
        return self.amino_acids.index(residue)


@lru_cache(maxsize=1)
def load_scan() -> Scan:
    """The precomputed scan: a flat local copy if CELLFM_LOCAL_WEIGHTS is set, else the model repo."""
    if LOCAL_WEIGHTS:
        path = os.path.join(LOCAL_WEIGHTS, os.path.basename(SCAN_FILE))
    else:
        from huggingface_hub import hf_hub_download

        path = hf_hub_download(repo_id=MODEL_REPO, filename=SCAN_FILE)

    saved = np.load(path)
    wt = {m: saved[f"wt_{m}"] for m in METRICS}
    delta = {m: saved[f"mut_{m}"] - wt[m] for m in METRICS}
    return Scan(
        amino_acids=tuple(str(a) for a in saved["amino_acids"]),
        idr_index=saved["idr_index"],
        wt_sequence=saved["wt_sequence"],
        wt=wt,
        delta=delta,
        mean={m: mean_matrix(delta[m]) for m in METRICS},
    )


def heatmap(scan: Scan, metric: str, selected=None, compared=None):
    """(image, layout): the matrix drawn as in the paper, as an RGB array, and where its parts sit
    in that image, so a click maps back to what it landed on (see hit).

    selected, a (from, to) index pair, is outlined. compared is ("row" | "col", [(index, colour)]):
    rows or columns being compared, each outlined, with its residue in bold and a mark in its
    series colour beside it -- the colour is the series', the text stays ink.
    """
    matrix = scan.mean[metric]
    n = len(scan.amino_acids)

    fig = Figure(figsize=(8, 7), dpi=100, facecolor=SURFACE)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot()
    im = ax.imshow(matrix, origin="upper", cmap="bwr",
                   norm=TwoSlopeNorm(vcenter=0, vmin=matrix.min(), vmax=matrix.max()))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=20)
    ax.set_xticks(np.arange(n), scan.amino_acids)
    ax.set_yticks(np.arange(n), scan.amino_acids)
    ax.tick_params(axis="x", labelsize=18, top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.tick_params(axis="y", labelsize=18)
    ax.set_xlabel("Mutate TO", fontsize=18)
    ax.xaxis.set_label_position("top")
    ax.set_ylabel("Mutate FROM", fontsize=18)
    ax.set_title(f"Mean Δ {LABELS[metric]}", fontsize=20)
    fig.tight_layout()

    # Everything below is drawn after the layout, so the matrix sits at the same pixels in every
    # render and a click always maps the same way. Marks go just outside it: rows beside the
    # right edge, before the colour bar; columns under the bottom edge, where no labels are.
    axis, marks = compared if compared is not None else ("row", [])
    for k, colour in marks:
        if axis == "row":
            ax.add_patch(Rectangle((-0.5, k - 0.5), n, 1, fill=False, edgecolor=INK, linewidth=1.6, clip_on=False))
            ax.add_patch(Rectangle((n - 0.38, k - 0.4), 0.34, 0.8, color=colour, clip_on=False))
            ax.get_yticklabels()[k].set_fontweight("bold")
        else:
            ax.add_patch(Rectangle((k - 0.5, -0.5), 1, n, fill=False, edgecolor=INK, linewidth=1.6, clip_on=False))
            ax.add_patch(Rectangle((k - 0.4, n - 0.38), 0.8, 0.34, color=colour, clip_on=False))
            ax.get_xticklabels()[k].set_fontweight("bold")
    if selected is not None:
        i, j = selected
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK,
                               linewidth=3, clip_on=False))

    canvas.draw()
    image = np.asarray(canvas.buffer_rgba())[..., :3].copy()

    # display coordinates count up from the bottom; image rows count down from the top
    height = image.shape[0]
    box = ax.get_window_extent()
    renderer = canvas.get_renderer()
    rows = [t.get_window_extent(renderer) for t in ax.get_yticklabels()]
    cols = [t.get_window_extent(renderer) for t in ax.get_xticklabels()]
    layout = {
        "matrix": (box.x0, height - box.y1, box.x1, height - box.y0),  # left, top, right, bottom
        "row_labels": min(b.x0 for b in rows) - 6,                     # left edge of the FROM labels
        "col_labels": height - max(b.y1 for b in cols) - 6,            # top edge of the TO labels
    }
    return image, layout


def hit(layout, x: float, y: float, n: int = 20):
    """What pixel (x, y) of a heatmap image is on: ("cell", i, j) in the matrix, ("row", i) on a
    Mutate FROM label, ("col", j) on a Mutate TO label, or None anywhere else."""
    left, top, right, bottom = layout["matrix"]
    i = int(np.floor((y - top) / (bottom - top) * n))
    j = int(np.floor((x - left) / (right - left) * n))
    on_row, on_col = 0 <= i < n, 0 <= j < n
    if on_row and on_col:
        return ("cell", i, j)
    if on_row and layout["row_labels"] <= x < left:
        return ("row", i)
    if on_col and layout["col_labels"] <= y < top:
        return ("col", j)
    return None


def assign_colours(residues, previous: dict) -> dict:
    """residue -> series colour for the compared residues. A residue keeps the colour it had, and a
    new one takes the first free slot, so adding or removing one never repaints the others."""
    colours = {r: c for r, c in previous.items() if r in residues}
    for residue in residues:
        if residue not in colours:
            colours[residue] = next(c for c in SERIES if c not in colours.values())
    return colours


def _profile(scan: Scan, metric: str, axis: str, k: int):
    """Mean change and its 95% confidence half-width along compared row/column k, NaN on the diagonal."""
    delta = scan.delta[metric][k] if axis == "row" else scan.delta[metric][:, k]  # (other, idp)
    scored = np.isfinite(delta)
    n = scored.sum(axis=1)
    mean = (scan.mean[metric][k] if axis == "row" else scan.mean[metric][:, k]).astype(float).copy()
    sd = np.array([np.std(d[s], ddof=1) if c > 1 else np.nan for d, s, c in zip(delta, scored, n)])
    half = 1.96 * sd / np.sqrt(np.maximum(n, 1))
    mean[k] = half[k] = np.nan
    return mean, half, n


def comparison(scan: Scan, axis: str, colours: dict):
    """The compared rows ("row": Mutate FROM) or columns ("col": Mutate TO) side by side: each one's
    mean change along the other axis, with 95% confidence intervals over the IDPs, AUC left and
    AAC right. Plotly, for the hover; a figure with a note when nothing is compared."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    aa = list(scan.amino_acids)
    other = "Mutate TO" if axis == "row" else "Mutate FROM"
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.09,
                        subplot_titles=[f"Mean Δ {LABELS[m]}" for m in METRICS])

    for col, metric in enumerate(METRICS, start=1):
        fig.add_hline(y=0, line=dict(color=MUTED, width=1, dash="dot"), row=1, col=col)
        for residue, colour in colours.items():
            k = scan.index(residue)
            mean, half, n = _profile(scan, metric, axis, k)
            name = f"FROM {residue}" if axis == "row" else f"TO {residue}"
            pairs = [f"{residue} → {a}" if axis == "row" else f"{a} → {residue}" for a in aa]
            fig.add_trace(go.Scatter(
                x=aa, y=mean, mode="lines+markers", name=name, legendgroup=name, showlegend=col == 1,
                line=dict(color=colour, width=2), marker=dict(color=colour, size=8, line=dict(color=SURFACE, width=2)),
                error_y=dict(type="data", array=half, color=colour, thickness=1, width=0),
                customdata=np.stack([pairs, half, n], axis=1),
                hovertemplate=(f"<b>%{{customdata[0]}}</b><br>mean Δ {LABELS[metric]} %{{y:+.4f}}"
                               " ± %{customdata[1]:.4f} (95% CI)<br>%{customdata[2]} IDPs<extra></extra>"),
            ), row=1, col=col)
            if len(colours) <= 4:  # a few series are named where they end, as well as in the legend
                last = int(np.flatnonzero(np.isfinite(mean))[-1])
                fig.add_annotation(x=aa[last], y=mean[last], text=residue, showarrow=False, xshift=12,
                                   font=dict(size=12, color=INK), row=1, col=col)
        fig.update_xaxes(title_text=other, categoryorder="array", categoryarray=aa, showgrid=False,
                         row=1, col=col)
        fig.update_yaxes(title_text=f"Mean Δ {LABELS[metric]} (mutant − wild type)", zeroline=False,
                         gridcolor="rgba(82,81,78,0.15)", row=1, col=col)

    if not colours:
        fig.add_annotation(text="Click a residue on either axis of a matrix, or pick some above",
                           xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                           font=dict(size=13, color=MUTED))
    fig.update_layout(
        height=430, margin=dict(l=10, r=10, t=40, b=10), paper_bgcolor=SURFACE, plot_bgcolor=SURFACE,
        font=dict(color=INK, size=12), legend=dict(orientation="h", y=1.14, x=0, font=dict(size=12)),
        hovermode="closest",
    )
    return fig


def comparison_table(scan: Scan, axis: str, colours: dict) -> pd.DataFrame:
    """The numbers the comparison plots: one row per residue of the other axis, a mean change and
    its confidence half-width per compared residue and metric."""
    table = {("Mutate TO" if axis == "row" else "Mutate FROM"): list(scan.amino_acids)}
    for residue in colours:
        k = scan.index(residue)
        name = f"FROM {residue}" if axis == "row" else f"TO {residue}"
        for metric in METRICS:
            mean, half, _ = _profile(scan, metric, axis, k)
            table[f"{name} Δ {LABELS[metric]}"] = mean
            table[f"{name} ± 95% CI ({LABELS[metric]})"] = half
    return pd.DataFrame(table).round(4)


def distribution(scan: Scan, i: int, j: int) -> Figure:
    """The per-IDP changes behind cell (i, j) of both matrices, side by side."""
    source, target = scan.amino_acids[i], scan.amino_acids[j]

    fig = Figure(figsize=(11, 3.9), dpi=100, facecolor=SURFACE)
    for ax, metric in zip(fig.subplots(1, 2), METRICS):
        values = scan.delta[metric][i, j]
        values = values[np.isfinite(values)]

        ax.set_facecolor(SURFACE)
        ax.hist(values, bins=20, color=COLORS[metric], alpha=0.8, edgecolor=SURFACE, linewidth=0.6)
        ax.axvline(0, color=MUTED, linewidth=1.0, linestyle=":")
        ax.axvline(values.mean(), color=INK, linewidth=1.8,
                   label=f"mean {values.mean():+.4f} (matrix cell)")
        ax.axvline(np.median(values), color=INK, linewidth=1.2, linestyle="--",
                   label=f"median {np.median(values):+.4f}")

        ax.set_title(f"Δ {LABELS[metric]}, {source} → {target}  ({len(values)} IDPs)",
                     fontsize=12, color=INK)
        ax.set_xlabel(f"Δ normalized {LABELS[metric]} (mutant − wild type)", fontsize=10, color=INK)
        ax.set_ylabel("IDPs", fontsize=10, color=INK)
        ax.tick_params(labelsize=9, colors=MUTED)
        ax.grid(True, axis="y", alpha=0.18, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
            ax.spines[side].set_linewidth(0.8)
        ax.legend(loc="upper left", frameon=False, fontsize=9, labelcolor=INK)

    fig.tight_layout()
    return fig


def idp_table(scan: Scan, i: int, j: int) -> pd.DataFrame:
    """One row per IDP scored for FROM -> TO: wild-type and mutant values and the change."""
    source, target = scan.amino_acids[i], scan.amino_acids[j]
    scored = np.isfinite(scan.delta["auc"][i, j])
    table = {
        "CondenSeq index": scan.idr_index,
        f"{source} residues mutated": [s.count(source) for s in scan.wt_sequence],
    }
    for metric in METRICS:
        label = LABELS[metric]
        table[f"WT {label}"] = scan.wt[metric]
        table[f"{source}→{target} {label}"] = scan.wt[metric] + scan.delta[metric][i, j]
        table[f"Δ {label}"] = scan.delta[metric][i, j]
    table["wild-type sequence"] = scan.wt_sequence
    return pd.DataFrame(table)[scored].round(4).reset_index(drop=True)
