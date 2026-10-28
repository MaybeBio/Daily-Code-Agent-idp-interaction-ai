"""Virtual OpenCell: CELL-FM's virtual staining of the OpenCell proteins, browsed like OpenCell.

The images are the dataset repo BoHuangLab/CELL-FM, one folder per protein:

    OpenCell/virtual_staining/<GENE>/<NNNN>.tif    (2, 256, 256) uint16

Channel 0 is the nucleus every sample was conditioned on: one real cell (ATG7, crop 1),
shared by all 1276 proteins, so differences between proteins are the model's and not the
cell's. Channel 1 is the protein CELL-FM generated from the sequence; a protein's samples
are independent draws for that same cell. They are fetched the first time a protein is
viewed and kept in memory.

The annotations shown beside the images (the catalog CSV, built by
notebooks/tools/build_virtual_opencell_assets.py) are OpenCell's measurements of the *real*
cell line -- localization grades, identifiers, abundance -- shown for comparison, not
predicted.

Images come from the Hub by default (HF_TOKEN is used if the repo is private). Set
VIRTUAL_OPENCELL_LOCAL_DIR to a local copy of the dataset, the folder holding OpenCell/, to
read them from disk instead.

The map of all samples is opencell/vs_umap_embedding.npz in the weights repo: one UMAP point
per sample (93,510, grouped by protein), embedded with opencell/vit.bin.

The 3D virtual staining (CELL-FM 3D) is in the same repo, 64 samples for every protein:

    OpenCell/virtual_staining_3d/<GENE>/<NNNN>.tif    (48, 2, 192, 192) uint16, z first
    OpenCell/virtual_staining_3d_thumbnails/<GENE>.jpg

A stack is about 7 MB, so a protein's 64 are fetched one at a time, when opened, and the sample
strip comes from the thumbnail sheet (notebooks/tools/build_virtual_opencell_3d_thumbnails.py).
Stacks are shown one z slice at a time, like OpenCell's "2D slice" view. There is no 3D map yet.
"""

import html
import os
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import tifffile
from PIL import Image

DATASET_REPO = os.environ.get("VIRTUAL_OPENCELL_REPO", "BoHuangLab/CELL-FM")
LOCAL_DIR = os.environ.get("VIRTUAL_OPENCELL_LOCAL_DIR", "")
IMAGE_DIR = "OpenCell/virtual_staining"
IMAGE_DIR_3D = "OpenCell/virtual_staining_3d"
THUMBNAIL_DIR_3D = "OpenCell/virtual_staining_3d_thumbnails"
SAMPLES_3D = 64
SLICES_3D = 48
# Crops are centred on the nucleus' 3D centroid, so its centre slice is the default view
DEFAULT_SLICE_3D = SLICES_3D // 2
THUMBNAIL_SIZE_3D = 96
THUMBNAIL_COLUMNS_3D = 8

MODEL_REPO = os.environ.get("CELLFM_MODEL_REPO", "BoHuangLab/CELL-FM")
UMAP_FILE = "opencell/vs_umap_embedding.npz"
# Backdrop points per protein. The figure is re-sent on every selection, so the map is drawn
# from an even subsample (~15k points, a few hundred kB); the selected protein shows all of its own
UMAP_POINTS_PER_PROTEIN = 12

DEFAULT_GENE = "TOMM20"
CHANNELS = ["Nucleus", "Target", "Both"]  # OpenCell's names for the 405 and 488 channels

# OpenCell draws the nucleus in blue and the tagged target in gray; colours add in "Both"
NUCLEUS_RGB = np.array([0.0, 0.45, 1.0], dtype=np.float32)
TARGET_RGB = np.array([1.0, 1.0, 1.0], dtype=np.float32)
# Auto-contrast: 0-100% of the intensity range maps this percentile window to black-white
AUTO_PERCENTILES = (0.1, 99.9)

THUMBNAIL_SIZE = 128
OPENCELL_URL = "https://opencell.sf.czbiohub.org"

# OpenCell's localization categories, as its pages print them
LOCALIZATION_NAMES = {
    "er": "ER",
    "golgi": "Golgi",
    "nucleolus_gc": "nucleolus GC",
    "nucleolus_fc_dfc": "nucleolus FC/DFC",
}

# The palette of the CELL-FM UMAP figures (vit_cls/embed_opencell.py), so the map reads the same
LOCATION_COLORS = {
    "Multilocalizing": "#999999", "big_aggregates": "#D50000", "cell_contact": "#018600",
    "centrosome": "#b400fe", "chromatin": "#04abc5", "cytoplasmic": "#96fe00",
    "cytoskeleton": "#fea42e", "er": "#fe8dc7", "focal_adhesions": "#78525d",
    "golgi": "#00fcce", "membrane": "#aea4fe", "mitochondria": "#92ab82",
    "negative": "#996900", "nuclear_membrane": "#356961", "nuclear_punctae": "#d2008b",
    "nucleolus_fc_dfc": "#fcf490", "nucleolus_gc": "#c76d66", "nucleoplasm": "#9de1fe",
    "vesicles": "#00c746", "NA": "#999999",
}


def load_catalog(path: str) -> pd.DataFrame:
    """The protein table, indexed by gene name."""
    return pd.read_csv(path, keep_default_na=False, na_values=[""]).set_index("gene_name")


def _sample_files(gene: str, n_samples: int) -> list:
    """Local paths of a protein's samples, downloading them first unless reading from disk.

    Files are numbered 0001.tif .. NNNN.tif without gaps, so they are requested by name;
    listing the 93k-file repo on every selection would be far slower.
    """
    names = [f"{IMAGE_DIR}/{gene}/{i:04d}.tif" for i in range(1, n_samples + 1)]
    if LOCAL_DIR:
        return [os.path.join(LOCAL_DIR, name) for name in names]

    from huggingface_hub import hf_hub_download

    def fetch(name):
        return hf_hub_download(repo_id=DATASET_REPO, filename=name, repo_type="dataset")

    with ThreadPoolExecutor(max_workers=16) as pool:
        return list(pool.map(fetch, names))


@lru_cache(maxsize=16)
def load_samples(gene: str, n_samples: int):
    """(stack, paths): stack is (N, 2, H, W) uint16, nucleus first."""
    paths = _sample_files(gene, n_samples)
    return np.stack([tifffile.imread(p) for p in paths]), tuple(paths)


def _download(name: str) -> str:
    if LOCAL_DIR:
        return os.path.join(LOCAL_DIR, name)
    from huggingface_hub import hf_hub_download

    return hf_hub_download(repo_id=DATASET_REPO, filename=name, repo_type="dataset")


def sample_file_3d(gene: str, index: int) -> str:
    """Local path of one 3D sample (0-based index), downloading it first unless reading from disk."""
    return _download(f"{IMAGE_DIR_3D}/{gene}/{index + 1:04d}.tif")


@lru_cache(maxsize=32)
def load_sample_3d(gene: str, index: int):
    """(stack, limits): stack is (Z, 2, H, W) uint16; limits are each channel's auto-contrast
    window over the whole stack, so a slice keeps the same contrast as the others."""
    stack = tifffile.imread(sample_file_3d(gene, index))
    limits = tuple(tuple(np.percentile(stack[:, c], AUTO_PERCENTILES)) for c in range(stack.shape[1]))
    return stack, limits


@lru_cache(maxsize=64)
def thumbnails_3d(gene: str) -> list:
    """The sample strip of a protein's 3D stacks: its centre slices, cut from the thumbnail sheet."""
    sheet = np.asarray(Image.open(_download(f"{THUMBNAIL_DIR_3D}/{gene}.jpg")).convert("RGB"))
    size = THUMBNAIL_SIZE_3D
    return [sheet[r:r + size, c:c + size] for r in range(0, sheet.shape[0], size)
            for c in range(0, sheet.shape[1], size)]


def _scale(channel: np.ndarray, low: float, high: float, gamma: float, limits=None) -> np.ndarray:
    """One channel to [0, 1]: auto-contrast (over `limits` if given, else this image's own
    percentiles), then OpenCell's intensity range (%) and gamma."""
    lo, hi = limits if limits is not None else np.percentile(channel, AUTO_PERCENTILES)
    x = (channel.astype(np.float32) - lo) / max(hi - lo, 1.0)
    low, high = low / 100.0, max(high, low + 1.0) / 100.0
    return np.clip((x - low) / (high - low), 0.0, 1.0) ** gamma


def render(sample: np.ndarray, channel: str = "Both",
           nucleus=(0, 100, 1.0), target=(0, 100, 1.0), limits=(None, None)) -> np.ndarray:
    """A (2, H, W) sample as an (H, W, 3) uint8 image; nucleus/target are (min %, max %, gamma),
    and limits each channel's auto-contrast window (None: from this image)."""
    rgb = np.zeros(sample.shape[1:] + (3,), dtype=np.float32)
    if channel in ("Nucleus", "Both"):
        rgb += _scale(sample[0], *nucleus, limits=limits[0])[..., None] * NUCLEUS_RGB
    if channel in ("Target", "Both"):
        rgb += _scale(sample[1], *target, limits=limits[1])[..., None] * TARGET_RGB
    return (np.clip(rgb, 0.0, 1.0) * 255).astype(np.uint8)


def thumbnails(stack: np.ndarray, size: int = THUMBNAIL_SIZE) -> list:
    """Both channels at default settings, downsized for the sample strip."""
    step = max(stack.shape[-1] // size, 1)
    return [render(sample)[::step, ::step] for sample in stack]


def _localization_name(category: str) -> str:
    return LOCALIZATION_NAMES.get(category, category.replace("_", " "))


@lru_cache(maxsize=1)
def umap_data():
    """(xy, labels, genes, backdrop): every sample's UMAP point, OpenCell label and gene, and
    the indices drawn behind the selected protein."""
    from huggingface_hub import hf_hub_download

    saved = np.load(hf_hub_download(repo_id=MODEL_REPO, filename=UMAP_FILE))
    xy, labels, genes = saved["embeddings"], saved["labels"], saved["gene_names"]
    # rows come grouped by protein; take evenly spaced samples from each block
    starts = np.flatnonzero(np.r_[True, genes[1:] != genes[:-1]])
    ends = np.r_[starts[1:], len(genes)]
    backdrop = np.concatenate([
        np.unique(np.linspace(s, e - 1, UMAP_POINTS_PER_PROTEIN).round().astype(int))
        for s, e in zip(starts, ends)
    ])
    return xy, labels, genes, backdrop


@lru_cache(maxsize=1)
def _protein_centres():
    """(names, centres): every protein, labelled at the mean of its points as in the paper figure."""
    xy, _, genes, _ = umap_data()
    starts = np.flatnonzero(np.r_[True, genes[1:] != genes[:-1]])
    centres = np.add.reduceat(xy, starts, axis=0) / np.diff(np.r_[starts, len(genes)])[:, None]
    return genes[starts], centres.round(2)


def umap_figure(gene: str, full: bool = False):
    """The map of every protein's samples, coloured by OpenCell localization, with `gene` outlined.

    The panel map draws the backdrop subsample; full=True draws all 93,510 points and names every
    protein, like umap_embedding_selected_all.svg (about 2 MB). An empty `gene` outlines nothing.
    """
    xy, labels, genes, backdrop = umap_data()
    points = np.arange(len(genes)) if full else backdrop
    selected = np.flatnonzero(genes == gene)
    fig = go.Figure()

    # multilocalizing proteins underneath, as in the paper figures; NA (unannotated) last
    names = sorted(set(labels[points]) - {"Multilocalizing", "NA"})
    for label in ["Multilocalizing"] + names + ["NA"]:
        idx = points[labels[points] == label]
        if not len(idx):
            continue
        name = _localization_name(label)
        size = (2 if full else 3) if label == "Multilocalizing" else (3 if full else 4)
        fig.add_trace(go.Scattergl(
            x=xy[idx, 0].round(2), y=xy[idx, 1].round(2), mode="markers", name=name,
            marker=dict(color=LOCATION_COLORS.get(label, "#999999"), size=size,
                        opacity=0.45 if len(selected) else 0.8),
            customdata=genes[idx], hovertemplate=f"<b>%{{customdata}}</b><br>{name}<extra></extra>",
        ))

    if full:
        # fixed-size labels: zooming in pulls them apart until each one is readable
        protein_names, centres = _protein_centres()
        fig.add_trace(go.Scatter(
            x=centres[:, 0], y=centres[:, 1], mode="text", text=protein_names, showlegend=False,
            textfont=dict(size=8, color="#000"), hoverinfo="skip",
        ))

    if len(selected):
        label = labels[selected[0]]
        fig.add_trace(go.Scatter(
            x=xy[selected, 0], y=xy[selected, 1], mode="markers", name=gene, showlegend=False,
            marker=dict(color=LOCATION_COLORS.get(label, "#999999"), size=8, line=dict(color="#000", width=1.2)),
            hovertemplate=f"<b>{html.escape(gene)}</b><br>{_localization_name(label)}<extra></extra>",
        ))
        x, y = np.median(xy[selected], axis=0)
        fig.add_annotation(x=x, y=y, text=f"<b>{html.escape(gene)}</b>", ax=34, ay=-34, arrowcolor="#000",
                           font=dict(size=14, color="#000"), bgcolor="rgba(255,255,255,0.85)")

    # no axes, and the UMAP1 / UMAP2 corner of the paper figures
    corner = dict(type="line", xref="paper", yref="paper", line=dict(color="#000", width=2))
    fig.update_layout(
        height=1000 if full else 600, margin=dict(l=0, r=0, t=0, b=0), paper_bgcolor="#fff", plot_bgcolor="#fff",
        xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor="x"),
        legend=dict(orientation="h", y=-0.02, yanchor="top", font=dict(size=11), itemsizing="constant"),
        shapes=[dict(corner, x0=0.02, x1=0.14, y0=0.02, y1=0.02), dict(corner, x0=0.02, x1=0.02, y0=0.02, y1=0.14)],
        annotations=list(fig.layout.annotations) + [
            dict(text="UMAP1", x=0.08, y=0.03, xref="paper", yref="paper", showarrow=False, yanchor="bottom", font=dict(size=11)),
            dict(text="UMAP2", x=0.03, y=0.08, xref="paper", yref="paper", showarrow=False, xanchor="left", textangle=-90, font=dict(size=11)),
        ],
    )
    return fig


def _metadata_item(label: str, value: str) -> str:
    return (f'<div class="vo-metadata-item"><div class="vo-metadata-label">{label}</div>'
            f'<div class="vo-metadata-value">{value}</div></div>')


def target_html(gene: str, row: pd.Series, three_d: bool = False) -> str:
    """The left-hand panel of an OpenCell target page: name, identifiers, localization."""
    esc = html.escape
    uniprot, ensg, cid = row["uniprot_id"], row["ensg_id"], row["opencell_cid"]

    items = [
        _metadata_item("UniProt ID", f'<a href="https://www.uniprot.org/uniprotkb/{esc(uniprot)}" '
                                     f'target="_blank">{esc(uniprot)}</a>'),
        _metadata_item("Ensembl ID", f'<a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={esc(ensg)}" '
                                     f'target="_blank">{esc(ensg)}</a>'),
        _metadata_item("Sequence length", f"{int(row['sequence_length'])} aa"),
    ]
    if pd.notna(row["hek_protein_conc_nm"]):
        items.append(_metadata_item("HEK293 abundance", f"{row['hek_protein_conc_nm']:,.0f} nM"))

    rows = []
    for grade in (3, 2, 1):
        value = row[f"grade_{grade}"]
        for category in ([] if pd.isna(value) else value.split(";")):
            rows.append(
                f'<div class="vo-localization-row"><div class="vo-grade-container">'
                f'<div class="vo-grade vo-grade-{grade}"></div></div>'
                f'<div>{esc(_localization_name(category))}</div></div>'
            )
    localization = "".join(rows) or '<div class="vo-muted">Not annotated by OpenCell</div>'

    name = row["protein_name"] if pd.notna(row["protein_name"]) else ""
    return f"""
<div class="vo-target">
  <div class="vo-target-name">{esc(gene)}</div>
  <div class="vo-protein-description">{esc(name)}</div>
  <div class="vo-status">CELL-FM {"3D " if three_d else ""}virtual staining &middot; {SAMPLES_3D if three_d else int(row['n_samples'])} samples</div>
  <div class="vo-metadata">{''.join(items)}</div>
  <div class="vo-section-header">Protein localization</div>
  <div class="vo-section-caption">OpenCell's annotation of the real cell line; bar height is the grade</div>
  {localization}
  <div class="vo-external-links">
    <a href="{OPENCELL_URL}/target/{esc(cid)}" target="_blank">Real images on OpenCell</a>
    <a href="https://www.uniprot.org/uniprotkb/{esc(uniprot)}" target="_blank">UniProt</a>
    <a href="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={esc(ensg)}" target="_blank">Ensembl</a>
  </div>
</div>
"""


def load_error(exc: Exception) -> str:
    """Why a protein's images could not be loaded, in words a visitor can act on."""
    # A private repo answers 404 to anyone without access, so "not found" means "no access" here
    if type(exc).__name__ in ("RepositoryNotFoundError", "GatedRepoError"):
        return f"the image dataset {DATASET_REPO} is not readable here (private, with no token that can read it)."
    return str(exc).splitlines()[0]


def error_html(gene: str, exc: Exception) -> str:
    return (f'<div class="vo-target"><div class="vo-target-name">{html.escape(gene)}</div>'
            f'<div class="vo-muted">Could not load the images: {html.escape(load_error(exc))}</div></div>')


UMAP_HTML = f"""
<div class="vo-target">
  <div class="vo-section-header">Virtual staining map</div>
  <div class="vo-section-caption">Each dot is one generated image, placed by UMAP of its embedding from CELL-FM's
  OpenCell ViT and coloured by OpenCell's localization of the real cell line. The background shows
  {UMAP_POINTS_PER_PROTEIN} images per protein; the selected protein shows all of its images, outlined in black.
  The full map, with every image and every protein named, is below.</div>
</div>
"""

NAVBAR_HTML = """
<div class="vo-navbar-brand">
  <span class="vo-navbar-title">Virtual OpenCell</span>
  <span class="vo-navbar-menu">CELL-FM virtual staining of 1,276 OpenCell proteins</span>
</div>
"""

FOOTER_HTML = f"""
<div class="vo-footer">
  Each image is CELL-FM's prediction from the protein sequence alone, drawn in one shared cell:
  the nucleus channel is the same real nucleus for every protein, so differences between
  proteins come from the sequence. Localization grades, identifiers and abundance are
  OpenCell's measurements of the real cell lines
  (<a href="{OPENCELL_URL}" target="_blank">OpenCell</a>; Cho et al., <i>Science</i> 2022).
</div>
"""

# Fonts OpenCell uses: Lato for names and titles, Nunito Sans for metadata labels
HEAD = """
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Lato:wght@300;400;600;700&family=Nunito+Sans:wght@300;400;700&display=swap" rel="stylesheet">
"""

# OpenCell's look (its stylesheet: navbar, clm-*, metadata-*, localization-grade-*,
# slice-viewer-*, roi-thumbnail-*), scoped to the tab so the other sections keep the theme
CSS = """
#virtual-opencell {
  --vo-blue: #01a1dd;
  --color-accent: #137cbd;
  --color-accent-soft: #e8f2fa;
  --slider-color: #137cbd;
  /* OpenCell is flat: no cards, small grey labels instead of the theme's badges */
  --block-background-fill: transparent;
  --block-border-width: 0px;
  --block-shadow: none;
  --block-title-background-fill: none;
  --block-title-border-width: 0px;
  --block-title-text-color: #999;
  --block-title-text-size: 12px;
  --block-title-text-weight: 400;
  --block-title-padding: 0;
  --block-label-background-fill: none;
  --block-label-text-color: #999;
  --panel-background-fill: transparent;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  color: #182026;
}
#virtual-opencell a { color: #333; text-decoration: underline; }
#virtual-opencell a:hover { color: #999; }

/* navbar: grey band, cyan Lato title, small-caps menu, bold uppercase search */
#virtual-opencell .vo-navbar {
  background: #eee; border-radius: 6px; padding: 8px 16px; align-items: center; margin-bottom: 8px;
}
#virtual-opencell .vo-navbar-title {
  color: var(--vo-blue); font-family: "Lato", sans-serif; font-size: 30px; font-weight: 600;
  letter-spacing: 1px; margin-right: 24px;
}
#virtual-opencell .vo-navbar-menu {
  color: #333; font-size: 20px; font-variant: all-small-caps;
}
#virtual-opencell .vo-search input {
  font-size: 22px !important; font-weight: bold; text-transform: uppercase; background: #fff;
}

/* target panel (OpenCell's clm-*) */
#virtual-opencell .vo-target { font-family: "Lato", sans-serif; padding-right: 12px; }
#virtual-opencell .vo-target-name { font-size: 60px; font-weight: 400; line-height: 1.1; }
#virtual-opencell .vo-protein-description { font-size: 20px; font-weight: 300; margin-top: 4px; }
#virtual-opencell .vo-status { font-size: 16px; font-weight: 300; font-style: italic; margin-top: 6px; color: #555; }
#virtual-opencell .vo-metadata { display: flex; flex-wrap: wrap; margin: 18px 0 8px 0; }
#virtual-opencell .vo-metadata-item { margin: 0 25px 12px 0; max-width: 200px; color: #333; font-size: 15px; }
#virtual-opencell .vo-metadata-label { color: #aaa; font-family: "Nunito Sans", sans-serif; font-size: 13px; }
#virtual-opencell .vo-section-header { font-size: 20px; font-weight: 400; margin-top: 14px; border-bottom: 1px solid #eee; }
#virtual-opencell .vo-section-caption { color: #999; font-size: 12px; margin: 4px 0 6px 0; }
#virtual-opencell .vo-localization-row { display: flex; align-items: center; padding: 3px 0; font-size: 16px; }
#virtual-opencell .vo-grade-container { width: 32px; height: 22px; padding: 0 3px; margin-right: 8px; display: flex; align-items: flex-end; }
#virtual-opencell .vo-grade { width: 100%; border-radius: 5px; }
#virtual-opencell .vo-grade-1 { height: 10px; background-color: #ccdeef; }
#virtual-opencell .vo-grade-2 { height: 14px; background-color: #51ade1cc; }
#virtual-opencell .vo-grade-3 { height: 18px; background-color: #00477bcc; }
#virtual-opencell .vo-external-links { margin-top: 18px; font-size: 14px; }
#virtual-opencell .vo-external-links a { margin-right: 16px; }
#virtual-opencell .vo-muted { color: #999; font-size: 14px; }

/* Gradio groups inputs in a grey .form box; OpenCell lays them straight on the page */
#virtual-opencell .form { background: transparent !important; border: none !important; box-shadow: none !important; }

/* viewer: black, rounded, the image filling it (slice-viewer-canvas-container) */
#virtual-opencell .vo-viewer { background: #000 !important; border-radius: 10px !important; overflow: hidden; }
#virtual-opencell .vo-viewer .image-container,
#virtual-opencell .vo-viewer .image-frame { width: 100%; height: 100%; }
#virtual-opencell .vo-viewer img { width: 100%; height: 100%; object-fit: contain; }

/* OpenCell's button groups: grey label, rounded simple buttons, bold when active */
#virtual-opencell .vo-buttons input[type="radio"] { display: none; }
#virtual-opencell .vo-buttons label {
  background: transparent !important; border: none !important; box-shadow: none !important;
  padding: 4px 10px; border-radius: 5px; color: #333; font-size: 15px;
}
#virtual-opencell .vo-buttons label:hover { background: #d2d2d2 !important; }
#virtual-opencell .vo-buttons label.selected { font-weight: bold; background: #d1d4d64d !important; }
#virtual-opencell .vo-download button { font-size: 13px; font-weight: 400; }
#virtual-opencell .vo-settings-label { font-size: 16px; color: #182026; margin: 6px 0 -6px 0; }

/* sample strip: white frame, cyan when active (roi-thumbnail-*) */
/* the gallery block clips to its height, so the grid must scroll inside it rather than overflow */
#virtual-opencell .vo-thumbnails .gallery-container { height: 100%; }
#virtual-opencell .vo-thumbnails .grid-wrap { height: 100%; overflow-y: auto; }
#virtual-opencell .vo-thumbnails .grid-container { grid-template-columns: repeat(var(--grid-cols), minmax(0, 1fr)) !important; }
#virtual-opencell .vo-thumbnails .thumbnail-item { border: 3px solid #fff; border-radius: 5px; box-shadow: none; }
#virtual-opencell .vo-thumbnails .thumbnail-item:hover { opacity: 0.7; }
#virtual-opencell .vo-thumbnails .thumbnail-item.selected { border-color: var(--vo-blue); box-shadow: none; }

#virtual-opencell .vo-footer {
  background: #eee; border-radius: 6px; padding: 10px 16px; margin-top: 10px; font-size: 13px; color: #333;
}
"""


def scoped_css(scope: str) -> str:
    """The same look for another container: every rule above is written against #virtual-opencell,
    so re-pointing them gives a second tab OpenCell's style without touching this one's."""
    return CSS.replace("#virtual-opencell", scope)
