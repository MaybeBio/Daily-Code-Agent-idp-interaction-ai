# -*- coding: utf-8 -*-
import os
import sys

import torch
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import umap
import plotly.graph_objects as go

from tqdm import tqdm
from torch.utils.data import DataLoader
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import (
    adjusted_rand_score, v_measure_score, fowlkes_mallows_score,
    normalized_mutual_info_score, adjusted_mutual_info_score,
    homogeneity_score, completeness_score, rand_score
)
from sklearn.metrics.cluster import contingency_matrix
from scipy.optimize import linear_sum_assignment

from cell_fm.data.opencell_3d_crop_data.dataset import OpenCell3DCropImageOnlyDataset
from cell_fm.models.vit_cls_3d.config import ViT3DConfig
from cell_fm.models.vit_cls_3d.model import ViT3DModel
from cell_fm.utils.cli_utils import cli


location_colors = {
    "Multilocalizing": "#999999",
    "big_aggregates": "#D50000",
    "cell_contact": "#018600",
    "centrosome": "#b400fe",
    "chromatin": "#04abc5",
    "cytoplasmic": "#96fe00",
    "cytoskeleton": "#fea42e",
    "er": "#fe8dc7",
    "focal_adhesions": "#78525d",
    "golgi": "#00fcce",
    "membrane": "#aea4fe",
    "mitochondria": "#92ab82",
    "negative": "#996900",
    "nuclear_membrane": "#356961",
    "nuclear_punctae": "#d2008b",
    "nucleolus_fc_dfc": "#fcf490",
    "nucleolus_gc": "#c76d66",
    "nucleoplasm": "#9de1fe",
    "vesicles": "#00c746",
    "NA": "#999999",
}


def purity_score(y_true, y_pred):
    cm = contingency_matrix(y_true, y_pred)
    return np.sum(np.max(cm, axis=0)) / np.sum(cm)


def b3_scores(y_true, y_pred):
    cm = contingency_matrix(y_true, y_pred).astype(float)
    a = cm.sum(axis=1)
    b = cm.sum(axis=0)
    n = cm.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        prec_num = (cm * (cm / b)).sum()
        rec_num  = (cm * (cm / a[:, None])).sum()
    prec = prec_num / n
    rec  = rec_num  / n
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return float(prec), float(rec), float(f1)


def hungarian_matched_accuracy(y_true, y_pred):
    cm = contingency_matrix(y_true, y_pred)
    r_ind, c_ind = linear_sum_assignment(-cm)
    return cm[r_ind, c_ind].sum() / cm.sum()


def variation_of_information(y_true, y_pred, *, base=np.e):
    cm = contingency_matrix(y_true, y_pred).astype(float)
    n = cm.sum()
    pi = cm.sum(axis=1) / n
    pj = cm.sum(axis=0) / n
    pij = cm / n

    def H(p):
        p = p[p > 0]
        return float(-np.sum(p * (np.log(p) / np.log(base))))

    H_C  = H(pi)
    H_K  = H(pj)
    H_CK = H(pij.flatten())
    return (H_CK - H_K) + (H_CK - H_C)


def _cluster_labels(X, method="kmeans", n_clusters=None, random_state=42, **kw):
    if method == "kmeans":
        return KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10, **kw).fit_predict(X)
    if method == "spectral":
        return SpectralClustering(n_clusters=n_clusters, affinity=kw.pop("affinity", "nearest_neighbors"),
                                  random_state=random_state, **kw).fit_predict(X)


def get_scores(embedding, labels, method="spectral"):
    labels = np.array(labels, dtype=object)
    mask = (labels != "Multilocalizing") & (labels != "NA")
    X, y = embedding[mask], labels[mask]

    unique_y = np.unique(y)
    if X.shape[0] < 2 or unique_y.size < 2:
        raise ValueError(f"Not enough samples/classes after filtering: n={X.shape[0]}, k={unique_y.size}")

    X_std = StandardScaler().fit_transform(X)
    n_clusters = unique_y.size
    y_pred = _cluster_labels(X_std, method=method, n_clusters=n_clusters)

    ari = adjusted_rand_score(y, y_pred)
    v   = v_measure_score(y, y_pred)
    fmi = fowlkes_mallows_score(y, y_pred)
    nmi = normalized_mutual_info_score(y, y_pred, average_method="arithmetic")
    ami = adjusted_mutual_info_score(y, y_pred)
    h   = homogeneity_score(y, y_pred)
    c   = completeness_score(y, y_pred)
    ri  = rand_score(y, y_pred)
    purity = purity_score(y, y_pred)
    b3_p, b3_r, b3_f1 = b3_scores(y, y_pred)
    hung_acc = hungarian_matched_accuracy(y, y_pred)
    vi = variation_of_information(y, y_pred, base=np.e)

    print(f"[{method}] n_samples={X.shape[0]}, n_clusters={n_clusters}")
    print(f"Adjusted Rand Index     : {ari:.4f}")
    print(f"V-measure               : {v:.4f}")
    print(f"Fowlkes-Mallows Index   : {fmi:.4f}")
    print(f"NMI                     : {nmi:.4f}")
    print(f"AMI                     : {ami:.4f}")
    print(f"Homogeneity             : {h:.4f}")
    print(f"Completeness            : {c:.4f}")
    print(f"Rand Index              : {ri:.4f}")
    print(f"Purity                  : {purity:.4f}")
    print(f"B3 P/R/F1               : {b3_p:.4f} / {b3_r:.4f} / {b3_f1:.4f}")
    print(f"Hungarian Accuracy      : {hung_acc:.4f}")
    print(f"Variation of Information: {vi:.4f}  (lower is better)")

    return {
        "ARI": ari, "V": v, "FMI": fmi, "NMI": nmi, "AMI": ami,
        "Homogeneity": h, "Completeness": c, "RandIndex": ri,
        "Purity": purity, "B3_P": b3_p, "B3_R": b3_r, "B3_F1": b3_f1,
        "Hungarian_Acc": hung_acc, "VI": vi,
        "n_samples": int(X.shape[0]), "n_clusters": int(n_clusters),
    }


_JS_HOVER = """
(function() {
    var gd = document.querySelector('.plotly-graph-div');
    var geneMap = {};
    var base = gd.data[0];
    base.customdata.forEach(function(cd, i) {
        var gene = cd[0];
        if (!geneMap[gene]) geneMap[gene] = {x: [], y: [], color: cd[2], loc: cd[1]};
        geneMap[gene].x.push(base.x[i]);
        geneMap[gene].y.push(base.y[i]);
    });
    var lastGene = null;
    gd.on('plotly_hover', function(data) {
        var gene = data.points[0].customdata[0];
        if (gene === lastGene) return;
        lastGene = gene;
        var entry = geneMap[gene];
        Plotly.restyle(gd, {
            x: [entry.x], y: [entry.y],
            'marker.color': [entry.color],
            'hovertemplate': ['<b>' + gene + '</b><br>' + entry.loc + '<extra></extra>']
        }, [1]);
    });
    gd.on('plotly_unhover', function() {
        lastGene = null;
        Plotly.restyle(gd, {x: [[null]], y: [[null]]}, [1]);
    });
})();
"""


def make_interactive_plot(embedding_2d, all_labels, all_genes, out_path):
    genes_arr = np.array(all_genes)
    labels_arr = np.array(all_labels)
    default_color = "#999999"

    gene_color = {
        g: location_colors.get(labels_arr[genes_arr == g][0], default_color)
        for g in np.unique(genes_arr)
    }

    customdata = [
        [g, l, gene_color.get(g, default_color)]
        for g, l in zip(all_genes, all_labels)
    ]

    base_trace = go.Scattergl(
        x=embedding_2d[:, 0], y=embedding_2d[:, 1],
        mode='markers',
        marker=dict(color=default_color, size=4, opacity=0.7),
        customdata=customdata,
        hovertemplate="<b>%{customdata[0]}</b><br>%{customdata[1]}<extra></extra>",
        name='',
    )
    highlight_trace = go.Scattergl(
        x=[None], y=[None], mode='markers',
        marker=dict(color=default_color, size=9, opacity=1.0,
                    line=dict(color='white', width=0.5)),
        hoverinfo='skip', name='', showlegend=False,
    )

    fig = go.Figure([base_trace, highlight_trace])
    fig.update_layout(
        width=800, height=800,
        paper_bgcolor='white', plot_bgcolor='white',
        xaxis=dict(visible=False, scaleanchor='y', scaleratio=1),
        yaxis=dict(visible=False),
        margin=dict(l=10, r=10, t=10, b=10),
        showlegend=False,
        hoverdistance=5,
    )

    html = fig.to_html(include_plotlyjs='cdn', post_script=_JS_HOVER)
    Path(out_path).write_text(html)
    print(f"Interactive plot saved to: {out_path}")


@cli(ViT3DConfig)
def main(args) -> None:
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    dataset = OpenCell3DCropImageOnlyDataset(args, split_key=args.split_key)

    # Build gene → location mapping from the metadata CSV
    meta = pd.read_csv(
        os.path.join(args.data_path, f'{args.split_key}_3d_192_crop_merged_meta_data.csv')
    )
    gene_to_location = {}
    for _, row in meta.iterrows():
        locs_str = str(row['locations']) if not pd.isna(row['locations']) else 'NA'
        loc_list = [l.strip() for l in locs_str.split(',')]
        gene_to_location[row['gene_name']] = loc_list[0] if len(loc_list) == 1 else 'Multilocalizing'

    config = ViT3DConfig(**vars(args))
    config.num_classes = len(dataset.gene_names)
    model = ViT3DModel(config=config)
    model.to(device)
    model.eval()

    loader = DataLoader(
        dataset, batch_size=args.per_device_eval_batch_size, shuffle=False,
        num_workers=args.dataloader_num_workers, collate_fn=dataset.collate,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_embeddings = []
    all_labels = []
    all_genes = []
    sample_idx = 0

    with torch.no_grad():
        for batch in tqdm(loader, desc="Embedding"):
            data = batch['batched_data']
            protein_img = data['protein_img'].to(device)   # (B, 1, D, H, W)
            nucleus_img = data['nucleus_img'].to(device)   # (B, 1, D, H, W)

            input_img = torch.cat([nucleus_img, protein_img], dim=1)  # (B, 2, D, H, W)
            embedding = model.embed(input_img)

            all_embeddings.append(embedding.detach().cpu().numpy())

            B = protein_img.shape[0]
            for b in range(B):
                gene = dataset.gene_name_labels[sample_idx]
                all_genes.append(gene)
                all_labels.append(gene_to_location.get(gene, 'NA'))
                sample_idx += 1

    embedding_matrix = np.concatenate(all_embeddings, axis=0)
    assert embedding_matrix.shape[0] == len(all_labels)

    np.savez_compressed(
        output_dir / "embeddings.npz",
        embeddings=embedding_matrix,
        labels=np.array(all_labels),
        gene_names=np.array(all_genes),
    )
    print(f"Embeddings saved to {output_dir / 'embeddings.npz'}")

    get_scores(embedding_matrix, all_labels, method="kmeans")
    get_scores(embedding_matrix, all_labels, method="spectral")

    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.3)
    embedding_2d = reducer.fit_transform(embedding_matrix)

    # save 2D embedding and metadata for plotting
    np.savez_compressed(
        output_dir / "umap_embeddings.npz",
        embedding_2d=embedding_2d,
        labels=np.array(all_labels),
        gene_names=np.array(all_genes),
    )
    print(f"2D embeddings saved to {output_dir / 'umap_embeddings.npz'}")

    labels_set = set(all_labels)
    middle_labels = sorted([l for l in labels_set if l not in ("Multilocalizing", "NA")])
    unique_labels = (["Multilocalizing"] if "Multilocalizing" in labels_set else []) + \
                    middle_labels + (["NA"] if "NA" in labels_set else [])

    default_color = "#999999"
    labels_arr = np.array(all_labels)
    genes_arr = np.array(all_genes)
    mask_multi = labels_arr == "Multilocalizing"

    fig = plt.figure(figsize=(8, 8), dpi=300)
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[1.0, 0.22], wspace=0.02)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax_leg.axis('off')

    ax.set_box_aspect(1)
    divider = make_axes_locatable(ax)
    ax_leg = divider.append_axes("right", size="5%", pad=0.02)
    ax_leg.axis("off")

    if mask_multi.any():
        ax.scatter(embedding_2d[mask_multi, 0], embedding_2d[mask_multi, 1],
                   c=[mcolors.to_rgba(location_colors.get("Multilocalizing", default_color))],
                   s=2, marker='o', edgecolors='white', linewidths=0.2, alpha=1.0, zorder=0)

    colors_other = [mcolors.to_rgba(location_colors.get(l, default_color)) for l in labels_arr[~mask_multi]]
    ax.scatter(embedding_2d[~mask_multi, 0], embedding_2d[~mask_multi, 1],
               c=colors_other, s=4, marker='o', edgecolors='white', linewidths=0.2, alpha=1.0, zorder=2)

    for gene in np.unique(genes_arr):
        pts = embedding_2d[genes_arr == gene]
        cx, cy = pts[:, 0].mean(), pts[:, 1].mean()
        ax.text(cx, cy, gene, fontsize=1, fontweight='bold', ha='center', va='center',
                color='black', bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.5))

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    x0, x1 = 0.06, 0.25
    y0, y1 = 0.06, 0.25
    ax.plot([x0, x1], [y0, y0], transform=ax.transAxes, color="#000000", lw=2)
    ax.plot([x0, x0], [y0, y1], transform=ax.transAxes, color="#000000", lw=2)
    ax.text((x0 + x1) / 2, y0 + 0.01, "UMAP1",
            transform=ax.transAxes, ha='center', va='bottom', fontsize=10, color="#000000")
    ax.text(x0 + 0.01, (y0 + y1) / 2, "UMAP2",
            transform=ax.transAxes, ha='left', va='center', fontsize=10, color="#000000", rotation=90)

    ax_leg.legend(
        handles=[plt.Line2D([0], [0], marker='o', color='w', label=lbl,
                            markerfacecolor=location_colors.get(lbl, default_color), markersize=6)
                 for lbl in unique_labels],
        loc='center left', bbox_to_anchor=(0.0, 0.5),
        frameon=False, ncol=2, fontsize='small',
    )

    fig.tight_layout()
    output_path = output_dir / "umap_embedding.png"
    fig.savefig(output_path, dpi=1200, bbox_inches='tight', pad_inches=0)
    output_path_svg = output_path.with_suffix(".svg")
    fig.savefig(output_path_svg, bbox_inches='tight', pad_inches=0)
    plt.close(fig)
    print(f"UMAP plot saved to {output_path} and {output_path_svg}")

    make_interactive_plot(embedding_2d, all_labels, all_genes,
                          output_dir / "umap_embedding_interactive.html")


if __name__ == "__main__":
    main()
