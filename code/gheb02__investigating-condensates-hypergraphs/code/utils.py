import xgi
import numpy as np
from nestedness.nestedness_calculator import NestednessCalculator
from scipy import sparse, stats
from scipy.stats import pearsonr, spearmanr
import pandas as pd
import matplotlib.pyplot as plt
import ast
import networkx as nx
from networkx.algorithms import bipartite
import seaborn as sns
import re


RANDOM_SEED = 42
# ---- Hypergraph construction ----
MIN_PROTEINS_HYPERGRAPH = 3

# ---- Disorder features ----
DISORDER_THRESHOLDS = (0.50, 0.60, 0.75)

# ---- Sequence-composition residue sets ----
STICKERS = set("YFR")
IDRS = set("SGQN")
REACTIVE = set("YFRSGQN")

# ---- Statistics ----
ALPHA = 0.05
CORRECTION_METHOD = "holm"          # one of: "fdr_bh", "holm", "bonferroni"
FAMILY_COLS = ("Scenario",)           # tests are corrected jointly within each group
METHODS = ("Pearson", "Spearman")
TARGETS = ("Hyperdegree", "Projection Degree")

CORRECTION_LABELS = {
    "fdr_bh": "Benjamini-Hochberg",
    "holm": "Holm",
    "bonferroni": "Bonferroni",
}
P_ADJ_COLUMNS = {
    "fdr_bh": "p_fdr_bh",
    "holm": "p_holm",
    "bonferroni": "p_bonferroni",
}

# ---- Scenario names ----
SCENARIO_ALL = "All proteins"
SCENARIO_NO_DEG1 = "Excluding hyper-degree-1 proteins"
SCENARIO_NO_LARGEST = "Excluding largest condensate"

# ---- Feature names (columns of the protein feature table) ----
IUPRED_FEATURES = (
    ["Mean IUPred Score", "IUPred Score std"]
    + [f"Num. Disordered Regions_{round(t * 100)}" for t in DISORDER_THRESHOLDS]
    + [f"Max Disordered Length_{round(t * 100)}" for t in DISORDER_THRESHOLDS]
)
SEQ_FEATURES = ["Sticker Fraction", "IDR Fraction", "Reactive Fraction"]
ALL_FEATURES = IUPRED_FEATURES + SEQ_FEATURES

# ---- Display labels (insertion order = row order in figures and tables) ----
HEATMAP_LABELS = {
    "Mean IUPred Score": "Mean IUPred",
    "IUPred Score std": "IUPred SD",
    "Num. Disordered Regions_50": "Disordered regions",
    "Max Disordered Length_50": "Max disorder length",
    "Sticker Fraction": "Sticker fraction",
    "IDR Fraction": "IDR fraction",
    "Reactive Fraction": "Reactive fraction",
}
LATEX_LABELS = {
    "Num. Disordered Regions_50": "Num. disordered regions",
    "Max Disordered Length_50": "Max disordered length",
    "Mean IUPred Score": "Mean IUPred Score",
    "IUPred Score std": "IUPred Score std",
    "Sticker Fraction": "Sticker fraction",
    "IDR Fraction": "IDR fraction",
    "Reactive Fraction": "Reactive fraction",
}
TARGET_LABELS = {"Hyperdegree": "Hypergraph", "Projection Degree": "Projection"}


def describe_hypergraph(hypergraph):

	connected = xgi.is_connected(hypergraph)
	print(f"Number of proteins (nodes): {hypergraph.num_nodes}")
	print(f"Number of biomolecular condensates (hyperedges): {hypergraph.num_edges}")

	print(f"\nHypergraph incidence density: {100*xgi.incidence_density(hypergraph):.2f}%")
	H_int = xgi.convert_labels_to_integers(hypergraph)
	clust_coeffs = xgi.local_clustering_coefficient(H_int)
	print(f"Global clustering coefficient: {100*np.mean(list(clust_coeffs.values())):.2f}%")

	if connected:
		print('The hypergraph is connected.')
	else:
		print('The hypergraph is NOT connected.')
		print(f"The graph is composed of {xgi.number_connected_components(hypergraph)} connected components.")


def nodf_rows_cols(mat):
    
	if hasattr(mat, "toarray"):
		mat = mat.toarray()
		
	mat = np.asarray(mat, dtype=int)

	row_deg = mat.sum(axis=1)
	col_deg = mat.sum(axis=0)

	row_order = np.argsort(-row_deg)
	col_order = np.argsort(-col_deg)
	mat = mat[row_order][:, col_order]

	calc = NestednessCalculator(mat)

	row_sum = calc.get_paired_nestedness(mat, rows=True)
	col_sum = calc.get_paired_nestedness(mat, rows=False)

	n_rows, n_cols = mat.shape
	row_norm = n_rows * (n_rows - 1) / 2
	col_norm = n_cols * (n_cols - 1) / 2

	nodf_rows = row_sum / row_norm if row_norm > 0 else np.nan
	nodf_cols = col_sum / col_norm if col_norm > 0 else np.nan
	nodf_total = (row_sum + col_sum) / (row_norm + col_norm) if (row_norm + col_norm) > 0 else np.nan

	return nodf_rows, nodf_cols, nodf_total

rng = np.random.default_rng(RANDOM_SEED)


def as_binary_csr(mat):
    """Return a binary CSR matrix from a dense or sparse incidence matrix."""
    if sparse.issparse(mat):
        out = mat.astype(bool).astype(np.int32).tocsr()
    else:
        out = sparse.csr_matrix(np.asarray(mat, dtype=bool).astype(np.int32))
    out.eliminate_zeros()
    return out


def get_incidence_objects(H):
    """Build a binary incidence matrix plus ordered node/edge ID mappings."""
    inc, node_dict, edge_dict = xgi.incidence_matrix(
        H, index=True, sparse=True
    )
    inc = as_binary_csr(inc)

    def ordered_ids(mapping, n):
        # Case 1: index -> ID, e.g. {0: "protein_A", 1: "protein_B"}
        if all(
            isinstance(k, (int, np.integer))
            for k in mapping.keys()
        ):
            return np.array(
                [mapping[i] for i in range(n)],
                dtype=object
            )

        # Case 2: ID -> index, e.g. {"protein_A": 0, "protein_B": 1}
        ids = np.empty(n, dtype=object)
        for node_id, idx in mapping.items():
            ids[int(idx)] = node_id
        return ids

    node_ids = ordered_ids(node_dict, inc.shape[0])
    edge_ids = ordered_ids(edge_dict, inc.shape[1])

    return inc, node_ids, edge_ids


def sample_edge_size_preserving(edge_sizes, n_nodes, random_state):
    """Null matrix preserving each hyperedge cardinality, randomizing protein memberships."""
    rows = []
    cols = []
    for j, size in enumerate(np.asarray(edge_sizes, dtype=int)):
        chosen = random_state.choice(n_nodes, size=size, replace=False)
        rows.extend(chosen)
        cols.extend([j] * size)
    return sparse.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_nodes, len(edge_sizes)))


def sample_node_degree_preserving(node_degrees, n_edges, random_state):
    """Null matrix preserving each protein hyperdegree, randomizing condensate memberships."""
    rows = []
    cols = []
    for i, degree in enumerate(np.asarray(node_degrees, dtype=int)):
        if degree == 0:
            continue
        chosen = random_state.choice(n_edges, size=degree, replace=False)
        rows.extend([i] * degree)
        cols.extend(chosen)
    return sparse.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(len(node_degrees), n_edges))


def sample_fill_only(n_nodes, n_edges, n_incidences, random_state):
    """Null matrix preserving only matrix dimensions and total number of incidences."""
    positions = random_state.choice(n_nodes * n_edges, size=n_incidences, replace=False)
    rows = positions // n_edges
    cols = positions % n_edges
    return sparse.csr_matrix((np.ones(n_incidences, dtype=np.int8), (rows, cols)), shape=(n_nodes, n_edges))


def sample_fixed_fixed_swap(mat, n_swaps, random_state, attempts_per_swap=100):
    """Approximate fixed-row/fixed-column null using bipartite double-edge swaps."""
    mat = as_binary_csr(mat)
    coo = mat.tocoo()
    rows = coo.row.copy()
    cols = coo.col.copy()
    row_sets = [set(mat.indices[mat.indptr[i]:mat.indptr[i + 1]]) for i in range(mat.shape[0])]
    n_ones = len(rows)

    if n_ones < 2:
        return mat.copy()

    accepted = 0
    max_attempts = max(attempts_per_swap * n_swaps, 1)
    attempts = 0

    while accepted < n_swaps and attempts < max_attempts:
        attempts += 1
        p, q = random_state.choice(n_ones, size=2, replace=False)
        r1, c1 = rows[p], cols[p]
        r2, c2 = rows[q], cols[q]

        if r1 == r2 or c1 == c2:
            continue
        if (c2 in row_sets[r1]) or (c1 in row_sets[r2]):
            continue

        row_sets[r1].remove(c1)
        row_sets[r2].remove(c2)
        row_sets[r1].add(c2)
        row_sets[r2].add(c1)
        cols[p], cols[q] = c2, c1
        accepted += 1

    if accepted < n_swaps:
        print(f"Warning: accepted {accepted:,} of {n_swaps:,} requested swaps; null may be less randomized than intended.")

    return sparse.csr_matrix((np.ones(n_ones, dtype=np.int8), (rows, cols)), shape=mat.shape)


def ccdf(values, support):
    values = np.asarray(values)
    return np.array([(values > x).mean() for x in support])


def distribution_null_ensemble(inc, n_null=200, random_state=None):
    """Generate distribution null ensembles for node hyperdegree and edge cardinality."""
    random_state = np.random.default_rng(random_state)
    inc = as_binary_csr(inc)
    n_nodes, n_edges = inc.shape
    node_degrees = np.asarray(inc.sum(axis=1)).ravel().astype(int)
    edge_sizes = np.asarray(inc.sum(axis=0)).ravel().astype(int)
    n_incidences = int(edge_sizes.sum())

    null_node_degree_edge_size = []
    null_edge_size_node_degree = []
    null_node_degree_fill = []
    null_edge_size_fill = []

    for _ in range(n_null):
        m_edge = sample_edge_size_preserving(edge_sizes, n_nodes, random_state)
        m_node = sample_node_degree_preserving(node_degrees, n_edges, random_state)
        m_fill = sample_fill_only(n_nodes, n_edges, n_incidences, random_state)

        null_node_degree_edge_size.append(np.asarray(m_edge.sum(axis=1)).ravel())
        null_edge_size_node_degree.append(np.asarray(m_node.sum(axis=0)).ravel())
        null_node_degree_fill.append(np.asarray(m_fill.sum(axis=1)).ravel())
        null_edge_size_fill.append(np.asarray(m_fill.sum(axis=0)).ravel())

    return {
        "observed_node_degree": node_degrees,
        "observed_edge_size": edge_sizes,
        "node_degree_edge_size_null": null_node_degree_edge_size,
        "edge_size_node_degree_null": null_edge_size_node_degree,
        "node_degree_fill_null": null_node_degree_fill,
        "edge_size_fill_null": null_edge_size_fill,
    }


def plot_distribution_nulls(dist_results, plot_fill_null = False, colors = ['steelblue', 'purple'], figsize=(10, 6)):
	"""Plot observed CCDFs against 95% pointwise null envelopes."""
	obs_deg = dist_results["observed_node_degree"]
	obs_size = dist_results["observed_edge_size"]

	deg_support = np.arange(1, max(obs_deg.max(), max(x.max() for x in dist_results["node_degree_edge_size_null"])) + 1)
	size_support = np.arange(1, max(obs_size.max(), max(x.max() for x in dist_results["edge_size_node_degree_null"])) + 1)

	deg_null = np.vstack([ccdf(x, deg_support) for x in dist_results["node_degree_edge_size_null"]])
	deg_fill = np.vstack([ccdf(x, deg_support) for x in dist_results["node_degree_fill_null"]])
	size_null = np.vstack([ccdf(x, size_support) for x in dist_results["edge_size_node_degree_null"]])
	size_fill = np.vstack([ccdf(x, size_support) for x in dist_results["edge_size_fill_null"]])

	fig, axes = plt.subplots(1, 2, figsize=figsize)

	axes[0].plot(deg_support, ccdf(obs_deg, deg_support), color="black", label="Observed")
	axes[0].fill_between(deg_support, *np.quantile(deg_null, [0.025, 0.975], axis=0), color=colors[0], alpha=0.25, label="Edge-size null 95% envelope")
	axes[0].plot(deg_support, deg_null.mean(axis=0), color=colors[0], ls="--", label="Edge-size null mean")
	if plot_fill_null:
		axes[0].fill_between(deg_support, *np.quantile(deg_fill, [0.025, 0.975], axis=0), color="tab:gray", alpha=0.18, label="Fill-only null 95% envelope")
	axes[0].set_xscale("log")
	axes[0].set_yscale("log")
	axes[0].set_xlabel("Node Degree")
	axes[0].set_ylabel("P($K \\geq k$)")
	axes[0].set_title("Node Degree vs Null Models")
	axes[0].legend(loc='lower left')

	axes[1].plot(size_support, ccdf(obs_size, size_support), color="black", label="Observed")
	axes[1].fill_between(size_support, *np.quantile(size_null, [0.025, 0.975], axis=0), color=colors[1], alpha=0.25, label="Node-degree null 95% envelope")
	axes[1].plot(size_support, size_null.mean(axis=0), color=colors[1], ls="--", label="Node-degree null mean")
	if plot_fill_null:
		axes[1].fill_between(size_support, *np.quantile(size_fill, [0.025, 0.975], axis=0), color="tab:gray", alpha=0.18, label="Fill-only null 95% envelope")
	axes[1].set_xscale("log")
	axes[1].set_yscale("log")
	axes[1].set_xlabel("Hyper-Edge Cardinality")
	axes[1].set_ylabel("P($S \\geq s$)")
	axes[1].set_title("Hyper-Edge Cardinality vs Null Models")
	axes[1].legend(loc='lower left')
	plt.tight_layout()

	return fig, axes

def summarize_distribution_nulls(dist_results):
    """Return compact observed-vs-null summaries for distributional shape."""
    def gini(x):
        x = np.sort(np.asarray(x, dtype=float))
        if len(x) == 0 or x.sum() == 0:
            return np.nan
        idx = np.arange(1, len(x) + 1)
        return (2 * np.sum(idx * x) / (len(x) * x.sum())) - ((len(x) + 1) / len(x))

    metrics = {
        "mean": np.mean,
        "std": np.std,
        "max": np.max,
        "gini": gini,
    }
    rows = []
    specs = [
        ("node_hyperdegree", dist_results["observed_node_degree"], dist_results["node_degree_edge_size_null"], "edge-size preserving"),
        ("node_hyperdegree", dist_results["observed_node_degree"], dist_results["node_degree_fill_null"], "fill only"),
        ("edge_cardinality", dist_results["observed_edge_size"], dist_results["edge_size_node_degree_null"], "node-degree preserving"),
        ("edge_cardinality", dist_results["observed_edge_size"], dist_results["edge_size_fill_null"], "fill only"),
    ]
    for target, observed, null_values, null_name in specs:
        for metric_name, fn in metrics.items():
            obs = fn(observed)
            null_metric = np.array([fn(x) for x in null_values])
            p_upper = (1 + np.sum(null_metric >= obs)) / (len(null_metric) + 1)
            p_lower = (1 + np.sum(null_metric <= obs)) / (len(null_metric) + 1)
            rows.append({
                "target": target,
                "null_model": null_name,
                "metric": metric_name,
                "observed": obs,
                "null_mean": null_metric.mean(),
                "null_q025": np.quantile(null_metric, 0.025),
                "null_q975": np.quantile(null_metric, 0.975),
                "empirical_p_two_sided": min(1.0, 2 * min(p_lower, p_upper)),
            })
    return pd.DataFrame(rows)


def nodf_sparse(mat):
    """Compute row, column, and total NODF fractions for a binary incidence matrix.

    Equal-degree row/column pairs contribute zero. Values are returned as fractions in [0, 1].
    """
    mat = as_binary_csr(mat)
    n_rows, n_cols = mat.shape

    def paired_nodf(block):
        deg = np.asarray(block.sum(axis=1)).ravel().astype(float)
        overlap = (block @ block.T).tocoo()
        mask = overlap.row < overlap.col
        r = overlap.row[mask]
        c = overlap.col[mask]
        data = overlap.data[mask].astype(float)
        deg_r = deg[r]
        deg_c = deg[c]
        comparable = (deg_r != deg_c) & (np.minimum(deg_r, deg_c) > 0)
        if comparable.sum() == 0:
            return 0.0
        return np.sum(data[comparable] / np.minimum(deg_r[comparable], deg_c[comparable]))

    row_sum = paired_nodf(mat)
    col_sum = paired_nodf(mat.T.tocsr())
    row_norm = n_rows * (n_rows - 1) / 2
    col_norm = n_cols * (n_cols - 1) / 2

    nodf_rows = row_sum / row_norm if row_norm > 0 else np.nan
    nodf_cols = col_sum / col_norm if col_norm > 0 else np.nan
    nodf_total = (row_sum + col_sum) / (row_norm + col_norm) if (row_norm + col_norm) > 0 else np.nan
    return nodf_rows, nodf_cols, nodf_total


def nodf_fixed_fixed_null(inc, n_null=100, swaps_per_incidence=10, random_state=None, attempts_per_swap=100):
    """Compare observed NODF with fixed-row/fixed-column swapped null matrices."""
    random_state = np.random.default_rng(random_state)
    inc = as_binary_csr(inc)
    n_swaps = int(swaps_per_incidence * inc.nnz)
    observed = np.array(nodf_sparse(inc))
    null_values = []
    for _ in range(n_null):
        null_mat = sample_fixed_fixed_swap(inc, n_swaps=n_swaps, random_state=random_state, attempts_per_swap=attempts_per_swap)
        null_values.append(nodf_sparse(null_mat))
    null_values = np.asarray(null_values)

    rows = []
    for i, label in enumerate(["rows", "columns", "total"]):
        p_upper = (1 + np.sum(null_values[:, i] >= observed[i])) / (n_null + 1)
        p_lower = (1 + np.sum(null_values[:, i] <= observed[i])) / (n_null + 1)
        null_mean = null_values[:, i].mean()
        null_sd = null_values[:, i].std(ddof=1)
        rows.append({
            "component": label,
            "observed": observed[i],
            "null_mean": null_mean,
            "null_sd": null_sd,
            "null_q025": np.quantile(null_values[:, i], 0.025),
            "null_q975": np.quantile(null_values[:, i], 0.975),
            "z_score": (observed[i] - null_mean) / null_sd if null_sd > 0 else np.nan,
            "empirical_p_two_sided": min(1.0, 2 * min(p_lower, p_upper)),
        })
    return pd.DataFrame(rows), observed, null_values


def plot_nodf_null(nodf_summary, nodf_null_values):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=False)
    labels = ["rows", "columns", "total"]
    titles = {"rows": "Proteins", "columns":"Condensates", "total":"total"}
    for i, label in enumerate(labels):
        obs = nodf_summary.loc[nodf_summary["component"] == label, "observed"].iloc[0]
        axes[i].hist(100 * nodf_null_values[:, i], bins=25, color="lightgray", edgecolor="white")
        axes[i].axvline(100 * obs, color="crimson", lw=2, label="Observed")
        axes[i].set_title(f"NODF {titles[label]}")
        axes[i].set_xlabel("NODF (%)")
        axes[i].set_ylabel("Null count")
        axes[i].legend()
    plt.tight_layout()
    plt.show()


def bh_fdr(p_values):
    """Benjamini-Hochberg FDR correction."""
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    n = len(p_values)
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for i in range(n - 1, -1, -1):
        cumulative = min(cumulative, ranked[i] * n / (i + 1))
        adjusted[order[i]] = cumulative
    return np.minimum(adjusted, 1.0)


def headline_stats(H_, inc_, fit_deg_, fit_card_):
    nodf_r_, nodf_c_, nodf_t_ = nodf_sparse(inc_)
    return {
        "proteins (nodes)": H_.num_nodes,
        "condensates (hyperedges)": H_.num_edges,
        "incidence density (%)": 100*xgi.incidence_density(H_),
        "mean clustering coefficient (%)": 100 * np.mean(list(xgi.clustering_coefficient(H_).values())),
        "connected components": xgi.number_connected_components(H_),
        "degree gamma": fit_deg_.power_law.alpha,
        "degree xmin": fit_deg_.power_law.xmin,
        "cardinality gamma": fit_card_.power_law.alpha,
        "cardinality xmin": fit_card_.power_law.xmin,
        "NODF rows (%)": 100 * nodf_r_,
        "NODF columns (%)": 100 * nodf_c_,
        "NODF total (%)": 100 * nodf_t_,
    }


def parse_scores(value):
    """IUPred scores are stored as stringified lists in the CSV; return a real list."""
    if isinstance(value, str):
        return ast.literal_eval(value)
    return value


def longest_run(arr, threshold=0.5):
    """Length of the longest stretch of consecutive values above `threshold`."""
    max_run = run = 0
    for v in arr:
        if v > threshold:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return max_run


def count_regions(arr, threshold=0.5):
    """Number of distinct stretches of consecutive values above `threshold`."""
    count = 0
    in_region = False
    for v in arr:
        if v > threshold and not in_region:
            count += 1
            in_region = True
        elif v <= threshold:
            in_region = False
    return count


def residue_fraction(seq, residues):
    """Fraction of residues in `seq` that belong to the set `residues` (NaN if no sequence)."""
    if not isinstance(seq, str) or not seq:
        return np.nan
    return sum(aa in residues for aa in seq) / len(seq)

def build_bipartite_graph(edges):
    """Bipartite protein-condensate graph from a table with `uniprot_id` and `condensate_id`.

    Proteins get bipartite=0 and condensates bipartite=1.
    Returns the graph and the list of protein nodes.
    """
    graph = nx.Graph()
    protein_nodes = list(edges["uniprot_id"].unique())
    graph.add_nodes_from(protein_nodes, bipartite=0)
    graph.add_nodes_from(edges["condensate_id"].unique(), bipartite=1)
    graph.add_edges_from(zip(edges["uniprot_id"], edges["condensate_id"]))
    return graph, protein_nodes


def compute_degrees(edges):
    """Build the bipartite graph and its protein projection from an edge table.

    Returns (degrees, bipartite_graph, projection), where `degrees` is a DataFrame
    indexed by uniprot_id with columns 'Projection Degree' and 'Hyperdegree'.
    """
    graph, protein_nodes = build_bipartite_graph(edges)
    projection = bipartite.projected_graph(graph, protein_nodes)
    degrees = pd.DataFrame({
        "Projection Degree": pd.Series(dict(projection.degree())),
        "Hyperdegree": pd.Series(dict(graph.degree(protein_nodes))),
    })
    degrees.index.name = "uniprot_id"
    return degrees, graph, projection


def describe_projected_graph(graph):
    """Print basic descriptive statistics of a projected graph."""
    print(f"Number of nodes in the projection: {graph.number_of_nodes()}")
    print(f"Number of edges in the projection: {graph.number_of_edges()}")
    print(f"\nGraph density: {nx.density(graph)}")
    print(f"Graph global clustering coefficient: {nx.average_clustering(graph)}")
    print(f"Graph degree assortativity: {nx.degree_assortativity_coefficient(graph)}")


def p_to_stars(p):
    """Significance stars for a (adjusted) p-value."""
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""



def correlate(x, y):
    """Pearson and Spearman correlation of two aligned Series (pairwise-complete cases).

    Returns {"n": n, "Pearson": (r, p), "Spearman": (r, p)}; NaNs if the
    correlation is undefined (fewer than 3 points or a constant input).
    """
    mask = x.notna() & y.notna()
    xv = x[mask].to_numpy(dtype=float)
    yv = y[mask].to_numpy(dtype=float)
    out = {"n": len(xv), "Pearson": (np.nan, np.nan), "Spearman": (np.nan, np.nan)}
    if len(xv) >= 3 and np.ptp(xv) > 0 and np.ptp(yv) > 0:
        pear_r, pear_p = pearsonr(xv, yv)
        spear_r, spear_p = spearmanr(xv, yv)
        out["Pearson"] = (pear_r, pear_p)
        out["Spearman"] = (spear_r, spear_p)
    return out


def adjust_pvalues(p_values, method):
    """Adjust p-values for multiple testing (NaNs are ignored and returned as NaN).

    method: "bonferroni" (FWER), "holm" (FWER, step-down) or "fdr_bh" (Benjamini-Hochberg FDR).
    """
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan)
    valid = ~np.isnan(p)
    pv = p[valid]
    m = len(pv)
    if m == 0:
        return adjusted

    order = np.argsort(pv)
    ranks = np.arange(1, m + 1)
    if method == "bonferroni":
        adj = np.minimum(pv * m, 1.0)
    elif method == "holm":
        adj_sorted = np.maximum.accumulate((m - ranks + 1) * pv[order])
        adj = np.empty(m)
        adj[order] = np.minimum(adj_sorted, 1.0)
    elif method == "fdr_bh":
        scaled = pv[order] * m / ranks
        adj_sorted = np.minimum.accumulate(scaled[::-1])[::-1]
        adj = np.empty(m)
        adj[order] = np.minimum(adj_sorted, 1.0)
    else:
        raise ValueError(f"Unknown correction method: {method}")

    adjusted[valid] = adj
    return adjusted


def add_adjusted_pvalues(results, family_cols=FAMILY_COLS, method=CORRECTION_METHOD, alpha=ALPHA):
    """Add Bonferroni, Holm and BH-FDR adjusted p-values to a long-format results table.

    Tests are corrected jointly within each group defined by `family_cols`.
    `p_adj` holds the chosen method and `significant` flags p_adj < alpha.
    """
    results = results.copy()
    for name, column in P_ADJ_COLUMNS.items():
        results[column] = results.groupby(list(family_cols))["p"].transform(
            lambda s, name=name: pd.Series(adjust_pvalues(s.to_numpy(), name), index=s.index)
        )
    results["p_adj"] = results[P_ADJ_COLUMNS[method]]
    results["significant"] = results["p_adj"] < alpha
    return results


def correlation_table(data, features, targets, scenario):
    """Long-format table of Pearson and Spearman correlations of each feature with each target."""
    rows = []
    for feature in features:
        for target in targets:
            res = correlate(data[feature], data[target])
            for method in METHODS:
                r, p = res[method]
                rows.append({
                    "Scenario": scenario,
                    "Feature": feature,
                    "Target": target,
                    "Method": method,
                    "r": r,
                    "p": p,
                    "n": res["n"],
                })
    return pd.DataFrame(rows)


def plot_correlation_heatmap(results, scenario, label_map=HEATMAP_LABELS):
    """Heatmap of Pearson and Spearman correlations for one scenario.

    Cell colours show r; stars show significance of the multiple-testing-adjusted p-value.
    """
    sub = results[(results["Scenario"] == scenario) & results["Feature"].isin(label_map)]
    columns = [(m, t) for m in METHODS for t in TARGETS]
    col_names = [f"{m}\n{TARGET_LABELS[t]}" for m, t in columns]

    r_mat = pd.DataFrame(index=list(label_map), columns=col_names, dtype=float)
    annot = pd.DataFrame(index=list(label_map), columns=col_names, dtype=object)
    for feature in label_map:
        for (method, target), name in zip(columns, col_names):
            row = sub[(sub["Feature"] == feature) & (sub["Method"] == method) & (sub["Target"] == target)]
            r = row["r"].iloc[0]
            r_mat.loc[feature, name] = r
            annot.loc[feature, name] = f"{r:.2f}{p_to_stars(row['p_adj'].iloc[0])}"
    r_mat.index = [label_map[f] for f in r_mat.index]
    annot.index = r_mat.index

    fig, ax = plt.subplots(figsize=(7.2, 5.5))
    sns.heatmap(
        r_mat,
        annot=annot,
        fmt="",
        cmap="coolwarm",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.5,
        linecolor="white",
        annot_kws={"size": 8},
        cbar_kws={"label": "Correlation", "shrink": 0.8},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=10)
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
    fig.suptitle(
        f"Sequence features vs. protein connectivity\n{scenario} ($n = {int(sub['n'].max())}$)",
        fontsize=12,
    )
    plt.tight_layout(rect=(0, 0.03, 1, 1))
    plt.show()


def scenario_slug(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def make_latex_table(results, scenario, label_map=LATEX_LABELS):
    """LaTeX table (Pearson and Spearman x projection degree and hyperdegree) for one scenario."""
    sub = results[results["Scenario"] == scenario]
    n = int(sub["n"].max())
    columns = [
        ("Pearson", "Projection Degree"),
        ("Pearson", "Hyperdegree"),
        ("Spearman", "Projection Degree"),
        ("Spearman", "Hyperdegree"),
    ]

    def cell(feature, method, target):
        row = sub[(sub["Feature"] == feature) & (sub["Method"] == method) & (sub["Target"] == target)].iloc[0]
        return f"{row['r']:.3f}{p_to_stars(row['p_adj'])}"

    caption = (
        "Correlation between sequence-derived features and node connectivity measures "
        f"({scenario[0].lower() + scenario[1:]}, $n = {n}$)."
    )
    note = (
        "\\multicolumn{5}{l}{\\footnotesize \\textit{Note:} *** $p_{\\mathrm{adj}} < 0.001$, "
        "** $p_{\\mathrm{adj}} < 0.01$, * $p_{\\mathrm{adj}} < 0.05$ "
        f"({CORRECTION_LABELS[CORRECTION_METHOD]} correction).}} \\\\"
    )
    lines = [
        "\\begin{table}",
        f"\\caption{{{caption}}}",
        f"\\label{{tab:correlations_{scenario_slug(scenario)}}}",
        "\\begin{tabular}{lcccc}",
        "\\toprule",
        "& Projection Degree & Hyperdegree & Projection Degree & Hyperdegree \\\\",
        "Feature & (Pearson) & (Pearson) & (Spearman) & (Spearman) \\\\",
        "\\midrule",
    ]
    for feature, label in label_map.items():
        cells = [cell(feature, m, t) for m, t in columns]
        lines.append(f"{label} & " + " & ".join(cells) + " \\\\")
    lines += ["\\bottomrule", note, "\\end{tabular}", "\\end{table}"]
    return "\n".join(lines)