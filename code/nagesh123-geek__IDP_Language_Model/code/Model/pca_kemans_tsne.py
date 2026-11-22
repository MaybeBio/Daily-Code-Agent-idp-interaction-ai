import os
import json
import numpy as np
import matplotlib.pyplot as plt
import umap

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score
)

# ============================================================
# INPUT
# ============================================================

LOCAL_EMB= "../seq_500/local_embeddings.json"
GLOBAL_EMB = "../seq_500/global_embeddings.json"
#NARDINI_EMB="features_len_le_250.json"

# ============================================================
# OUTPUT
# ============================================================

SAVE_DIR = "embedding_analysis_500_Local"
os.makedirs(SAVE_DIR, exist_ok=True)

# ============================================================
# LOAD EMBEDDINGS
# ============================================================

with open(LOCAL_EMB) as f:
    X = np.array(json.load(f),dtype=np.float32)

print("Original shape of local embedding:", X.shape)

with open(GLOBAL_EMB) as f:
    Y = np.array(json.load(f),dtype=np.float32)
print("Original shape of global embedding:",Y.shape)

#with open(NARDINI_EMB) as f :
#    N = np.array(json.load(f),dtype=np.float32)
x = X
#x = np.concatenate([X, Y], axis=1)
print("Original shape of concatenated embedding:",x.shape)


# ============================================================
# STANDARDIZE
# ============================================================

print("\nStandardizing embeddings...")

scaler = StandardScaler()

X_scaled = scaler.fit_transform(x)

# ============================================================
# PCA
# ============================================================

print("\nRunning PCA...")

pca = PCA(n_components=0.95,random_state=42)

X_pca = pca.fit_transform(X_scaled)

explained = np.sum(pca.explained_variance_ratio_)

print("PCA shape:",X_pca.shape)

print("Explained variance:",explained)

np.save(os.path.join(SAVE_DIR,"pca_embeddings.npy"),X_pca)

# ============================================================
# KMEANS
# ============================================================

print("\nRunning KMeans...")

kmeans = KMeans(n_clusters=20,random_state=42,n_init=20)

labels = kmeans.fit_predict(X_pca)

np.save(os.path.join(SAVE_DIR,"kmeans_labels.npy"),labels)

np.savetxt(os.path.join(SAVE_DIR,"kmeans_labels.txt"),labels,fmt="%d")

np.save(os.path.join(SAVE_DIR,"cluster_centers.npy"),kmeans.cluster_centers_)

print("KMeans done.")

# ============================================================
# CLUSTER SCORES
# ============================================================

print("\nComputing clustering metrics...")

sil = silhouette_score(X_pca,labels)

ch = calinski_harabasz_score(X_pca,labels)

db = davies_bouldin_score(X_pca,labels)

print(f"Silhouette       : {sil:.6f}")
print(f"Calinski-Harabasz: {ch:.6f}")
print(f"Davies-Bouldin   : {db:.6f}")

with open(os.path.join(SAVE_DIR,"clustering_scores.txt"),"w") as f:

    f.write(f"PCA dimensions     : {X_pca.shape[1]}\n")

    f.write(f"Explained variance : {explained:.6f}\n")

    f.write(f"K                  : 20\n")

    f.write(f"Silhouette         : {sil:.6f}\n")

    f.write(f"Calinski-Harabasz  : {ch:.6f}\n")

    f.write(f"Davies-Bouldin     : {db:.6f}\n")

# ============================================================
# UMAP
# ============================================================

print("\nRunning UMAP...")

umap_model = umap.UMAP(
    n_components=2,
    n_neighbors=30,
    min_dist=0.0,
    metric="cosine",
    random_state=42
)

X_umap = umap_model.fit_transform(X_pca)

np.save(os.path.join(SAVE_DIR,"umap_2d.npy"),X_umap)

print("UMAP done.")

# ============================================================
# UMAP PLOT
# ============================================================

plt.figure(figsize=(8, 7))

plt.scatter(X_umap[:, 0],X_umap[:, 1],c=labels,cmap="tab20",s=6,alpha=0.7)

plt.title("PCA → KMeans → UMAP")

plt.xlabel("UMAP 1")

plt.ylabel("UMAP 2")

plt.grid(True)

plt.tight_layout()

plt.savefig(
    os.path.join(SAVE_DIR,"umap_kmeans.png"),dpi=300)

plt.close()

print("UMAP plot saved.")

# ============================================================
# TSNE
# ============================================================

print("\nRunning t-SNE...")

tsne = TSNE(
    n_components=2,
    perplexity=30,
    init="pca",
    learning_rate="auto",
    random_state=42
)

X_tsne = tsne.fit_transform(X_pca)

np.save(os.path.join(SAVE_DIR,"tsne_2d.npy"),X_tsne)

print("t-SNE done.")

# ============================================================
# TSNE PLOT
# ============================================================

plt.figure(figsize=(8, 7))

plt.scatter( X_tsne[:, 0],X_tsne[:, 1],c=labels,cmap="tab20",s=6,alpha=0.7)

plt.title("PCA → KMeans → t-SNE")

plt.xlabel("t-SNE 1")

plt.ylabel("t-SNE 2")

plt.grid(True)

plt.tight_layout()

plt.savefig(os.path.join(SAVE_DIR,"tsne_kmeans.png"),dpi=300)

plt.close()

print("t-SNE plot saved.")

# ============================================================
# SAVE SUMMARY
# ============================================================

summary_file = os.path.join(SAVE_DIR,"summary.txt")

with open(summary_file, "w") as f:

    f.write(f"Original shape      : {X.shape}\n")

    f.write(f"PCA shape           : {X_pca.shape}\n")

    f.write(f"Explained variance  : {explained:.6f}\n")

    f.write(f"Clusters            : 20\n")

    f.write(f"Silhouette          : {sil:.6f}\n")

    f.write(f"Calinski-Harabasz   : {ch:.6f}\n")

    f.write(f"Davies-Bouldin      : {db:.6f}\n")

print("\n========================================")
print("Analysis complete.")
print("Results saved to:", SAVE_DIR)
print("========================================")
