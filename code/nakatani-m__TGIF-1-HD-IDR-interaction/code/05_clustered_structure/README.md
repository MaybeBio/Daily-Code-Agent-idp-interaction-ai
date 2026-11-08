# 05_clustered_structure

This directory contains representative structures extracted from the major conformational clusters identified by DBSCAN clustering in PCA space.

These structures were used for structural analysis, visualization, and discussion of HD-IDR interactions.

---

## Structure Selection Procedure

Representative structures were selected using the following workflow:

```text
mD-VcMD trajectories
        ↓
PCA
        ↓
DBSCAN clustering
        ↓
Identification of cluster centroids
        ↓
Selection of structures nearest to the centroid
        ↓
Representative structure analysis
```

---

## Directory Structure

```text
05_clustered_structure/
├── pIDR-HD/
└── non-pIDR-HD/
```

---

## pIDR-HD/

Representative structures from the phosphorylated model.

Included clusters:

### C1_centroid_nearest10/

Contains the 10 structures closest to the centroid of cluster C1.

### C2_centroid_nearest10/

Contains the 10 structures closest to the centroid of cluster C2.

These clusters were selected because they represent the major conformational states observed in the phosphorylated system.

---

## non-pIDR-HD/

Representative structures from the non-phosphorylated model.

Included clusters:

### C1_centroid_nearest10/

Contains the 10 structures closest to the centroid of cluster C1.

This cluster was selected because it represents the dominant conformational state observed in the non-phosphorylated system.

---

## Representative Structures

For each cluster:

```text
Cluster
    ↓
Cluster centroid
    ↓
Distance calculation
    ↓
10 nearest structures
```

The selected structures serve as representative conformations of the corresponding cluster.

These structures were used for:

- structural comparison between pIDR-HD and non-pIDR-HD systems
- visualization of HD-IDR interactions
- identification of characteristic intermolecular contacts
- interpretation of cluster-specific conformational features

---

## Relationship with DBSCAN Analysis

The clusters provided in this directory were generated using the DBSCAN analyses described in:

```text
04_analysis/DBSCAN/
```

Cluster labels correspond to those assigned during DBSCAN analysis.

For example:

- C1: largest cluster
- C2: second largest cluster

---

## Notes

Each cluster directory contains the ten structures closest to the cluster centroid rather than all structures assigned to the cluster.

These representative structures were selected to facilitate structural interpretation of the major conformational states identified during the simulations.

Only clusters discussed in this study are included in this repository.
