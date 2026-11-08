# Data Pipeline

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/data/__init__.py](src/data/__init__.py)
- [src/data/balancing.py](src/data/balancing.py)
- [src/data/dataset.py](src/data/dataset.py)

</details>



## Purpose and Scope

The data pipeline in SurPass is responsible for transforming raw protein structural data, protein Language Model (pLM) embeddings, and PPI cluster annotations into model-ready batches. It manages the ingestion of paired chain features, spatial or contiguous cropping for memory efficiency, distance binning for distogram prediction, negative sampling for contrastive/binary learning, dataset balancing across multiple data sources, and collation into variable-length concatenated pair tensors.

This page provides a high-level parent overview of the data pipeline sub-systems. For deep technical implementation details, refer to the individual child sections:
- [PepoTrainDataset and Feature Loading](#2.1)
- [Cropping, Labels and Collation](#2.2)
- [Multi-Source Balancing and Splits](#2.3)
- [Dataset Preparation Scripts](#2.4)

Sources: [src/data/dataset.py:26-167](), [src/data/balancing.py:1-106]()

---

## Data Pipeline Architecture Overview

The pipeline integrates multiple stages from raw dataset preparation scripts to batched training inputs. The diagram below bridges the natural language conceptual flow to the core codebase entities involved in dataset construction and batching.

```mermaid
graph TD
    sub-DataIngestion["Data Ingestion & Loading"]
      Prep["Scripts & Preparation (scripts/)"] --> TSV["Cluster TSV Metadata (src/data/dataset.py)"]
      TSV --> PTD["PepoTrainDataset (src/data/dataset.py:26-40)"]
    end

    sub-PairSampling["Pair Sampling & Cropping"]
      PTD -->|Positive Pair| Pos["_positive_sample_from_row (src/data/dataset.py:108-134)"]
      PTD -->|Negative Pair| Neg["get_negative_item / get_static_negative_item (src/data/dataset.py:135-149)"]
      Pos --> Crop["crop_positive_pair / crop_negative_pair (src/data/cropping.py)"]
      Neg --> Crop
    end

    sub-LabelsAndBatching["Labels & Collation"]
      Crop --> Labels["distance_to_bins / build_negative_pair_labels (src/data/labels.py)"]
      Labels --> Collate["collate_concatenated_pairs (src/data/collate.py)"]
      Collate --> Sampler["DistributedWeightedSampler (src/data/balancing.py:131-182)"]
    end

    Collate --> Model["Model Training / Inference (src/models/)"]

    classDef default font-family:sans-serif;
```
*Figure 2.1: High-level mapping from dataset preparation and ingestion to collation and model consumption.*

Sources: [src/data/dataset.py:26-167](), [src/data/balancing.py:131-182](), [src/data/__init__.py:1-25]()

---

## 2.1 PepoTrainDataset and Feature Loading

The core dataset class `PepoTrainDataset` (`[src/data/dataset.py:26]`) indexes protein-protein interaction (PPI) clusters from a metadata TSV file (`[src/data/dataset.py:62-76]`). It coordinates positive sampling from true interacting pairs within the same cluster and negative sampling across distinct clusters (`[src/data/dataset.py:81-149]`).

Chain features (such as pLM embeddings and residue positions) are loaded from `.pt` files and preprocessed via `_process_complex` before pairing. Distance matrices and labels are dynamically fetched or calculated to prepare training items.

For full implementation details, see [PepoTrainDataset and Feature Loading](#2.1).

Sources: [src/data/dataset.py:26-149]()

---

## 2.2 Cropping, Labels and Collation

Because protein complexes vary dramatically in sequence length, the pipeline applies spatial or contiguous cropping via `src/data/cropping.py` to bound memory consumption (`[src/data/dataset.py:114-123]`). Ground-truth inter-chain and intra-chain distances are converted into discrete distance bins (`[src/data/dataset.py:127-132]`) defined in `src/data/labels.py`.

Finally, variable-sized residue pairs are batched using specialized collation routines like `collate_concatenated_pairs` (`[src/data/collate.py]`), handling padding and mask generation for PyTorch data loaders.

For full implementation details, see [Cropping, Labels and Collation](#2.2).

Sources: [src/data/dataset.py:114-133](), [src/data/collate.py:1-9]()

---

## 2.3 Multi-Source Balancing and Splits

To handle multiple data sources and prevent class imbalance between positive interactions and negative non-interacting pairs, SurPass provides dataset wrappers and samplers (`[src/data/balancing.py]`). 

`BalancedClusterDataset` (`[src/data/balancing.py:33]`) pairs positive items with a configurable `negative_ratio` (`[src/data/balancing.py:39]`). Additionally, `DistributedWeightedSampler` (`[src/data/balancing.py:131]`) facilitates distributed training (DDP) by drawing replacement samples according to source-specific dataset weights across ranks (`[src/data/balancing.py:172-181]`).

```mermaid
graph TD
    sub-BalancingAndSplits["Balancing & Splitting Module (src/data/balancing.py)"]
      Split["split_dataset (src/data/balancing.py:13-30)"] --> TrainSet["Train Subset"]
      Split --> TestSet["Test Subset"]
      TrainSet --> Bal["BalancedClusterDataset (src/data/balancing.py:33-64)"]
      TestSet --> BalTest["BalancedClusterDataset (Test)"]
      Bal --> Sampler["DistributedWeightedSampler (src/data/balancing.py:131-182)"]
    end

    classDef default font-family:sans-serif;
```
*Figure 2.2: Data balancing and distributed sampling structure.*

For full implementation details, see [Multi-Source Balancing and Splits](#2.3).

Sources: [src/data/balancing.py:13-182]()

---

## 2.4 Dataset Preparation Scripts

Prior to training, raw databases must be processed through preparation scripts located in `scripts/`. These pipelines handle PDB entry parsing, AlphaFold Database (AFDB) domain verification, MMSeqs2 sequence clustering and redundancy filtering, pLM embedding generation (e.g., via ESMC), and UniProt sequence fetching.

For full implementation details, see [Dataset Preparation Scripts](#2.4).

Sources: [src/data/dataset.py:41-43]()

---