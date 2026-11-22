# Hierarchical Transformer-Based IDP-PLM

Welcome to the repository of **Hierarchical Transformer-Based Intrinsically Disordered Protein Language Model (IDP-PLM)**.

This repository contains the implementation of a hierarchical Transformer-based protein language model for learning protein sequence representations, extracting hierarchical embeddings, and predicting properties of intrinsically disordered proteins (IDPs).

The repository also includes multiple regression frameworks for IDP property prediction, including LSTM, PLM-based regression, and Transformer-based models.

---

## Repository Structure

The main branch contains three directories:

```text
.
├── data/
├── model/
└── Albatross/
```

---

## 1. Data

The `data/` directory contains protein sequences and datasets used for model training, validation, and property prediction.

| File | Description |
|---|---|
| `sequences_len_le_500.txt` | Protein sequences of IDPs with sequence lengths ≤ 500. |
| `combined_all_prop.csv` | Combined dataset containing all protein properties used for Albatross regression. |
| `train_all_prop.csv` | IID training dataset generated from the combined dataset. |
| `test_all_prop.csv` | IID test dataset generated from the combined dataset for regression evaluation. |
| `train_all_prop_minmax.csv` | Min-max normalized training dataset. |
| `test_all_prop_minmax.csv` | Min-max normalized test dataset. |

---

## 2. Model

The `model/` directory contains the hierarchical Transformer-based protein language model, embedding extraction scripts, and clustering utilities.

### Hierarchical Transformer PLM

| Script | Description |
|---|---|
| `hieratchical_transformer.py` | Hierarchical Transformer-based protein language model, including its training implementation. |
| `extract_token_representation.py` | Inference script for extracting token-level, local-level, and global-level protein sequence representations from a trained model. |
| `pc_kmeans_tsne.py` | Clustering and visualization of extracted protein embeddings using PCA, K-means, and t-SNE. |

The model provides hierarchical representations at three levels:

- **Token representations:** Representations associated with individual sequence tokens.
- **Local representations:** Representations capturing local sequence-level information.
- **Global representations:** Representations capturing global protein sequence information.

### Nardini Regression

The `model/NardiniRegression/` directory contains scripts for fine-tuning and evaluating the learned PLM representations against Nardini-derived protein features.

| Script | Description |
|---|---|
| `Fine_tune_on_NARDINI.py` | Regression model mapping PLM embeddings to Nardini protein features. |
| `EXTRACT_EMBEDDING.py` | Extracts embeddings from the trained regression model. |
| `pca_kmeans_tsne.py` | Clustering and visualization of extracted embeddings using PCA, K-means, and t-SNE. |

---

## 3. Albatross

The `Albatross/` directory contains three different approaches for predicting IDP properties from protein sequences or learned protein representations:

1. LSTM
2. PLM
3. Transformer

### 3.1 LSTM

The `Albatross/LSTM/` directory implements bidirectional LSTM-based regression for predicting IDP properties directly from protein sequences.

| Script | Description |
|---|---|
| `lstm.py` | Bidirectional LSTM model for predicting IDP features from protein sequences. |
| `train.py` | Training script for sequence-to-feature regression using the BiLSTM model. |
| `infere.py` | Inference script for predicting IDP properties using the trained LSTM model on the provided dataset. |

### 3.2 PLM

The `Albatross/PLM/` directory implements regression using embeddings extracted from the pretrained protein language model.

| Script | Description |
|---|---|
| `model_2.py` | Regression model mapping PLM embeddings to IDP properties. |
| `train_from_embeddings.py` | Training script for embedding-to-feature regression. |
| `predict_from_embeddings.py` | Inference script for predicting IDP properties from protein embeddings. |

### 3.3 Transformer

The `Albatross/Transformer/` directory implements Transformer-based regression directly from protein sequences.

| Script | Description |
|---|---|
| `protein_transformer_regression.py` | Transformer regression model for predicting IDP features directly from protein sequences. |
| `train.py` | Training script for sequence-to-feature regression using the Transformer model. |
| `predict.py` | Prediction script using the trained Transformer regression model. |

---

## 4. Model Overview

The repository supports two complementary modeling tasks.

### Task 1: Protein Language Modeling

The hierarchical Transformer learns protein sequence representations at multiple levels:

- Token-level representations
- Local sequence representations
- Global sequence representations

These representations can subsequently be used for downstream analysis and regression tasks.

### Task 2: IDP Property Prediction

Three approaches are implemented for predicting IDP properties:

| Approach | Input | Model |
|---|---|---|
| LSTM | Protein sequences | Bidirectional LSTM |
| PLM | Precomputed PLM embeddings | Regression network |
| Transformer | Protein sequences | Transformer regression |

This enables comparison between conventional sequence-based models and models that leverage pretrained protein language model representations.

---

## 5. Workflow

The general workflow is:

1. Prepare the IDP protein sequences and property datasets.
2. Train the hierarchical Transformer-based protein language model.
3. Extract token-level, local-level, and global-level representations.
4. Perform embedding analysis and clustering.
5. Fine-tune PLM representations using Nardini-derived features.
6. Train and evaluate the LSTM, PLM-based, and Transformer-based regression models for IDP property prediction.

---

## 6. Citation

If you use this repository or build upon the implementations, please cite the associated research work.(coming soon)

