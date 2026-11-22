# FlexBind: Multi-Scale Complex Learning for Imbalance-Robust IDR and Binding Residue Prediction

FlexBind is a sequence-only deep learning framework that **jointly** predicts, at the residue level:
- Intrinsically Disordered Regions (IDRs)
- Protein-binding residues
- RNA-binding residues
- DNA-binding residues

It needs neither structural information nor pre-annotated IDRs, so it can be used as a
high-throughput computational filter on raw protein sequences.

---

## Architecture

```
protein sequence
   │  ProtT5 (frozen)                                   H ∈ R^{B×L×1024}
   ▼
┌──────────────── shared multi-scale complex representation ────────────────┐
│ CWT with complex Morlet wavelet, 7 scales {1,3,4,5,8,16,32}, kernel 5      │
│   → Re / Im streams → learnable 1×1 conv (1024 → 128)   ∈ R^{B×128×S×L}     │
└───────────────────────────────────────────────────────────────────────────┘
   │ (one independently parameterised branch per task: IDR / protein / RNA / DNA)
   ▼
 Residue-wise Scale Attention (RWSA)  – query from H_l, keys/values from each scale
 Complex Axis Alignment (CAA)         – channel-wise learnable complex scalar (a_k + i b_k)
 Dual-stream integration              – C̄ = W_F[axis1 ‖ axis2],  H_fused = H + W_O C̄
 3 × Mamba block (d = 1024)           – padding-aware selective state-space model
 MLP + Sigmoid                        – residue-level probability
```

Training objective: task-specific **asymmetric focal loss** (γ⁺ = 0; γ⁻ = 2 for IDR/protein,
3 for RNA/DNA) combined with **uncertainty-based task weighting** (σ_t initialised to 1).

| Component | File |
|---|---|
| CWT, RWSA, CAA, dual-stream integration, Mamba, task branch | `models/layers.py` |
| Full model | `models/flexbind.py` |
| Asymmetric focal loss, uncertainty weighting | `models/losses.py` |
| Data loading (full-length, dynamic padding) | `scripts/dataset.py` |
| Training (early stopping, seeds, thresholds) | `scripts/train.py` |
| Metrics / checkpoint evaluation | `scripts/evaluate.py` |
| Protein-level paired bootstrap test | `scripts/bootstrap_test.py` |
| ProtT5 embedding extraction | `scripts/extract_embeddings.py` |
| Inference from FASTA | `scripts/predict.py` |

---

## 📌 Reproducibility
- **Preprocessing & split commands**: `data/preprocess.sh`
- **Configuration & random seeds**: `configs/default_config.yaml` (all hyperparameters, 10 seeds)
- **Pretrained weights**: `weights/` (base ProtT5 weights via `weights/download_weights.sh`)
- **Baseline protocol**: `baseline/baseline_instructions.md`
- **Environment**: `environment.yaml`

---

## Installation

```bash
conda env create -f environment.yaml
conda activate flexbind
# optional, recommended for speed (fused CUDA selective scan; a pure-PyTorch fallback is built in)
pip install mamba-ssm
```

The parameterisation of the Mamba blocks is identical with and without `mamba-ssm`, so
checkpoints can be used in either setting.

## Pretrained Weights
Download the base ProtT5 encoder weights:

```bash
cd weights
bash download_weights.sh
cd ..
```

## Data format
- Embeddings: `.npz` with `ids` and `per_residue` (object array of `(L, 1024)` ProtT5 embeddings).
  Create one from a FASTA file with
  `python scripts/extract_embeddings.py --fasta seqs.fasta --out data/processed/X.npz`.
- Labels: one TSV per task, each line `<id>\t<0/1 string of length L>`.

Full-length sequences are always used (no truncation or chunking); mini-batches are padded to the
longest sequence in the batch and padded positions are masked in the CWT, the attention, the Mamba
state updates and the losses.

## Training and evaluation

```bash
python scripts/train.py --config configs/default_config.yaml --name DP81 \
    --npz_path data/processed/DP81/train_val.npz \
    --disorder_txt data/processed/DP81/train_disorder.txt \
    --protein_txt  data/processed/DP81/train_protein.txt \
    --rna_txt      data/processed/DP81/train_rna.txt \
    --dna_txt      data/processed/DP81/train_dna.txt \
    --test_npz_path data/processed/DP81/test.npz \
    --test_disorder_txt data/processed/DP81/test_disorder.txt \
    --test_protein_txt  data/processed/DP81/test_protein.txt \
    --test_rna_txt      data/processed/DP81/test_rna.txt \
    --test_dna_txt      data/processed/DP81/test_dna.txt
```
(Adjust the label-file paths to your processed data.)

- Use `--val_npz_path` and `--val_{disorder,protein,rna,dna}_txt` to supply the benchmark's
  validation set. Without them, a fixed random split of `--npz_path` (`val_ratio`, `split_seed`)
  is used, identical for every seed.
- AdamW (lr 1e-4, weight decay 0.01), batch size 32; checkpoint selection and early stopping on
  validation AUPR (mean over the four tasks) with patience 10.
- The decision threshold of each task is the F1-maximising threshold on the validation set; it is
  stored in the checkpoint and applied unchanged to the test set.
- All 10 seeds in the config are run by default (`--seeds 42` for a single run). Outputs:
  `checkpoints/<name>/seed<k>/{best.pt, results.json, test_predictions.pkl}` and
  `checkpoints/<name>/summary.json` (mean ± SD).

Evaluate a trained model on another test set (e.g. DP94 with the DP93 model):

```bash
python scripts/evaluate.py --checkpoint checkpoints/DP93/seed42/best.pt \
    --npz_path data/processed/DP94/test.npz \
    --disorder_txt ... --protein_txt ... --rna_txt ... --dna_txt ... \
    --out_dir results/DP94_seed42
```

Metrics (AUPR, F1max, AUROC, Precision, Recall, MCC) are micro-averaged over all test residues.

## Statistical significance

```bash
python scripts/bootstrap_test.py --pred_a checkpoints/DP81/seed42/test_predictions.pkl \
    --pred_b <second_best_method>.pkl --n_iter 1000
```
Proteins are resampled with replacement (paired between methods), residues are pooled, and
micro-averaged metrics are compared; `*` marks p < 0.05.

## Inference

```bash
python scripts/predict.py --checkpoint checkpoints/DP81/seed42/best.pt \
    --fasta my_proteins.fasta --out predictions.tsv
```
