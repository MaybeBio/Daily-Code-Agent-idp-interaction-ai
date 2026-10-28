#!/usr/bin/env bash
# Full reproduction of the analysis (≈5–6 h on 2 CPU cores; no GPU required).
set -euo pipefail
cd "$(dirname "$0")"
bash scripts/download_data.sh                      # public inputs -> data/raw
PSPRED=ext/pspred                                  # von Buelow et al. 2025 predictor (GPL-3.0)

python scripts/00_fetch_me_sequences.py data/raw
python scripts/01_build_isoforms.py data/raw data/processed
python scripts/05_neural_psi.py data/raw results/neural_dpsi.tsv
python - <<'EOF'
import pandas as pd
d = pd.read_csv('data/processed/isoform_pairs.tsv', sep='\t')
d = d[(d.mapped == True) & (d.in_frame == True)].drop_duplicates(subset=['transcript', 'aa_pos', 'me_peptide'])
d = d.sample(frac=1, random_state=1)
for i in range(2):
    d.iloc[i::2].to_csv(f'data/processed/shard{i}.tsv', sep='\t', index=False)
EOF
python scripts/02_llps_features.py data/processed/shard0.tsv $PSPRED results/features_shard0.tsv 50 2026 &
python scripts/02_llps_features.py data/processed/shard1.tsv $PSPRED results/features_shard1.tsv 50 2026 &
python scripts/03_alphafold_plddt.py data/processed/isoform_pairs.tsv data/raw/Homo_sapiens.GRCh38.115.uniprot.tsv.gz results/alphafold_plddt.tsv &
wait
python scripts/04_benchmark_predictors.py $PSPRED data/raw results/benchmark 3000 7
python scripts/08_orthogonal_llphyscore.py ext/llphyscore/standalone_package/LLPhyScore 40
python scripts/06_analysis.py $PSPRED                 # switch calls (table S2)
python scripts/09_switch_permutation.py $PSPRED 1000   # 1,000-null tests for switches
python scripts/06_analysis.py $PSPRED                 # final statistics, tables, figures
