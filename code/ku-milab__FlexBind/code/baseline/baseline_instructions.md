# Baseline evaluation protocol

All baselines are evaluated on the same DP81, DP93 and DP94 test sets as FlexBind, with identical
test proteins, label definitions and metric computation (`scripts/evaluate.py`).

## Compared methods

| Purpose | Methods |
|---|---|
| Protein/RNA/DNA-binding residues | DeepDISOBind, DisoRDPbind, MTDsite, ProNA2020, DisoFLAG, flDPnn |
| IDR identification | DisoFLAG, flDPnn, SPOT-Disorder2, SPOT-Disorder-Single, ESpritz, AUCpreD, IUPred3, DisoMine, MobiDB-lite |

## Protocol

1. **Methods with public training code** are re-trained from scratch on our training sets
   (same splits as FlexBind) and evaluated over the same 10 random seeds
   (`configs/default_config.yaml → random_seeds`); results are reported as mean ± SD.
2. **Methods without public training code** are run through their official web servers /
   released predictors and reported as single-run results.
3. Convert every baseline's residue-level scores to the FlexBind prediction format so that the
   same metric and significance code is used:

   ```python
   {
     "tasks": ["disorder", "protein", "rna", "dna"],
     "proteins": {"<id>": {"prob": np.ndarray (L, 4),     # NaN for tasks the method does not predict
                           "label": np.ndarray (L, 4)}},  # same labels as the test set
     "thresholds": {"disorder": t1, "protein": t2, "rna": t3, "dna": t4},
     "method": "<name>"
   }
   ```
   Thresholds are the F1-maximising thresholds on the validation set (for re-trained methods), or
   the method's default operating point when no validation predictions are available.
4. Significance against the second-best method:

   ```bash
   python scripts/bootstrap_test.py --pred_a <flexbind.pkl> --pred_b <baseline.pkl> --n_iter 1000
   ```
