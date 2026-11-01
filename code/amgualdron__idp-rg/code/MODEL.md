# Model card: idprg 1.0

## What it predicts

The model gives the equilibrium distribution P(R_g) of the radius of gyration of a single intrinsically disordered chain.
The input is the amino-acid sequence alone: the 20 standard residues, N ≥ 20.

The predictions reproduce HPS1 coarse-grained molecular dynamics (Dignon et al. 2018) at 298 K and 150 mM salt. They do
not reproduce experiment, other force fields, or other temperatures and salt concentrations.

There are two predictors with the same output: a network and a 21-parameter closed form.

## The output distribution

With s = R_g / S_m:

    p(s) ds ∝ exp(−A [s^δ/δ + s^−g/g]) ds,   δ = 1/(1−ν),   g = 3/(3ν−1),   ν = 1/3 + (2/3) expit(u)

This is Victor and Lhuillier's two-exponent form (J. Chem. Phys. 92, 1362, 1990) in d = 3.
- S_m is the mode.
- A sets the width.
- ν sets the asymmetry. The excluded-volume wall at small R_g falls as s^−g and the stretched tail as s^δ.

The normaliser, moments, CDF and quantiles are computed by quadrature on t = ln s ∈ [−3, 2] with 4,001 points
(`idprg/distribution.py`). The density integrates to 1 within 2·10⁻¹³.

## Network

    ln S_m = b_S + ν_F ln N + o₁(x),   u = u₀ + o₂(x),   ln A = c · [1, ν, ν², 1/N, ⟨λ⟩, s₊]

- **MLP.** o = MLP(x): 9 → 512 → 512 → 2 with GELU. The last layer starts at zero, so training starts from the
  population-average distribution.
- **Fixed constants.** ν_F = 0.57874 is the Flory exponent of the training proteins' ⟨R_g⟩. b_S and u₀ are fixed
  offsets: the means of the per-protein fits with ln A on the law (`data/md_fits.csv`).
- **Amplitude law.** The six law coefficients c are global and learned with the MLP. They start from least squares on
  the per-protein maximum-likelihood fits of the training proteins.
- **Parameters: 268,808**, of which 268,802 are the MLP and 6 the law. The input means and standard deviations and the
  two offsets are fixed from the training split and not trained.

**Inputs x.** The nine inputs are standardized with the training-split mean and sd (`model/network.json`), in this
order:

| input | definition |
|---|---|
| `log_seq_length` | ln N |
| `scd_pos`, `scd_neg` | s = asinh(5·SCD/N), split at zero: s₊ = max(s, 0) means like charges dominate; s₋ = min(s, 0) means opposite charges are segregated into blocks. SCD = Σ_{i<j} q_i q_j \|i−j\|^½ / N |
| `shd` | sequence hydropathy decoration, Σ_{i<j} (λ_i + λ_j) \|i−j\|⁻¹ / N |
| `contact_energy` | loop-entropy-weighted contact score of residue windows from the HPS1 pair energies (`idprg/descriptors.py`, `paper/closed_form.tex`) |
| `fcr` | fraction of charged residues |
| `patch_neg`, `patch_neu` | fraction of adjacent residue pairs that are both −1, or both neutral |
| `frac_sticky` | fraction of residues with λ ≥ 0.85 (I L F P W Y V) |

The amplitude law also uses ⟨λ⟩, the mean HPS1 hydropathy, and s₊, both unstandardized. q and λ are the HPS1 residue
charges (histidine +0.5) and hydropathies.

**Training** (`idprg train`). Training uses the 3,223 training proteins.
- **Loss:** the mean −ln p(R_g) per frame over 1,024 stored MD frames per protein, redrawn every batch. No other terms.
- **Optimisation:** AdamW with learning rate 10⁻³ and cosine decay to 0, weight decay 10⁻⁴, gradient-norm clip 1.
  Batch 64, 40 epochs, float64, seed 42.
- **Checkpoint:** the last epoch is kept; there is no checkpoint selection.
- **Cost:** about 1 minute on a CPU.
- **Reproducibility:** retraining reproduces `model/network.eqx` bit for bit (SHA-256 in `model/network.json`).

## Closed form

The closed form uses the same distribution, with the network replaced by:
- ln S_m (9 coefficients) linear in s₊, s₋, ⟨λ⟩, h and ln N, plus s₋ ln N, ⟨λ⟩ ln N and (ln N)²;
- u (6) linear in s₊, s₊/N, s₋, ⟨λ⟩ and h;
- ln A from the law [1, ν, 1/N, ⟨λ⟩] (4).

h = max(CE − c(N), 0) is the excess of the contact energy over the CE threshold c(N) = c₀ + 0.25 ln N. In total the
closed form has 19 coefficients and 2 threshold constants, 21 parameters.

It was fitted by maximum likelihood on all 30,000 stored frames of each training protein, with a weak ridge penalty on
the slopes (`idprg fit-closed-form`, about 10 minutes on a CPU). `paper/closed_form.tex` gives the formula with every
coefficient, the descriptors in full, a worked example and its scores, including three-fold cross-validation on the
training split (`idprg fit-closed-form --fold 0|1|2`, pooled in `results/closed_form_cv.csv`).

## Data

The data are HPS1 simulations of 4,609 IDP sequences, each with 30,000 stored R_g frames. They are split by 30%-identity
sequence clusters into 3,223 training, 695 validation and 691 test proteins. `data/README.md` covers the simulations,
the cleaning rule, the split and the tables.

## Metrics

Each metric is computed per protein against its MD ensemble and then averaged over proteins (`idprg/metrics.py`).

- **NLL:** −ln p(R_g) per frame over all 30,000 stored frames, density in nm⁻¹. It grows with N, so it is comparable
  only on the same proteins.
- **⟨R_g⟩ error and σ error:** mean |relative error| against the full-trajectory MD moments. **Bias** is the signed
  mean.
- **Skew:** mean |Δ skewness|.
- **W₁:** exact Wasserstein-1 distance to the stored frames, in % of the MD ⟨R_g⟩; the median over proteins.
- **Compact chains:** MD ⟨R_g⟩ < 0.8 × exp(−1.56603) N^0.57874, which is about 5% of proteins.

## Results

The tables between `<!-- numbers -->` markers are written by `paper/make_numbers.py` from `results/`; do not edit them
by hand.

### Validation (695 proteins)

<!-- numbers:val -->
|  | NLL | ⟨R_g⟩ err | bias | σ err | skew err | W₁ median [p90] | compact ⟨R_g⟩ err |
|---|---|---|---|---|---|---|---|
| network | 1.1525 | 1.84% | +0.35% | 4.73% | 0.093 | 1.47% [4.24] | 4.90% |
| closed form | 1.1620 | 2.53% | +0.64% | 5.81% | 0.098 | 1.79% [5.59] | 12.6% |
<!-- /numbers -->

By chain length (⟨R_g⟩ error / σ error; number of proteins in brackets):

<!-- numbers:val-length -->
| N | ≤ 50 (51) | 51–100 (105) | 101–200 (223) | 201–300 (192) | > 300 (124) |
|---|---|---|---|---|---|
| network | 1.51% / 4.5% | 1.43% / 3.8% | 1.90% / 5.0% | 1.94% / 5.3% | 2.04% / 4.2% |
| closed form | 1.62% / 4.8% | 1.63% / 3.7% | 2.27% / 5.5% | 3.53% / 8.0% | 2.58% / 5.2% |
<!-- /numbers -->

By contact energy, relative to the training quantiles (⟨R_g⟩ error / σ error):

<!-- numbers:val-ce -->
| CE | ≤ q75 = 0.333 (509) | q75–q90 (91) | ≥ q90 = 1.035 (95) |
|---|---|---|---|
| network | 1.40% / 3.5% | 2.20% / 4.0% | 3.82% / 11.7% |
| closed form | 1.79% / 3.7% | 2.56% / 4.8% | 6.47% / 18.0% |
<!-- /numbers -->

### Test (691 proteins)

The test split was scored once for the network, on 2026-10-27, and once for the closed form, on 2026-10-31, each after
that model was fixed (`data/README.md`, test-split log). The closed form was revised after the 2026-10-27 opening of its
earlier 24-parameter version, so its scores are not a blind estimate.

<!-- numbers:test -->
|  | NLL | ⟨R_g⟩ err | bias | σ err | skew err | W₁ median [p90] | compact ⟨R_g⟩ err |
|---|---|---|---|---|---|---|---|
| network | 1.1379 | 2.02% | −0.18% | 4.69% | 0.104 | 1.53% [4.34] | 7.02% |
| closed form | 1.1462 | 2.53% | −0.18% | 5.09% | 0.112 | 1.82% [5.12] | 12.6% |
<!-- /numbers -->

By chain length and by contact energy (⟨R_g⟩ error / σ error):

<!-- numbers:test-length -->
| N | ≤ 50 (69) | 51–100 (102) | 101–200 (209) | 201–300 (152) | > 300 (159) |
|---|---|---|---|---|---|
| network | 1.90% / 4.3% | 1.36% / 4.5% | 1.93% / 4.9% | 2.41% / 4.6% | 2.23% / 4.8% |
| closed form | 1.84% / 3.7% | 1.61% / 4.3% | 2.27% / 5.3% | 2.99% / 5.4% | 3.30% / 5.5% |
<!-- /numbers -->

<!-- numbers:test-ce -->
| CE | ≤ q75 = 0.333 (485) | q75–q90 (127) | ≥ q90 = 1.035 (79) |
|---|---|---|---|
| network | 1.41% / 3.6% | 2.41% / 4.1% | 5.14% / 12.5% |
| closed form | 1.71% / 3.5% | 3.00% / 3.7% | 6.76% / 17.0% |
<!-- /numbers -->

The test scores are close to the validation scores, and the weak regimes are the same: sticky chains (CE ≥ q90) and
compact chains. The test split has 22 compact chains, against 37 in validation.

### Extrapolation

Each model was retrained, or refitted, on part of the training split and scored on training proteins outside that
part. Sequence-cluster mates of the scored proteins were removed from the fit set.
- **Length:** fit on N ≤ 200, scored on N ≥ 300. Scored chains are up to 2.45× longer than any fitted chain.
- **Sticky:** fit on CE ≤ 0.333, scored on CE ≥ 1.035.

<!-- numbers:extrapolation -->
|  | NLL | ⟨R_g⟩ err | bias | σ err | skew err | W₁ median [p90] | compact ⟨R_g⟩ err |
|---|---|---|---|---|---|---|---|
| **length** (602 scored, 35 compact) |  |  |  |  |  |  |  |
| network | 1.6510 | 2.38% | −0.07% | 4.85% | 0.075 | 1.69% [5.02] | 6.74% |
| closed form | 1.6628 | 3.26% | −0.98% | 5.56% | 0.100 | 2.48% [6.48] | 10.0% |
| **sticky** (323 scored, 82 compact) |  |  |  |  |  |  |  |
| network | 1.3195 | 5.56% | +1.86% | 15.6% | 0.171 | 3.94% [13.2] | 10.9% |
| closed form | 1.3305 | 5.83% | +2.42% | 19.3% | 0.246 | 4.64% [12.6] | 11.0% |
<!-- /numbers -->

The network models trained this way are rows `reference`, seed 42, of `results/ablations.csv`; the scores are in
`results/metrics.json` under `length` and `sticky`.

### Experimental benchmark (122 proteins)

The benchmark proteins have a SAXS R_g and an HPS1 trajectory. They were not used for any fitting or selection.

<!-- numbers:benchmark -->
| ⟨R_g⟩ | against HPS1: error / bias | against SAXS: error / bias / MAE |
|---|---|---|
| network | 1.73% / +0.61% | 13.6% / +7.60% / 0.416 nm |
| closed form | 1.98% / +0.93% | 13.7% / +7.92% / 0.419 nm |
| HPS1 simulation itself | — | 12.5% / +6.81% / 0.386 nm |
<!-- /numbers -->

Both predictors reproduce the simulations they were trained on. Their deviation from experiment is that of the HPS1
force field.

## Ablations

`idprg ablations` retrains the network, or refits the closed form, with one change from the released configuration
(`idprg/ablations.py`; `results/ablations.csv`). Network rows use seed 42; the reference row gives the spread over
seeds 42–46. Rows marked *pending* have not been run yet: `idprg ablations` resumes the grid where it stopped, and
`python paper/make_numbers.py` then fills them in.

<!-- numbers:ablations-network -->
| network, seed 42 | parameters | val NLL | val ⟨R_g⟩ err | val σ err | length ⟨R_g⟩ err | sticky ⟨R_g⟩ err |
|---|---|---|---|---|---|---|
| reference (seeds 42–46, mean ± sd) | 268,808 | 1.1523 ± 0.0001 | 1.82 ± 0.01% | 4.72 ± 0.01% | 2.32 ± 0.08% | 5.58 ± 0.19% |
| `-scd_pos` | 268,296 | 1.1686 | 2.97% | 5.58% | 5.38% | 10.0% |
| `-scd_neg` | 268,296 | 1.1563 | 2.05% | 4.83% | 2.80% | 6.74% |
| `-shd` | 268,296 | 1.1537 | 1.99% | 4.87% | 2.46% | 5.87% |
| `-contact_energy` | 268,296 | 1.1616 | 2.67% | 5.46% | 2.92% | 13.9% |
| `-fcr` | 268,296 | 1.1524 | 1.83% | 4.72% | 2.27% | 5.73% |
| `-patch_neg` | 268,296 | 1.1533 | 1.92% | 4.79% | 2.29% | 5.93% |
| `-patch_neu` | 268,296 | 1.1527 | 1.84% | 4.72% | 2.37% | 5.73% |
| `-frac_sticky` | 268,296 | 1.1527 | 1.85% | 4.73% | 2.39% | 5.69% |
| `ln N only` | 264,712 | 1.3044 | 8.56% | 17.3% | 12.3% | 22.3% |
| width × depth 128x2 | 18,056 | 1.1533 | 1.91% | 4.75% |  |  |
| width × depth 1024x2 | 1,061,896 | 1.1527 | 1.86% | 4.72% |  |  |
| width × depth 512x1 | 6,152 | 1.1547 | 1.99% | 4.84% |  |  |
| width × depth 512x3 | 531,464 | 1.1521 | 1.81% | 4.69% |  |  |
<!-- /numbers -->

<!-- numbers:ablations-closed-form -->
| closed form | parameters | val NLL | val ⟨R_g⟩ err | val σ err | val W₁ median |
|---|---|---|---|---|---|
| reference | 21 | 1.1620 | 2.53% | 5.81% | 1.79% |
| `-scd_pos` | 18 | 1.2320 | 6.18% | 12.6% | 4.65% |
| `-scd_neg` | 18 | 1.1644 | 2.67% | 6.33% | 1.86% |
| `-mean_hydro` | 17 | 1.1673 | 3.10% | 6.08% | 2.27% |
| `-h` | 17 | 1.1949 | 3.97% | 7.03% | 2.27% |
| no ln N interactions | 19 | 1.1632 | 2.61% | 5.99% | 1.83% |
| constant CE threshold | 20 | 1.1623 | 2.56% | 5.82% | 1.80% |
| free ln A | 23 | 1.1619 | 2.55% | 6.22% | 1.82% |
<!-- /numbers -->

## Domain of validity

**Training ranges.** The table gives the minimum and maximum over the training proteins, and their 1st and 99th
percentiles.

<!-- numbers:domain -->
|  | range | 1–99% |
|---|---|---|
| N | 20 to 490 | 24 to 427 |
| s₊ | 0 to 3.80 | 0 to 2.25 |
| s₋ | −0.984 to 0 | −0.255 to 0 |
| SHD | 1.96 to 9.50 | 3.13 to 7.99 |
| CE | −4.43 to 8.05 | −1.79 to 2.75 |
| FCR | 0 to 0.981 | 0.0447 to 0.787 |
| patch_neg | 0 to 0.966 | 0 to 0.389 |
| patch_neu | 0 to 1 | 0.0405 to 0.917 |
| frac_sticky | 0 to 0.778 | 0.00883 to 0.489 |
| ⟨λ⟩ | 0.267 to 0.912 | 0.404 to 0.766 |
<!-- /numbers -->

- **Confident.** Within these ranges, with CE below about 1 (the training q90) and a chain that is not strongly compact,
  the network's errors are at validation level: ⟨R_g⟩ 1.4–2%, σ 3.5–5%, W₁ about 1.4%.
- **Long chains.** Trained up to N = 200, the network predicted N = 300–490 almost as well as in distribution. Moderate
  length extrapolation is supported. Beyond N ≈ 500 nothing was tested.
- **Sticky chains** are the weak regime. This holds in distribution (CE ≥ q90: σ error 12%) and more so beyond the
  fitted range.
- **Compact chains** (MD ⟨R_g⟩ < 0.8 × the Flory baseline, about 5%) are the least reliable:
  - network: ⟨R_g⟩ over-predicted by about 3% and σ off by 18%;
  - closed form: over-predicted by 11.5%, with σ off by 36%.
- **Physics.** HPS1 at 298 K and 150 mM only. Nothing in the model depends on temperature or salt.
