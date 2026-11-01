# Data

## Simulations

Every protein was simulated with the HPS1 coarse-grained model (Dignon et al., *PLoS Comput. Biol.* 14, e1005941, 2018):
one bead per residue, Ashbaugh–Hatch short-range interactions with the residue hydropathies λ, and screened
Coulomb interactions. All runs used the same conditions: 298 K, 150 mM salt (Debye screening κ = 1.274 nm⁻¹, Bjerrum
length 0.7147 nm). The radius of gyration R_g was recorded at every saved frame.

- **4,743 sequences** were simulated. Their trajectories were run for about 0.5 N^2.2 frames.
- **The frame store** holds 30,000 R_g values per protein, rounded to 0.001 nm. They are a random subsample of the
  trajectory, in random order. Shorter runs were padded by resampling; 1,890 of the proteins used here were run for
  fewer than 30,000 frames.
- **The MD moments** in `proteins.csv` (mean, standard deviation, skewness, excess kurtosis) are computed from the full
  trajectory, not the subsample.

## Cleaning: 4,609 proteins

The 4,743 runs were filtered in three steps, leaving 4,610 proteins.

1. **Integrity.** The coefficient of variation of R_g must be at most 0.5. This removes exactly the 89 trajectories
   with frames beyond the rigid-rod limit, which are simulation blow-ups.
2. **Statistical precision.** The relative standard error of ⟨R_g⟩ must be at most 3%. It is computed as
   CV / √n_eff, with an effective sample count n_eff = 170 × (run length / median run length). Run length is measured
   in units of N^2.2 frames. The 170 is calibrated on replicate simulations of identical sequences. This step removes
   runs that stopped well short of their design length.
3. **N ≥ 20.**

Protein 6029 is then dropped. Its stored frames are corrupt: 44% are a single value, although its full-trajectory
moments are fine. That leaves **4,609 proteins**.

## Split

Sequences were clustered with MMseqs2: `easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0 -s 7.5 --mask 0`, giving
3,603 clusters among the 4,609 proteins. Whole clusters were assigned to the splits, so no validation or test protein
shares a cluster with a training protein.

| Split | Proteins |
|---|---|
| training | 3,223 |
| validation | 695 |
| test | 691 |

- **Assignment.** Clusters were stratified by the tercile of their median N and by the majority sign of the SCD.
  Every cluster holding more than 1% of the proteins went to training. The random seed was 2026.
- **What each split was used for.** All fitting used the training split only. Model choices were made by
  cross-validation inside the training split and on the length and sticky extrapolation splits (`idprg.data.split`).
  The validation split was used to compare final candidates.
- **The test split** was scored once for each released model, after that model was fixed (see the log below). It had
  also been scored for earlier versions, so it is not entirely unseen. It was never used to fit the released models.

## Tables

### `proteins.csv`: one row per protein, 4,609 rows

| column | content |
|---|---|
| `id` | simulation id |
| `split` | `train`, `val` or `test` |
| `cluster` | sequence cluster (30% identity) |
| `N` | chain length |
| `frame_row` | row of the protein's frames in the frame store |
| `n_frames` | length of the full trajectory, in frames |
| `rg_mean`, `rg_std`, `rg_skew`, `rg_exkurt` | MD moments of R_g over the full trajectory (nm, nm, –, excess) |
| `log_seq_length` … `mean_hydro` | the network's inputs and ⟨λ⟩, from `idprg.descriptors`; recomputing them from `sequence` agrees to 1e-15 |
| `sequence` | one-letter amino-acid sequence |

### `md_fits.csv`: maximum-likelihood fits of the two-exponent distribution to each protein's frames

There is one row per training and validation protein (3,918 rows). Regenerate the table with `idprg fit-md`.

| column | content |
|---|---|
| `lnA`, `lnSm`, `u`, `nll` | free fit: all three parameters, and the NLL per frame |
| `lnSm_<split>`, `u_<split>` | fit with ln A fixed to the amplitude law of that training split (`standard`, `length`, `sticky`) |

How training and fitting use these values:
- The amplitude law's starting coefficients are least squares on the free fits of the training proteins.
- The network's two fixed offsets are means of the law-constrained fits.
- The closed form's starting point is ridge regressions on the free fits.

### `benchmark.csv`: 122 proteins with an experimental R_g

Each row is a unique sequence. Entries that share a sequence are merged and their SAXS values averaged.

| column | content |
|---|---|
| `name` | protein name |
| `family` | sequence family: proteins linked by sequence similarity ≥ 0.6 or the same composition |
| `aliases` | other names of the merged literature entries (`;`-separated), empty if none |
| `n_entries` | number of literature entries merged into this row |
| `N` | chain length |
| `saxs_rg_nm` | SAXS R_g from the literature |
| `hps1_rg_nm` | ⟨R_g⟩ of the protein's HPS1 simulation |
| `trajectories` | trajectory files of that simulation (`;`-separated) |
| `n_frames` | total frames in those trajectories |
| `sequence` | one-letter amino-acid sequence |

The benchmark was not used for any fitting or selection. For hNHE1cdt, the compiled literature table gives 130 residues, without
the first M that `sequence` has.

### `benchmark_methods.csv`: ⟨R_g⟩ of the same 122 proteins from other methods

The rows are those of `benchmark.csv`, in the same order. All values are in nm, converted from the Å of the compiled
literature table. They are used only in the paper's methods comparison (`agree.*` and `signtest.*` in
`paper/numbers.tex`).

| column | content |
|---|---|
| `name` | protein name, as in `benchmark.csv` |
| `albatross_rg_nm` | ALBATROSS, a network trained on Mpipi-GG simulations |
| `idpgan_rg_nm` | idpGAN, a generative model of conformational ensembles |
| `zheng_rg_nm` | the sequence-based predictor of Zheng et al. (2020) |
| `afrc_rg_nm` | AFRC, the analytical Flory random coil (chain length only) |
| `mpipi_gg_rg_nm` | Mpipi-GG simulations |
| `calvados2_rg_nm` | CALVADOS 2 simulations |

## Files not in this repository

| file | content |
|---|---|
| `idp_rg_frames.npy` | the frame store: float32, shape (4743, 30000), R_g in nm; 543 MiB; SHA-256 `eef04250701014b3b311dcf535bddd6af92cb9c4c28707d8fc33844dac5dda4a` |
| `<trajectory>_pooled_timeseries.csv` | the benchmark proteins' HPS1 R_g trajectories (column `rg`, nm); used by `paper/make_figures.py` only |

Commands that need the frame store take `--frames PATH`, or read the `IDPRG_FRAMES` environment variable.

## Test-split log

The test split was opened five times in total. The last two openings concern the models in this repository; the
earlier ones scored development versions.

| date | model scored |
|---|---|
| 2026-10-21 | an earlier network |
| 2026-10-24 | an earlier 8-input network and its 30-parameter closed form |
| 2026-10-25 | an earlier 24-parameter two-exponent closed form |
| 2026-10-27 | **the released network**: `model/network.eqx` (SHA-256 `bd757251…`), and the 24-parameter closed form (`f1ea4ca1…`), scored once with `idprg evaluate --set test` |
| 2026-10-31 | **the released closed form**: the 21-parameter `model/closed_form.json` (`3f7f7061…`), frozen beforehand and scored once with `idprg evaluate --model closed-form --set test`; the network was not rescored |

The released models' scores are in `results/metrics.json` under `test` and, per protein, in
`results/network_test.csv` and `results/closed_form_test.csv`. The 21-parameter form was chosen on cross-validation and the
extrapolation sets after the 24-parameter form's test score was known, so it is not entirely blind. Nothing was fitted or
selected after the last opening, and the test split is not to be scored again.
