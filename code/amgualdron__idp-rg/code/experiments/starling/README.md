# STARLING on the experimental benchmark

`run_starling.py` samples STARLING ensembles for `data/benchmark.csv` and records R_g for every conformation. The aim is
the distributions: mean R_g from Mpipi-GG simulations is already in `data/benchmark_methods.csv`. The script does not
touch `model/` or `results/`. Output goes to `runs/starling/`, which is tracked, unlike the rest of `runs/`. The settings below were
chosen from small tests: HMPV, A1_LCD_10R and Tau_ht23 for step counts and speed, and six proteins for the comparison
with HPS1.

## The full run

Finished on 2026-10-31 at the default settings, in `runs/starling/benchmark/` (4.6 MB): 121 proteins with 2000
conformations each, and Tau_ht40 recorded as `too_long`. It took 2 h 15 min. No protein failed or needed a smaller batch,
and no conformation failed STARLING's distance check. Each histogram holds all 2000 samples and integrates to 1.

| mean R_g, 121 proteins | MAE (nm) | bias (nm) | mean absolute error, % | Pearson r |
|---|---:|---:|---:|---:|
| STARLING vs Mpipi-GG | 0.053 | -0.043 | 1.6 | 0.999 |
| STARLING vs SAXS | 0.327 | +0.173 | 9.7 | 0.940 |
| Mpipi-GG vs SAXS | 0.342 | +0.216 | 9.9 | 0.941 |
| HPS1 vs SAXS | 0.385 | +0.193 | 12.6 | 0.884 |

STARLING stays within 6% of Mpipi-GG on every protein. The worst cases are TRF2_BR (+5.9%), RpII215_portz (-4.4%),
stath_NTD (+4.2%) and Tau_K23 (-4.1%). The shape of all 121 is scored under "Shape: network, STARLING and HPS1".

## Setup and running

    mamba env create -f experiments/starling/environment.yml

    # check the input and print the plan, without using the GPU
    mamba run -n starling python experiments/starling/run_starling.py --dry-run

    # the set: 2000 conformations per protein, 45 DDIM steps, float16 denoiser, compiled; about 2.3 h
    mamba run -n starling python experiments/starling/run_starling.py --out runs/starling/benchmark

The first run downloads 1.9 GB of weights to `~/.cache/torch/hub/checkpoints` and verifies them against the release's
SHA-256 digests. The script's docstring describes the output files. A run can be interrupted and restarted with the
same command: it continues from the last finished protein. Use a separate `--out` for each set of settings; the script
refuses to mix them.

## Why 2000 conformations

STARLING's samples are independent, because each starts from its own noise. MD frames are not, so the number to match
is HPS1's independent samples, not its frames. The trajectories of all 122 benchmark proteins (`runs/adopt/traj`, links
into the idp-ml runs) have a median of 3750 frames and an integrated autocorrelation time of about 2 frames, which leaves
these independent samples per protein:

| | median | 10th to 90th percentile | range |
|---|---:|---:|---:|
| all 122 | 2065 | 1258 to 2456 | 308 to 6631 |
| N up to 50 (6 proteins) | 788 | | 308 to 1079 |
| N over 300 (15 proteins) | 2406 | | 1795 to 2552 |

The training set's rule of about 170 independent samples per run (`data/README.md`) describes other runs and does not
apply to these. `data/timeseries/` holds only 27 of the trajectories, and its hNHE1cdt has 3000 of the 6750 frames.

So 2000 conformations matches the HPS1 benchmark trajectories at the median. 69 of the 122 proteins have more than 2000
independent samples, most only slightly: the highest are alpha_syn and redAFP, which pool several trajectories.

| samples | error of mean | error of SD | error of skew |
|---:|---:|---:|---:|
| HPS1, 2065 independent | 0.44% | 1.6% | 0.054 |
| STARLING, 2000 | 0.45% | 1.6% | 0.055 |

The percentages assume SD = 0.2 mean, as the test proteins roughly have. Fewer conformations are a real option, because
a comparison with HPS1 also carries HPS1's noise. 1500 conformations adds 8% to that noise and takes 1.75 h; 1000 adds 23%
and takes 1.2 h.

## Why 45 steps

The table compares each step count with 200 steps, on the same 500 starting noises, so the differences are measured
tightly. The ranges cover HMPV, A1_LCD_10R and Tau_ht23. Times are plain sampling, 2000 conformations.

| steps | model evaluations | SD | mean | W1, % of mean | 95th percentile | time |
|---:|---:|---:|---:|---:|---:|---:|
| 30 (STARLING's default) | 31 | -2.1 to -3.4% | -0.8 to +0.1% | 0.4 to 0.8 | -0.9 to -1.5% | 2.7 h |
| **45** | **46** | **-0.6 to -1.3%** | **-0.6 to +0.4%** | **0.2 to 0.6** | **0 to -0.8%** | **3.7 h** |
| 50 | 50 | +1.5 to +2.6% | -1.2 to +1.4% | 0.8 to 1.4 | 0 to +1.6% | 3.9 h |
| 60 | 63 | -0.5 to -0.8% | -0.4 to +0.2% | 0.2 to 0.4 | 0 to -0.5% | 4.8 h |
| 66 | 67 | 0 to +0.6% | -0.5 to +0.5% | 0.2 to 0.5 | -0.4 to +0.4% | 5.1 h |
| 75 | 77 | +0.9 to +1.6% | -0.6 to +0.7% | 0.5 to 0.7 | -0.2 to +1.0% | 5.8 h |
| 90 | 91 | +0.6 to +1.2% | -0.5 to +0.5% | 0.4 to 0.5 | -0.1 to +0.8% | 6.7 h |
| 100 | 100 | +0.8 to +1.4% | -0.5 to +0.6% | 0.4 to 0.6 | -0.1 to +0.9% | 7.3 h |
| 200 | 200 | reference | reference | | | 14 h |

- From 45 steps up, except 50, the SD stays within HPS1's sampling error (1.6%), with mean differences about that error
  (0.44%). 45 is the cheapest of these, and 100 is no closer to 200 than it is.
- Accuracy does not rise steadily with steps, because the requested count also moves where the schedule starts. The
  spacing is `1000 // steps`, so the highest timestep is 981 for 50 steps, 989 to 993 for the others, and 996 for 200.
  The outlier at 50 fits this, but the cause was not tested further.
- 200 steps is not proven to be the truth, only the most finely resolved run that works. DDPM, the full-schedule reference,
  gave SDs of 0.557 and 1.382 nm at 500 conformations (A1_LCD_10R, Tau_ht23), known to about 8% when compared. That is
  consistent with every working step count from 30 up and cannot separate them.

## Speed

The denoiser takes 92% of the time at 100 steps; decoding takes 8 ms per conformation and everything else under 1 ms. Two
changes to the denoiser are on by default on the GPU, and `--plain` turns them off:

| at 45 steps, 2000 conformations | ms per conformation (N = 120 to 352) | full set |
|---|---:|---:|
| plain | 54 to 64 | 3.7 h |
| float16 denoiser | 40 to 52 | 2.8 h |
| float16 denoiser and `torch.compile` | 34 to 42 | **2.3 h** |

Neither changes the distributions by a meaningful amount. With the same noise and 500 conformations of Tau_ht23, the final
pipeline matches plain 45-step sampling to W1 of about 0.002 nm (the error of its mean is 0.06 nm). At 100 steps, on all
three proteins, the worst single-conformation difference was 0.013 nm. On six proteins at 2000 conformations, the final
pipeline against plain 100 steps gives SD ratios of 0.977 to 0.985, as expected from the table above.

Decoding stays in float32, because squared distances could overflow float16. Compilation uses `mode="default"`, because
`reduce-overhead` is left unusable for the whole process by a GPU out-of-memory error, so the script's retry would fail.
It needs a C compiler and may print harmless gcc warnings. A new sequence length adds no measurable time, and memory stayed
flat over 18 lengths.

Rejected:

| option | why |
|---|---|
| PLMS sampler, 20 to 50 steps | means off by up to 5.5%, SDs by up to 13% |
| batch 50 | about 6% faster, with less room on the GPU |
| `reduce-overhead` compilation | about the same speed; fails after an out-of-memory error |
| DDPM sampler | 1.1 to 1.4 s per conformation; out of memory at batch 25 |

## Resources

Measured on an RTX 3060 Ti (8 GB), Python 3.11, torch 2.6.0+cu124, STARLING 2.0.2, at batch size 25. With the defaults,
sampling takes 27.5 + 0.043 N ms per conformation (fitted to six proteins, with the first one paying about 20 s for
compilation, so slightly pessimistic). The totals apply this to the 121 lengths, plus about 1 s per protein.

| resource | use |
|---|---|
| GPU memory | 5.5 GB in `nvidia-smi`, including 0.7 GB for the desktop |
| GPU | 98% utilisation, 185 W average, 73 C at most, measured over 7 min; about 0.4 kWh for the run; the desktop may feel slow |
| CPU, RAM | one core; 5.1 GiB at the peak, on Tau_ht23 with 2000 conformations |
| disk, R_g only | under 10 MB |
| disk, `--save-ensembles` | about 6 GB at 2000 conformations, plus about 80 min of compression (extrapolated from Tau_ht23 as N^2) |
| already installed | environment 6.0 GB, weights 1.8 GB |

## Check against HPS1

Six proteins at the default settings, chosen to span length and charge rather than at random. `plot_vs_hps1.py` makes
the figure and `vs_hps1.csv` in `runs/starling/check_hps1_fast/`:

    mamba run -n starling python experiments/starling/plot_vs_hps1.py runs/starling/check_hps1_fast

| protein | N | mean vs HPS1 | mean vs Mpipi-GG | SD ratio to HPS1 | shape W1 | noise, 95th percentile | HPS1 vs other HPS1 proteins |
|---|---:|---:|---:|---:|---:|---:|---:|
| stath_NTD | 15 | +5.4% | +4.2% | 1.13 +- 0.05 | 0.010 | 0.012 | 0.051 |
| DSS1 | 71 | -1.5% | -0.1% | 1.04 +- 0.03 | 0.011 | 0.011 | 0.021 |
| P1_100 | 100 | +6.8% | -0.3% | 0.99 +- 0.02 | 0.012 | 0.010 | 0.015 |
| p150L_342_475 | 144 | -2.2% | -1.3% | 0.98 +- 0.02 | 0.004 | 0.010 | 0.014 |
| DeltaN_FATZ_1 | 209 | -0.5% | -2.3% | 0.98 +- 0.02 | 0.006 | 0.010 | 0.014 |
| Tau_ht23 | 352 | +14.0% | -2.2% | 1.05 +- 0.02 | 0.016 | 0.010 | 0.024 |

Shape W1 is between R_g / <R_g> distributions. Noise is W1 between resamples of HPS1's own shape at 2000 and at the
protein's independent sample count. The last column is the median shape W1 to the other five proteins.

- **The size follows Mpipi-GG, not HPS1,** as it should for a model trained on it. The large offsets from HPS1 (P1_100,
  Tau_ht23) are the offsets between the two force fields. Use these histograms as Mpipi-GG-like distributions.
- **The shape agrees closely,** at or near the noise for five of the six; Tau_ht23 is 1.6 times the noise. The SD is within
  5% for N of at least 71, and 13% wide for the 15-residue stath_NTD.
- **This is weak evidence of sequence-specific accuracy.** Different proteins differ only 1.3 to 5 times as much as
  STARLING differs from HPS1, so much of the agreement is a shape that disordered chains share. Stronger tests would
  compare proteins with unusual shapes, and Mpipi-GG distributions if they can be found.

## Shape: network, STARLING and HPS1

`score_shape.py` divides each distribution by its own mean and compares the released network (`model/network.eqx`),
STARLING and HPS1 by Wasserstein-1, so W1 reads as a percentage of <R_g>. It writes `shape_w1.png`, `shape_w1.csv` and,
when redirected as here, the report `shape_w1.txt` to the run directory. The network's size-included W1 reproduces the
values in `paper/figure_metrics.json`, and no benchmark sequence is in the network's MD data.

    mamba run -n idp-proteome python experiments/starling/score_shape.py runs/starling/benchmark > runs/starling/benchmark/shape_w1.txt

| 121 proteins, shape W1 (%) | median | 95% interval, resampling families | 90th percentile | noise floor |
|---|---:|---:|---:|---:|
| network vs HPS1 | 0.51 | 0.49 to 0.55 | 1.11 | 0.31 |
| network as 2000 samples vs HPS1 | 0.66 | | 1.16 | 0.53 |
| network vs STARLING | 0.92 | 0.70 to 1.06 | 1.61 | 0.42 |
| STARLING vs HPS1 | 1.05 | 0.74 to 1.27 | 1.80 | 0.53 |
| one shared shape vs HPS1 | 0.56 | 0.48 to 1.02 | 1.92 | 0.31 |

- **Compare STARLING with the network at equal noise.** The network is an exact curve and STARLING is 2000 samples. As 2000
  samples, the network is closer to HPS1 on 82 of 121 proteins, by a median of 0.25 (interval -0.003 to 0.51). The exact
  curve's 101 of 121 overstates the gap.
- **Most of these proteins have nearly the same shape.** A single shape from other families' HPS1 trajectories matches
  the network at the median. Among the 66 proteins measurably away from it, the network is closer on 53, or 51 as 2000
  samples, and STARLING on 42. That selection favours both predictors somewhat.
- **STARLING tracks HPS1's protein-to-protein differences less well.** Spearman correlation with HPS1 for N of at least 30:
  coefficient of variation 0.74 for the network and 0.38 for STARLING, skewness 0.63 and 0.21.
- **The force field is part of the gap.** The network is trained on HPS1 and STARLING on Mpipi-GG. STARLING vs HPS1 grows
  from 0.72 to 1.43 as the size offset goes from under 5% to over 20%, while the network stays near 0.5. Separating force
  field from model error needs Mpipi-GG distributions.
- **p27_v78 and p27_v56 are real misses.** Their HPS1 trajectories look stationary and are narrower than both predictors
  (coefficient of variation 0.13 and 0.16; network 0.22; STARLING 0.18 and 0.22).

The noise floor is calibrated on HPS1 split-half W1, because n_eff from R_g's autocorrelation overstates W1's noise: the
factor is 0.74 on W1, and STARLING's independent samples need none (ratio 1.00). Two proteins, stath_NTD and RS-peptide,
have N below the network's stated minimum of 20 and are drawn open in the figure.

`fit_nu.py` fits the two-exponent distribution to each ensemble, over its own mean and with ln A free (the fit paper fig24 uses
for MD), into `nu_fits.csv`, with the sd of nu over 20 resamples. The paper text quotes it against the network's nu (`sec:starling`).

## Things to know before using the numbers

- **Tau_ht40 cannot be run.** N = 441 is over STARLING's limit of 380, so the set is 121 proteins. It stays in
  `summary.csv` with status `too_long`. Sampling a fragment would give an R_g that cannot be compared.
- **STARLING 2.0.2's DDIM silently breaks from 334 steps up, so the script refuses those values.** The schedule then
  includes t = 999, where alpha_bar = 2.4e-10 and the estimate of the clean sample is amplified 64,000 times. At 400 steps,
  Tau_ht23's mean went from 6.7 to 9.2 nm and A1_LCD_10R's from 2.6 to 1.9 nm, with no error from STARLING. Its check for
  impossible distances flagged only 2 to 6 of 500 conformations (none at 50 or 200 steps), so that check does not catch this.
- **STARLING 2.0.2 samples in training mode, and this script corrects that.** In that mode `SequenceEncoder.forward`
  zeroes the ionic-strength embedding at random. This hit 3 of 20 batches in the test (the code's rate is 20%). The script
  calls `.eval()` on the models, so results differ from a plain `starling.generate()`. Always masking moved Tau_ht23 by
  0.12 nm at 20 mM (from 6.45 to 6.57), while at 150 mM all three test proteins stayed within 0.021 nm, about the
  sampling error. Other STARLING numbers obtained with this version away from 150 mM may carry this bias.
- **STARLING is trained on Mpipi-GG, not HPS1.** Compare it with the rest of this repository as a different force field.
  `summary.csv` carries `saxs_rg_nm`, `hps1_rg_nm` and `mpipi_gg_rg_nm`. SAXS gives only a mean, so it cannot test the histograms.
- **Salt is 150 mM for every protein,** the same as the HPS1 simulations, not each experiment's own condition. STARLING
  was trained at 20, 150 and 300 mM only, so the script rejects anything else.
- **Sampling is reproducible given the seed, device, batch size, step count and `--plain`.** Each protein's seed depends
  only on `--seed` and its name, and rerunning in a new process and another order gave bit-identical R_g. An out-of-memory
  retry lowers the batch size, and `summary.csv` records the batch size used.
- No conformation in the full run failed STARLING's check for impossible distances, so `rg_mean_clean_nm` equals
  `rg_mean_nm`. The 15-residue stath_NTD runs without complaint, but its accuracy is the least certain.
- `sequence` is used as it stands, including the 131-residue hNHE1cdt. The run read an older copy of `benchmark.csv` with
  the same names, sequences, SAXS and HPS1 values as the current table, so its results apply unchanged.
- `psutil` is pinned in `environment.yml` because the `starling` command imports it without declaring it.
- The documentation lists a `dpmpp` sampler, but 2.0.2 has `ddim`, `ddpm` and `plms`.
