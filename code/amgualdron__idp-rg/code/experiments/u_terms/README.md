# Closed-form u-term tests

Can a few extra terms remove the closed form's hard edge in the nu-vs-size plots (short chains with excess like charge,
and h) without giving up parameter economy? Nothing here changes `model/` or `results/`. Started 2026-10-30, paused
before these fits ran, to continue on a larger machine. This directory is working material, not part of the release.

**The test split is sealed.** `fit.py` cannot reach it and `compare.py` exits if it sees its rows. Choose on
cross-validation; validation is scored once, for one finalist (it has already compared s+/N and k=6).

## Already known (fits in `reference/`; all against the published 24-parameter form)

| Variant | Parameters | Result |
|---|---:|---|
| u + s+/N (slope about 21) | 25 | CV -0.00022 NLL (3.3 SE). Validation -0.00037 (2.8 SE), nu r 0.787 to 0.816; s+ > 0.3, N <= 100: sigma error 7.33 to 4.61%. Length set -0.0027 (passes). **Contact-energy set +0.081, fails the c1 rule's +0.001 limit**: protein 761 goes from 2.47 to 25.7, and it is still +0.0087 without that protein, median +0.0022. |
| No ln S_m products; u + s+/N, s-^2, s+ h, (ln N)^2, s+ <lambda>, s+ ln N, s- h | 27 | Validation -0.00012 (0.1 SE), nu r falls to 0.742, nu reaches 0.953. Rejected. |
| u + s+ N, 1/N or s+ ln N | 25 | CV -0.00006, +0.00000, -0.00009: all weaker than s+/N. |

Three things learned the hard way:
- Least squares on the free-fit u is a poor guide: it promised R2 0.63 for the 27-parameter form, and the likelihood fit
  did not follow. Judge on likelihood fits only.
- Extra terms can make L-BFGS-B report false convergence at iteration 1 (first step sends the objective to about 1e298).
  `fit.py` rescales the first step, restarts, and aborts at 20 iterations or fewer. `idprg/closed_form.py`'s `fit` has
  no such guard.
- Training residuals of s+/N show u too low by 0.14-0.19 at h > 2 at every length, and in training only 2 chains with
  N <= 50 have h > 0.5. So h/N may be the wrong shape, and Delta is left free to take either sign. A Delta near zero
  would point to curvature in h (a saturating term), which is not tested here.

## The tests

Terms are `x (N/100)^-Delta`, so Delta and the coefficient stay separate. In every model ln A and the CE threshold are
unchanged.

| Model | u relative to the published form | Parameters |
|---|---|---:|
| `check` | + s+/N | 25 |
| `dS` | + s+ N^-Delta, starting at 1 | 26 |
| `dH` | + s+/N, plain h replaced by h N^-Delta, starting at 0 (the plain term) | 26 |
| `hN` | + s+/N + h/N | 26 |
| `dSH` | both Deltas fitted, plain h replaced | 27 |

Saturation tests, added after the fits above: s+/N overshoots at high s+ on short chains. All on the 23-parameter form
(s+/N with s+ ln N and h ln N removed), jobs `MODEL:standard:FOLD:scd_pos,h`. Both limits nest s+/N.

| Model | Term in u | Parameters |
|---|---|---:|
| `satN` | s+ / (1 + N/N0) | 24 |
| `satS` | s0 tanh(s+/s0) / N | 24 |
| `satSN` | s0 tanh(s+/s0) / (1 + N/N0) | 25 |

## Outcome (2026-10-31)

Adopted model `a1`, 21 parameters: s+/N added; s+ ln N and h ln N removed from ln S_m; ln N removed from u; s+ removed from
ln A. Removing s+ from ln A is what fixed the contact-energy set (+0.024 to -0.0098; protein 761 no longer stands out).
The next removal, u's s+ (`a2`), failed at +2.1 SE. Saturation in s+ was not wanted (s0 ran to hundreds or more); saturation in N
did not pass. c1 = 0.25 was reselected on the new form. `fit.py` keeps the published form as its base, so its results stay
reproducible; `reference/published/` holds that form and its scores. The validation split was scored for the 23-parameter
form before `a1` was chosen, so validation is not a clean check of `a1` either.

## Running

```bash
export IDPRG_FRAMES=/path/to/idp_rg_frames.npy      # not in the repo: ml-idps data/idp_rg_frames.npy
export JOBS=4                                       # fits at once; a fold needs about 4 GB, a full-training fit about 6 GB
R=experiments/u_terms/run.sh

$R smoke          # under a minute; outputs in runs/u_terms_smoke are meaningless, it only checks the plumbing
$R u              # runs `check` first and stops if it misses the reference; then 12 fits
python experiments/u_terms/compare.py
```

**`check` has never completed**: the first attempt was stopped during it. It must reproduce objective 1.145896 on fold 0
(within 5e-5). If it fails, no other result from `fit.py` counts. Each fit took about 5-10 minutes on 12 cores.

Then, keeping the best u by the rule below (call it MODEL):

```bash
$R drop MODEL scd_pos                      # step 2: remove one product, score, and only then remove the next
$R drop MODEL scd_pos,h
$R drop MODEL scd_pos,h,mean_hydro
$R drop MODEL scd_pos,h,mean_hydro,scd_neg
python experiments/u_terms/compare.py --against MODEL
$R final MODEL scd_pos                     # step 3: length and contact-energy sets, with the products removed so far
$R validation MODEL scd_pos                # step 4: only if step 3 passes, and only once
```

The removal order is weakest first, from the spread of each product's contribution over the three s+/N folds (log units):
s+ ln N 0.001-0.025, h ln N 0.010-0.031, <lambda> ln N 0.024-0.042, s- ln N 0.044-0.072. s- ln N is not small. Stop at the first
removal that fails.

## Rule

Adapted from the c1 selection in `paper/closed_form.tex`:
- **Adding a term:** pooled out-of-fold NLL lower than its simpler parent by more than 2 standard errors, paired per
  protein with one sequence cluster as one unit.
- **Removing a product:** out-of-fold NLL not higher than its parent by more than 2 standard errors.
- **Extrapolation:** NLL on the length and contact-energy sets not higher than the published form by more than 0.001.
  s+/N fails this on the contact-energy set, so that is the check most likely to decide the outcome.

Also check Delta across the three folds: if it is unstable, the term is not well determined. These terms were chosen
after seeing training residuals, so small gains are slightly optimistic.

## If a form wins

Adopting it means changing `model/closed_form.json` and `idprg/closed_form.py`, and the paper's numbers via
`paper/make_numbers.py`. The closed form's test scores in the paper would then describe the old formula, and with the test
split sealed the new one would be reported on cross-validation and validation only. That is a decision to make
deliberately, not a by-product of a good fit.

Outputs go to `runs/u_terms/`, which git ignores: copy them back to keep them.
