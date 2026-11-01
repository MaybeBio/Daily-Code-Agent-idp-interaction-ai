# idprg: the radius-of-gyration distribution of disordered proteins from sequence

`idprg` predicts the full equilibrium distribution P(R_g) of the radius of gyration of an intrinsically disordered
protein (IDP) from its amino-acid sequence. Mean, width, skewness and tails are all included. The predictions
reproduce HPS1 coarse-grained molecular dynamics at 298 K and 150 mM salt.

Two predictors share one output distribution:

- **The network** (default): an MLP on nine sequence descriptors, trained by maximum likelihood on the MD frames of
  3,223 proteins.
- **The closed form**: 21 parameters, linear in four descriptors, ln N and s₊/N, fitted the same way. It is written out in
  full in `paper/closed_form.tex`.

The output is the two-exponent distribution of Victor and Lhuillier. With s = R_g / S_m:

    p(s) ∝ exp(−A [s^δ/δ + s^−g/g]),   δ = 1/(1−ν),   g = 3/(3ν−1)

S_m is the most probable R_g, ν ∈ (1/3, 1) is an effective scaling exponent that sets the asymmetry, and A sets the
width.

## Install

```bash
pip install -e .            # numpy, scipy, pandas, jax, equinox, optax
pip install -e ".[paper]"   # + matplotlib, for paper/make_figures.py
```

## Use

```bash
idprg predict MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA
idprg predict --fasta proteins.fasta --out predictions.csv
idprg predict --fasta proteins.fasta --closed-form --density densities.csv
```

```
name,N,mean_rg,std_rg,skew,exkurt,nu,S_m,lnA
seq1,140,3.3515,0.72294,0.73944,0.63618,0.57867,3.0401,1.1484
```

- `mean_rg`, `std_rg` and `S_m` are in nm. `S_m` is the most probable R_g.
- `--density` also writes P(R_g) (nm⁻¹) on a 200-point grid for each sequence.
- Only the 20 standard amino acids are accepted.

From Python:

```python
import idprg
from idprg import distribution

df = idprg.predict(["MDVFMKGLSK...", "GSHMSSG..."])            # closed_form=True for the formula
p = distribution.pdf(rg_nm, df.lnA[0], np.log(df.S_m[0]), df.nu[0])
```

## Accuracy

All errors are measured against the HPS1 simulations of proteins that share no 30%-identity sequence cluster with a
training protein. ⟨R_g⟩ and σ errors are mean absolute relative errors. W₁ is the median Wasserstein-1 distance
between the predicted and simulated distributions, as a percentage of ⟨R_g⟩.

<!-- numbers:accuracy -->
|  | proteins | ⟨R_g⟩ error | σ error | W₁ | NLL per frame |
|---|---|---|---|---|---|
| network, validation | 695 | 1.84% | 4.73% | 1.47% | 1.1525 |
| closed form, validation | 695 | 2.53% | 5.81% | 1.79% | 1.1620 |
| network, test | 691 | 2.02% | 4.69% | 1.53% | 1.1379 |
| closed form, test | 691 | 2.53% | 5.09% | 1.82% | 1.1462 |
<!-- /numbers -->

The test split was scored once for each released model, after that model was fixed (`data/README.md`). Further results are in `MODEL.md`:
- length and stickiness extrapolation;
- an experimental benchmark of 122 proteins with SAXS data;
- ablations;
- the domain of validity.

The model predicts what HPS1 gives, not experiment. Against SAXS its error is essentially that of the force field.

## Reproduce

Everything below needs the MD frame store (`data/README.md`): pass `--frames PATH` or set `IDPRG_FRAMES`.

```bash
idprg evaluate --model network --set val              # validation scores -> results/
idprg train --split standard --seed 42 --out runs/net # retrains the released network (bit-identical weights)
idprg fit-closed-form --split standard --out runs/cf  # refits the closed form
idprg fit-md --out runs/md_fits.csv                   # per-protein maximum-likelihood fits (data/md_fits.csv)
idprg ablations                                       # results/ablations.csv (resumes where it stopped)
python paper/make_figures.py --frames FRAMES --trajectories DIR
python paper/make_numbers.py                          # paper/numbers.tex and the tables in MODEL.md, README.md
python -m pytest tests/                               # or python tests/test_reproduce.py; no frames needed
```

Every number in the write-ups comes from `results/`. After rerunning anything there, run `python paper/make_numbers.py`;
the tests fail while a generated file is out of date.

On a CPU:
- training takes about 1 minute;
- a closed-form fit takes about 10 minutes;
- the ablation grid takes about 1.5 hours with the default two parallel jobs.

## Repository

| path | content |
|---|---|
| `idprg/` | the package: `descriptors`, `distribution`, `network`, `closed_form`, `data` (tables, splits), `metrics`, `mdfit`, `ablations`, `cli` |
| `model/` | the released network (`network.eqx`, `network.json`) and closed form (`closed_form.json`) |
| `data/` | protein table, per-protein MD fits, benchmark; `data/README.md` describes the simulations, cleaning and split |
| `results/` | `metrics.json` (all summary scores), per-protein scores `<model>_<set>.csv` (sets: `val`, `test`, `length`, `sticky`, `benchmark`, and `cv` for the closed form), `ablations.csv`, `ce_threshold_cv.csv` |
| `paper/` | the write-ups: `summary.tex` (the model and its results, end to end), `formulas.tex` (every formula and constant), `closed_form.tex` (the closed form on its own), `captions.tex` and `figures/`; `numbers.tex` holds every number they quote, written by `make_numbers.py` |
| `MODEL.md` | model card |
