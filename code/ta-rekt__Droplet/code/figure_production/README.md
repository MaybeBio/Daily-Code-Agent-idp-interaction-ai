# Paper figures

One folder per figure, each holding its notebook, the code that feeds it, and
the PDF it produces. `processed_data/` holds the reduced data of Figs. 1, 2, 3 and 5.

```
figure_production/
├── README.md             this file
├── FIGURE_STYLE.md       typography and geometry rules -- read before changing a figure
├── processed_data/       small CSV/JSON, committed, enough to rebuild Figs. 1, 2, 3 and most of 5
├── fig1/   fig1.ipynb · prepare_fig1_data.py · art_{a,b,c,d}.{svg,pdf} · art_bc_labels.json · extract_art.py · fig1.pdf
├── fig2/   fig2.ipynb · prepare_fig2_data.py · fig2.pdf
├── fig3/   fig3.ipynb · prepare_fig3_data.py · fig3.pdf
├── fig4/   fig4.ipynb · run_fig4_solver.sh · fig4_cartoon.svg · fig4.pdf
├── fig5/   fig5.ipynb · prepare_fig5_data.py · run_fig5_rouse.sh · rouse_toy.jl ·
            art_*.svg · extract_art.py · fig5.pdf
└── figs1/  figs1.ipynb · figs1.pdf
```

The notebooks are **hand-maintained**: edit them directly. Nothing generates
them.

Each notebook finds its own way around the repo — it locates this directory (or
the repo root, for Fig. 4) by looking for marker files, then `chdir`s there. So
it does not matter whether the kernel starts in the notebook's folder, at the
repo root, or at some ancestor. If the search ever fails it raises with
instructions; each has an override constant near the top of its first cell
(`DIR_OVERRIDE`, or `REPO_OVERRIDE` in Fig. 4).

## Requirements

Python with `numpy`, `scipy`, `pandas` and `matplotlib`, and a LaTeX
installation — the figures use `text.usetex`, so `latex`, `dvips`, Ghostscript
and **`cm-super`** must be present. Without `cm-super` the notebooks print a
warning and silently fall back to mathtext: the numbers stay right, the
typography does not.

Figs. 1, 4 and 5 additionally need `pypdf` (`pip install pypdf`) and an SVG converter,
either `cairosvg` or Inkscape, to composite their artwork. See below.

## Fig. 1 — lattice Monte Carlo model

The remake of `figs/fig1.svg` of the paper canvas, at `\textwidth`: 481.9 ×
486.1 pt. **Include it with `width=\linewidth`** — the old composite was placed at
`0.95\linewidth`, which would now shrink the 9 pt type.

**B and C are pasted** from the source canvas — the published matplotlib
rendering of the two 3D configurations (c = 0.003 run 13, c = 0.017 run 2 of the
20260423 sweep). Re-plotting them depends on the matplotlib version (3D depth
shading made the ligands vanish on a newer one). Only their tick labels are
re-set, at the published positions (`art_bc_labels.json`).

**Data** (`processed_data/fig1_*`, committed) for E come from one run of the
20260423 proximal-promoter sweep (the raw output Fig. 2 uses), and for F from
the source SVG:

```
python figure_production/fig1/prepare_fig1_data.py --data-root ~/Archive/Work/droplet_data
```

| panel | run | notes |
|---|---|---|
| E | c = 0.006, `output_cpu3` | frames 0–40000, transcription smoothed over 4000 frames (cell 155) |
| F | — | a hand-drawn sketch; its step traces, shaded intervals and markers are read back from the SVG's path data |

Cell 155 says `run = 12`, but the published trace is `output_cpu3` (its activity
peaks at 0.58 near t = 1.71 × 10⁵ inside the long distance dip; no other run
does). Its smoothing window is 4000 frames ≈ 3.2 × 10⁴ MC samples, not the
"500-MCS" of the old caption.

**Artwork.** `extract_art.py` cut the drawing of panels A–D out of the
source SVG by element id, and converts each piece to PDF with Inkscape once:
the notebook stamps the committed `art_*.pdf`, because SVG renderers disagree on
this file's gradients (cairosvg drew D's glows wrongly on another machine). D's
clouds and glows are also made opaque (pre-composited over what lies under
them): macOS Preview fills the transparent gradients of an Inkscape/cairo PDF
as solid discs (all lettering there is glyph outlines, and is dropped).
The pieces keep the whole source canvas, and the notebook places every label in
the same canvas units, so text and drawing stay registered. The source SVG had
stale `sodipodi:type="arc"` attributes on several circles; Inkscape redrew those
from the stale centres (D's clouds vanished off canvas, a stray dot appeared by
A's "Off"), so the extraction strips all `sodipodi:*` attributes. It also re-centres D's three
small contact glows, which every SVG renderer draws about a unit away from the
published Illustrator PDF (offsets measured, in `GLOW_SHIFT`). Source files now
live in `Tarek_Droplets/fig_sources/`.

**Changes from the source:** E is as tall as F and its time axis starts at 0;
E has a full box; F's rows moved down a few points to give its headers two rows;
Anchor/Probe and D's open-polymer labels moved apart; the F letter now aligns
with C and D.

## Fig. 2 — concentration sensitivity

The committed data is enough:

```
jupyter lab fig2/fig2.ipynb      # run all cells -> fig2/fig2.pdf
```

To rebuild `processed_data/fig2_*` from the raw lattice Monte Carlo runs
(not in the repo — pass the folder holding `20260423/`):

```
python figure_production/fig2/prepare_fig2_data.py --data-root ~/Archive/Work/droplet_data
```

It writes to `figure_production/processed_data/` by default (`--out` overrides).
Values are written at full precision (`%.17g`) and read back with
`float_precision="round_trip"`, which is load-bearing: the smoothing parameter
in the shape fits is chosen by leave-one-out CV on an ill-conditioned grid, and
a one-ULP difference in the input selects a visibly different fit.

Goes in a `figure*` (full width).

## Fig. 3 — sensitivity of a single probe

Same pattern:

```
jupyter lab fig3/fig3.ipynb      # -> fig3/fig3.pdf
python  figure_production/fig3/prepare_fig3_data.py --data-root ~/Archive/Work/droplet_data   # optional, rebuilds the data
```

**Fig. 3 is single column.** It must go in a `figure`, not a `figure*`: it is
built 231.95 pt wide, and a `figure*` scales it ~2x, so its 9 pt type prints at
nearly 18 pt.

## Fig. 4 — surface condensation on a spherical anchor

This one cannot be run from committed data: panels B–G read the output of eight
Julia BVP solves, which are too large to commit. Four steps, in order.

**1. Install Julia 1.11.1.** `Manifest.toml` was resolved with exactly that
version, so pin it rather than tracking `release`:

```
curl -fsSL https://install.julialang.org | sh    # then open a NEW shell
juliaup add 1.11.1
juliaup default 1.11.1
julia --version                                  # -> julia version 1.11.1
```

**2. Run the solver**, from the repo root (the directory with `Project.toml`):

```
bash figure_production/fig4/run_fig4_solver.sh
```

It instantiates the Julia project, runs the eight sweeps (~10–20 s each), and
then verifies the output: it prints each run's radius or surface attraction and
fails if any is missing or out of order. Results land in `data/20261004/` and
`data/20261007/`, which are gitignored.

**3. Install the two Python extras** panel A's artwork needs:

```
pip install pypdf          # pure Python, always required
pip install cairosvg       # needs the cairo C library: brew install cairo
```

If `cairosvg` will not install, Inkscape works instead — the notebook calls its
command line as a fallback and finds it inside the macOS app bundle. Without
either, the figure still builds, with panel A as text only and a printed note.

**4. Run the notebook:**

```
jupyter lab fig4/fig4.ipynb      # -> fig4/fig4.pdf
```

Panel A is assembled in two layers. The notebook draws all of its *text* in
Computer Modern at the caption size; the drawing itself is `fig4_cartoon.svg`,
converted to PDF and stamped **underneath** the matplotlib layer so the labels
sit on top of it. Labels are positioned through a `cart(u, v)` helper that maps
the SVG's own viewBox coordinates into the panel, so they follow the artwork if
it is resized or moved.

Goes in a `figure*` (full width).

## Fig. 5 — Rouse polymer with surface attraction

The remake of the old Fig. 4 of Tarek's draft (`figs/fig4.svg`,
`main-fig4.pdf`), with the model that figure used: the **June toy model**,
dimensionless, from `git show 9f78af2:spec/rouse_mfpt_attractive.jl` (Tarek,
2026-06-15) and `viz-rouse.ipynb` cells 23 and 37. Panels A (plot), C and E,
and the theory curves of D, are analytic. The histograms of panel D need eight
Rouse simulations, which are not committed. Without them the notebook still
runs, and panel D shows the theory alone with a red "simulation pending" note.

**The model.** Spheres of radius R = 2.5 (contact at d = R), box edge L = 10;
nb = 1, χ_b = 2.5, h₀ = ω_b = −0.5, κ = 1. λ = √(f_b''(φ∞)/κ),
U_ads(d) = −4πR³h₀² e^{−λ(d−R)} / (κ(λ + 1/R) d), and
U_tot = U_ads − 2 log d + 3d²/(2 s b²) with s = 99, b = 1. The sixteen
concentrations are φ∞ = linspace(0.003, 0.15, 16); low and high c are k = 0
and k = 9 (0.003 and 0.0912, as in cell 37), panel A's three potentials
k = 3, 9, 15. The ridge d*(c) is the first maximum of U_tot out from contact;
d^o = 1.2 max_c d*(c) and d^x = contact, as in the old caption.

**1. Julia 1.11.1**, as for Fig. 4.

**2. Run the Rouse simulations**, from the repo root:

```
bash figure_production/fig5/run_fig5_rouse.sh        # or: ... run_fig5_rouse.sh 4  (number of processes)
```

It instantiates the Julia project and runs `fig5/rouse_toy.jl` — the June
spec, changed only to run the two plotted concentrations with four seeds each
(43, 1043, 2043, 3043 and 52, 1052, …; the June seeds were 43 and 52) and to
store the pair potential too. N = 100 beads, dt = 5e-4, 5e7 steps, distance
kept every 500 steps. About 2 min per run and 0.4 GB RAM per process; the
number of processes is picked from cores and memory (override with the
argument). Output goes to `data/<today>/fig5_rouse_toy/` (gitignored); the
script ends by checking that all eight runs are there.

**3. Process:**

```
python figure_production/fig5/prepare_fig5_data.py
```

Computes the potentials and landmarks and writes `processed_data/fig5_*`. With
Rouse output present (the latest `data/*/fig5_rouse_toy/`, or `--rouse DIR`) it
adds `fig5_hist.csv`: 100 bins on [R, L], density over the in-range samples, no
burn-in, as cell 37. About 40 % (low c) and 20 % (high c) of the samples lie
beyond L — the chain's size is √99 ≈ 10 — and are left out, as in the original.
It checks that the potential each run felt is the plotted U_tot (to 1e-9) and
prints P(d < d*) from the simulation next to the theory.

**4. Run the notebook:**

```
jupyter lab fig5/fig5.ipynb      # -> fig5/fig5.pdf
```

The artwork (droplets, capillary bridge, Rouse chain, the inset sketch of
panel C, the polymer of panel E) was cut out of Tarek's SVG into
`fig5/art_*.svg` by `extract_art.py`, which is kept only as a record of how
(rerun it only if the source drawing changes). The pieces keep the source
canvas's coordinates in mm; the notebook's `GROUPS` table maps each onto the
page (scale 0.92, panel C's inset 1.45×), and all labels are placed in the
same canvas coordinates, so they move with the artwork. As in Fig. 4 the
artwork is stamped underneath the matplotlib layer. The capillary bridge's
glow, a filter effect in the source, is stored as a blurred bitmap inside
`art_a_bridge.svg`.

At high c the contact well (U_tot = −8.2 at R) is clipped below the frame in
C, as in the original, so the barrier stays visible. Panel E is redrawn on the
original's frame: x from 1.48 to 14.24 and y from −4.14 to −1.25 put d*, d^o
and the broad minimum where the original has them, so its curve runs past L
(`processed_data/fig5_panel_e.csv`, to d = 15) and climbs on the right; the
landmark lines take the original's light tints, d* stops at the barrier, and
the exit / re-entry arrows are cubic Béziers measured on the original. The schematic
panels mark contact as 2R, as the original did, while in the model contact is
at d = R (D's axis starts at 2.5).

Goes in a `figure*` (full width), 481.9 × 304 pt.

## Fig. S1 — mechanical analogy

The radial equation read as a unit mass rolling in the inverted potential
−f(φ)/κ, with "time" r and drag 2/r. Three panels, c/c_sat = 0.2, 0.5, 0.8
(nearest points of the Fig. 4D sweep), at the parameters of Fig. 4D. Each
panel shows the starting position φ(R) for R = 20, 50 and 3000 nm (the
`small_R`, `medium_R` and `large_R` runs, colours as in Fig. 4D), the initial "speed"
−h₀/κ as an arrow of fixed length (not to scale; it is the same for every R),
and the end point φ∞ at the top of the curve, enlarged in an inset because
φ∞ ~ 1e-5 sits on the y axis at this scale.

It reads the same solver output as Fig. 4, so run
`bash figure_production/fig4/run_fig4_solver.sh` first if `data/20261004/` is
missing, then:

```
jupyter lab figs1/figs1.ipynb    # -> figs1/figs1.pdf
```

Built at 481.9 pt: goes in a `figure*`.

## Changing a figure

Read `FIGURE_STYLE.md` first. The short version: every figure is built at the
width it is included at, so the scale factor is 1 and the sizes in the code are
the sizes on the page — 9 pt text, 8 pt tick labels, Arial 12 pt panel letters,
0.5 pt strokes throughout. Never save with `bbox_inches="tight"`: it re-crops
the page, LaTeX then rescales it, and the type sizing is undone.

After changing anything, check the output rather than the settings:

```
python - <<'PY'
import pdfplumber, collections
p = pdfplumber.open("fig4/fig4.pdf").pages[0]
scale = 481.9 / p.width                      # what LaTeX will apply
for s, n in sorted(collections.Counter(round(c["size"], 2) for c in p.chars).items()):
    print(f"{s:6.2f} pt in file -> {s*scale:5.2f} pt printed  (n={n})")
PY
```
