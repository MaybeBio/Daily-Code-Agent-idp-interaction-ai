# Figure style rules

How figures in `figure_production/` are built so that type in the figure matches
type in the paper. Written for `Tarek_Droplets`, `\documentclass[twocolumn]{revtex4}`
on A4.

## 1. Measure the target, never assume it

Everything below follows from four numbers, all measured from the compiled
`main.pdf` rather than guessed:

| quantity | value | how it was measured |
|---|---|---|
| body text | 10 pt | most common glyph size on a text page |
| **caption** | **9 pt** (CMR9 / CMBX9) | glyph size of a caption line |
| `\textwidth` (a `figure*`) | **481.9 pt** = 6.693 in | x-extent of body text on a page |
| `\columnsep` | 18 pt | widest interior gap with no glyphs, on a two-column text page |
| `\columnwidth` (a `figure`) | **231.95 pt** = 3.222 in | (`\textwidth` − `\columnsep`)/2 |

Re-measure after any class or geometry change:

```python
import pdfplumber, collections
p = pdfplumber.open("main.pdf").pages[4]
print(collections.Counter(round(c["size"], 2) for c in p.chars).most_common(5))
body = [c for c in p.chars if abs(c["size"] - 9.96) < 0.2]
print(min(c["x0"] for c in body), max(c["x1"] for c in body))
```

## 2. Sizes are what the reader sees

- **Everything is 9 pt printed** — axis labels, legends, colorbar labels, annotations.
- **Tick labels may be smaller**: 8 pt printed.
- **Panel letters: Arial, 12 pt printed**, upper case, no LaTeX.

"Printed" is the operative word. A font size in the code is only meaningful
together with the scale factor LaTeX applies on inclusion.

## 3. Control the scale factor explicitly

`\includegraphics[width=\linewidth]` scales the PDF by `target_width / pdf_width`.
Two acceptable ways to get 9 pt on the page:

**(a) Build at the final width** — scale factor 1, sizes in the code are printed sizes.

```python
W = 481.9                      # or 231.95 for a single-column figure
fig = plt.figure(figsize=(W / 72, H / 72))
RC = {"font.size": 9, "xtick.labelsize": 8, ...}
```

**(b) Keep an existing design width and compensate** — use when a layout already
exists and only the type should change.

```python
FIG_W_IN = 11.0                        # the figure's own width
TEXTWIDTH_PT = 481.9                   # what it is included at
SCALE = FIG_W_IN * 72 / TEXTWIDTH_PT   # 1.6435
RC = {"font.size": 9 * SCALE, "xtick.labelsize": 8 * SCALE, ...}
```

Option (b) also preserves panel aspect ratios, which (a) does not: enlarging the
type at fixed figure width forces `constrained_layout` to take margin from the
axes, reshaping the panels.

## 4. Never `bbox_inches="tight"` on a print-sized figure

It re-crops the page to the drawn content, so the PDF is no longer the width you
designed for, and LaTeX silently rescales it — undoing the font sizing.
`constrained_layout` already trims margins. Save with plain
`fig.savefig(path, format="pdf")`.

## 5. Match the float environment to the figure width

This is the easiest way to get it badly wrong, and it is invisible in the notebook:

| PDF width | environment | scale | 9 pt prints as |
|---|---|---|---|
| 481.9 pt | `figure*`, `width=\linewidth` | 1.00 | 9.0 pt |
| 231.95 pt | `figure`, `width=\linewidth` | 1.00 | 9.0 pt |
| 231.95 pt | `figure*`, `width=0.95\linewidth` | 1.97 | **17.8 pt** |

A single-column figure placed in a `figure*` is scaled ~2× and its type is
twice the caption size.

## 6. Fonts

```python
plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "text.latex.preamble": r"\usepackage{amsmath}"})
```

- `usetex` gives Computer Modern, matching the paper.
- **`amsmath` is required** for `\dfrac`; matplotlib's preamble does not load it
  and plain LaTeX has no `\dfrac`. Pre-flight the renderer with a string that
  actually uses the macros, so a missing package fails loudly rather than
  silently falling back.
- A fallback to mathtext should print a warning: the numbers stay right, the
  typography does not.
- Panel letters are Arial, so they must bypass LaTeX: `usetex=False`,
  `fontfamily="Arial"`.
- Consider `"pdf.fonttype": 42` so non-LaTeX fonts (Arial) embed as TrueType
  rather than Type 3, which publisher PDF checks often flag.

## 7. Panel letters

- Upper case, Arial, 12 pt printed, flush left, each on its own line above its
  panel where the layout allows.
- Place them from known coordinates, not from a renderer-dependent bounding box
  where possible.
- `constrained_layout` sizes axes to fill the figure, so there may be no room
  above the top row. Reserve it with an empty first gridspec row and make the
  figure taller by that band — then compensate the panel heights, because the
  band row also claims padding.

## 8. matplotlib version traps

The notebooks must run on whatever kernel the authors use (3.5.2 has been seen):

- `ax.get_tightbbox()` requires a `renderer` argument before 3.8.
- `fig.canvas.get_renderer()` does not exist on a plain `FigureCanvasBase`; fall
  back to `FigureCanvasAgg(fig).get_renderer()`.
- `fig.get_layout_engine()` is 3.6+.
- `constrained_layout` sizes colorbar columns differently across versions, so a
  layout that depends on it is not reproducible between machines. Prefer
  explicit `add_axes` rectangles when the exact geometry matters.

## 9. Verify the output, don't trust the settings

After building, check the page size and the embedded sizes:

```python
import pdfplumber, collections
p = pdfplumber.open("fig2/fig2.pdf").pages[0]
scale = 481.9 / p.width                       # the factor LaTeX will apply
for s, n in sorted(collections.Counter(round(c["size"], 2) for c in p.chars).items()):
    print(f"{s:6.2f} pt in file -> {s*scale:5.2f} pt printed  (n={n})")
```

`pdffonts figure.pdf` shows the embedded font types (Type 1 for Computer Modern,
and ideally TrueType rather than Type 3 for Arial).

## 10. Data precision, where fits are involved

`%.17g` on write and `float_precision="round_trip"` on read. The shape fits
behind Fig. 2 are ill-conditioned at the small end of their λ grid: a one-ULP
difference in the input selects a different smoothing weight. pandas' default
CSV float parser is not correctly rounded.

## 11. Compositing hand-drawn artwork, as vector

matplotlib cannot read SVG, so Illustrator artwork is converted to PDF and
stamped into the saved figure with `pypdf` (see the last cell of `fig4.ipynb`).
Three things make this reliable:

**Never assume the converted page box.** `cairosvg` writes the PDF at 0.75x the
viewBox (96 -> 72 dpi) and ignores `output_width` for PDF output; other
converters differ again. Read `page.mediabox` back and derive the factor:

```python
f = float(art.mediabox.width) / VB_W        # viewBox units -> that page's points
```

**Size on a landmark both files share, not on the bounding box.** A reference
assembly exported from Illustrator has its *text converted to outlines*, so a
vector bbox there mixes type with drawing and is useless. The artwork is
usually the only coloured content, so threshold on saturation in both and
compare those boxes — if the two aspect ratios agree, it is the same drawing:

```python
import numpy as np, io, cairosvg
from PIL import Image
png = cairosvg.svg2png(url=svg, output_width=int(VB_W*8), output_height=int(VB_H*8))
im = Image.open(io.BytesIO(png)).convert("RGBA")
rgba = np.asarray(im).astype(int)
comp = (rgba[..., :3]*rgba[..., 3:] + 255*(255-rgba[..., 3:]))//255   # over white
ys, xs = np.where((comp.max(2) - comp.min(2)) > 40)
print(xs.min()/8, xs.max()/8, ys.min()/8, ys.max()/8)                 # viewBox units
```

An alpha-only bbox is no good here: a gradient fading to transparent reaches far
past the visible edge, which inflates the box and changes the aspect ratio.

**Anchor the labels to the drawing.** Put the type in matplotlib, not in the
SVG, and place it through one mapping from viewBox coordinates to axes
coordinates (`cart(u, v)` in `fig4.ipynb`). Positions read off the drawing then
survive any later change to its size or position; fractions of the panel do not.
