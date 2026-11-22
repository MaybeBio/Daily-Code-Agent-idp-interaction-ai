# Figure data

These are the tables that the figure scripts in `runnable/code/` read: 77 CSV and TSV files, about 176 MB. They are kept in the project's folder layout (`BENDER-MT/...`, `BENDER-BIO/...`), so the runnable scripts find them once `CS_DATA` points here:

```bash
export CS_DATA="$(pwd)/figure_data"
source runnable/env.sh
python runnable/code/fig_apps.py
```

`MANIFEST.csv` gives the size, SHA-256 and figure scripts of each file.

The list comes from running every figure script with each file it opened recorded. No figure was saved and nothing was written, so it holds exactly the tables the figures use.

## What this folder does not contain

- **Per-frame caches and image assets.** Six scripts also read binary per-frame caches (`.npz`) or rendered glyphs (`.png`), and these are not included. So the CSV and TSV files here are enough to redraw only `fig_extra.py`, `fig_apps.py`, `fig_supplementary.py`, `fig_supp_redraw.py`, `fig_acetyl.py`, `fig_lineage.py` and `fig_phosphosite.py`.

  | Script | Not included |
  |---|---|
  | `fig_heldout.py` | per-frame caches of the held-out runs |
  | `fig_manuscript2.py` | per-frame caches |
  | `fig_forcing.py` | per-frame caches of the calibration runs |
  | `fig_page4.py` | per-frame caches of the PAGE4 runs |
  | `fig_overview.py` | an example trajectory cache and glyph images |
  | `fig_ogt_lineup.py` | glyph images |

- **Trajectories.** The 300 K trajectories are on Hugging Face (see the main README).
