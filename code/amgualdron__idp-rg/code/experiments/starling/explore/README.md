# Exploration behind paper fig26

Not part of the release. Run from the repo root with `PYTHONPATH=.:paper`; each needs `paper/cache/collected.pkl`.

- `nu_check.py`: why nu disagrees while W1 is small (normalization, fit quality, A convention, W1 per difference in nu). Writes `nu_check.csv`.
- `worst_nu.py DIR`: histograms of the proteins with the worst nu, `worst_nu.png`. Reads `DIR/nu_check.csv`.
- `variants.py DIR`: Bowley skewness and tail probabilities as replacements for the nu panel: `variants_table.csv`, `variants_values.csv`, `variant_A_bowley.png`, `variant_B_tails.png`. Variant B became fig26.
