### Table 1. Sparse-autoencoder quality.

| Model | Layer | Dict. | L0 | Dead % | #F1>0.5 | #F1>0.3 | Families | Median r (n) |
|---|---|---|---|---|---|---|---|---|
| ESM-2-8M | 4 | 10,240 | 64 | 0.01 | 228 | 680 | 45 | 0.71 (192) |
| SeqDance | 10 | 15,360 | 64 | 0.02 | 44 | 611 | 35 | 0.55 (206) |

### Table 2. Best-feature test F1 by disorder axis.

| Axis | ESM-2 | SeqDance |
|---|---|---|
| AlphaFold disorder | 0.83 | 0.60 |
| Curated disorder | 0.43 | 0.28 |
| Phase separation | 0.43 | 0.37 |
| Proline-rich | 0.55 | 0.51 |
| Extended coil | 0.40 | 0.31 |
| Fuzzy (D-to-D) | 0.32 | 0.18 |
| Ordered (control) | 0.18 | 0.10 |
| Compact (ν proxy) | 0.18 | 0.24 |
| MoRF (MobiDB) | 0.11 | 0.09 |
| MoRF (experimental) | 0.18 | 0.15 |
| Molten globule | 0.09 | 0.30 |
| Pre-molten globule | 0.12 | 0.33 |

### Table 3. Representative validated ESM-2 features spanning compositional families.

| Feature | Identity (from validated description) | Axis | F1 | Held-out r |
|---|---|---|---|---|
| 8306 | Acidic low-complexity (poly-D/E) | AlphaFold disorder | 0.83 | 0.72 |
| 2474 | Polar / charged homopolymer repeats | AlphaFold disorder | 0.75 | 0.69 |
| 4236 | Arginine-rich / NLS-like | AlphaFold disorder | 0.47 | 0.91 |
| 7063 | Proline-rich (polyproline-II) | AlphaFold disorder | 0.74 | 0.89 |
| 3350 | Glycine-rich linker | AlphaFold disorder | 0.72 | 0.91 |
| 3401 | Ser/Thr-rich (phospho-acceptor) | AlphaFold disorder | 0.71 | 0.67 |
| 8708 | Histidine-rich | AlphaFold disorder | 0.74 | 0.50 |
