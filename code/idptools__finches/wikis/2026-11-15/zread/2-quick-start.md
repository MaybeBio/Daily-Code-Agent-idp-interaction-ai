

Get from zero to your first interaction prediction in under five minutes. This guide walks you through installing FINCHES, initializing a frontend object, and computing your first **epsilon** interaction value — the core metric that quantifies the chemical specificity between two intrinsically disordered regions (IDRs).

## Prerequisites

FINCHES requires **Python 3.9+** (3.11 or 3.12 recommended for performance) and depends on `numpy`, `scipy`, `cython`, `pytorch`, `mdtraj`, and `metapredict`. The table below summarizes the supported installation paths:

| Method | Best For | Command |
|--------|----------|---------|
| pip from GitHub | Most users | `pip install git+https://git@github.com/idptools/finches.git` |
| Local editable install | Contributors / developers | `git clone git@github.com:idptools/finches.git` → `cd finches` → `pip install -e .` |
| Google Colab | Zero-install exploration | [finches-colab notebooks](https://github.com/idptools/finches-colab) |
| Web server | Quick lookups, no code | [finches-online.com](https://www.finches-online.com/) |

Sources: [README.md](/README.md#L75-L155), [pyproject.toml](/pyproject.toml#L13-L31)

## Installation

### Step 1 — Create and activate a conda environment

```bash
conda create -n finches python=3.12 -y
conda activate finches
```

### Step 2 — Install core dependencies via conda

This ensures consistent NumPy/PyTorch builds (critical on macOS where mixing conda and PyPI builds causes failures):

```bash
conda install numpy pytorch scipy cython matplotlib jupyter -c pytorch
conda install mdtraj
pip install metapredict
```

### Step 3 — Install FINCHES

```bash
pip install git+https://git@github.com/idptools/finches.git
```

### Step 4 — Verify the installation

```bash
python -c "from finches.frontend.mpipi_frontend import Mpipi_frontend; mf = Mpipi_frontend(); print('Success')"
```

> **Note:** The first import takes 60–90 seconds while Cython-compiled modules are loaded. Subsequent imports are fast.

Sources: [README.md](/README.md#L83-L158), [docs/getting_started.rst](/docs/getting_started.rst#L22-L40)

## Your First Calculation

The diagram below shows the minimal workflow from sequence to prediction. FINCHES follows a **frontend pattern**: you instantiate a forcefield-specific frontend object, then call its methods with amino acid sequences.

```mermaid
flowchart LR
    A["Import frontend"] --> B["Instantiate object<br/>(choose forcefield)"]
    B --> C["Define sequence(s)<br/>as Python strings"]
    C --> D["Call .epsilon()"]
    D --> E["Get interaction value<br/>(float)"]
    style A fill:#e8f5e9,stroke:#2e7d32
    style B fill:#e3f2fd,stroke:#1565c0
    style C fill:#fff3e0,stroke:#ef6c00
    style D fill:#f3e5f5,stroke:#7b1fa2
    style E fill:#fce4ec,stroke:#c62828
```

### Initialize a frontend

FINCHES exposes two frontend classes — **`Mpipi_frontend`** and **`CALVADOS_frontend`** — that provide identical interfaces backed by different forcefield parameterizations:

```python
from finches import Mpipi_frontend, CALVADOS_frontend

# Mpipi-GG forcefield (supports protein AND RNA via 'U' residues)
mf = Mpipi_frontend(salt=0.150, dielectric=80.0)

# CALVADOS2 forcefield (protein only; pH- and temperature-sensitive)
cf = CALVADOS_frontend(salt=0.150, pH=7.4, temp=288)
```

Both constructors default to roughly physiological conditions. The key difference: `Mpipi_frontend` accepts `salt` and `dielectric`, while `CALVADOS_frontend` accepts `salt`, `pH`, and `temp` (in Kelvin).

Sources: [finches/__init__.py](/finches/__init__.py#L16-L17), [finches/frontend/mpipi_frontend.py](/finches/frontend/mpipi_frontend.py#L12-L21), [finches/frontend/calvados_frontend.py](/finches/frontend/calvados_frontend.py#L34-L43)

### Compute an epsilon value

The **epsilon** value is a scalar proxy for the interaction free energy between two IDR sequences — more negative means more attractive, more positive means more repulsive, and zero corresponds to a neutral GS-linker reference:

```python
# DDX4 N-terminal domain (a well-studied phase-separating IDR)
ddx4_ntd = 'MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSRFEDGDSSGFWRESSNDCEDNPTRNRGFSKRGGYRDGNNSEASGPYRRGGRGSFRGCRGGFGLGSPNNDLDPDECMQRTGGLFGSRRPVLSGTGNGDTSQSRSGSGSERGGYKGLNEEVITGSGKNSWKSEAEGGES'

# Homotypic interaction (DDX4 with itself)
print(mf.epsilon(ddx4_ntd, ddx4_ntd))
```

For **heterotypic** interactions, simply pass two different sequences:

```python
seq_A = 'GGGGGRRRRR'   # a short arginine-rich motif
seq_B = 'GGGGGYYYYY'   # a short tyrosine-rich motif

print(mf.epsilon(seq_A, seq_B))  # attractive if π–cation contacts dominate
```

> [!TIP]
> The epsilon value is **order-dependent** (extrinsic): `epsilon(A, B) ≠ epsilon(B, A)` in general. Think of the first sequence as "bathed in" the second. To recover an intrinsic (order-independent) value, divide by the length of the first sequence: `epsilon(A, B) / len(A) == epsilon(B, A) / len(B)`.

Sources: [finches/frontend/frontend_base.py](/finches/frontend/frontend_base.py#L206-L241), [README.md](/README.md#L36-L43), [demo/overview_uses/basic_uses.ipynb](/demo/overview_uses/basic_uses.ipynb#L140-L175)

### Generate an interaction map

Beyond a single scalar, the **interaction figure** produces a 2D heatmap showing position-resolved interaction strength across both sequences — ideal for identifying which subregions drive binding:

```python
mf.interaction_figure(ddx4_ntd, ddx4_ntd)
```

This computes a sliding-window interaction matrix and renders it with parallel disorder prediction tracks. The returned tuple contains the matplotlib figure and axes objects for further customization.

Sources: [finches/frontend/frontend_base.py](/finches/frontend/frontend_base.py#L288-L416), [finches/frontend/mpipi_frontend.py](/finches/frontend/mpipi_frontend.py#L128-L149)

## Core Frontend API at a Glance

Both `Mpipi_frontend` and `CALVADOS_frontend` expose the same set of high-level methods. The table below summarizes the primary entry points you'll use most often:

| Method | Returns | Purpose |
|--------|---------|---------|
| `epsilon(seq1, seq2)` | `float` | Scalar interaction energy between two full sequences |
| `epsilon_vectors(seq1, seq2)` | tuple of arrays | Per-residue attractive and repulsive interaction vectors |
| `intermolecular_idr_matrix(seq1, seq2)` | `(matrix_tuple, disorder_1, disorder_2)` | Raw sliding-window interaction matrix + disorder profiles |
| `interaction_figure(seq1, seq2)` | `(fig, im, ax_main, ...)` | Rendered 2D interaction heatmap with disorder tracks |
| `per_residue_attractive_vector(seq1, seq2)` | `(indices, values)` | Per-residue average of attractive interactions (identifies "stickers") |
| `per_residue_repulsive_vector(seq1, seq2)` | `(indices, values)` | Per-residue average of repulsive interactions (identifies "spacers") |

All methods accept `use_aliphatic_weighting` and `use_charge_weighting` flags (default `True`) that apply correction terms for hydrophobic and charge-context effects.

Sources: [finches/frontend/frontend_base.py](/finches/frontend/frontend_base.py#L206-L282), [finches/frontend/frontend_base.py](/finches/frontend/frontend_base.py#L570-L696)

## Choosing a Forcefield

| Feature | `Mpipi_frontend` | `CALVADOS_frontend` |
|---------|-------------------|----------------------|
| Underlying model | Mpipi-GG v1 | CALVADOS2 |
| Constructor params | `salt`, `dielectric` | `salt`, `pH`, `temp` |
| RNA support (`U` residue) | ✅ Yes | ❌ No (raises `ValueError`) |
| Best for | Protein–protein and protein–RNA | Protein–protein with pH/temperature sensitivity |
| Primary citation | Joseph et al. (2021) + Lotthammer et al. (2024) | Tesei et al. (2022) |

For most protein-only use cases, either forcefield works well. Use `Mpipi_frontend` when working with RNA-containing sequences (use `U` for uracil); use `CALVADOS_frontend` when you need pH- or temperature-dependent predictions.

Sources: [finches/frontend/mpipi_frontend.py](/finches/frontend/mpipi_frontend.py#L12-L21), [finches/frontend/calvados_frontend.py](/finches/frontend/calvados_frontend.py#L16-L25), [finches/frontend/calvados_frontend.py](/finches/frontend/calvados_frontend.py#L28-L43)

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `ModuleNotFoundError: No module named 'finches'` | FINCHES not installed | Run `pip install git+https://git@github.com/idptools/finches.git` |
| First import takes >60s | Cython compilation on first load | Normal — subsequent imports are fast |
| `ValueError: CALVADOS2 cannot handle RNA ('U')` | `U` residue in CALVADOS input | Switch to `Mpipi_frontend` which supports RNA |
| `ImportError` for `metapredict` | Missing pip-only dependency | Run `pip install metapredict` |
| macOS NumPy/PyTorch crash | Mixed conda + PyPI builds | Install both via conda: `conda install numpy pytorch -c pytorch` |

Sources: [README.md](/README.md#L94-L109), [finches/frontend/calvados_frontend.py](/finches/frontend/calvados_frontend.py#L16-L25)

## Where to Go Next

Now that you can compute interaction energies and generate interaction maps, explore FINCHES' deeper capabilities:

1. **[Installation Guide](3-installation-guide)** — Advanced setup: editable installs, Cython compilation details, `uv` package manager
2. **[Mpipi and CALVADOS Frontends](5-mpipi-and-calvados-frontends)** — Full parameter reference for every frontend method
3. **[Interaction Maps and Figures](6-interaction-maps-and-figures)** — Customizing heatmap appearance, domain overlays, and null-shuffle normalization
4. **[Epsilon Calculation and Weighting](9-epsilon-calculation-and-weighting)** — Understanding aliphatic and charge weighting corrections
5. **[Flory-Huggins Phase Diagrams](14-flory-huggins-phase-diagrams)** — Going from epsilon values to phase behavior predictions
