# MSM: How Small Molecules Bind Disordered Proteins

Code for **building Markov State Models (MSMs) of small molecules moving over intrinsically disordered proteins (IDPs)**.

Disordered proteins have no fixed binding pocket, so docking does not apply. A ligand slides and hops across the protein surface instead. This code takes an MD simulation of a protein + ligand and turns it into a **map of binding states and the paths between them**. It shows where the ligand prefers to sit, how it moves between sites, and how tightly it binds overall.

  <p align="center">
    <img src="https://github.com/louetadelie1/MSM/raw/HEAD/images/simu_time.png" width="90%">
  </p>


> **Paper:** Adélie Louet, Gerhard Hummer, Michele Vendruscolo.
> *Binding without Docking: Small-Molecule Binding Paths on Disordered Proteins.*
> J. Am. Chem. Soc. 2026. https://doi.org/10.1021/jacs.6c05891

---

## Pipeline

```
 MD trajectory                Binding states                 Network of states
 (protein + ligand)           (which residues touch          (MSM + communities)
                               the ligand?)
                                                              ●──●
  frame 1  ──►  residues {12, 15, 16, 19, 20}  = state A      │╲ │   ← ligand "hops"
  frame 2  ──►  residues {12, 15, 16, 19, 20}  = state A      ●──●     between groups
  frame 3  ──►  residues {15, 16, 19, 20, 23}  = state B       ╲      of related sites
  frame 4  ──►  no contacts                    = unbound        ●──●──●
     ...                                                        (communities = binding regions)
```

1. **Define a state for every frame.** In each frame, find the *n* residues closest to the ligand (default *n* = 5). That set of residues, called an "n-uplet", is the binding state. Frames with too few contacts count as **unbound**.
2. **Count transitions.** Count how often the ligand moves from one state to the next between frames. This gives the **transition matrix** (the MSM).
3. **Equilibrium populations.** The leading eigenvector of the matrix gives the probability of finding the ligand in each state.
4. **Group states into binding regions.** Louvain community detection on the state network merges microstates into **macrostates**: regions of the protein the ligand explores together.
5. **Binding affinity.** The bound/unbound populations and binding/unbinding rates give an estimate of **K<sub>d</sub>**.

---

## Repository layout

```
MSM/
├── community_identification/   ← the core MSM pipeline (start here)
├── chapman_kolmogorov_test/    ← checks that the MSM is valid
├── VAMPNET/                    ← deep-learning MSM used for comparison
├── matrix_handoff/             ← residue-to-residue "handoff" analysis
└── visualization/              ← export states for viewing in VMD
```

### `community_identification/`: the core pipeline

| File | What it does |
|---|---|
| **`run_auto.py`** | **Main script.** Builds the MSM for each system, finds communities, and saves the network figures and results. |
| `run_auto_with_kd_ts.py` | Extended version that also tracks K<sub>d</sub> and timescales across settings. |
| `population_equilibrium.py` | Turns distances into states, builds the transition matrix, and computes equilibrium populations. Optionally reweights frames from biased (metadynamics) simulations. |
| `clustering_uplets.py` | Builds the microstate network and runs Louvain community detection to get macrostates. |
| `kd_calculation.py` | Estimates K<sub>d</sub> two ways: from bound/unbound populations and from on/off rates. |
| `distance_cal.py` | Pre-processing: computes ligand–residue distances for every frame from the trajectory. |
| `additional_analysis.py` | Extra analyses: linear binding motifs, and **gliding vs hopping** between states (Jaccard similarity). |

### Validation and comparison

| File | What it does |
|---|---|
| `chapman_kolmogorov_test/enumeration_ck_test.py` | **Chapman–Kolmogorov test.** Compares MSM-predicted transition probabilities with those observed in the trajectory over increasing lag times. |
| `VAMPNET/vampnet_2_iterations.py` | Builds an alternative MSM with **VAMPNet** (deeptime), using two rounds of clustering. |
| `VAMPNET/ami_score.py` | Compares VAMPNet clusters with the community clusters (adjusted mutual information). |

### Other analyses

| File | What it does |
|---|---|
| `matrix_handoff/matrix_handoff.py` | Tracks how the ligand is "handed off" from one residue to another over time and plots it as a heatmap. |
| `visualization/visualizing_clusters.py` | Writes the residues of each cluster as Tcl scripts, so binding states can be highlighted in **VMD**. |
| `visualization/visualizing_clusters_vmd.ipynb` | Notebook version of the VMD visualisation. |

---

## Getting started

### 1. Install dependencies

```bash
conda create -n msm -c conda-forge python=3.11 \
    numpy scipy pandas matplotlib seaborn scikit-learn networkx \
    python-louvain mdanalysis mdtraj tqdm
pip install smplotlib

# only needed for VAMPNet
pip install torch deeptime
```

### 2. Compute distances from your trajectory

You need a topology (`.pdb`) and a trajectory (`.xtc`) of the protein + one ligand:

```bash
python community_identification/distance_cal.py <sim_dir> <protein_name> <ligand_name> <trajectory.xtc>
```

This saves two arrays of shape `(n_frames, n_residues)`:

- `d_24_t_com_avg.pkl`: average ligand–residue distance
- `d_24_t_closest.pkl`: closest ligand–residue distance

### 3. Build the MSM

Open `community_identification/run_auto.py` and add a block for your system with:

| Setting | Meaning |
|---|---|
| `pdb` | Topology file |
| `distances_com`, `distances_closest` | The two pickles from step 2 |
| `w_com`, `w_closest` | How much to weight each distance type (should sum to 1) |
| `combined_threshold` | Contact cutoff in nm, or `True` to derive it from the weights |
| `uplet_type` | Number of residues that define a state (default `5`) |
| `w_eq` | Per-frame weights (`1.0` for unbiased simulations) |

Then run it **from the repository root** (the scripts import `community_identification` as a package):

```bash
python community_identification/run_auto.py
```

### Outputs

- Network figures of microstates and macrostates (PDF/PNG)
- Pickled results: equilibrium populations (`p_eq_keys.pckl`), transitions (`dictionary_transitions.pckl`), community assignments (`inv_map.pckl`)
- Printed K<sub>d</sub> estimates and the fraction of time the ligand is bound

---

## Note

File paths are currently hard-coded to the original analysis machine. Update them to your own data before running.

## Citation

If you use this code, please cite:

```bibtex
@article{louet2026binding,
  title   = {Binding without Docking: Small-Molecule Binding Paths on Disordered Proteins},
  author  = {Louet, Ad{\'e}lie and Hummer, Gerhard and Vendruscolo, Michele},
  journal = {Journal of the American Chemical Society},
  year    = {2026},
  doi     = {10.1021/jacs.6c05891}
}
```
