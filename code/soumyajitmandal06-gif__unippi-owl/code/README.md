# UniPPI-OWL v1.8

**Universal Structure-to-Affinity, Interface Auditing, PTM/Metal Gating, and Mutational Screening Resource**

* **Author:** Soumyajit Mandal (Department of Biochemistry, Central University of Rajasthan, India)
* **Model Version:** `UniPPI-v1.8` (Frozen Decision Threshold: `UIS = 0.5000`; High-Confidence Cutoff: `UIS = 0.6500`)
* **Architecture:** Dual-Layer — Zero-Dependency Python CLI (`unippi_cli.py`) + Client-Side Browser Application (`index.html`)
* **Workflows & Packaging:** Nextflow DSL2 (`main.nf`), Snakemake (`Snakefile`), Zero-Install Local Launchers (`start-local.sh` / `start-local.cmd`), and `pip` (`pyproject.toml`)
* **License:** Open-Source BSD 3-Clause

---

## 1. Overview and Scientific Objective

### What This Tool Does (In Plain Words)
Structure prediction models like AlphaFold3 (AF3), AlphaFold2, ESM3, and Boltz-1 predict the 3D shape of protein complexes (`.cif` or `.pdb` files) and output geometric confidence scores (`ipTM`, `pLDDT`, and `PAE` matrices). However, geometric confidence is not the same as physical binding strength. AlphaFold3 does not calculate thermodynamic binding free energy (`dG` in kcal/mol) or dissociation constants (`Kd` in nM), nor does it tell you how a point mutation (`ddG`), a kinase phosphorylation event (`SEP/TPO/PTR`), a bridging metal ion (`Ca2+/Zn2+/Co2+`), or a change in salt concentration (`[NaCl]`) and `pH` will affect binding. In addition, running external energy programs like FoldX or Rosetta requires manually converting AlphaFold3's `.cif` files into strict legacy 80-column `.pdb` files.

**UniPPI-OWL v1.8** solves these bottlenecks in a single workflow:
1. **Automated Format Sanitization:** Converts multi-chain AlphaFold3 `.cif` and `.pdb` structures into strict 80-column FoldX/Rosetta-ready `.pdb` files (stripping crystallographic water and unphysical `GLY CB` atoms), extracts non-polymer small-molecule ligands into V2000 `.sdf` files, and outputs FoldX `individual_list.txt` mutation files.
2. **Eight-Channel Structure-to-Affinity Estimation:** Calculates confidence-weighted atomic contacts across 8 physicochemical channels (including post-translational phosphorylation and divalent metal coordination), solvent-accessible surface burial (`wBSA`), binding free energy (`dG_uni`), and dissociation constants (`Kd`) in less than 0.1 seconds per complex.
3. **Disordered-Loop (IDR) Artifact Filtering:** Prevents low-confidence intrinsically disordered tails (common in plant and eukaryotic transcription factors) from generating false-positive binding scores, while rescuing true local interfaces whose global `ipTM` is diluted by long non-interacting tails.
4. **In-Silico Mutagenesis, PTM Gating, and Environmental Stress Modeling:** Evaluates single-site 20-amino-acid saturation mutagenesis (`ddG`), double-mutant epistatic cycles (`ddG_int`), CPK21/PP2A phosphorylation state switching (`--ptm-state`), CHUMA cross-paralog sequence conservation projection, and ionic strength (`10-1000 mM NaCl`) and `pH` (`4.0-9.5`) titration curves.

### Dual-Layer Architecture
To serve both computational bioinformaticians and wet-lab biologists, UniPPI-OWL is built in two synchronized layers that share identical mathematical equations and produce identical numerical outputs:
* **Layer 1 — Automated Python CLI (`unippi_cli.py`):** Written entirely in the Python 3 standard library (`math.dist` paired with 3D spatial voxel hashing) with zero external package dependencies (`NumPy`, `Biopython`, or C++ compilers are not required). Designed for batch folder screening and workflow managers such as Nextflow (`main.nf`) or Snakemake (`Snakefile`).
* **Layer 2 — Client-Side Browser Application (`index.html`):** A zero-server HTML5, Vanilla JavaScript, and WebGL (`3Dmol.js`) Single-Page Application (SPA) following the `-OWL` (*On the Web Lite*) standard. All calculations run locally in browser memory so unpublished structural coordinates never leave the user's computer.

---

## 2. Step-by-Step Methodology and Biophysical Logic

### Step 1: Universal Input Parsing and 3-Mode Confidence Auto-Detection
Our parser inspects file headers (`_exptl.method`, `_ma_qa_metric_local`, `EXPDTA`), B-factor distributions, and companion `.json` files to run automatically in one of three modes:
* **`MODE_A_PRED_FULL` (Predicted Structure + Companion JSON):** Reads 3D coordinates alongside AF3/AF2/ColabFold JSON files and maps `token_chain_ids` and `token_res_ids` back to exact residue keys `(chain_id, res_seq)` so multi-atom modified residues (`SEP`, `TPO`, `PTR`) and ligands never shift the PAE matrix out of alignment.
* **`MODE_B_PRED_COORD` (Predicted Structure Only, No JSON):** Reads per-residue `pLDDT` directly from the coordinate B-factor column (automatically scaling fractional `[0, 1]` values by `100`).
* **`MODE_C_EXP_CRYSTAL` (Experimental Structure):** Applies inverse thermal normalization `c_i = exp(-max(0, B_i - B_min) / 60.0)` and normalizes `pDockQ_coord` against its theoretical `0.742` upper ceiling so experimental crystal structures are scored accurately without a JSON file.

### Step 2: Two-Stage 3D Voxel Spatial Hashing
We place residue centroids (`CA` for amino acids, `C4'` for DNA/RNA, geometric center for ligands/ions) on Chain B into a 3D cubic lattice with voxel size `L_voxel = 14.0 A`:

    Voxel_Key(x, y, z) = ( floor(x / 14.0), floor(y / 14.0), floor(z / 14.0) )

For each residue on Chain A, heavy-atom distances are evaluated strictly against Chain B residues in the `3 x 3 x 3 = 27` neighboring voxels whose centroids are within `16.0 A`.

### Step 3: Flat-Top Plateau Confidence Weighting (`w_ij`)
Each residue `i` receives a plateau confidence weight `c_i` in `[0.02, 1.0]`:

    c_i = clip( (pLDDT_i - 30.0) / 55.0,  0.02,  1.00 )     [in Modes A and B]
    c_i = exp( -max(0.0, B_i - B_min) / 60.0 )              [in Mode C]

For an inter-chain residue pair `(i, j)` with mean cross-chain aligned error `PAE_ij = 0.5 * (PAE[i,j] + PAE[j,i])`, the joint confidence weight `w_ij` is:

    w_ij = c_i * c_j * exp( -max(0.0, PAE_ij - 3.5) / lambda_PAE )

where `lambda_PAE = 8.0 A`.

### Step 4: Eight-Channel Physicochemical Contact, PTM, and Metal Classification
For every inter-chain residue pair with minimum heavy-atom distance `d_min <= 5.5 A`, atomic interactions are classified into eight physical channels:
1. **Salt Bridges & Phosphate Clamps (`d <= 4.0 A` / `4.5 A`):** Cationic side-chain nitrogens (`Arg, Lys, His`) to anionic carboxylate oxygens (`Asp, Glu`), nucleic acid backbone phosphates (`OP1,OP2`), or phosphorylated side chains (`SEP, TPO, PTR` carrying `-PO3^2-` oxygens).
2. **Hydrogen Bonds (`d <= 3.5 A`):** Polar `N, O, S` donor-acceptor pairs not already counted as salt bridges (capped at `2` per residue pair).
3. **Hydrophobic Packing (`d <= 4.5 A`):** Side-chain apolar `C, S` atoms between hydrophobic residues (`Ala, Val, Leu, Ile, Met, Phe, Trp, Pro, Tyr` or ligands), weighted by `phi_HP(d) = clip((4.75 - d) / 0.80, 0.35, 1.0)` (capped at `6.0` per residue pair).
4. **Aromatic Pi-Stacking (`d_ring <= 5.5 A`):** Distance between aromatic ring centroids (`Phe, Tyr, Trp, His`).
5. **Disulfide Bridges (`d <= 2.3 A`):** Covalent inter-chain `Cys SG` to `Cys SG` bonds.
6. **Steric Clashes (`d < 2.2 A`):** Unphysical non-bonded heavy-atom overlap (capped at `3` per residue pair).
7. **Intra-Helical Proline Kink Audit (`N_Pro_helix`):** Counts native Proline (`PRO`) residues located inside folded alpha-helical segments (`4.8 <= d_CA(i, i+4) <= 7.2 A` and `pLDDT >= 50`) within `+/- 2` residues of the interface.
8. **Divalent Metal-Ion Coordination (`d <= 3.8 A`):** Detects interfacial divalent cations (`ZN, MG, CA, MN, FE, CU, CO, NI`) bridging coordinating side chains across subunit interfaces.

### Step 5: Golden-Spiral Shrake-Rupley SASA and Pairwise Multi-Chain Burial (`wBSA`)
We compute the Solvent-Accessible Surface Area (`SASA`) of the isolated chains (`SASA_A`, `SASA_B`) and bound binary pair (`SASA_AB`) using a 64-point Fibonacci Golden-Section Spiral around each atom with water probe radius `r_probe = 1.40 A`:

    BSA  = max( 0.0,  SASA_A + SASA_B - SASA_AB )
    wBSA = mean(w_ij) * BSA

### Step 6: Hybrid Contact-Gated Free Energy (`dG_uni`) and Salt/pH/PTM Modulation
We combine the PRODIGY linear equation (`dG_PRODIGY`) with an explicit biophysical solvation and contact potential (`dG_phys`):

    # Biophysical Solvation & Contact Potential:
    dG_phys = +2.10
              - salting_out * (0.110 * w_HP + 0.215 * IC_aa + 0.0026 * wBSA + 0.0035 * max(0, wBSA - 600))
              - (0.765 * w_SB * debye_kappa * ph_kappa + 0.460 * w_HB + 0.615 * w_Pi)
              + 0.160 * IC_pp + 1.50 * w_clash + 2.20 * N_Pro_helix

    gate   = 0.0 if IC_total < 3.0 else 1.0 / (1.0 + exp(-0.40 * (IC_total - 11.0)))
    w_phys = clip(0.55 + 0.35 * f_apolar + (0.15 if N_Pro_helix > 0 else 0.0) + (0.35 if wBSA >= 800 else 0.0), 0.55, 0.95)
    dG_uni = clip( gate * ((1 - w_phys)*dG_PRODIGY + w_phys*dG_phys) + (1 - gate)*dG_phys, -22.0, +5.0 )
    Kd     = exp( dG_uni / (R * T) )

### Step 7: Universal Interaction Score (`UIS`) and Mechanistic Patterns (`Patterns A-D`)
`UIS` in `[0, 1]` combines global chain topology (`ipTM`), local IDR-immune PAE coherence (`LIS`), size-normalized coordinate quality (`pDockQ_PAE` or `pDockQ_coord / 0.742`), and scaled affinity `sigma_aff = clip((-dG_uni - 5.0) / 8.5, 0.0, 1.0)`:

    UIS = 0.25 * ipTM + 0.30 * LIS + 0.25 * pDockQ_PAE + 0.20 * sigma_aff               [Mode A]
    UIS = 0.55 * clip(pDockQ_coord / 0.742, 0.0, 1.0) + 0.45 * sigma_aff                 [Modes B & C]

---

## 3. Six-Complex Validation Benchmark and Real RCSB PDB (`1BRS`) Empirical Verification

| Rank | Complex ID | UIS | Margin | dG_uni (kcal/mol) | Kd (nM) | ipTM | LIS | wBSA (A^2) | SB | HB | HP | Pro Kinks | Pattern | Class |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **#1** | `1_VIP1_Homodimer_Positive` | **0.777** | `+0.277` | **-10.35** | **25.89** | 0.90 | 0.808 | 479.4 | 4 | 0 | 48 | 0 | **Pattern A** | `HIGH_CONFIDENCE` |
| **#2** | `3_Barnase_Barstar_Electrostatic_Pos` | **0.771** | `+0.271` | **-10.09** | **39.87** | 0.89 | 0.817 | 474.7 | 16 | 4 | 0 | 0 | **Pattern A** | `HIGH_CONFIDENCE` |
| **#3** | `2_VIP1_bZIP61_GroupE_Positive` | **0.732** | `+0.232` | **-10.10** | **39.68** | 0.83 | 0.742 | 406.8 | 4 | 0 | 48 | 0 | **Pattern A** | `HIGH_CONFIDENCE` |
| **#4** | `4_VIP1_GRF6_14-3-3_IDR_Diluted_Pos` | **0.609** | `+0.109` | **-7.61** | **2644.23** | 0.52 | 0.783 | 473.2 | 4 | 2 | 36 | 0 | **Pattern B** | `MODERATE` |
| **#5** | `5_bZIP34_Proline_Homodimer_Negative` | **0.120** | `-0.380` | **+5.00** | `>10mM (NB)` | 0.38 | 0.000 | 16.8 | 0 | 0 | 5 | 6 | **Pattern C** | `LOW_OR_DECOY` |
| **#6** | `6_VIP1_IDR_Hallucination_Decoy` | **0.061** | `-0.439` | **+2.10** | `>10mM (NB)` | 0.24 | 0.000 | 0.1 | 0 | 0 | 0 | 0 | **Pattern D** | `LOW_OR_DECOY` |

* **Empirical Crystal Structure Validation (`RCSB PDB: 1BRS`):** Tested directly on the 6-chain X-ray crystal structure of Barnase-Barstar (`1BRS.cif`) in `MODE_C_EXP_CRYSTAL`, the engine strips crystal waters, isolates pairwise subunit `SASA`, scores the biological dimer (`Chains A-D`) at `UIS = 0.882` (`dG_uni = -12.33 kcal/mol`, `BSA = 1546.7 A^2`, `Pattern A`), and recovers the experimental hotspots (`D:ASP39`, `A:HIS102`, `A:ARG83`, `D:TYR29`, `A:ARG59`, `A:ARG87`, `D:ASP35`).

---

## 4. Automated 16-Gate Regression Suite (`test_unippi_regression.py`)

Run all 16 automated verification gates with:

    python3 test_unippi_regression.py

* **Gates 1–5:** Frozen 6-complex benchmark parity, Mode B fractional `[0,1]` pLDDT auto-scaling, strict 80-column FoldX/Rosetta PDB sanitization, 800 mM NaCl salt-stress divergence, and real 6-chain RCSB PDB `1BRS` crystal structure evaluation.
* **Gates 6–16:** Layer 1 CLI vs. Layer 2 Browser JS constant parity, multi-modal Protein-DNA (`OP1/OP2`) & V2000 `.sdf` ligand extraction, PTM phosphoserine (`SEP`) & CPK21/PP2A phosphorylation switch (`-7.61 -> -9.81 kcal/mol` with Gibbs-van 't Hoff `dG <-> Kd` coupling), Channel 8 divalent metal (`Co2+/Ca2+`) coordination, and Phase 5 release packaging hygiene.

---

## 5. Quick-Start Usage

### Option A: Client-Side Web Application (`index.html` or Local Launchers)
Open **https://soumyajitmandal06-gif.github.io/unippi-owl/** or run `./start-local.sh` (Linux/macOS) / `start-local.cmd` (Windows) or open `index.html` directly in any browser.

### Option B: Python CLI (`unippi_cli.py`), Nextflow (`main.nf`), and Snakemake (`Snakefile`)

    # 1. Generate and run the 6-complex benchmark suite with custom point mutations:
    python3 unippi_cli.py --generate-demo ./demo_structures -o ./demo_results --mutate A:LEU5PRO,A:PRO5LEU,A:ASP2ALA,A:TYR5ALA

    # 2. Screen a directory under 400 mM NaCl salt stress and CPK21 phosphorylation:
    python3 unippi_cli.py -i ./my_af3_outputs -o ./screening_results --temp 298.15 --ionic-strength 400.0 --ph 7.2 --ptm-state phospho

    # 3. Run via Nextflow DSL2 or Snakemake:
    nextflow run main.nf --input_dir ./demo_structures --output_dir ./demo_results --ptm_state native
    snakemake --cores 1


## Empirical Benchmark & Ablation Summary (`benchmark_conference_suite.py`)
- **24-Complex Affinity Benchmark:** Pearson `r = 0.984` (Binders `RMSE = 0.30 kcal/mol`), `ROC-AUC = 1.000`, `MCC = 1.00` (`./demo_results/EMPIRICAL_AFFINITY_BENCHMARK.tsv`)
- **32-Mutation SKEMPI 2.0 + VIP1 Benchmark:** Pearson `r = 0.997`, Spearman `rho = 0.995`, `RMSE = 0.22 kcal/mol` (`./demo_results/SKEMPI2_MUTATION_BENCHMARK.tsv`)
- **5-Variant Ablation & 14.0 A Voxel Scaling:** `53.8x` pair-pruning speedup (`./demo_results/ABLATION_STUDY_REPORT.tsv`, `./demo_results/RUNTIME_SCALING_BENCHMARK.tsv`)
