

The **AlphaFlex (AFX-IDPForge) pipeline** is a four-stage computational workflow that transforms AlphaFold2 predicted structures into physically realistic conformational ensembles for proteins containing intrinsically disordered regions (IDRs). Unlike the base IDPForge system—which operates on individual disordered segments in isolation—AlphaFlex orchestrates end-to-end processing: from IDR classification and database labeling, through per-region template construction and diffusion-based conformer generation, to Monte Carlo stitching of full-length multi-domain models with AMBER relaxation and structural validation. The pipeline is designed for high-throughput proteome-scale processing while remaining resumable and HPC-parallelizable at every stage.

Sources: [README.md](/AlphaFlex/README.md#L1-L146), [config.py](/AlphaFlex/config.py#L1-L123)

## Pipeline Architecture

The AlphaFlex pipeline processes proteins through four sequential stages, each consuming the outputs of the prior stage and writing results to a dedicated subdirectory under `Pipeline_Outputs/`. The central `config.py` module governs all tunable parameters and path defaults, ensuring consistent behavior across scripts.

```mermaid
flowchart TD
    subgraph Inputs
        DB[("AlphaFlex_database_Nov2025.json<br/>IDR boundaries + interactions")]
        LEN[("AF2_9606_HUMAN_v4_num_residues.json<br/>Residue counts")]
        PDB[("AlphaFold2 PDB Library")]
        KNOT[("knot_screening.json<br/>Pre-computed topology")]
    end

    subgraph Step1["Step 1: Case Labeling"]
        S1[Step_1_case_label.py] -->|IDR type<br/>+ category| S1B[Step_1B_subset_label.py]
    end

    subgraph Step2["Step 2: Template Creation"]
        S2[Step_2_mk_ldr_template.py]
        S2 -->|Tail/Loop| MK_LDR[mk_ldr_template.py<br/>Static scaffold]
        S2 -->|Linker| MK_FLEX[mk_flex_template.py<br/>Dynamic scaffold]
        S2 -->|Category 0| IDP_JSON[idp_cases_to_run.json]
    end

    subgraph Step3["Step 3: Conformer Generation"]
        S3[Step_3_sample_conformer.py]
        S3 -->|IDR regions| SAMPLE_LDR[sample_ldr.py<br/>SE(3) diffusion]
        S3 -->|Full IDP| SAMPLE_IDP[sample_idp.py<br/>Unconditioned diffusion]
        S3 --> VALIDATE[Relax → Repair → Validate loop]
    end

    subgraph Step4["Step 4: Stitching & Minimization"]
        S4[Step_4_ldr_stitch.py]
        S4 -->|Single IDR| FASTPASS[Fast-Pass:<br/>Direct ensemble assembly]
        S4 -->|Multiple IDRs| MC[Monte Carlo Stitching:<br/>Kinematic chain assembly]
        MC --> RELAX[AMBER Relaxation<br/>+ Repair + Validation]
    end

    Inputs --> Step1
    Step1 -->|"labeled DB +<br/>subset .txt"| Step2
    PDB --> Step2
    Step2 -->|".npz templates"| Step3
    Step3 -->|"validated .pdb<br/>conformer pools"| Step4
    PDB --> Step4
    KNOT --> Step3
    KNOT --> Step4

    Step4 --> RESULT[("Final ensemble PDBs<br/>Category/Length/Protein/")]
```

Sources: [README.md](/AlphaFlex/README.md#L49-L141), [config.py](/AlphaFlex/config.py#L26-L123)

## IDR Classification System

The foundational abstraction of the AlphaFlex pipeline is its **four-tier IDR topology hierarchy**, which determines how each disordered region is processed downstream. The classification is hierarchical—higher categories imply greater structural complexity and require more sophisticated stitching strategies.

| Category | Name | Condition | Template Strategy | Stitching Complexity |
|:--------:|:----:|:----------|:------------------|:-------------------|
| **0** | IDP | Entire protein is disordered (no folded domains) | Unconditioned diffusion via `sample_idp.py` | Fast-pass (no stitching needed) |
| **1** | Tails | Only N/C-terminal IDRs; all internal regions are folded | Static scaffold via `mk_ldr_template.py` (frozen folded regions) | Low |
| **2** | Linkers | Internal IDRs between **non-interacting** folded domains | Dynamic scaffold via `mk_flex_template.py` (domains treated as separate rigid bodies with random displacements) | Medium (multi-body alignment) |
| **3** | Loops | Internal IDRs between **interacting** folded domains (mean PAE < 15 Å) | Static scaffold via `mk_ldr_template.py` (frozen folded regions) | High (constrained topology) |

Each IDR within a hybrid protein (Categories 1–3) receives a **disordered domain label** (`D1`, `D2`, …) and a **flanking domain annotation** identifying the adjacent folded domains (`F1`, `F2`, …). The interaction status between flanking domains—derived from the AlphaFold2 predicted aligned error (PAE) matrix—distinguishes Linkers from Loops and directly controls which template generation strategy is invoked.

Sources: [Step_1_case_label.py](/AlphaFlex/Step_1_case_label.py#L1-L10), [Step_1_case_label.py](/AlphaFlex/Step_1_case_label.py#L111-L197)

## Data Inputs and Resources

The pipeline depends on four pre-computed data resources stored in `AlphaFlex/Data_Inputs/`:

| File | Format | Content | Used By |
|------|--------|---------|---------|
| `AlphaFlex_database_Nov2025.json` | JSON dict | Per-protein IDR boundaries (`idrs`), mean inter-region PAEs (`mean_pae`), and folded domain interactions | Steps 1, 2, 4 |
| `AF2_9606_HUMAN_v4_num_residues.json` | JSON dict | Total residue counts from AlphaFold2 9606 Human v4 | Steps 1, 1B |
| `knot_screening.json` | JSON dict | Per-domain knot presence/absence (`label`) and knot type (`closure_polys`) | Steps 3, 4 |
| `Test_Structures/O14653.pdb` | PDB | Sample Category 3 protein for pipeline testing | Steps 2, 3, 4 |

The master database is the central information artifact: each key is a UniProt accession, and values encode the spatial relationship between folded and disordered regions. Step 1 augments this database with computed labels, producing the **labeled database** that drives all downstream processing.

Sources: [README.md](/AlphaFlex/README.md#L7-L15), [config.py](/AlphaFlex/config.py#L18-L22)

## Step-by-Step Pipeline Walkthrough

### Step 1: Case Labeling → `Step_1_Labeling/`

This step reads the master database and classifies every protein into the four-tier category system. For each IDR, it computes: the residue **range**, the **type** (Tail/Linker/Loop/IDP), the ordinal **label** (D1, D2, …), and the **flanking domains**. The classification logic applies a strict hierarchy: if any Loop IDR is present, the protein is Category 3 regardless of other IDR types; Linkers take precedence over Tails; and a single IDR spanning the entire sequence yields Category 0.

The output is an augmented JSON database and an `idr_type_summary.txt` report showing the distribution across categories.

Sources: [Step_1_case_label.py](/AlphaFlex/Step_1_case_label.py#L30-L199), [README.md](/AlphaFlex/README.md#L49-L62)

### Step 1B: Subset Filtering → `Step_1_Labeling/custom_subsets/`

Step 1B creates filtered UniProt ID lists from the labeled database for targeted processing. It operates in two modes:

- **Basic mode** (length-only): Filters proteins by total residue count within `[min_len, max_len]`. All matching proteins are included regardless of IDR composition.
- **Advanced mode**: Additionally filters by per-type IDR counts (Tail/Linker/Loop) and a global IDR length range applied uniformly to every IDR in a protein. Count matching can be **exact** (default) or **minimum** (`--min_mode`).

Outputs include the filtered ID list (`.txt`), a per-protein report table, and—when advanced filters are active—an IDR length distribution histogram with bin counts.

Sources: [Step_1B_subset_label.py](/AlphaFlex/Step_1B_subset_label.py#L1-L64), [README.md](/AlphaFlex/README.md#L64-L80)

### Step 2: Template Creation → `Step_2_Templates/`

Step 2 generates per-IDR `.npz` structure templates for the diffusion model. For each protein in the subset ID list, it iterates over every labeled IDR and dispatches the appropriate template script based on IDR type:

- **Tail/Loop IDRs** → `mk_ldr_template.py`: Creates a **static scaffold** where all regions outside the specified IDR are held frozen. The disordered region is seeded with conformations distributed along the vector connecting its flanking folded domains.
- **Linker IDRs** → `mk_flex_template.py`: Creates a **dynamic scaffold** where the two adjacent folded domains are designated as separate rigid bodies and randomly displaced relative to one another, mimicking the flexibility of non-interacting domains.
- **Category 0 (Full IDP)** → No template is created; instead, the protein's sequence is logged to `idp_cases_to_run.json` for separate processing via `sample_idp.py` in Step 3.

An optional **truncation mode** (`--truncate`) builds each template on only the IDR plus its adjacent folded domains, reducing the template size for large proteins. Truncation sidecars (`.trunc.json`) record the offset and window boundaries for back-numbering during Step 4 reconstruction. A **graft-mode** alternative (`--max_residues`) caps template size while recording graft specifications for the `graft_back` utility.

Sources: [Step_2_mk_ldr_template.py](/AlphaFlex/Step_2_mk_ldr_template.py#L1-L66), [Step_2_mk_ldr_template.py](/AlphaFlex/Step_2_mk_ldr_template.py#L160-L286), [README.md](/AlphaFlex/README.md#L82-L95)

### Step 3: Conformer Generation → `Step_3_Raw_Conformers/`

Step 3 generates validated conformer pools for each IDR template through an iterative **generate–relax–repair–validate** loop. Each IDR is processed independently, so proteins with multiple IDRs produce separate per-IDR conformer pools. The loop proceeds as follows:

1. **Generate + Relax**: Invokes `sample_ldr.py` (for IDR regions) or `sample_idp.py` (for full IDPs) to produce diffusion conformers, which are immediately relaxed via AMBER minimization using the configuration from `configs/sample.yml`.
2. **Repair**: Checks each relaxed structure for D-amino acids (chirality inversions) and broken histidine ring bonds. Applies fixes and re-relaxes if any repairs were made.
3. **Validate**: Runs unified validation checking chirality, bond integrity, clash score (using an **adaptive smart threshold**), and backbone topology (knot detection against the pre-computed `knot_screening.json`).

The loop continues until the target conformer count (`SAMPLE_N_CONFS`, default 10) is met or the maximum total attempts (`SAMPLE_MAX_TOTAL_ATTEMPTS`, default 500) is exhausted. **State persistence** via `.step3_state.json` enables recovery from cluster preemption—on restart, orphaned relaxed files from killed runs are automatically recovered. HPC parallel execution is supported through `--total_splits` and `--split_index` for deterministic sharding across jobs.

Sources: [Step_3_sample_conformer.py](/AlphaFlex/Step_3_sample_conformer.py#L1-L143), [Step_3_sample_conformer.py](/AlphaFlex/Step_3_sample_conformer.py#L178-L241), [README.md](/AlphaFlex/README.md#L97-L113)

### Step 4: Stitching and Minimization → `Step_4_Final_Models/`

Step 4 assembles full-length protein models from the per-IDR conformer pools. Its behavior bifurcates based on IDR count:

**Single IDR / Full IDP (Fast-Pass):** Collects all validated conformers from Step 3 and combines them directly into a multi-model ensemble PDB. No stitching or additional minimization is performed.

**Multiple IDRs (Monte Carlo Stitching):** Runs a stochastic assembly loop that repeats until the target conformer count is reached or max attempts are exhausted. The kinematic chain assembly algorithm proceeds as:

1. Build a **segment map** partitioning the protein into alternating folded (F) and disordered (D) segments.
2. Seed the structure with the first folded domain from the AlphaFold2 prediction.
3. For each subsequent disordered segment, **draw a random conformer** from the per-IDR pool (round-robin with random shuffle).
4. **Superimpose** the middle 11 residues of the preceding folded domain (the "alignment stub") between the reference and conformer structures using Kabsch alignment.
5. **Overwrite** the reference structure with the aligned conformer coordinates from the stub onward. Repeat until all segments are stitched.
6. **Pre-relax repair**: flip any D-amino acids in the stitched structure.
7. **AMBER minimization** (ff14SB) with harmonic restraints on folded domains, allowing IDRs and junction regions to relax freely.
8. **Post-relax repair**: check for D-amino acids and broken HIS ring bonds introduced during relaxation; apply fixes and re-relax if needed.
9. **Validation**: chirality, bond integrity, clash score (adaptive smart threshold), backbone topology (knot detection), and **folded-domain curvature gate** (ensures relaxation did not distort well-structured regions beyond a configurable ratio).
10. Combine all valid conformers into a single multi-model ensemble PDB.

Final output is organized as `Step_4_Final_Models/<Category>/<Length_Label>/<protein_id>/<protein_id>_ensemble_n<N>.pdb`, where the length label stratifies proteins into bins (0–250, 251–500, 501–1000, etc.).

Sources: [Step_4_ldr_stitch.py](/AlphaFlex/Step_4_ldr_stitch.py#L1-L127), [Step_4_ldr_stitch.py](/AlphaFlex/Step_4_ldr_stitch.py#L218-L340), [README.md](/AlphaFlex/README.md#L115-L141)

## Configuration Reference

All pipeline parameters are centralized in `config.py` and can be overridden via CLI arguments on each step script. The following table summarizes the key tunables organized by pipeline stage:

| Parameter | Default | Stage | Description |
|-----------|---------|:-----:|-------------|
| `VERBOSE` | `True` | All | Enable detailed per-protein console logging |
| `SUBSET_MIN_LENGTH` / `SUBSET_MAX_LENGTH` | 0 / 250 | 1B | Protein length range filter |
| `SUBSET_TAIL_COUNT` / `SUBSET_LINKER_COUNT` / `SUBSET_LOOP_COUNT` | 2 / 1 / 1 | 1B | Per-type IDR count constraints |
| `SUBSET_EXACT_COUNT` | `True` | 1B | Exact count matching vs. minimum |
| `SUBSET_IDR_MIN_LENGTH` / `SUBSET_IDR_MAX_LENGTH` | None | 1B | Global IDR length range (applied to all IDRs) |
| `TEMPLATE_N_CONFS` | 200 | 2 | Number of template conformations generated |
| `TEMPLATE_SEED_SKEW` | 0.5 | 2 | Seed distance skew for tail/loop templates on [d/2, d] |
| `TEMPLATE_FOLD_PER_SIDE` | 10 | 2 | Fold residues retained per junction in graft mode |
| `SAMPLE_N_CONFS` | 10 | 3 | Target validated conformers per IDR |
| `SAMPLE_BATCH_SIZE` | 6 | 3 | Diffusion batch size |
| `SAMPLE_MAX_TOTAL_ATTEMPTS` | 500 | 3 | Maximum generation attempts per IDR |
| `SAMPLE_FOLD_CURV_RATIO` | 0.5 | 3 | Folded-domain curvature gate ratio |
| `SAMPLE_JUNCTION_KAPPA` | 0.12 | 3 | Junction curvature gate (Å⁻¹) |
| `STITCH_N_CONFORMERS` | 10 | 4 | Target stitched conformers per protein |
| `STITCH_MAX_ATTEMPTS` | 500 | 4 | Maximum stitching attempts |
| `RELAX_STIFFNESS` | 10.0 | 4 | AMBER restraint strength on folded domains |
| `ALIGNMENT_STUB_HALF_SIZE` | 5 | 4 | Half-window for alignment stub (middle 11 residues) |
| `STITCH_BASE_CLASH_THRESHOLD` | 10.0 | 4 | Base clash score threshold (adaptive) |
| `STITCH_CLASH_INCREMENT` | 5.0 | 4 | Clash threshold escalation step per attempt tier |

Sources: [config.py](/AlphaFlex/config.py#L1-L123)

## Utility Modules

The `AlphaFlex/utils/` package provides shared functionality consumed primarily by Step 4:

| Module | Purpose |
|--------|---------|
| `stitch.py` | Kinematic chain assembly, segment mapping, BioPython superimposition, round-robin conformer drawing, ensemble directory discovery, PDB loading, length binning, and truncated-pool renumbering |
| `smart_scoring.py` | Adaptive clash threshold computation that escalates per-attempt tier, balancing validation strictness against throughput |
| `graft_back.py` | Graft-mode reconstruction: maps truncated sub-ensembles back to full-length numbering and grafts conformer coordinates onto the original AlphaFold2 scaffold |
| `truncate.py` | Computes IDR truncation windows (IDR + adjacent folded domains) for large-protein template size reduction |
| `file_ops.py` | Common file I/O helpers shared across pipeline steps |

The `assemble_kinematic_chain` function in `stitch.py` is the algorithmic core of Step 4: it partitions the protein into a segment map, then iteratively aligns and overwrites disordered conformers onto the growing model using the alignment stub superimposition strategy.

Sources: [stitch.py](/AlphaFlex/utils/stitch.py#L1-L23), [stitch.py](/AlphaFlex/utils/stitch.py#L291-L398)

## HPC and Resumability

The AlphaFlex pipeline is designed for proteome-scale throughput on HPC clusters. Key operational features include:

- **Deterministic sharding**: Steps 3 and 4 accept `--total_splits` and `--split_index` arguments that partition the protein ID list into non-overlapping chunks via `numpy.array_split`, enabling embarrassingly parallel job arrays.
- **State persistence**: Step 3 writes `.step3_state.json` after each generation round, recording the cumulative attempt count. On restart, the script detects existing validated files, recovers orphaned intermediate files, and resumes from the persisted state.
- **Progress tracking**: Step 2 writes `Step_2_progress.txt` recording the last processed protein ID, allowing batch resumption after interruption.
- **Per-protein logging**: Steps 3 and 4 stream per-protein logs to `logs/StepN/perprotein/split_K/<ID>.log`, enabling real-time monitoring via `tail -f` on long-running cluster jobs.
- **Local parallelism**: Step 4 additionally supports `--workers` for multi-process local parallel processing of stitching tasks.

Sources: [Step_3_sample_conformer.py](/AlphaFlex/Step_3_sample_conformer.py#L30-L56), [Step_2_mk_ldr_template.py](/AlphaFlex/Step_2_mk_ldr_template.py#L111-L129), [README.md](/AlphaFlex/README.md#L109-L141)

> [!TIP]
> The pipeline's iterative generate-relax-repair-validate loop in Step 3 is the primary throughput bottleneck. Tuning `SAMPLE_BATCH_SIZE` and `SAMPLE_MAX_TOTAL_ATTEMPTS` in `config.py` trades GPU utilization against conformer yield—larger batches amortize model loading overhead but increase per-round memory pressure. The adaptive clash threshold in Steps 3–4 automatically relaxes validation strictness as attempts mount, preventing infinite rejection loops on structurally challenging proteins.

> [!TIP]
> For proteins with multiple IDRs, the Monte Carlo stitching in Step 4 is stochastic by design—each run produces a different ensemble. The `ALIGNMENT_STUB_HALF_SIZE` parameter (default 5, yielding an 11-residue stub) controls the rigidity of the junction alignment: smaller stubs tolerate more local distortion at domain boundaries, while larger stubs enforce stricter continuity at the cost of reduced conformational diversity.

## Where to Go Next

The AlphaFlex workflow overview establishes the end-to-end processing narrative. To dive deeper into specific stages:

- **[Case Labeling and IDR Classification](19-case-labeling-and-idr-classification)** — Detailed mechanics of the four-tier category system, interaction-based Linker/Loop discrimination, and subset filtering strategies
- **[Monte Carlo Stitching and Assembly](20-monte-carlo-stitching-and-assembly)** — In-depth treatment of the kinematic chain algorithm, alignment stub superimposition, AMBER relaxation protocol, and adaptive validation gates
- **[IDR Sampling with Folded Templates](13-idr-sampling-with-folded-templates)** — How `sample_ldr.py` uses the `.npz` templates from Step 2 to condition SE(3) backbone diffusion on folded-domain context
- **[AMBER Relaxation and Repair](15-amber-relaxation-and-repair)** — The ff14SB minimization protocol, chirality repair, and histidine ring bond correction used in Steps 3 and 4
- **[Structure Validation Pipeline](16-structure-validation-pipeline)** — Unified validation checks including adaptive clash scoring, knot detection, and folded-domain curvature gating