# Sampling Full IDPs (sample_idp.py)

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [configs/sample.yml](configs/sample.yml)
- [sample_idp.py](sample_idp.py)

</details>



The `sample_idp.py` script serves as the primary entrypoint for generating all-atom structural ensembles of Intrinsically Disordered Proteins (IDPs). It leverages a diffusion-based transformer model to reverse-diffuse random noise into physically plausible protein conformers, optionally guided by experimental constraints.

### Sources:
- [sample_idp.py:1-20]()
- [README.md:1-4]()

---

## CLI Arguments and Configuration

The script uses a combination of command-line arguments and a YAML configuration file to define the sampling environment.

| Argument | Description |
| :--- | :--- |
| `seq` | The amino acid sequence of the IDP to be sampled. |
| `ckpt_path` | Path to the trained IDPForge model weights (`.ckpt` or `.pt`). |
| `output_dir` | Directory where generated PDB files will be stored. |
| `sample_cfg` | Path to the YAML configuration file (e.g., `configs/sample.yml`). |
| `--batch` | Number of conformers to generate in a single GPU/CPU pass (default: 32). |
| `--nconf` | Total target number of conformers for the ensemble (default: 100). |
| `--cuda` | Flag to enable GPU acceleration. |
| `--ss_db` | Path to a secondary structure database (overrides config `data_path`). |

### Sources:
- [sample_idp.py:172-183]()
- [configs/sample.yml:1-56]()

---

## Initialization and Data Flow

The sampling process begins by initializing the `Diffuser` and `Denoiser` objects, which manage the noise schedules and the iterative reverse-diffusion steps. The `IDPForge` model is then loaded with weights, specifically extracting the Exponential Moving Average (EMA) parameters.

### System Entity Mapping: Initialization
The following diagram illustrates the relationship between the CLI entrypoint and the core classes defined in the `idpforge` module.

```mermaid
graph TD
    subgraph "CLI Entrypoint (sample_idp.py)"
        MAIN["main()"]
        ARGS["argparse"]
    end

    subgraph "Core Entities (idpforge/)"
        MODEL["IDPForge Class (model.py)"]
        DIFF["Diffuser Class (utils/diff_utils.py)"]
        DENOISE["Denoiser Class (utils/diff_utils.py)"]
        OUT["output_to_pdb() (misc.py)"]
    end

    ARGS --> MAIN
    MAIN --> DIFF
    MAIN --> DENOISE
    MAIN --> MODEL
    MODEL --> DENOISE
    MAIN --> OUT
```

### Sources:
- [sample_idp.py:31-60]()
- [idpforge/model.py:16-16]()
- [idpforge/utils/diff_utils.py:17-17]()

---

## Secondary Structure Database Lookup

IDPForge utilizes a secondary structure (SS) database to provide structural context during sampling. The script follows a fallback logic to assign SS strings to the target sequence:

1.  **Direct Path**: If `sec_path` is provided in the config, it reads SS strings directly from that file [sample_idp.py:105-109]().
2.  **Database Lookup**: It uses `fetch_sec_from_seq` to find similar sequences in a pickle-based database (`example_data.pkl`) and retrieves their associated SS [sample_idp.py:90-100]().
3.  **Fallback**: If lookup fails, it defaults to an "all-coil" (`C`) assignment [sample_idp.py:101-103]().

### Sources:
- [sample_idp.py:89-109]()
- [idpforge/utils/prep_sec.py:19-19]()

---

## Generation Loop and Batching

The generation process is wrapped in a `while` loop that persists until the `nsample` target is met. This design allows for **resume-capability**: the script checks the output directory for existing files using `count_done()` and only generates the remaining required conformers [sample_idp.py:122-140]().

### Sampling Logic Flow

```mermaid
sequenceDiagram
    participant S as sample_idp.py
    participant D as Denoiser
    participant M as IDPForge Model
    participant P as output_to_pdb

    loop While current_count < nsample
        S->>D: init_samples(seq_list)
        D-->>S: xt_list, tor_list (Gaussian/IGSO3 noise)
        S->>M: model.sample(denoiser, ss_list, xt_list, ...)
        Note over M: Iterative T-step Reverse Diffusion
        M-->>S: outputs (Final Coordinates)
        S->>P: Save batch to PDB
        S->>S: count_done()
    end
```

### Sources:
- [sample_idp.py:142-170]()
- [sample_idp.py:122-136]()

---

## Experimental Guidance Integration

If enabled in the configuration (`potential: true`), the script initializes guidance potentials that influence the diffusion gradient. This allows the generated ensemble to satisfy experimental constraints such as Paramagnetic Relaxation Enhancement (PRE) or Radius of Gyration (Rg).

*   **PRE/NOE**: Converted into a contact map potential using `get_contact_map` [sample_idp.py:70-77]().
*   **Rg**: Targets a specific ensemble average radius of gyration [sample_idp.py:78-81]().
*   **Hyperparameters**: Guidance is controlled by `timescale` (when guidance starts) and `grad_clip` (maximum gradient magnitude) [sample_idp.py:64-68]().

### Sources:
- [sample_idp.py:64-87]()
- [idpforge/utils/np_utils.py:71-71]()

---

## Output Naming and Relaxation

Generated conformers are saved using the `output_to_pdb` function. The naming convention depends on whether AMBER relaxation is performed:

1.  **Raw Outputs**: Named `{idx}_raw.pdb` if `no_relax` is true [sample_idp.py:111-113]().
2.  **Validated/Relaxed**: Named `{idx}_validated.pdb` if the full relaxation and validation pipeline is active [sample_idp.py:114-117]().

The `next_available_idx()` function ensures that new batches do not overwrite existing files by finding the lowest unused integer prefix in the output directory [sample_idp.py:125-136]().

### Sources:
- [sample_idp.py:111-120]()
- [sample_idp.py:158-160]()
- [idpforge/misc.py:18-18]()

---