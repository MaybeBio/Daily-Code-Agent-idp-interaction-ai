# Docker Container

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.dockerignore](.dockerignore)
- [Dockerfile](Dockerfile)

</details>



The IDP-o pipeline is distributed as a containerized environment to ensure reproducible execution across different GPU-accelerated infrastructures. The container encapsulates the complex intersection of JAX-based structural biology tools, NVIDIA CUDA libraries, and high-performance data structures required for fragment-based ensemble generation.

## Base Image and Environment

The container is built upon the official NVIDIA JAX image: `nvcr.io/nvidia/jax:24.10-py3` [Dockerfile:1-1](). This base provides a pre-configured environment with `JAX`, `CUDA`, and `cuDNN` optimized for NVIDIA hardware, which is critical for the vectorized operations in the assembly and extraction stages.

The `PYTHONPATH` is explicitly set to `/IDP-o/` to ensure that the orchestration scripts can resolve internal module dependencies when executed from the root [Dockerfile:6-6]().

### Build Exclusions
To maintain a lean image size and prevent local development artifacts from polluting the container, the `.dockerignore` file excludes the `.docker-cache` directory [.dockerignore:1-1]().

## Dependency Stack

The container installs several specialized libraries via `pip` to handle structural data, high-speed indexing, and coordinate geometry [Dockerfile:3-3]().

| Dependency | Purpose in IDP-o |
| :--- | :--- |
| `nerfax` | A JAX-based implementation of the Natural Extension Reference Frame (NeRF). Used for backbone reconstruction and hydrogen placement [Dockerfile:3-3](). |
| `foldcomp` | Provides the binary decompression logic for `.fcmp` files to retrieve discretized torsion angles [Dockerfile:3-3](). |
| `cupy-cuda12x` | Enables GPU-accelerated byte-stream matching for sequence fragment searches [Dockerfile:3-3](). |
| `hirola` | High-performance hash maps used for efficient indexing and deduplication of fragment hits [Dockerfile:3-3](). |
| `tables` (PyTables) | Manages the HDF5 storage backend for intermediate structural fragments and final ensembles [Dockerfile:3-3](). |
| `pynmrstar` | Provides parsing capabilities for NMR-STAR files, supporting experimental data integration [Dockerfile:3-3](). |
| `joblib` | Facilitates CPU-bound parallelization for tasks such as FASTA indexing and batch file I/O [Dockerfile:3-3](). |

**Sources:**
- [Dockerfile:1-6]()
- [.dockerignore:1-1]()

## Execution Model and Entrypoint

The container is designed as an executable tool rather than a persistent service. The `ENTRYPOINT` is configured to trigger the single-sequence ensemble generation pipeline [Dockerfile:8-8]().

### Code Entity Space to Natural Language Space

The following diagram illustrates how the Docker configuration maps to the execution of the Python source code within the container filesystem.

**Container Execution Mapping**
```mermaid
graph TD
    subgraph "Docker Host"
        A["docker run idp-o --sequence ATGC..."]
    end

    subgraph "Container Image"
        B["ENTRYPOINT"] -- "executes" --> C["/IDP-o/scripts/build_ensemble.py"]
        D["ENV PYTHONPATH"] -- "points to" --> E["/IDP-o/"]
        
        subgraph "Internal Modules"
            E --> F["scripts/fasta_search_in_foldcomp_database.py"]
            E --> G["scripts/extract_structures_from_foldcomp_database.py"]
            E --> H["scripts/join_fragments.py"]
        end
    end

    A --> B
    C -.-> F
    C -.-> G
    C -.-> H
```

**Sources:**
- [Dockerfile:5-8]()

## Data Flow and Dependency Interaction

The containerized environment facilitates the flow of data from raw FoldComp databases to refined ensembles by providing the necessary C++ extensions (via `foldcomp`) and GPU kernels (via `cupy` and `jax`).

**Dependency Interaction Diagram**
```mermaid
graph LR
    subgraph "External Input"
        DB[".fcmp / .fasta DB"]
        SEQ["Input Sequence"]
    end

    subgraph "IDP-o Container Space"
        subgraph "Search Phase"
            CP["cupy-cuda12x"]
            HI["hirola"]
        end

        subgraph "Extraction Phase"
            FC["foldcomp"]
            JX["JAX"]
        end

        subgraph "Assembly Phase"
            NX["nerfax"]
            TB["tables (PyTables)"]
        end
    end

    DB --> CP
    SEQ --> CP
    CP --> HI
    HI --> FC
    FC --> JX
    JX --> NX
    NX --> TB
    TB --> OUT["Final Ensemble (.h5)"]
```

**Sources:**
- [Dockerfile:3-3]()
- [Dockerfile:8-8]()

---