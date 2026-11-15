# Infrastructure and Deployment

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.dockerignore](.dockerignore)
- [.gitignore](.gitignore)
- [Dockerfile](Dockerfile)

</details>



This page provides a high-level overview of the containerized deployment model and the development environment configuration for the IDP-o project. IDP-o is designed to leverage GPU acceleration for both sequence searching and structural assembly, necessitating a specific stack of high-performance computing dependencies.

### Deployment Strategy

IDP-o is deployed primarily as a containerized application to ensure reproducibility across different high-performance computing (HPC) environments. The deployment model focuses on providing a stable environment for JAX-based computations and GPU-accelerated operations.

The core orchestration logic is encapsulated in the `build_ensemble.py` script, which serves as the primary entrypoint for the containerized environment.

For detailed information on the container configuration, see [Docker Container](#4.1).

### Dependency Stack

The IDP-o infrastructure relies on a specialized stack of libraries to handle large-scale structural data and accelerated geometric computations:

*   **GPU Acceleration:** Utilizes `nvcr.io/nvidia/jax` as the base image to provide a robust JAX environment [Dockerfile:1-1]().
*   **Structural Data Handling:** Uses `foldcomp` for efficient binary structure storage and `nerfax` for JAX-native Natural Extension Reference Frame (NeRF) coordinate reconstruction [Dockerfile:3-3]().
*   **High-Performance Storage:** Employs `tables` (PyTables) for HDF5 ensemble storage and `hirola` for fast hash-table operations [Dockerfile:3-3]().
*   **Vectorized Searching:** Leverages `cupy-cuda12x` for GPU-accelerated byte-stream matching during the fragment search phase [Dockerfile:3-3]().

### Development and Repository Layout

The repository is structured to separate core logic from execution scripts and assets. It includes configurations for standard Python tooling to maintain code quality and environment consistency.

Key aspects of the development environment include:
*   **Source Organization:** A dedicated `scripts/` directory for execution logic and an `assets/` directory for static resources.
*   **Environment Management:** Support for standard `venv` or `conda` environments, with a comprehensive `.gitignore` to prevent local artifacts and large datasets (e.g., `IDR-30`, `IDRome`) from being committed [Dockerfile:5-5](), [.gitignore:122-130](), [.gitignore:168-172]().
*   **Tooling:** Integration with `pytest` for testing, `mypy` for static type checking, and `Jupyter` for interactive development [.gitignore:51-51](), [.gitignore:142-142](), [.gitignore:79-79]().

For details on the repository structure and local setup, see [Development Environment and Repository Layout](#4.2).

### System Mapping

The following diagrams illustrate the relationship between the conceptual infrastructure components and the specific files or configurations that define them.

#### Infrastructure to Code Mapping
This diagram maps the high-level deployment components to their respective definitions in the codebase.

```mermaid
graph TD
    subgraph "Deployment Entities"
        DOCKER["Docker Container"]
        ENV["Environment Variables"]
        EP["System Entrypoint"]
    end

    subgraph "Code Implementation"
        DFILE["Dockerfile"]
        PTPATH["PYTHONPATH"]
        B_ENS["scripts/build_ensemble.py"]
    end

    DOCKER -->|defined by| DFILE
    ENV -->|sets| PTPATH
    EP -->|executes| B_ENS

    DFILE -- "line 1" --> DOCKER
    DFILE -- "line 6" --> PTPATH
    DFILE -- "line 8" --> B_ENS
```
**Sources:** [Dockerfile:1-8]()

#### Dependency and Asset Mapping
This diagram shows how external dependencies and local file exclusions are managed within the repository structure.

```mermaid
graph LR
    subgraph "External Dependencies"
        JAX["nvcr.io/nvidia/jax"]
        CUPY["cupy-cuda12x"]
        NERF["nerfax"]
    end

    subgraph "Repository Config"
        GITIGN[".gitignore"]
        DOCKIGN[".dockerignore"]
    end

    subgraph "Excluded Assets"
        D_CACHE[".docker-cache"]
        IDR_DATA["IDR-30 / IDRome"]
    end

    JAX -->|Base Image| DFILE["Dockerfile"]
    CUPY -->|pip install| DFILE
    NERF -->|git+https| DFILE
    
    D_CACHE -.->|ignored by| DOCKIGN
    IDR_DATA -.->|ignored by| GITIGN
```
**Sources:** [Dockerfile:1-3](), [.dockerignore:1-1](), [.gitignore:168-172]()

---
**Sources:**
- [Dockerfile:1-8]()
- [.gitignore:1-172]()
- [.dockerignore:1-1]()

---