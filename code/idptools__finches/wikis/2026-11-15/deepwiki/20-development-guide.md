# Development Guide

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [MANIFEST.in](MANIFEST.in)
- [finches/tests/test_data/test_sequences.py](finches/tests/test_data/test_sequences.py)
- [finches/tests/test_data/update_test_data.py](finches/tests/test_data/update_test_data.py)
- [finches/tests/test_epsilon_calculation.py](finches/tests/test_epsilon_calculation.py)
- [pyproject.toml](pyproject.toml)
- [setup.cfg](setup.cfg)
- [setup.py](setup.py)

</details>



This document provides essential information for developers contributing to the FINCHES codebase. It covers the build system, testing framework, packaging configuration, and development workflow. For information about using FINCHES as an end user, see [User Guide](#3). For detailed API documentation, see [API Reference](#5).

## Development Environment Overview

FINCHES uses a modern Python packaging approach with `pyproject.toml` configuration and Cython extensions for performance-critical components. The development workflow centers around pytest for testing and setuptools for building.

### Build System Architecture

```mermaid
graph TB
    subgraph "Build Configuration"
        PYPROJECT["pyproject.toml<br/>Modern packaging config"]
        SETUP["setup.py<br/>Cython compilation"]
        MANIFEST["MANIFEST.in<br/>File inclusion rules"]
    end
    
    subgraph "Dependencies"
        BUILD_DEPS["Build Dependencies<br/>setuptools, wheel, cython, numpy"]
        RUNTIME_DEPS["Runtime Dependencies<br/>numpy, scipy, pandas, metapredict"]
        TEST_DEPS["Test Dependencies<br/>pytest, pytest-runner"]
    end
    
    subgraph "Cython Compilation"
        PYX["matrix_manipulation.pyx<br/>Performance-critical code"]
        COMPILED["finches.utils.matrix_manipulation<br/>Compiled extension"]
    end
    
    subgraph "Version Management"
        VERSIONINGIT["versioningit<br/>Git-based versioning"]
        VERSION_FILE["finches/_version.py<br/>Generated version"]
    end
    
    PYPROJECT --> BUILD_DEPS
    PYPROJECT --> RUNTIME_DEPS  
    PYPROJECT --> TEST_DEPS
    PYPROJECT --> VERSIONINGIT
    
    SETUP --> PYX
    PYX --> COMPILED
    
    VERSIONINGIT --> VERSION_FILE
    
    MANIFEST --> PYX
    MANIFEST --> "finches/data/*"
```

The build system follows modern Python packaging standards while accommodating legacy Cython compilation requirements.

**Sources:** [pyproject.toml:1-76](), [setup.py:1-32](), [MANIFEST.in:1-17]()

### Package Configuration Details

| Configuration File | Purpose | Key Components |
|-------------------|---------|----------------|
| `pyproject.toml` | Primary package metadata | Dependencies, build backend, version config |
| `setup.py` | Cython extension compilation | Extension definition, numpy integration |
| `MANIFEST.in` | Package file inclusion | Data files, Cython sources, exclusions |
| `setup.cfg` | Tool configurations | Coverage, code formatting, test aliases |

The package uses `versioningit` for automatic version management based on Git tags and commits, writing version information to `finches/_version.py`.

**Sources:** [pyproject.toml:60-76](), [setup.cfg:1-24]()

## Testing Framework Architecture

```mermaid
graph TB
    subgraph "Test Organization"
        TEST_DIR["finches/tests/"]
        TEST_MAIN["test_epsilon_calculation.py<br/>Main test suite"]
        TEST_DATA_DIR["test_data/<br/>Test fixtures"]
    end
    
    subgraph "Test Data Management"
        TEST_SEQS["test_sequences.py<br/>Test sequence definitions"]
        UPDATE_DATA["update_test_data.py<br/>Reference data generation"]
        NPZ_FILES["*.npz files<br/>Precomputed reference results"]
    end
    
    subgraph "Test Categories"
        MATRIX_TESTS["Matrix Constructor Tests<br/>Homotypic, heterotypic matrices"]
        MANIPULATION_TESTS["Matrix Manipulation Tests<br/>Masking, flattening operations"]
        EPSILON_TESTS["Epsilon Calculation Tests<br/>Vector calculations, weighting"]
        FH_TESTS["Flory-Huggins Tests<br/>Phase diagram generation"]
    end
    
    subgraph "Test Models"
        MPIPI_MODEL["mpipi_model('mPiPi_GGv1')<br/>Primary test model"]
        CALVADOS_MODEL["calvados_model<br/>Secondary test model"]
        IMC["Interaction_Matrix_Constructor<br/>Test harness"]
    end
    
    TEST_MAIN --> MATRIX_TESTS
    TEST_MAIN --> MANIPULATION_TESTS
    TEST_MAIN --> EPSILON_TESTS
    
    TEST_SEQS --> MATRIX_TESTS
    UPDATE_DATA --> NPZ_FILES
    NPZ_FILES --> MATRIX_TESTS
    
    MPIPI_MODEL --> IMC
    IMC --> MATRIX_TESTS
```

The testing framework uses regression testing against precomputed reference data stored in NumPy archives.

**Sources:** [finches/tests/test_epsilon_calculation.py:1-187](), [finches/tests/test_data/test_sequences.py:1-21]()

### Test Data Management

The test suite relies on a sophisticated data management system for regression testing:

```mermaid
flowchart TD
    subgraph "Test Sequence Definitions"
        TEST_SEQS["test_sequences.py<br/>7 predefined sequences"]
        CONDITIONS["test_condition_dict<br/>Salt, pH, dielectric values"]
    end
    
    subgraph "Reference Data Generation"
        UPDATE_FUNCS["update_test_data.py functions"]
        HOMOTYPIC["write_test_pairwise_homotypic_matrix()"]
        HETEROTYPIC["write_test_pairwise_heterotypic_matrix()"]
        WEIGHTED["write_test_weighted_matrix()"]
        MANIPULATION["write_test_matrix_manipulation()"]
        EPSILON["write_test_seq_epsilon_and_vectors()"]
        FH["write_FH_out_data()"]
    end
    
    subgraph "Stored Reference Data"
        NPZ_HOMO["test_mPiPi_GGv1_homotypic_matrix.npz"]
        NPZ_HETERO["test_mPiPi_GGv1_heterotypic_matrix.npz"]
        NPZ_WEIGHTED["test_mPiPi_GGv1_weighted_matrix.npz"]
        NPZ_MANIP["test_matrix_manipulation.npz"]
        NPZ_EPSILON["mPiPi_GGv1_seq_epsilon_and_vectors.npz"]
        NPZ_FH["mPiPi_GGv1_FH_outdata.npz"]
    end
    
    TEST_SEQS --> HOMOTYPIC
    CONDITIONS --> FH
    
    HOMOTYPIC --> NPZ_HOMO
    HETEROTYPIC --> NPZ_HETERO
    WEIGHTED --> NPZ_WEIGHTED
    MANIPULATION --> NPZ_MANIP
    EPSILON --> NPZ_EPSILON
    FH --> NPZ_FH
```

**Critical Note:** The `update_test_data.py` functions can overwrite reference data. Use with extreme caution to avoid invalidating the test suite.

**Sources:** [finches/tests/test_data/update_test_data.py:1-204](), [finches/tests/test_data/test_sequences.py:11-20]()

## Core Test Categories

### Matrix Construction Tests

Tests for the `Interaction_Matrix_Constructor` class verify fundamental matrix operations:

- **Homotypic matrices:** Self-interaction calculations using `calculate_pairwise_homotypic_matrix()`
- **Heterotypic matrices:** Cross-interaction calculations using `calculate_pairwise_heterotypic_matrix()`  
- **Weighted matrices:** Charge and aliphatic weighting using `calculate_weighted_pairwise_matrix()`

### Matrix Manipulation Tests

Tests for core matrix processing functions in `epsilon_calculation` module:

- **Attractive/repulsive separation:** `get_attractive_repulsive_matrixes()` with threshold `-0.15`
- **Matrix masking:** `mask_matrix()` with binary masks
- **Vector flattening:** `flatten_matrix_to_vector()` with orientation parameters

### Epsilon Calculation Tests

Tests for sequence-level epsilon value calculations:

- **Vector generation:** `get_sequence_epsilon_vectors()` with various weighting options
- **Prefactor handling:** Custom prefactor and baseline parameters
- **Weighting variations:** Charge weighting and aliphatic weighting toggles

**Sources:** [finches/tests/test_epsilon_calculation.py:15-187]()

## Development Workflow

### Setting Up Development Environment

1. **Clone and install in development mode:**
```bash
git clone https://github.com/idptools/finches
cd finches
pip install -e ".[test]"
```

2. **Run the test suite:**
```bash
pytest finches/tests/
```

3. **Build Cython extensions explicitly:**
```bash
python setup.py build_ext --inplace
```

### Key Development Files

| File Path | Purpose | When to Modify |
|-----------|---------|----------------|
| `pyproject.toml` | Package metadata, dependencies | Adding dependencies, changing build config |
| `setup.py` | Cython compilation | Modifying Cython extensions |
| `finches/tests/test_epsilon_calculation.py` | Main test suite | Adding new functionality tests |
| `finches/tests/test_data/test_sequences.py` | Test sequences | Adding new test cases |

### Adding New Tests

When adding new functionality:

1. **Add test sequences** to `test_sequences.py` if needed
2. **Generate reference data** using appropriate `update_test_data.py` functions
3. **Write test functions** in `test_epsilon_calculation.py` following existing patterns
4. **Verify** tests pass with `pytest finches/tests/test_your_new_test.py`

### Modifying Cython Code

For changes to `finches/utils/matrix_manipulation.pyx`:

1. **Modify the `.pyx` file**
2. **Recompile:** `python setup.py build_ext --inplace`
3. **Run performance-critical tests** to verify functionality
4. **Update reference data** if computational results change

**Sources:** [setup.py:12-23](), [MANIFEST.in:6-8]()

## Package Distribution

The package uses modern Python packaging with backward compatibility for Cython:

- **Build backend:** `setuptools.build_meta` 
- **Version management:** Git-based via `versioningit`
- **License:** CC BY-NC 4.0 (non-commercial)
- **Python compatibility:** >=3.7

The `MANIFEST.in` ensures all necessary files are included in distributions, including Cython sources, pickle parameter files, and reference data.

**Sources:** [pyproject.toml:1-31](), [MANIFEST.in:6-16]()

---