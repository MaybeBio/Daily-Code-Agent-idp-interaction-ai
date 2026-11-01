# Continuous Integration

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.github/workflows/soursop-ci.yml](.github/workflows/soursop-ci.yml)
- [soursop/ssutils.py](soursop/ssutils.py)
- [soursop/tests/test_ssutils.py](soursop/tests/test_ssutils.py)

</details>



This page documents SOURSOP's continuous integration and continuous deployment (CI/CD) infrastructure. The CI system ensures code quality by automatically running the full test suite across multiple platforms and Python versions on every code change.

**Scope**: This page covers the GitHub Actions workflow configuration, test matrix, coverage reporting, and cross-platform execution strategy. For information about the test suite itself, see [Testing Infrastructure](#9.2). For documentation building automation, see [Documentation Building](#9.4). For release automation, see [Packaging and Release](#9.5).

## Overview

SOURSOP uses GitHub Actions to provide automated testing on every push and pull request. The CI pipeline validates code across 9 different configurations (3 operating systems × 3 Python versions) and enforces code coverage standards through Codecov integration.

**Sources**: [.github/workflows/soursop-ci.yml:1-70]()

---

## GitHub Actions Workflow

The CI configuration is defined in a single workflow file that orchestrates all testing activities. The workflow is triggered automatically on repository events and uses a matrix strategy to parallelize testing across multiple environments.

### Workflow Definition

| Property | Value |
|----------|-------|
| **Workflow Name** | `Soursop CI` |
| **Trigger Events** | `push`, `pull_request` |
| **Configuration File** | `.github/workflows/soursop-ci.yml` |
| **Job Name** | `build` |
| **Parallelization** | Matrix-based (9 configurations) |

**Sources**: [.github/workflows/soursop-ci.yml:1-4]()

### Trigger Conditions

The workflow executes on two event types:

- **Push events**: Any push to any branch triggers the full test suite
- **Pull request events**: Opening, updating, or reopening a pull request triggers testing

This ensures that both direct commits and proposed changes are validated before integration.

**Sources**: [.github/workflows/soursop-ci.yml:2]()

---

## Test Matrix Configuration

The CI system tests SOURSOP across a comprehensive matrix of operating systems and Python versions to ensure broad compatibility.

```mermaid
graph TB
    subgraph "Test Matrix: 9 Configurations"
        subgraph Ubuntu["ubuntu-latest"]
            U37["Python 3.7"]
            U38["Python 3.8"]
            U39["Python 3.9"]
        end
        
        subgraph MacOS["macos-latest"]
            M37["Python 3.7"]
            M38["Python 3.8"]
            M39["Python 3.9"]
        end
        
        subgraph Windows["windows-latest"]
            W37["Python 3.7"]
            W38["Python 3.8"]
            W39["Python 3.9"]
        end
    end
    
    Workflow["soursop-ci.yml<br/>matrix strategy"]
    
    Workflow --> Ubuntu
    Workflow --> MacOS
    Workflow --> Windows
```

**Test Matrix: Platform and Python Version Combinations**

### Supported Platforms

| Operating System | GitHub Runner | Purpose |
|-----------------|---------------|---------|
| `ubuntu-latest` | Ubuntu Linux | Primary development/deployment platform |
| `macos-latest` | macOS | Apple Silicon/Intel Mac support |
| `windows-latest` | Windows Server | Windows user support |

### Supported Python Versions

The CI system tests against Python 3.7, 3.8, and 3.9 based on active support status referenced from https://endoflife.date/python. Python 3.10+ support is noted as requiring additional Conda compatibility work.

**Sources**: [.github/workflows/soursop-ci.yml:6-10]()

---

## Pipeline Stages

The CI pipeline executes in six sequential stages, each building on the previous stage's setup. The workflow uses explicit shell configuration to ensure consistent behavior across operating systems.

```mermaid
flowchart TD
    Trigger["Trigger Event<br/>(push or pull_request)"]
    
    Checkout["Stage 1: Checkout Code<br/>actions/checkout@v3"]
    
    CondaSetup["Stage 2: Setup Miniconda<br/>conda-incubator/setup-miniconda@v2<br/>miniconda-version: latest<br/>python-version: matrix.python-version"]
    
    CondaInit["Stage 3: Initialize Conda<br/>conda init bash"]
    
    InstallDeps["Stage 4: Install Dependencies<br/>shell: bash -l {0}<br/>conda activate test<br/>conda install --file anaconda_requirements.txt"]
    
    RunTests["Stage 5: Run Test Suite<br/>shell: bash -l {0}<br/>conda activate test<br/>python -m pytest --capture=sys"]
    
    Coverage["Stage 6: Generate Coverage<br/>shell: bash -l {0}<br/>conda activate test<br/>python -m pytest --cov=soursop --cov-report=xml -n auto"]
    
    Upload["Stage 7: Upload to Codecov<br/>codecov/codecov-action@v3<br/>fail_ci_if_error: true"]
    
    Trigger --> Checkout
    Checkout --> CondaSetup
    CondaSetup --> CondaInit
    CondaInit --> InstallDeps
    InstallDeps --> RunTests
    RunTests --> Coverage
    Coverage --> Upload
```

**CI Pipeline Stages and Data Flow**

### Stage 1: Code Checkout

The pipeline begins by checking out the repository code using `actions/checkout@v3`. This action is critical for Conda-based workflows, as noted in the configuration comments.

**Sources**: [.github/workflows/soursop-ci.yml:13-15]()

### Stage 2: Conda Environment Setup

The workflow uses `conda-incubator/setup-miniconda@v2` to establish a Conda environment:

- **Miniconda version**: `latest` (explicitly specified to ensure availability)
- **Auto-update**: Enabled via `auto-update-conda: true`
- **Python version**: Dynamically selected from matrix configuration

This approach was updated in July 2024 to ensure Miniconda installation reliability.

**Sources**: [.github/workflows/soursop-ci.yml:17-23]()

### Stage 3: Conda Initialization

The workflow explicitly runs `conda init bash` to configure the shell environment. This step is necessary because subsequent steps use bash login shells (`bash -l {0}`).

**Sources**: [.github/workflows/soursop-ci.yml:25-27]()

### Stage 4: Dependency Installation

Dependencies are installed using `anaconda_requirements.txt` rather than a traditional `requirements.txt`. This approach allows the Conda solver to determine appropriate package versions based on Python version and library compatibility:

```bash
conda activate test
conda install --file anaconda_requirements.txt \
  --channel default \
  --channel anaconda \
  --channel conda-forge
```

The workflow searches three Conda channels in order: `default`, `anaconda`, `conda-forge`.

**Sources**: [.github/workflows/soursop-ci.yml:29-44]()

### Stage 5: Test Execution

The test suite is executed with basic output capture:

```bash
conda activate test
python -m pytest --capture=sys
```

The `--capture=sys` flag redirects stdout/stderr to pytest's capture mechanism for cleaner test output.

**Sources**: [.github/workflows/soursop-ci.yml:46-50]()

### Stage 6: Coverage Generation

A second pytest run generates code coverage data with parallel execution:

```bash
conda activate test
python -m pytest --cov=soursop --cov-report=xml -n auto
```

| Flag | Purpose |
|------|---------|
| `--cov=soursop` | Measure coverage for the `soursop` package |
| `--cov-report=xml` | Generate XML coverage report for Codecov |
| `-n auto` | Run tests in parallel using pytest-xdist (auto-detect CPU count) |

The XML report is written to `./coverage.xml` in the repository root.

**Sources**: [.github/workflows/soursop-ci.yml:52-57]()

### Stage 7: Codecov Upload

Coverage results are uploaded to Codecov using `codecov/codecov-action@v3`:

| Configuration | Value | Purpose |
|--------------|-------|---------|
| `token` | `secrets.CODECOV_TOKEN` | Authentication for Codecov API |
| `name` | `codecov-umbrella` | Identifier for this upload |
| `flags` | `unittests` | Tag coverage as from unit tests |
| `env_vars` | `OS,PYTHON` | Include OS and Python version in report |
| `fail_ci_if_error` | `true` | **Fail build if upload fails** |
| `files` | `./coverage.xml` | Coverage data file |
| `verbose` | `true` | Detailed logging |
| `path_to_write_report` | `codecov_report.txt` | Local report artifact |

The `fail_ci_if_error: true` setting is critical—it ensures that coverage reporting failures cause CI failure, preventing merges with incomplete coverage data.

**Sources**: [.github/workflows/soursop-ci.yml:59-69]()

---

## Cross-Platform Shell Configuration

SOURSOP's CI workflow uses a specialized shell configuration strategy to ensure consistent behavior across Linux, macOS, and Windows platforms.

### The Login Shell Pattern

All steps after Conda initialization use the shell specification:

```yaml
shell: bash -l {0}
```

This configuration solves a cross-platform compatibility issue:

1. **Problem**: `conda init bash` updates different files on different operating systems:
   - Linux: Updates `~/.bashrc`
   - macOS: Updates `~/.bash_profile`
   - Windows: Updates appropriate configuration for bash emulation

2. **Solution**: The `-l` flag forces a **login shell**, which automatically sources the appropriate Conda initialization script regardless of OS.

3. **Benefit**: Eliminates need for OS-specific sourcing logic in the workflow.

**Sources**: [.github/workflows/soursop-ci.yml:33-39](), [.github/workflows/soursop-ci.yml:41](), [.github/workflows/soursop-ci.yml:47](), [.github/workflows/soursop-ci.yml:54]()

### Why Bash on All Platforms

The workflow uses `bash` even on Windows for consistency:

- **Windows compatibility**: Modern GitHub Actions Windows runners support bash
- **Unified commands**: Same commands work on all platforms
- **Conda compatibility**: Conda works reliably with bash across platforms

**Sources**: [.github/workflows/soursop-ci.yml:33-39]()

---

## Code Coverage Requirements

### Coverage Metrics

The CI system tracks coverage for the entire `soursop` package and enforces quality standards through Codecov integration.

### Coverage Reporting Configuration

```mermaid
graph LR
    PytestCov["pytest --cov=soursop<br/>--cov-report=xml"]
    
    XMLFile["coverage.xml<br/>(generated in repo root)"]
    
    Codecov["codecov/codecov-action@v3<br/>fail_ci_if_error: true"]
    
    CodecovAPI["Codecov Service<br/>(codecov.io)"]
    
    BuildStatus["CI Build Status<br/>(pass/fail)"]
    
    PytestCov --> XMLFile
    XMLFile --> Codecov
    Codecov --> CodecovAPI
    Codecov -->|"Upload failure"| BuildStatus
    CodecovAPI -->|"Coverage analysis"| BuildStatus
```

**Coverage Data Flow from Test Execution to Build Status**

### Failure Conditions

The CI build fails under these conditions:

1. **Test failures**: Any test in the suite fails (standard pytest behavior)
2. **Coverage upload failure**: The `fail_ci_if_error: true` setting causes build failure if coverage cannot be uploaded to Codecov
3. **Coverage degradation**: Codecov may be configured (via Codecov dashboard settings, not shown in workflow) to fail builds when coverage decreases

### Coverage Report Artifacts

The workflow generates multiple coverage artifacts:

| Artifact | Location | Purpose |
|----------|----------|---------|
| `coverage.xml` | Repository root | Machine-readable coverage data for Codecov |
| `codecov_report.txt` | Repository root | Human-readable upload report |

**Sources**: [.github/workflows/soursop-ci.yml:53-69]()

---

## Testing Thread Control

SOURSOP includes functionality for controlling NumPy threading behavior, which is tested as part of the CI pipeline. The `ssutils.set_numpy_threads()` function manages BLAS library threading across different platforms.

### Thread Control Testing

The test suite validates thread control on all platforms:

```python
# From test_ssutils.py
def test_set_numpy_threads():
    num_threads = 2
    set_threads, blas_library = ssutils.set_numpy_threads(num_threads)
    assert blas_library != 'unknown'
    assert set_threads == num_threads
```

This test ensures:
- Thread control works on the CI platform
- The BLAS library is properly detected (MKL or OpenBLAS)
- Thread limits can be successfully set

### Platform-Specific Behavior

The `set_numpy_threads()` function handles platform differences:

| Platform | BLAS Library | Detection Method |
|----------|-------------|------------------|
| Windows | MKL | Uses `mkl` package directly |
| Linux | OpenBLAS or MKL | Searches virtual environment for library files |
| macOS | OpenBLAS or MKL | Searches virtual environment (`.dylib` or `.so`) |

**Sources**: [soursop/ssutils.py:137-151](), [soursop/tests/test_ssutils.py:15-19]()

---

## CI Workflow Architecture Summary

```mermaid
graph TB
    subgraph GitHubEvents["GitHub Events"]
        Push["git push<br/>(any branch)"]
        PR["Pull Request<br/>(open/update)"]
    end
    
    subgraph WorkflowFile["Workflow Configuration"]
        YAML["soursop-ci.yml<br/>matrix: 3 OS × 3 Python"]
    end
    
    subgraph Runners["GitHub Actions Runners"]
        RunnerU["ubuntu-latest runner"]
        RunnerM["macos-latest runner"]
        RunnerW["windows-latest runner"]
    end
    
    subgraph Execution["Execution Stages"]
        Setup["1. checkout@v3<br/>2. setup-miniconda@v2<br/>3. conda init bash<br/>4. conda install deps"]
        Test["5. pytest --capture=sys"]
        Cov["6. pytest --cov=soursop<br/>--cov-report=xml -n auto"]
        Upload["7. codecov-action@v3<br/>fail_ci_if_error: true"]
    end
    
    subgraph External["External Services"]
        CodecovSvc["Codecov Service<br/>(codecov.io)<br/>Coverage analysis"]
    end
    
    subgraph Results["Build Results"]
        Pass["✓ Build Passes<br/>All tests pass<br/>Coverage uploaded"]
        Fail["✗ Build Fails<br/>Test failures or<br/>Coverage upload error"]
    end
    
    Push --> YAML
    PR --> YAML
    
    YAML --> RunnerU
    YAML --> RunnerM
    YAML --> RunnerW
    
    RunnerU --> Setup
    RunnerM --> Setup
    RunnerW --> Setup
    
    Setup --> Test
    Test --> Cov
    Cov --> Upload
    
    Upload --> CodecovSvc
    Upload --> Pass
    Upload --> Fail
    
    CodecovSvc --> Pass
    CodecovSvc --> Fail
```

**Complete CI Workflow Architecture: From Trigger to Result**

**Sources**: [.github/workflows/soursop-ci.yml:1-70]()

---