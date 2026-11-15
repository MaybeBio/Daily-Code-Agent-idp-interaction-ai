# Build and Packaging

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.readthedocs.yml](.readthedocs.yml)
- [MANIFEST.in](MANIFEST.in)
- [pyproject.toml](pyproject.toml)
- [readthedocs.yml](readthedocs.yml)
- [setup.cfg](setup.cfg)
- [setup.py](setup.py)

</details>



This document covers the build system, packaging configuration, and distribution mechanisms for FINCHES. It explains how the project is structured for development, built into distributable packages, and deployed to users.

For information about testing and development workflows, see page [7.1](#7.1). For details about the core computational modules being packaged, see page [4.3](#4.3).

## Build System Architecture

FINCHES uses a hybrid build system that combines modern Python packaging standards with legacy components required for Cython compilation. The system is configured to handle both pure Python modules and compiled Cython extensions.

### Build System Components

```mermaid
graph TB
    subgraph "Modern Build System"
        PT["pyproject.toml<br/>PEP 518 configuration"]
        SB["setuptools build backend"]
        VI["versioningit<br/>dynamic versioning"]
    end
    
    subgraph "Legacy Build Components"
        SP["setup.py<br/>Cython compilation"]
        CE["Cython Extensions<br/>matrix_manipulation.pyx"]
        NI["NumPy includes"]
    end
    
    subgraph "Package Configuration"
        MI["MANIFEST.in<br/>file inclusion rules"]
        SC["setup.cfg<br/>tool configuration"]
        PD["Package data<br/>pickle files, data"]
    end
    
    subgraph "Build Output"
        WH["Python wheel"]
        SD["Source distribution"]
        CY["Compiled extensions<br/>.so/.pyd files"]
    end
    
    PT --> SB
    SB --> VI
    SP --> CE
    CE --> NI
    MI --> SB
    SC --> SB
    
    SB --> WH
    SB --> SD
    SP --> CY
    
    PT -.-> SP
```

**Build System Flow**
The build process starts with `pyproject.toml` defining build requirements and metadata, while `setup.py` handles the specialized Cython compilation step.

Sources: [pyproject.toml:1-8](), [setup.py:6-30]()

## Modern Python Packaging Configuration

The primary build configuration uses `pyproject.toml` following PEP 518 standards for declarative project metadata and build system specification.

### Project Metadata and Dependencies

| Configuration Section | Purpose | Key Elements |
|----------------------|---------|--------------|
| `[build-system]` | Build tool requirements | setuptools, wheel, versioningit, cython, numpy |
| `[project]` | Package metadata | name, description, authors, license, dependencies |
| `[project.optional-dependencies]` | Optional extras | test dependencies (pytest) |
| `[tool.setuptools]` | Package discovery | include patterns, package data |

The build system requires Cython and NumPy at build time due to the compiled extensions:

```python
requires = ["setuptools>=61.0", "wheel", "versioningit~=2.0", "cython", "numpy"]
```

Runtime dependencies include scientific computing libraries and domain-specific tools:

```python
dependencies = ["numpy", "afrc>=0.3.4", "scipy", "soursop>=0.2.4", "pandas", "metapredict", "ipython"]
```

Sources: [pyproject.toml:1-31](), [pyproject.toml:38-42]()

### Package Discovery and Data Inclusion

The package discovery configuration automatically finds all subpackages within the `finches` namespace:

```python
[tool.setuptools.packages.find]
include = ["finches", "finches.*"]
```

Package data inclusion is configured to include type information:

```python
[tool.setuptools.package-data]
finches = ["py.typed"]
```

Sources: [pyproject.toml:48-58]()

## Cython Extension Compilation

FINCHES includes performance-critical code written in Cython that must be compiled during the build process. The `setup.py` file handles this compilation step.

### Extension Module Configuration

```mermaid
graph LR
    subgraph "Source Files"
        PYX["matrix_manipulation.pyx<br/>Cython source"]
        NPY["numpy headers<br/>C API"]
    end
    
    subgraph "Build Process"
        EXT["Extension definition<br/>setup.py"]
        CYT["cythonize()<br/>compilation"]
        GCC["C compiler<br/>platform specific"]
    end
    
    subgraph "Output"
        SO["matrix_manipulation.so<br/>Linux/macOS"]
        PYD["matrix_manipulation.pyd<br/>Windows"]
    end
    
    PYX --> EXT
    NPY --> EXT
    EXT --> CYT
    CYT --> GCC
    GCC --> SO
    GCC --> PYD
```

**Cython Extension Build Flow**
The Cython source file is compiled into a platform-specific binary extension that provides optimized matrix operations.

The extension is defined in `setup.py` with NumPy includes for array operations:

```python
extensions = [
    Extension(
        "finches.utils.matrix_manipulation",
        [cython_file],
        include_dirs=[numpy.get_include()], 
    )]
```

Sources: [setup.py:12-23](), [setup.py:25-30]()

## Version Management

FINCHES uses `versioningit` for automatic version generation from Git tags and commits. This provides consistent versioning across development and release cycles.

### Version Configuration

The version management system is configured in `pyproject.toml`:

```python
[tool.versioningit]
default-version = "1+unknown"

[tool.versioningit.format]
distance = "{base_version}+{distance}.{vcs}{rev}"
dirty = "{base_version}+{distance}.{vcs}{rev}.dirty"

[tool.versioningit.vcs]
method = "git"
match = ["*"]
default-tag = "1.0.0"
```

The version information is written to `finches/_version.py` at build time, making it available to the package at runtime.

Sources: [pyproject.toml:60-76]()

## Package Data and Distribution

The build system carefully controls which files are included in the distribution through `MANIFEST.in` and package data configuration.

### File Inclusion Rules

```mermaid
graph TB
    subgraph "Source Repository"
        SRC["Source code<br/>.py files"]
        CY_SRC["Cython source<br/>.pyx files"]
        DATA["Package data<br/>.pickle files"]
        DOCS["Documentation<br/>LICENSE, README"]
    end
    
    subgraph "MANIFEST.in Rules"
        INC["include rules<br/>recursive-include"]
        EXC["exclude rules<br/>global-exclude"]
        GRAFT["graft finches<br/>include all"]
    end
    
    subgraph "Distribution Package"
        WHEEL["Binary wheel<br/>.whl file"]
        TARBALL["Source distribution<br/>.tar.gz file"]
    end
    
    SRC --> INC
    CY_SRC --> INC
    DATA --> INC
    DOCS --> INC
    
    INC --> EXC
    EXC --> GRAFT
    
    GRAFT --> WHEEL
    GRAFT --> TARBALL
```

**Package Data Inclusion Flow**
Files are selectively included based on `MANIFEST.in` rules, ensuring only necessary files are distributed.

Key inclusion rules:
- All `.pyx` files are included recursively
- All `.pickle` files containing model parameters are included  
- All data under `finches/data/` is included
- Temporary files (`.pyc`, `__pycache__`, `.so`) are excluded

Sources: [MANIFEST.in:6-16]()

## Documentation Building

FINCHES documentation is built automatically using ReadTheDocs with Sphinx. The configuration supports the project's scientific computing dependencies.

### ReadTheDocs Configuration

```mermaid
graph LR
    subgraph "ReadTheDocs Config"
        RTD[".readthedocs.yml<br/>build environment"]
        REQS["docs/requirements.txt<br/>doc dependencies"]
        CONF["docs/conf.py<br/>Sphinx config"]
    end
    
    subgraph "Build Environment"
        UBUNTU["ubuntu-22.04<br/>build OS"]
        PYTHON["Python 3.9/3.11<br/>interpreter"]
        SPHINX["Sphinx<br/>doc generator"]
    end
    
    subgraph "Documentation Output"
        HTML["HTML documentation<br/>readthedocs.io"]
        API["API reference<br/>auto-generated"]
    end
    
    RTD --> UBUNTU
    RTD --> PYTHON
    REQS --> SPHINX
    CONF --> SPHINX
    
    SPHINX --> HTML
    SPHINX --> API
```

**Documentation Build Pipeline**
ReadTheDocs automatically builds documentation from the repository, handling scientific dependencies and Sphinx configuration.

The build environment is configured for Ubuntu 22.04 with Python 3.9, and documentation dependencies are specified separately from the main package requirements.

Sources: [.readthedocs.yml:1-18](), [readthedocs.yml:1-15]()

## Development Configuration

Additional build and development tool configuration is provided through `setup.cfg` for consistency across development environments.

### Tool Configuration

The configuration includes settings for:
- **Coverage reporting**: Omits test files and generated version files
- **Code formatting**: YAPF configuration with 119 character line limit
- **Linting**: Flake8 configuration matching the formatting settings
- **Test execution**: Pytest integration with setup.py

Key configuration example:
```python
[yapf]
COLUMN_LIMIT = 119
INDENT_WIDTH = 4
USE_TABS = False
```

Sources: [setup.cfg:1-24]()

## Build Commands and Workflow

The build system supports standard Python packaging commands for development and distribution:

| Command | Purpose | Output |
|---------|---------|---------|
| `python -m build` | Create distribution packages | wheel and source distribution |
| `pip install -e .` | Development installation | Editable install with Cython compilation |
| `python setup.py build_ext --inplace` | Compile Cython only | In-place compiled extensions |
| `pip install .` | Standard installation | Compiled package installation |

The development workflow typically involves editable installation for active development, with the build system automatically handling Cython compilation when source files change.

Sources: [pyproject.toml:1-76](), [setup.py:1-32]()