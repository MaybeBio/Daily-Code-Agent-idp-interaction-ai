# Development

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.gitattributes](.gitattributes)
- [.github/workflows/python_format.yaml](.github/workflows/python_format.yaml)
- [.github/workflows/python_lint.yaml](.github/workflows/python_lint.yaml)
- [.gitignore](.gitignore)
- [requirements.txt](requirements.txt)
- [setup.py](setup.py)

</details>



This page provides essential information for developers who want to contribute to the PaCS-Toolkit codebase. It covers the development environment setup, code organization, contribution workflow, and coding standards. For information about using the toolkit, see [Overview](#1) and [Quick Start Guide](#3).

## Development Environment Setup

Setting up a development environment for PaCS-Toolkit involves a few straightforward steps:

1. **Clone the repository**:
   ```bash
   git clone https://github.com/Kitaolab/PaCS-Toolkit.git
   cd PaCS-Toolkit
   ```

2. **Install development dependencies**:
   ```bash
   pip install -r requirements.txt
   pip install -r requirements_dev.txt
   ```

3. **Install the package in development mode**:
   ```bash
   pip install -e .
   ```

### Dependencies Overview

PaCS-Toolkit has a modular dependency structure that allows users to install only what they need.

```mermaid
graph TD
    subgraph "Core Dependencies"
        core["Core Requirements"]
        core --> numpy["numpy<2.0.0"]
        core --> tomli["tomli>=2.0.1"]
    end
    
    subgraph "Optional Extensions"
        mdtraj["MDTraj Extension"]
        mdtraj --> mdtrajLib["mdtraj==1.9.9"]
        mdtraj --> scipy["scipy>=1.7.3"]
        mdtraj --> sklearn["scikit-learn>=1.0.2"]
        
        msm["MSM Extension"]
        msm --> mdtrajLib
        msm --> scipy
        msm --> sklearn
        msm --> matplotlib["matplotlib>=3.5.3"]
        msm --> pandas["pandas>=1.1.5"]
        msm --> deeptime["deeptime==0.4.4"]
        msm --> ipykernel["ipykernel>=6.16.2"]
    end
    
    dev["Development Requirements"]
    dev --> linting["Linting Tools"]
    dev --> testing["Testing Tools"]
```

Sources: [setup.py:32-62](), [requirements.txt:1-10]()

## Project Structure

The PaCS-Toolkit follows a modular architecture that separates different components of the system.

```mermaid
graph TD
    root["PaCS-Toolkit Root"]
    root --> pacs["pacs/ (Main Package)"]
    root --> MSM["MSM/ (MSM Analysis)"]
    root --> docs["docs/ (Documentation)"]
    root --> examples["examples/ (Example Files)"]
    root --> tests["tests/ (Test Files)"]
    
    pacs --> core["core/ (Core Framework)"]
    pacs --> analyzers["analyzers/ (Analyzer Modules)"]
    pacs --> utils["utils/ (Utility Functions)"]
    pacs --> cli["__main__.py (CLI Entry Point)"]
    pacs --> version["_version.py (Version Info)"]
    
    core --> cycle["cycle.py (Simulation Cycle)"]
    core --> mdsettings["mdsettings.py (Settings Manager)"]
    core --> simulator["simulator/ (MD Engine Interfaces)"]
    
    analyzers --> super["superanalyzer.py (Base Class)"]
    analyzers --> rmsd["rmsd.py (RMSD Analyzer)"]
    analyzers --> target["target.py (Target Analyzer)"]
    analyzers --> other["Other Analyzers..."]
    
    utils --> fileutils["fileutils.py (File Management)"]
    utils --> genfeature["genfeature.py (Feature Generation)"]
    utils --> other_utils["Other Utilities..."]
```

Sources: The project structure is inferred from the system architecture overview and patterns typical of Python projects

## Contribution Workflow

### Branch Management and Pull Requests

```mermaid
flowchart TD
    fork["Fork Repository"]
    clone["Clone Fork"]
    branch["Create Feature Branch"]
    develop["Make Changes"]
    test["Run Tests"]
    lint["Run Linting"]
    push["Push Changes"]
    pr["Create Pull Request"]
    review["Code Review"]
    merge["Merge to Main"]
    
    fork --> clone
    clone --> branch
    branch --> develop
    develop --> test
    test --> lint
    lint -- "Pass" --> push
    lint -- "Fail" --> develop
    push --> pr
    pr --> review
    review -- "Approved" --> merge
    review -- "Changes Requested" --> develop
```

Sources: Common GitHub workflow practices

### Automatic Code Formatting and Linting

PaCS-Toolkit uses GitHub Actions to automatically enforce code quality standards:

1. **Automatic Formatting**: The `python_format.yaml` workflow runs on the main branch to format code.
2. **Linting Checks**: The `python_lint.yaml` workflow runs on push and pull requests to verify code quality.

```mermaid
flowchart LR
    push["Push to GitHub"]
    
    subgraph "GitHub Actions"
        format["Format Check (pysen format)"]
        lint["Lint Check (pysen lint)"]
    end
    
    push --> format
    push --> lint
    
    format -- "Auto-commit" --> fixed["Formatted Code"]
    lint -- "Pass/Fail" --> pr["Pull Request Status"]
```

Sources: [.github/workflows/python_format.yaml:1-36](), [.github/workflows/python_lint.yaml:1-34]()

## Coding Standards

### Code Style

PaCS-Toolkit follows standard Python coding practices:

- **PEP 8** for code style
- **Type hints** for function and method signatures
- **Docstrings** for all public functions, classes, and methods

### Tools Used for Code Quality

| Tool | Purpose | Configuration |
|------|---------|--------------|
| pysen | Orchestrates code quality tools | pyproject.toml |
| flake8 | Style guide enforcement | Via pysen |
| mypy | Static type checking | Via pysen |
| black | Code formatting | Via pysen |
| isort | Import sorting | Via pysen |

Sources: [.github/workflows/python_format.yaml:1-36](), [.github/workflows/python_lint.yaml:1-34]()

## Extending PaCS-Toolkit

### Adding a New Analyzer

One common extension point is adding a new analyzer. Analyzers inherit from the `SuperAnalyzer` class and implement specific analysis methods.

```mermaid
classDiagram
    class SuperAnalyzer {
        +__init__(mdsetting)
        +analysis(outputs)
        +ranking(results)
        +check_threshold(results)
    }
    
    class RMSDAnalyzer {
        +__init__(mdsetting)
        +analysis(outputs)
        +ranking(results)
        +check_threshold(results)
    }
    
    class YourCustomAnalyzer {
        +__init__(mdsetting)
        +analysis(outputs)
        +ranking(results)
        +check_threshold(results)
    }
    
    SuperAnalyzer <|-- RMSDAnalyzer
    SuperAnalyzer <|-- YourCustomAnalyzer
```

To create a new analyzer:

1. Create a new Python file in the `pacs/analyzers/` directory
2. Define a class that inherits from `SuperAnalyzer`
3. Implement the required methods: `analysis()`, `ranking()`, and `check_threshold()`
4. Register your analyzer in the appropriate configuration systems

Sources: Based on system architecture diagram showing analyzer inheritance

### Adding a New MD Engine Interface

Similarly, you can extend PaCS-Toolkit to support a new MD engine:

```mermaid
classDiagram
    class BaseSimulator {
        +__init__(mdsetting)
        +simulate(cycle, rank)
        +preprocess()
        +postprocess()
    }
    
    class GromacsSimulator {
        +__init__(mdsetting)
        +simulate(cycle, rank)
        +preprocess()
        +postprocess()
    }
    
    class AmberSimulator {
        +__init__(mdsetting)
        +simulate(cycle, rank)
        +preprocess()
        +postprocess()
    }
    
    class NewEngineSimulator {
        +__init__(mdsetting)
        +simulate(cycle, rank)
        +preprocess()
        +postprocess()
    }
    
    BaseSimulator <|-- GromacsSimulator
    BaseSimulator <|-- AmberSimulator
    BaseSimulator <|-- NewEngineSimulator
```

Sources: Based on system architecture showing simulator integration

## Package Management

PaCS-Toolkit uses setuptools for package management. The package has a modular structure with optional components.

```mermaid
flowchart TD
    setup["setup.py"]
    
    setup --> metadata["Package Metadata"]
    setup --> dependencies["Dependencies"]
    setup --> entry_points["Entry Points"]
    
    dependencies --> core["Core Dependencies"]
    dependencies --> extras["Extra Dependencies"]
    
    extras --> mdtraj["MDTraj Extension"]
    extras --> msm["MSM Extension"]
    extras --> all["All Extensions"]
    
    entry_points --> cli["pacs CLI Command"]
```

The package uses extras_require to define optional dependency sets:

| Extra | Purpose | Dependencies |
|-------|---------|-------------|
| mdtraj | Enables analyzers using MDTraj | mdtraj, scipy, scikit-learn |
| msm | Enables MSM analysis | mdtraj, matplotlib, scipy, pandas, deeptime, ipykernel, scikit-learn |
| all | All optional features | All dependencies from mdtraj and msm |

Sources: [setup.py:8-69]()

## Testing and Documentation

### Testing Guidelines

When contributing code to PaCS-Toolkit, it's important to:

1. Write tests for any new functionality
2. Ensure all existing tests pass before submitting a pull request
3. Test with different configurations (core, mdtraj, msm, all)

### Documentation Guidelines

Documentation in PaCS-Toolkit comes in several forms:

1. **Code comments**: Explain complex algorithms and non-obvious decisions
2. **Docstrings**: Document the purpose, parameters, and return values of functions and classes
3. **Wiki pages**: Provide high-level explanations of system components
4. **Examples**: Demonstrate how to use various features

When adding new code, ensure that you update or create appropriate documentation.

## Version Control Guidelines

### Git Best Practices

```mermaid
flowchart LR
    subgraph "File Categories"
        code["Python Code"]
        config["Configuration Files"]
        docs["Documentation"]
        tests["Test Files"]
        temp["Temporary Files"]
    end
    
    subgraph "Version Control"
        include["Include in Git"]
        exclude["Exclude from Git"]
    end
    
    code --> include
    config --> include
    docs --> include
    tests --> include
    temp --> exclude
```

The project uses `.gitignore` to exclude common temporary files and build artifacts.

Sources: [.gitattributes:1-2](), [.gitignore:1-168]()

## Release Process

The release process for PaCS-Toolkit involves:

1. Updating the version number in `pacs/_version.py`
2. Running final tests and quality checks
3. Creating a release tag in Git
4. Building and uploading the package to PyPI

Sources: [setup.py:6-10]()

---