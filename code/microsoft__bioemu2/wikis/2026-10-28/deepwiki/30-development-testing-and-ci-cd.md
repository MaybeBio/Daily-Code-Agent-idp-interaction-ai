# Development, Testing and CI/CD

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.github/workflows/ci.yml](.github/workflows/ci.yml)
- [CONTRIBUTING.md](CONTRIBUTING.md)
- [DEVELOPER.md](DEVELOPER.md)

</details>



This section provides a high-level overview of the contributor workflow, test suite organization, code quality tooling, CI/CD pipelines, and governance policies for BioEmu2. The codebase utilizes modern Python packaging via `uv` [DEVELOPER.md:3-8](), automated linting and formatting via pre-commit hooks [.github/workflows/ci.yml:35-37](), and multi-version matrix testing through GitHub Actions [.github/workflows/ci.yml:17-23]().

For detailed information, refer to the child pages:
- [Test Suite Layout](#8.1)
- [Code Quality Tooling](#8.2)
- [CI Workflows and Release Publishing](#8.3)
- [Contributing and Governance](#8.4)

---

## 8.1 Test Suite Layout

The test suite in BioEmu2 is structured under the `tests/` directory, mirroring the modular layout of the core library under `src/bioemu2/`. It encompasses unit tests, integration tests, and model regression checks for components such as score models, context parallelism, denoising/steering routines, structure generation, and MSA features. 

For detailed layouts of test directories, custom fixtures, and regression data assets, see [Test Suite Layout](#8.1).

Sources: `[DEVELOPER.md:25-29]()`

---

## 8.2 Code Quality Tooling

BioEmu2 enforces strict code formatting, linting, and type checking using tools managed through `uv` and `pre-commit` [.github/workflows/ci.yml:35-37](). Developers can initialize their environment with CUDA support and install pre-commit hooks to ensure compliance before pushing changes:

```bash
uv sync --extra cuda
uv run pre-commit install
```

```mermaid
title "Developer Workflow and Tooling"
graph TD
    A["DeveloperWorkspace"] --> B["uv sync --extra cuda" > DEVELOPER.md:7]
    B --> C["pre-commit run --all-files" > .github/workflows/ci.yml:37]
    C --> D["pytest tests" > DEVELOPER.md:28]
    D --> E["LocalVerificationPassed"]
```
*Sources: [.github/workflows/ci.yml:35-37](), [DEVELOPER.md:1-38]()*

For configuration details regarding `ruff`, `mypy`, and documentation generation, see [Code Quality Tooling](#8.2).

---

## 8.3 CI Workflows and Release Publishing

Continuous integration is handled via GitHub Actions workflows (such as `[.github/workflows/ci.yml:1-39]()`). The CI pipeline validates code across multiple Python versions using matrix strategies, synchronizes dependencies with lockfiles via `uv`, runs style checks, and executes the complete `pytest` suite [.github/workflows/ci.yml:17-39]().

```mermaid
title "CI Pipeline Architecture"
graph TD
    A["GitHubPushOrPR" > .github/workflows/ci.yml:3-8] --> B["MatrixJobPython312to314" > .github/workflows/ci.yml:20-23]
    B --> C["InstallUvAndSync" > .github/workflows/ci.yml:31-34]
    C --> D["PreCommitCheck" > .github/workflows/ci.yml:35-37]
    C --> E["PytestSuiteExecution" > .github/workflows/ci.yml:38-39]
```
*Sources: [.github/workflows/ci.yml:1-39]()*

For details on automated PyPI publishing, trusted publishers, and dependency update bots, see [CI Workflows and Release Publishing](#8.3).

---

## 8.4 Contributing and Governance

Contributions to BioEmu2 are managed through pull requests on GitHub [CONTRIBUTING.md:3-4](). Code contributions require agreement to the Microsoft Open Source Contributor License Agreement (CLA) to grant necessary usage permissions [CONTRIBUTING.md:6-11](). 

For full details on the pull request workflow, code of conduct, and security vulnerability reporting, see [Contributing and Governance](#8.4).

Sources: `[CONTRIBUTING.md:1-11]()`

---