import argparse
import os
from pathlib import Path

import tempfile
import shutil

import nox

try:
    import tomllib
except ModuleNotFoundError or ImportError:
    import tomli as tomllib


def get_python_version():
    """Get the python version from .python-version file unless it is github actions."""
    if "GITHUB_ACTIONS" in os.environ:
        return Path(os.environ["Python_ROOT_DIR"]).parent.name
    return Path(".python-version").read_text().strip()


def get_package_name():
    """Get the package name from pyproject.toml."""
    pyproject_path = Path("pyproject.toml")
    pyproject_data = tomllib.loads(pyproject_path.read_text())
    return pyproject_data["project"]["name"].replace("-", "_")


PYTHON_VERSION = get_python_version()
PACKAGE_NAME = get_package_name()
PACKAGE_DIR_PATH = Path(PACKAGE_NAME).resolve()
SINGULARITY_DIR: Path | None = None
INTERPRO_VERSION = "5.77-108.0"

nox.options.reuse_existing_virtualenvs = True
nox.options.default_venv_backend = "uv"


def uv(session, *args):
    """Run uv commands inside a nox session."""
    session.run("uv", *args, "--active", external=True)


@nox.session(python=PYTHON_VERSION)
def install(session):
    """
    Equivalent of:
      make install
    """
    uv(session, "sync", "--group", "dev")
    # git clone the main of the domain-annotation-pipeline into FolderSeekStrucAnnoFlow/external/
    external_dir = PACKAGE_DIR_PATH.joinpath("external")
    if not external_dir.exists():
        # make the external directory if it doesn't exist
        external_dir.mkdir()
    session.run(
        "git",
        "clone",
        "https://github.com/UCLOrengoGroup/domain-annotation-pipeline.git",
        str(external_dir.joinpath("domain-annotation-pipeline")),
        external=True,
    )


#
# Formatting / chores
#


@nox.session(python=PYTHON_VERSION)
def ruff_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "check", ".", "--fix")


@nox.session(python=PYTHON_VERSION)
def black_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "format", ".")


@nox.session(python=PYTHON_VERSION)
def tomlsort_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run(
        "tombi",
        "format",
        # Only format toml files that are not in .venv directories
        *[str(p) for p in Path(".").rglob("*.toml") if ".venv" not in p.parts and ".nox" not in p.parts],
        external=True,
    )


@nox.session(python=PYTHON_VERSION)
def isort_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run("isort", PACKAGE_NAME, "tests")


@nox.session(python=PYTHON_VERSION)
def chores(session):
    """
    Equivalent of:
      make chores

    This runs all formatting and linting fixes but also mypy checks.
    """
    session.notify("isort_fixes")
    session.notify("ruff_fixes")
    session.notify("black_fixes")
    session.notify("tomlsort_fixes")
    session.notify("mypy_check")
    session.notify("nextflow_check")


@nox.session(python=PYTHON_VERSION)
def github_actions_chores(session):
    """
    Same as chores but with more verbose pytest output for better debugging in CI.

    and no cloud nextflow check as this repo relies on non-downloaded modules.

    """

    session.notify("isort_fixes")
    session.notify("ruff_fixes")
    session.notify("black_fixes")
    session.notify("tomlsort_fixes")
    session.notify("mypy_check")


#
# Checks / linting
#


@nox.session(python=PYTHON_VERSION)
def ruff_check(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "check")


@nox.session(python=PYTHON_VERSION)
def black_check(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "format", ".", "--check")


@nox.session(python=PYTHON_VERSION)
def mypy_check(session):
    uv(session, "sync", "--group", "dev")
    session.run("mypy", PACKAGE_NAME)


@nox.session(python=PYTHON_VERSION)
def tomlsort_check(session):
    uv(session, "sync", "--group", "dev")
    tomls = [str(p) for p in Path(".").rglob("*.toml") if ".venv" not in p.parts and ".nox" not in p.parts]
    session.run("tombi", "lint", *tomls, external=True)
    session.run("tombi", "format", *tomls, "--check", external=True)


@nox.session(python=PYTHON_VERSION)
def isort_check(session):
    uv(session, "sync", "--group", "dev")
    session.run("isort", PACKAGE_NAME, "tests", "--check-only")


@nox.session(python=PYTHON_VERSION)
def nextflow_check(session):
    """
    Lint Nextflow files using built in Nextflow linting.
    """
    session.run(
        "nextflow",
        "lint",
        *[str(p) for p in PACKAGE_DIR_PATH.rglob("*.nf")]
        + [str(p) for p in PACKAGE_DIR_PATH.rglob("*.config")],
        "-tabs",
        "-sort-declarations",
        "-harshil-alignment",
        external=True,
    )  

#
# Testing
#


@nox.session(python=PYTHON_VERSION)
def tests(session):
    """
    Equivalent of:
      make tests
    """

    session.notify("ruff_check")
    session.notify("black_check")
    session.notify("mypy_check")
    session.notify("tomlsort_check")
    # session.notify("isort_check")
    session.notify("nextflow_check")
    session.notify("pytest")


@nox.session(python=PYTHON_VERSION)
def github_actions_tests(session):
    """
    Same as tests but with more verbose pytest output for better debugging in CI.

    and no cloud nextflow check as this repo relies on non-downloaded modules.

    """

    session.notify("ruff_check")
    session.notify("black_check")
    session.notify("mypy_check")
    session.notify("tomlsort_check")
    # session.notify("isort_check")
    session.notify("pytest_loud")


@nox.session(python=PYTHON_VERSION)
def precommit(session):
    """
    Equivalent of:
      make precommit
    """

    session.notify("ruff_check")
    session.notify("black_check")
    session.notify("mypy_check")
    session.notify("tomlsort_check")
    # session.notify("isort_check")


@nox.session(python=PYTHON_VERSION)
def pytest(session):
    uv(session, "sync", "--group", "dev")
    session.run(
        "pytest",
        f"--cov={PACKAGE_NAME}",
        "--cov-report=term-missing",
        "tests",
    )


@nox.session(python=PYTHON_VERSION)
def pytest_loud(session):
    uv(session, "sync", "--group", "dev")
    session.run(
        "pytest",
        "--log-cli-level=DEBUG",
        "-log_cli=true",
        f"--cov=./{PACKAGE_NAME}",
        "--cov-report=term-missing",
        "tests",
    )


#
# Packaging
#


def parse_args(session: nox.Session) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Docker and Apptainer containers")
    parser.add_argument(
        "--output",
        help=f"Directory to save .sif files", required=True, type=Path
    )
    return parser.parse_args(session.posargs)
 

@nox.session(name="build_apptainer", python=False)
def build_apptainer(session: nox.Session) -> None:
    """Build all Docker containers locally and convert to Apptainer .sif files.

    Usage:
        nox -s build_apptainer -- --output /path/to/dir
    """
    args = parse_args(session)
    sif_dir: Path = args.output

    base_dir = Path(__file__).parent.resolve()
    services = [d for d in base_dir.joinpath("docker").iterdir() if d.is_dir()]

    sif_dir.mkdir(parents=True, exist_ok=True)
    session.log(f"Apptainer .sif files will be saved to: {sif_dir}")

    docker_failed = []
    apptainer_failed = []

    for service in services:
        sif_name = f"GAQA_{service.name}.sif"
        sif_path = sif_dir / sif_name

        if sif_path.exists():
            session.log(f"[{service.name}] Exists. Do you want to overwrite? [y/N]")
            answer = input().strip().lower()
            if answer != "y":
                session.log(f"[{service.name}] Skipping...")
                continue

        dockerfile = service / "Dockerfile"
        tag = f"gaqa_{service.name.lower()}:latest"

        session.log("")
        session.log(f"[{service.name}] Building Docker image...")

        try:
            cmd = [
                "docker",
                "build",
                "-t",
                tag,
                "-f",
                str(dockerfile),
            ]
            if service.name == "interproscan":
                cmd += ["--build-arg", f"VERSION={INTERPRO_VERSION}"]
       
            cmd += [str(service)]

            if service.name == "python":
                # temporarily copy the pyproject.toml to the docker/python directory for the build context
                temp_pyproject_path = service / "pyproject.toml"
                temp_pyproject_path.write_text((base_dir / "pyproject.toml").read_text())
                session.log(f"[{service.name}] Copied pyproject.toml to {temp_pyproject_path} for Docker build context")
            

            session.run(
                *cmd,
                external=True,
            )

            if service.name == "python":
                temp_pyproject_path.unlink()
                session.log(f"[{service.name}] Removed temporary pyproject.toml from {temp_pyproject_path}")
       
            session.log(f"[{service.name}] ✓ Docker image built as {tag}")
        except Exception:
            session.log(f"[{service.name}] ✗ Docker build failed — skipping Apptainer conversion")
            docker_failed.append(service)
            continue

        session.log(f"[{service.name}] Converting to Apptainer .sif -> {sif_path}")
        try:
            session.run(
                "apptainer",
                "build",
                "--force",
                str(sif_path),
                f"docker-daemon:{tag}",
                external=True,
            )
            session.log(f"[{service.name}] ✓ Apptainer .sif saved to {sif_path}")
        except Exception:
            session.log(f"[{service.name}] ✗ Apptainer conversion failed")
            apptainer_failed.append(service)

    session.log("")
    session.log("=" * 50)
    session.log("Build Summary")
    session.log("=" * 50)

    for service in services:
        if service in docker_failed:
            session.log(f"  ✗ {service.name} (docker build failed)")
        elif service in apptainer_failed:
            session.log(f"  ~ {service.name} (docker ok, apptainer failed)")
        else:
            session.log(f"  ✓ {service.name} ({sif_dir / f'GAQA_{service.name}.sif'})")

    all_failed = docker_failed + apptainer_failed
    if all_failed:
        session.error(f"\n{len(all_failed)} container(s) had errors: {', '.join([s.name for s in all_failed])}")


@nox.session(python=PYTHON_VERSION)
def build(session):
    """
    Equivalent of:
      make build
    """
    uv(session, "sync", "--group", "dev")
    session.run("python", "-m", "build")


#
# Nextflow
#


@nox.session()
def nextflow(session):
    """
    Run Nextflow pipeline with a specified profile (default: local).
    Usage:
      nox -s nextflow -- [profile] [nextflow args...]
      e.g. nox -s nextflow -- slurm --input_dir /path/to/inputs
    """
    args = list(session.posargs)

    session.run(
        "nextflow",
        "run",
        PACKAGE_DIR_PATH.joinpath("main.nf"),
        *args,
        external=True,
        env={
            "NXF_LOG_FILE": str(PACKAGE_DIR_PATH.joinpath("work", ".nextflow.log")),
            **os.environ,
        },
    )


@nox.session()
def nextflow_clean(session):
    """
    Clean all Nextflow work directories for this project.
    """

    session.run(
        "nextflow",
        "clean",
        "-f",
        external=True,
    )


@nox.session()
def nextflow_obliterate(session):
    """
    Obliterate all Nextflow work directories for this project.
    WARNING: This will delete all work directories and cannot be undone.
    """
    session.chdir(PACKAGE_DIR_PATH)

    shutil.rmtree(PACKAGE_DIR_PATH.joinpath("work"), ignore_errors=True)
