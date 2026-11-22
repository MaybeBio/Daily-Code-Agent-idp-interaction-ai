"""OpenDNA from Python, for the notebook this work usually happens in.

    >>> import opendna
    >>> run = opendna.predict("disorder", {"ubiquitin": "MQIFVKTLTGK..."})
    >>> run.frame.head()
    >>> run.manifest["head"]["build"]

This drives the `opendna` command-line tool rather than binding to the Rust
library. That is a deliberate first step, not a shortcut around a harder one:
the command-line path is the one with the exit-code contract, the model
verification and the run manifest already tested around it, and shelling out
inherits all three. A compiled extension would mean a wheel per platform and a
second code path to keep honest, for the same four predictions.

What it adds over calling the binary yourself is that **the provenance comes
back with the numbers**. A `DataFrame` that has lost which build produced it is
how a prediction turns into a fact in a spreadsheet; every result here carries
its manifest and its limitations, and `Run.cite()` prints what a methods
section needs.

Nothing here estimates anything. A missing head raises, exactly as the
command-line tool exits 3, because the alternative is a plausible number from
no model.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "Run",
    "OpenDnaError",
    "ToolMissing",
    "HeadMissing",
    "InputRejected",
    "predict",
    "models",
    "selftest",
    "binary",
    "TASKS",
]

#: The four heads. `stability` additionally needs `mutations`.
TASKS = ("disorder", "burial", "solubility", "stability")

# The command-line tool's contract, which this mirrors rather than reinvents.
_EXIT_USAGE = 2
_EXIT_INSTALL = 3
_EXIT_INPUT = 4
_EXIT_CHECKS = 5


class OpenDnaError(RuntimeError):
    """Anything OpenDNA refused to do."""


class ToolMissing(OpenDnaError):
    """The `opendna` binary could not be found."""


class HeadMissing(OpenDnaError):
    """A trained head is not installed, or does not verify."""


class InputRejected(OpenDnaError):
    """The sequences or mutations were not acceptable."""


def binary() -> str:
    """Where the `opendna` executable is.

    `OPENDNA_BIN` wins, then whatever is on the path. Raises rather than
    guessing, because a wrong binary would produce numbers from the wrong
    build.
    """
    explicit = os.environ.get("OPENDNA_BIN")
    if explicit:
        if not Path(explicit).is_file():
            raise ToolMissing(f"OPENDNA_BIN points at {explicit}, which is not a file")
        return explicit
    found = shutil.which("opendna")
    if found:
        return found
    raise ToolMissing(
        "no `opendna` on PATH. Take a build from "
        "https://github.com/dyber-pqc/OpenDNA/releases, or set OPENDNA_BIN to one."
    )


@dataclass
class Run:
    """One prediction, with everything needed to say where it came from."""

    #: Which property was predicted.
    task: str
    #: The rows, as dictionaries. Keys are the command-line tool's columns.
    rows: list[dict[str, Any]]
    #: The run manifest: input and output digests, head build, device, licence.
    #: Empty only if the manifest was suppressed.
    manifest: dict[str, Any] = field(default_factory=dict)
    #: The model card's limitations, carried so they cannot be left behind.
    limitations: list[str] = field(default_factory=list)
    #: The head's own version string.
    head_version: str = ""
    #: The embedding model the head reads.
    embedding_model: str = ""
    #: What it ran on.
    device: str = ""

    @property
    def frame(self):
        """The rows as a pandas `DataFrame`.

        Raises if pandas is not installed, rather than being a hard dependency
        of a package whose job is to run a binary.
        """
        try:
            import pandas
        except ImportError as error:  # pragma: no cover - environment-dependent
            raise OpenDnaError(
                "pandas is not installed. `pip install pandas`, or use `.rows`."
            ) from error
        return pandas.DataFrame(self.rows)

    @property
    def build(self) -> str | None:
        """The head's build id, which identifies the exact weights."""
        head = self.manifest.get("head") or {}
        return head.get("build")

    def cite(self) -> str:
        """What a methods section has to say about this run."""
        parts = [f"{self.task} predicted with OpenDNA head version {self.head_version}"]
        if self.build:
            parts.append(f"build {self.build}")
        if self.embedding_model:
            parts.append(f"on frozen {self.embedding_model} embeddings")
        out = ", ".join(parts)
        if self.device:
            out += f" ({self.device})"
        if self.limitations:
            out += "\n\nLimitations stated by the model card:\n"
            out += "\n".join(f"  - {limitation}" for limitation in self.limitations)
        return out

    def __len__(self) -> int:
        return len(self.rows)


def _fasta(sequences: Mapping[str, str] | Iterable[tuple[str, str]]) -> str:
    items = sequences.items() if isinstance(sequences, Mapping) else sequences
    out = []
    for name, sequence in items:
        if not sequence:
            raise InputRejected(f"{name}: empty sequence")
        if "\n" in name or ">" in name:
            raise InputRejected(f"a sequence name cannot contain a newline or '>': {name!r}")
        out.append(f">{name}\n{sequence}\n")
    if not out:
        raise InputRejected("no sequences given")
    return "".join(out)


def _refusal(output: str) -> str:
    """The part of the tool's output that says what went wrong.

    It writes progress to standard error before it fails, so the first line is
    usually "N sequences from ..." and the refusal is further down. A caller
    reading `str(exception)` gets the first line, so the first line has to be
    the reason.
    """
    lines = output.strip().splitlines()
    for index, line in enumerate(lines):
        if line.startswith("Error:"):
            return "\n".join(lines[index:]).strip()
    return "\n".join(lines).strip()


def _run(args: Sequence[str], *, models: str | os.PathLike[str] | None = None) -> str:
    command = [binary()]
    if models is not None:
        command += ["--models", str(models)]
    command += list(args)
    finished = subprocess.run(command, capture_output=True, text=True, check=False)
    if finished.returncode == 0:
        return finished.stdout
    message = _refusal(finished.stderr or finished.stdout)
    if finished.returncode == _EXIT_INSTALL:
        raise HeadMissing(message)
    if finished.returncode in (_EXIT_INPUT, _EXIT_USAGE):
        raise InputRejected(message)
    if finished.returncode == _EXIT_CHECKS:
        raise OpenDnaError(f"this installation does not reproduce its recorded results:\n{message}")
    raise OpenDnaError(f"opendna exited {finished.returncode}:\n{message}")


def predict(
    task: str,
    sequences: Mapping[str, str] | Iterable[tuple[str, str]],
    *,
    mutations: Iterable[tuple[str, str]] | None = None,
    models: str | os.PathLike[str] | None = None,
    min_confidence: float | None = None,
    batch_tokens: int | None = None,
    device: str | None = None,
    offline: bool = False,
) -> Run:
    """Predict one property for some sequences.

    `sequences` maps a name to a sequence. For `stability`, `mutations` is an
    iterable of `(name, label)` pairs such as `("ubiquitin", "L8P")`, and the
    names must match.

    Raises `HeadMissing` if the head is not installed: nothing here falls back
    to an estimate.
    """
    if task not in TASKS:
        raise InputRejected(f"unknown task {task!r}; expected one of {', '.join(TASKS)}")
    if task == "stability" and not mutations:
        raise InputRejected(
            "stability scores named substitutions, so it needs `mutations`, "
            "for example [('ubiquitin', 'L8P')]"
        )
    if mutations and task != "stability":
        raise InputRejected(f"{task} takes no mutations; only stability does")

    # A directory of its own, so the manifest the tool writes can be read back
    # and nothing is left behind in the notebook's working directory.
    with tempfile.TemporaryDirectory(prefix="opendna-") as work:
        work_dir = Path(work)
        fasta = work_dir / "input.fasta"
        fasta.write_text(_fasta(sequences), encoding="utf-8")
        output = work_dir / "out.json"

        args: list[str] = []
        if offline:
            args.append("--offline")
        if batch_tokens is not None:
            args += ["--batch-tokens", str(batch_tokens)]
        if device is not None:
            args += ["--device", device]
        args += ["predict", "--task", task, "--format", "json", "-i", str(fasta), "-o", str(output)]
        if min_confidence is not None:
            args += ["--min-confidence", str(min_confidence)]
        if mutations:
            rows = ["name,mutation"] + [f"{name},{label}" for name, label in mutations]
            mutation_file = work_dir / "mutations.csv"
            mutation_file.write_text("\n".join(rows) + "\n", encoding="utf-8")
            args += ["--mutations", str(mutation_file)]

        _run(args, models=models)

        document = json.loads(output.read_text(encoding="utf-8"))
        manifest_path = output.with_suffix("").with_suffix(".manifest.json")
        if not manifest_path.is_file():
            # `out.json` -> `out.manifest.json`, which is how the tool names it.
            manifest_path = output.parent / (output.stem + ".manifest.json")
        manifest = (
            json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
        )

    return Run(
        task=document.get("task", task),
        rows=document.get("results", []),
        manifest=manifest,
        limitations=document.get("limitations", []),
        head_version=document.get("head_version", ""),
        embedding_model=document.get("embedding_model", ""),
        device=document.get("device", ""),
    )


def models(*, models_dir: str | os.PathLike[str] | None = None) -> str:
    """What is installed, as the command-line tool reports it.

    Returns its report verbatim rather than parsing it, because the report is
    prose meant for a person and parsing prose into a promise is how a check
    stops checking.
    """
    return _run(["models"], models=models_dir)


def selftest(*, models_dir: str | os.PathLike[str] | None = None) -> str:
    """Check that this machine reproduces the heads' recorded values.

    Needs all four heads installed. Raises if any check fails.
    """
    return _run(["selftest"], models=models_dir)
