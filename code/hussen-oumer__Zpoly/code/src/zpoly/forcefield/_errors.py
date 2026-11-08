"""Exceptions raised by the force-field layer"""

from __future__ import annotations

from collections.abc import Sequence


# rq-2802a0d0
class ForceFieldError(Exception):
    """Base class for force-field errors"""


class UnknownForceFieldError(ForceFieldError, KeyError):
    def __init__(self, name: str, available: Sequence[str]) -> None:
        self.name = name
        self.available = tuple(available)
        super().__init__(f"Unknown force field '{name}'. Available: {', '.join(self.available) or 'none'}.")

    def __str__(self) -> str:
        return str(self.args[0])


class ForceFieldUnavailableError(ForceFieldError):
    def __init__(self, name: str, reason: str, hint: str) -> None:
        self.name = name
        self.hint = hint
        super().__init__(f"Force field '{name}' is not available: {reason}. {hint}")


class ParametrizationError(ForceFieldError):
    def __init__(self, tool: str, returncode: int, log: str, message: str | None = None) -> None:
        self.tool = tool
        self.returncode = returncode
        self.log = "\n".join(log.splitlines()[-40:])
        super().__init__(message or f"{tool} failed (exit code {returncode}):\n{self.log}")


class MissingParameterError(ForceFieldError):
    """Bonded terms without parameters"""

    def __init__(self, missing, template: str = "") -> None:
        if isinstance(missing, str):  # backwards-compatible MissingParameterError(kind, types)
            missing = {missing: tuple(template)}
            template = ""
        self.missing = {k: tuple(v) for k, v in missing.items() if v}
        self.kind = next(iter(self.missing))
        self.types = self.missing[self.kind]
        self.template = template
        parts = []
        for kind, types in self.missing.items():
            shown = ", ".join("-".join(t) for t in types[:10])
            more = f" (and {len(types) - 10} more)" if len(types) > 10 else ""
            parts.append(f"No {kind} parameters for atom types {shown}{more}.")
        hint = ("\nAdd them with a parameter file (PARAMS = FILE in the input file). Fill in the values, or use "
                "same_as: [...] to copy an existing parameter:\n" + template) if template else ""
        super().__init__(" ".join(parts) + hint)


class ChargeFileError(ForceFieldError):
    def __init__(self, message: str, path: str | None = None, residue: str | None = None,
                 atom: str | None = None) -> None:
        self.path = path
        self.residue = residue
        self.atom = atom
        where = f"{path}: " if path else ""
        super().__init__(f"{where}{message}")


class TypingError(ForceFieldError):
    """OPLS-AA typing failed: a bad rule, or atoms no rule types"""

    def __init__(self, message: str, rule: str | None = None, atoms: Sequence[str] = ()) -> None:
        self.rule = rule
        self.atoms = tuple(atoms)
        super().__init__(message)


class ParameterFileError(ForceFieldError):
    """A user parameter file is invalid"""

    def __init__(self, message: str, path: str | None = None, entry: str | None = None) -> None:
        self.path = path
        self.entry = entry
        where = f"{path}: " if path else ""
        super().__init__(f"{where}{message}")
