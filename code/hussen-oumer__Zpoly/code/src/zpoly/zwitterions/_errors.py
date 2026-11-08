"""Exceptions raised by the zwitterion library"""

from __future__ import annotations

import os
from collections.abc import Sequence


# rq-ef81db1c rq-687220b8
class ZwitterionError(Exception):
    """Base class for all zwitterion-library errors"""


class UnknownEntryError(ZwitterionError, KeyError):
    """No library entry has this name or alias"""

    def __init__(self, name: str, suggestions: Sequence[str] = ()) -> None:
        self.name = name
        self.suggestions = tuple(suggestions)
        hint = f" Did you mean: {', '.join(self.suggestions)}?" if self.suggestions else ""
        super().__init__(f"Unknown zwitterion '{name}'.{hint}")

    def __str__(self) -> str:  # KeyError would otherwise repr() the message
        return str(self.args[0])


class UnknownBlockError(ZwitterionError, KeyError):
    """No building block has this key"""

    def __init__(self, key: str) -> None:
        self.key = key
        super().__init__(f"Unknown building block '{key}'.")

    def __str__(self) -> str:
        return str(self.args[0])


class InvalidStructureError(ZwitterionError, ValueError):
    """A SMILES string does not describe a valid structure"""

    def __init__(self, smiles: str, reason: str) -> None:
        self.smiles = smiles
        self.reason = reason
        super().__init__(f"Invalid structure '{smiles}': {reason}")


class ChargeError(InvalidStructureError):
    """The structure is not a net-neutral zwitterion"""

    def __init__(self, smiles: str, net_charge: int, reason: str) -> None:
        self.net_charge = net_charge
        super().__init__(smiles, reason)


class IncompatibleBlocksError(ZwitterionError, ValueError):
    """Block attachment points do not form a chain"""

    def __init__(self, block_keys: Sequence[str], position: int, reason: str) -> None:
        self.block_keys = tuple(block_keys)
        self.position = position
        super().__init__(
            f"Blocks {list(self.block_keys)} cannot be joined at position {position}: {reason}"
        )


class NoPolymerizableGroupError(ZwitterionError, ValueError):
    """No polymerizable group was found where one is required"""


class AmbiguousPolymerizableGroupError(ZwitterionError, ValueError):
    """Several C=C groups could polymerize and none was chosen"""

    def __init__(self, smiles: str, candidates: Sequence[tuple[int, int]]) -> None:
        self.smiles = smiles
        self.candidates = tuple(candidates)
        super().__init__(
            f"'{smiles}' has several candidate polymerizable C=C groups {list(self.candidates)}; "
            "pass vinyl_atoms to choose one."
        )


class DuplicateEntryError(ZwitterionError, ValueError):
    """A name, alias or block key is defined twice"""

    def __init__(self, name: str, sources: Sequence[str]) -> None:
        self.name = name
        self.sources = tuple(sources)
        super().__init__(f"'{name}' is defined more than once (in {', '.join(self.sources)}).")


class LibraryFileError(ZwitterionError):
    """A library file is missing, malformed or contains an invalid entry"""

    def __init__(self, path: str | os.PathLike[str], reason: str) -> None:
        self.path = os.fspath(path)
        self.reason = reason
        super().__init__(f"{self.path}: {reason}")


class UnknownInitiatorError(UnknownEntryError):
    """No initiator has this name or alias"""

    def __init__(self, name: str, suggestions: Sequence[str] = ()) -> None:
        super().__init__(name, suggestions)
        hint = f" Did you mean: {', '.join(self.suggestions)}?" if self.suggestions else ""
        self.args = (f"Unknown initiator '{name}'.{hint}",)


class ChainError(ZwitterionError, ValueError):
    """A chain cannot be built from the given units"""


class EndGroupError(ChainError):
    """The requested end groups are not available"""

    def __init__(self, ends: str, initiator: str | None, reason: str) -> None:
        self.ends = ends
        self.initiator = initiator
        super().__init__(f"Cannot cap ends '{ends}' with initiator {initiator!r}: {reason}")


class VariantError(ZwitterionError, ValueError):
    """The entry cannot take the requested variant option"""

    def __init__(self, name: str, option: str, reason: str) -> None:
        self.name = name
        self.option = option
        super().__init__(f"Cannot apply the {option} option to '{name}': {reason}")
