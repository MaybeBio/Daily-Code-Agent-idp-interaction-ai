"""Exceptions raised by the polymer chain builder"""

from __future__ import annotations


# rq-280c4e6d
class StructureError(Exception):
    """Base class for chain-builder errors"""


class SequenceError(StructureError, ValueError):
    """The sequence description is invalid"""


class BuildError(StructureError):
    """No acceptable 3D structure could be built"""

    def __init__(self, message: str, min_distance: float | None = None) -> None:
        self.min_distance = min_distance
        super().__init__(message)


class PdbLimitError(StructureError):
    """The structure exceeds the PDB format's field limits"""
