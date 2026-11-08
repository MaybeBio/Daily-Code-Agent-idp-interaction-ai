"""Errors of the system builder"""

from __future__ import annotations

from ..structure import StructureError


class SystemBuildError(StructureError):
    """A system cannot be built."""


class SystemAnalysisError(StructureError):
    """Results cannot be analysed."""
