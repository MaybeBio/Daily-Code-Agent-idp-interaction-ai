"""Shared localCIDER accessor for physicochemical sequence features.

Composition, charge, charge patterning, and hydropathy are all read from localCIDER's
SequenceParameters (Holehouse et al., 2017), the field-standard reference, so every value
matches it exactly instead of a local reimplementation. One SequenceParameters object is built
per unique segment and cached, so the charge, composition, and hydropathy calls for a segment
share it.
"""

from __future__ import annotations

from functools import lru_cache

from .disorder import standardize_sequence


@lru_cache(maxsize=8192)
def params(seq: str):
    """Cached localCIDER SequenceParameters, built on the standardized sequence."""
    from localcider.sequenceParameters import SequenceParameters

    return SequenceParameters(standardize_sequence(seq))
