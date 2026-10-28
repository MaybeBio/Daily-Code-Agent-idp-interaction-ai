"""The Virtual OpenCell app: browse CELL-FM's virtual staining of the OpenCell proteins.

Kept in the tracked package, like the CondenSeq app, so the logic behind the Space is on
GitHub; the Space imports the same module from its vendored copy of cell_fm.
"""

from . import viewer  # noqa: F401
