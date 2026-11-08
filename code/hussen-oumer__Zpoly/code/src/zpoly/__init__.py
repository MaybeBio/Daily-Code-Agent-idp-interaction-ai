"""Zwitterionic polymers for LLPS simulations"""

from ._version import __version__
from .inputfile import parse_input, run_input  # noqa: E402,F401

# rq-b39965b3
from .simple import BuildResult, build, dft_inputs  # noqa: E402,F401
from .structure import alternating, block, homopolymer, pattern, random_copolymer  # noqa: E402,F401
from .system import analyze_slab, analyze_tg, bulk_system, fit_tg, slab_system  # noqa: E402,F401
