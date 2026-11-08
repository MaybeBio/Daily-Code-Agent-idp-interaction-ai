"""Multi-chain systems"""

# rq-4a8e3ff0
from ._analysis import SlabResult, TgResult, analyze_slab, analyze_tg, fit_tg, phase_concentrations
from ._bulk import bulk_system
from ._errors import SystemAnalysisError, SystemBuildError
from ._slab import slab_system

__all__ = ["SlabResult", "SystemAnalysisError", "SystemBuildError", "TgResult", "analyze_slab", "analyze_tg",
           "bulk_system", "fit_tg", "phase_concentrations", "slab_system"]
