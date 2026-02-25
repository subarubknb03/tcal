from .tcal import Tcal, PairAnalysis

try:
    from .tcal_pyscf import TcalPySCF
    __all__ = ["Tcal", "PairAnalysis", "TcalPySCF"]
except ImportError:
    __all__ = ["Tcal", "PairAnalysis"]
