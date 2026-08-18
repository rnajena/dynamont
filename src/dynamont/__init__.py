from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dynamont")
except PackageNotFoundError:
    __version__ = "unknown"

from dynamont._dynamont import Aligner, PoreType
__all__ = [
    "Aligner",
    "PoreType",
    "__version__"
]