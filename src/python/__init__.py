from importlib.metadata import PackageNotFoundError, version

try:
  __version__ = version("dynamont")
except PackageNotFoundError:
  __version__ = "unknown"

try:
  from python._dynamont import Aligner, PoreType
except ImportError:
  __all__ = ["__version__"]
else:
  __all__ = ["Aligner", "PoreType", "__version__"]