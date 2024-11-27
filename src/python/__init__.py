import os
__build__ = os.path.join(os.path.dirname(__file__), '..', '..', 'build')
try:
  from src.python._version import version as __version__
except ImportError:
  __version__ = "unknown"
