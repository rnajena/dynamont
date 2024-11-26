from . import _version
__version__ = _version.get_versions()['version']

import os
__build__ = os.path.join(os.path.dirname(__file__), '..', '..', 'build')