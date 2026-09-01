"""
Package initialization.

The current main purpose of this is to provide package-level globals
that can be imported by submodules.
"""

# Set version
from .pkg.version import version
__version__ = version

# Start the log
import logging
from .pkg.logger import get_logger
log = get_logger(level=logging.INFO)

# Do a global import of ktl so that it is only done once.
from .pkg.ktl import ktl
