"""
Package initialization.

The current main purpose of this is to provide package-level globals
that can be imported by submodules.
"""

from .pkg.version import version
from .pkg.ktl import ktl

# Set version
__version__ = version

# Start the log
import logging
from .pkg.logger import get_logger
log = get_logger(level=logging.DEBUG)
