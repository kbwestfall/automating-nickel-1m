"""
Test configuration for GUI-specific tests.

Skips the entire `tests/gui/` subtree if PySide6 isn't installed, so the
rest of the test suite (Qt-free) is unaffected either way. Also forces
the offscreen Qt platform plugin, so these tests never need a real
display, and adds the repository root to :obj:`sys.path` so ``gui`` can
be imported as a top-level package (mirroring what the top-level
`tests/conftest.py` already does for `scripts/`).
"""
import os
import sys
from pathlib import Path

import pytest

pytest.importorskip('PySide6')

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))


@pytest.fixture(scope='session')
def qapp():
    """The (session-wide, singleton) QApplication needed to create any Qt widget."""
    from gui.qt import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
