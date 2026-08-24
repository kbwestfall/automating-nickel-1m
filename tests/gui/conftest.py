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


@pytest.fixture(autouse=True)
def isolate_qsettings(monkeypatch, tmp_path):
    """
    Redirect `MainWindow`'s settings persistence (GUI_DESIGN.md §9 phase
    5) to an isolated, per-test ini file instead of the real application
    preferences store. Without this, every test that constructs --  and
    especially every test that closes -- a `MainWindow` would read from
    and write to the developer's actual saved GUI settings on disk.
    """
    from gui.qt import QtCore
    from gui.views.main_window import MainWindow
    settings_path = str(tmp_path / 'test_settings.ini')
    monkeypatch.setattr(
        MainWindow, '_settings',
        lambda self: QtCore.QSettings(settings_path, QtCore.QSettings.Format.IniFormat))
