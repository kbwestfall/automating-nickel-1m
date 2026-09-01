"""
Test configuration for GUI-specific tests.

Skips the entire `tests/gui/` subtree if PySide6 isn't installed, so the
rest of the test suite (Qt-free) is unaffected either way. Also forces
the offscreen Qt platform plugin, so these tests never need a real
display.
"""
import os

import pytest

pytest.importorskip('PySide6')

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from nickel_focus import log
from nickel_focus.gui.log_handler import QtLogHandler


@pytest.fixture(scope='session')
def qapp():
    """The (session-wide, singleton) QApplication needed to create any Qt widget."""
    from nickel_focus.gui.qt import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(autouse=True)
def _remove_qt_log_handlers():
    """
    Remove any `QtLogHandler` a test's `Controller` added to the
    `nickel_focus.log` singleton once that test finishes, so a handler
    pointing at a destroyed window never lingers into a later test module.
    """
    yield
    log.remove_handlers_of_type(QtLogHandler)


@pytest.fixture(autouse=True)
def isolate_qsettings(monkeypatch, tmp_path):
    """
    Redirect `MainWindow`'s settings persistence (GUI_DESIGN.md §9 phase
    5) to an isolated, per-test ini file instead of the real application
    preferences store. Without this, every test that constructs --  and
    especially every test that closes -- a `MainWindow` would read from
    and write to the developer's actual saved GUI settings on disk.
    """
    from nickel_focus.gui.qt import QtCore
    from nickel_focus.gui.views.main_window import MainWindow
    settings_path = str(tmp_path / 'test_settings.ini')
    monkeypatch.setattr(
        MainWindow, '_settings',
        lambda self: QtCore.QSettings(settings_path, QtCore.QSettings.Format.IniFormat))
