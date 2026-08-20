"""Tests for :mod:`gui.qt`, the single PySide6 import point."""
import sys

import pytest


def test_qt_shim_exposes_expected_names(qapp):
    import gui.qt
    assert hasattr(gui.qt, 'QtCore'), 'gui.qt should expose QtCore'
    assert hasattr(gui.qt, 'QtGui'), 'gui.qt should expose QtGui'
    assert hasattr(gui.qt, 'QtWidgets'), 'gui.qt should expose QtWidgets'


def test_qt_shim_raises_clear_error_without_pyside6(monkeypatch):
    # Simulate PySide6 being uninstalled: setting a module to None in
    # sys.modules makes any subsequent `import <name>` raise ImportError
    # immediately, without touching the real (already-installed) package.
    monkeypatch.setitem(sys.modules, 'PySide6', None)
    sys.modules.pop('gui.qt', None)

    with pytest.raises(ImportError, match='PySide6'):
        import gui.qt
