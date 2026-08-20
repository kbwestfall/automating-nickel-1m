"""
GUI entry point.

Run with ``python -m gui.main`` from an environment with PySide6
installed (see :mod:`gui.qt` and `requirements-gui.txt`).

This is currently just scaffolding: it opens an empty main window to
prove the package structure and the optional-Qt import boundary both
work, before any of the real panels (`ImagePanel`, `FocusCurvePanel`,
`FocusControlPanel`) exist.
"""
import sys

from gui.qt import QtWidgets


def build_app():
    """Return the (singleton) :class:`~PySide6.QtWidgets.QApplication`."""
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)


def build_window():
    """Construct the (currently empty) main window."""
    window = QtWidgets.QMainWindow()
    window.setWindowTitle('Nickel Focus GUI')
    window.resize(800, 600)
    return window


def main():
    app = build_app()
    window = build_window()
    window.show()
    app.exec()


if __name__ == '__main__':
    main()
