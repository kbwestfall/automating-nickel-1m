"""
GUI entry point.

Run with ``python -m gui.main`` from an environment with PySide6
installed (see :mod:`gui.qt` and `requirements-gui.txt`).
"""
import sys

from gui.qt import QtWidgets
from gui.views.main_window import MainWindow
from gui.controller import Controller


def build_app():
    """Return the (singleton) :class:`~PySide6.QtWidgets.QApplication`."""
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)


def build_window():
    """Construct the main window (panels only; unwired -- see :func:`main`)."""
    return MainWindow()


def main():
    app = build_app()
    window = build_window()
    controller = Controller(window)  # noqa: F841 (kept alive for the life of main())
    window.show()
    app.exec()


if __name__ == '__main__':
    main()
