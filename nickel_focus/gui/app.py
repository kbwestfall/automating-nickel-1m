"""
GUI entry point.

Run with ``python -m gui.main`` from an environment with PySide6
installed (see :mod:`gui.qt` and `requirements-gui.txt`).
"""
from nickel_focus.gui.qt import QtWidgets
from nickel_focus.gui.views.main_window import MainWindow


def build_app():
    """Return the (singleton) :class:`~PySide6.QtWidgets.QApplication`."""
    return QtWidgets.QApplication.instance() # or QtWidgets.QApplication(sys.argv)


def build_window():
    """Construct the main window (panels only; unwired -- see :func:`main`)."""
    return MainWindow()
