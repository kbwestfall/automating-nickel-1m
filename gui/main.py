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
    """
    Build the window and its `Controller`, show the window, and run the
    Qt event loop until the user quits.
    """
    app = build_app()
    window = build_window()
    # `controller` looks unused (hence the noqa), but it must stay alive
    # for as long as the window does: PySide6 only keeps a signal
    # connection working while *something* in Python still holds a
    # reference to the connected object. Once `main()`'s local variable
    # `controller` would normally go out of scope, nothing in Python
    # would reference this Controller instance any more (Qt's C++ side
    # doesn't parent it to the window automatically), and it -- along
    # with every signal connection it made in `Controller.__init__` --
    # could be garbage-collected out from under the running window.
    controller = Controller(window)  # noqa: F841 (kept alive for the life of main())
    window.show()
    # QApplication.exec() starts Qt's event loop: it blocks here,
    # dispatching mouse clicks, key presses, timers, and cross-thread
    # signal deliveries (e.g. from SequenceWorker) to the right widgets
    # one at a time, until the application quits (e.g. the last window
    # closes). Nothing in this GUI runs unless this loop is running --
    # that's true of every Qt application, not just this one.
    app.exec()


if __name__ == '__main__':
    main()
