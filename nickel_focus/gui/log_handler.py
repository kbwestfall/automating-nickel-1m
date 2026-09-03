"""
A :class:`logging.Handler` that safely delivers log records into a GUI widget.
"""
import logging

from nickel_focus.gui.qt import QtCore


class QtLogHandler(QtCore.QObject, logging.Handler):
    """
    Forward formatted log records to any connected Qt slot via
    :attr:`~nickel_focus.gui.log_handler.QtLogHandler.record_logged`, instead of
    writing them anywhere itself. Log calls
    can happen on a background thread, so this handler never touches a widget
    directly -- Qt automatically queues a :class:`PySide6.QtCore.Signal`
    emission onto the receiving slot's own thread, which is what keeps this
    safe.

    Parameters
    ----------
    parent : :class:`PySide6.QtCore.QObject`, optional
        Passed through to the base :class:`PySide6.QtCore.QObject` constructor.

    Attributes
    ----------
    record_logged : :class:`PySide6.QtCore.Signal`
        Emitted with one formatted line of text per log record accepted by
        this handler (i.e., past its level filter).
    """
    record_logged = QtCore.Signal(str)

    def __init__(self, parent=None):
        QtCore.QObject.__init__(self, parent)
        logging.Handler.__init__(self)

    def emit(self, record):
        """
        Format ``record`` and emit it via
        :attr:`~nickel_focus.gui.log_handler.QtLogHandler.record_logged`.
        """
        self.record_logged.emit(self.format(record))
