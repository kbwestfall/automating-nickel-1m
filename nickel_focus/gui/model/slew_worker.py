"""
The `QThread` engine that drives a `slew.NickelTelescopePointing.slew_to`
call off the GUI thread.

See `gui.model.focus_worker.FocusWorker`'s module docstring: every
hardware call here is blocking, and `slew_to` can wait up to five minutes
for the telescope to reach its target, so running it directly on the
GUI's event-loop thread would freeze the whole interface.
"""
from nickel_focus.gui.qt import QtCore


class SlewWorker(QtCore.QThread):
    """
    Drive one `slew.NickelTelescopePointing.slew_to` call on a background
    thread.

    A new `SlewWorker` is expected to be created for each slew, matching
    `~PySide6.QtCore.QThread`'s own expectation that `start` is called at
    most once per instance.

    Parameters
    ----------
    telescope : `slew.NickelTelescopePointing`
        The telescope to command. Must already be constructed; this
        class only calls `~slew.NickelTelescopePointing.slew_to` on it.
    ra
        Target right ascension, passed through unchanged to
        `~slew.NickelTelescopePointing.slew_to`.
    dec
        Target declination, passed through unchanged to
        `~slew.NickelTelescopePointing.slew_to`.

    Attributes
    ----------
    slewFinished : `~PySide6.QtCore.Signal`
        Emitted once the telescope has reached its target successfully.
    slewFailed : `~PySide6.QtCore.Signal`
        Emitted once, with a human-readable message, if
        `~slew.NickelTelescopePointing.slew_to` raises -- movement
        disabled, tracking disabled, not ready to move, or a timeout
        waiting to reach the target.
    """
    slewFinished = QtCore.Signal()
    slewFailed = QtCore.Signal(str)

    def __init__(self, telescope, ra, dec, parent=None):
        super().__init__(parent)
        self.telescope = telescope
        self.ra = ra
        self.dec = dec

    def run(self):
        """
        `QThread`'s entry point: runs ``telescope.slew_to(ra, dec)`` on
        the background thread started by `start`. See
        `~gui.model.focus_worker.FocusWorker.run` for why emitting
        signals from here is still safe to connect to GUI-thread slots.
        """
        try:
            self.telescope.slew_to(self.ra, self.dec)
        except ValueError as e:
            self.slewFailed.emit(f'Could not move to target: {e}')
            return
        self.slewFinished.emit()
