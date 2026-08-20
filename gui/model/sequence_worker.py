"""
The `QThread` engine that drives a `focus.FocusSequence` off the GUI
thread.

See GUI_DESIGN.md §4.3: every hardware call in `scripts/focus.py` is
blocking, and a real exposure takes seconds to tens of seconds, so
running a sequence directly on the GUI's event-loop thread would freeze
the whole interface. `SequenceWorker` runs it here instead and reports
progress back via Qt signals.
"""
from gui.qt import QtCore


class SequenceWorker(QtCore.QThread):
    """
    Drive one :class:`focus.FocusSequence` run on a background thread.

    A new :class:`SequenceWorker` is expected to be created for each run
    (start a sequence, reanalyze one) rather than reused, matching
    :class:`~PySide6.QtCore.QThread`'s own expectation that :func:`start`
    is called at most once per instance.

    Parameters
    ----------
    sequence : :class:`focus.FocusSequence`
        The sequence to drive. Must already be constructed (e.g., a
        :class:`focus.GridFocusSequence` or
        :class:`focus.ArchiveFocusSequence`); this class only calls
        :func:`~focus.FocusSequence.step`,
        :func:`~focus.FocusSequence.reanalyze`, and
        :func:`~focus.FocusSequence.fit_best_focus` on it.
    method : :obj:`str`, :obj:`tuple`, optional
        The photometry method to use; see
        :func:`photometry.image_quality`.
    mode : :obj:`str`, optional
        ``'step'`` to advance the sequence (:func:`~focus.FocusSequence.step`,
        taking new exposures or replaying archived ones), or
        ``'reanalyze'`` to re-run photometry on exposures already
        collected (:func:`~focus.FocusSequence.reanalyze`), without
        taking any new ones.

    Attributes
    ----------
    stepComplete : :class:`~PySide6.QtCore.Signal`
        Emitted with one :class:`focus.StepResult` each time the driven
        generator yields.
    sequenceFinished : :class:`~PySide6.QtCore.Signal`
        Emitted once, with ``(best_focus, best_fwhm)``, after the driven
        generator is exhausted (or stopped) and a quadratic fit to the
        results collected so far succeeds.
    sequenceFailed : :class:`~PySide6.QtCore.Signal`
        Emitted once, with a human-readable message, if the sequence
        raises while stepping/reanalyzing, or if too few points remain
        for :func:`~focus.FocusSequence.fit_best_focus` to fit (e.g.,
        after an early :func:`request_stop`). Mutually exclusive with
        ``sequenceFinished``: exactly one of the two fires per run.
    """
    stepComplete = QtCore.Signal(object)
    sequenceFinished = QtCore.Signal(float, float)
    sequenceFailed = QtCore.Signal(str)

    def __init__(self, sequence, method='brightest', mode='step', parent=None):
        super().__init__(parent)
        if mode not in ('step', 'reanalyze'):
            raise ValueError(f"mode must be 'step' or 'reanalyze', got {mode!r}")
        self.sequence = sequence
        self.method = method
        self.mode = mode
        self.stop_requested = False

    def request_stop(self):
        """
        Ask the sequence to stop after the exposure currently in progress
        finishes -- never mid-exposure (GUI_DESIGN.md §4.3). Has no effect
        once the run has already finished.
        """
        self.stop_requested = True

    def run(self):
        self.stop_requested = False
        generator = (
            self.sequence.step(method=self.method) if self.mode == 'step'
            else self.sequence.reanalyze(method=self.method)
        )

        try:
            for result in generator:
                self.stepComplete.emit(result)
                if self.stop_requested:
                    break
        except Exception as e:
            self.sequenceFailed.emit(f'Sequence failed: {e}')
            return

        try:
            best_focus, best_fwhm = self.sequence.fit_best_focus(
                self.sequence.observed_focus, self.sequence.img_quality)
        except ValueError as e:
            self.sequenceFailed.emit(f'Stopped early -- not enough points for a focus fit: {e}')
            return

        self.sequenceFinished.emit(best_focus, best_fwhm)
