"""
The `QThread` engine that drives a `focus.FocusSequence` off the GUI
thread.

See GUI_DESIGN.md §4.3: every hardware call in `scripts/focus.py` is
blocking, and a real exposure takes seconds to tens of seconds, so
running a sequence directly on the GUI's event-loop thread would freeze
the whole interface. `SequenceWorker` runs it here instead and reports
progress back via Qt signals.
"""
from nickel_focus.gui.qt import QtCore


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
        taking new exposures or replaying archived ones), ``'reanalyze'``
        to re-run photometry on exposures already collected
        (:func:`~focus.FocusSequence.reanalyze`), without taking any new
        ones, or ``'single'`` to take one confirmation exposure at
        ``focus_value`` (:func:`~focus.FocusSequence.take_single_exposure`)
        -- used for "Move to Best Focus" (GUI_DESIGN.md §5.4) and the
        standalone single-exposure workflow (§5.5).
    exp_kwargs : :obj:`dict`, optional
        Exposure settings (``record``, ``speed``, ``binning``,
        ``exptime``) applied via
        :func:`~focus.ExposureConfig.configure` before stepping, matching
        what the old CLI-only :func:`~focus.FocusSequence.execute` did.
        Meaningful for ``mode='step'``/``'single'`` against a sequence
        with real (or fake) exposure hardware; ignored for
        ``mode='reanalyze'`` (no new exposures are taken) and harmless
        for archive/replay sequences (which have no exposure hardware to
        configure at all).
    focus_value : :obj:`int`, :obj:`float`, optional
        The focus value to move to before exposing. Required for
        ``mode='single'``; ignored otherwise.

    Attributes
    ----------
    stepComplete : :class:`~PySide6.QtCore.Signal`
        Emitted with one :class:`focus.StepResult` each time the driven
        generator yields. Not emitted for ``mode='single'``, which has no
        generator to drive -- see ``singleExposureFinished``.
    sequenceFinished : :class:`~PySide6.QtCore.Signal`
        Emitted once, with ``(best_focus, best_fwhm)``, after the driven
        generator is exhausted (or stopped) and a quadratic fit to the
        results collected so far succeeds. Not emitted for ``mode='single'``.
    sequenceFailed : :class:`~PySide6.QtCore.Signal`
        Emitted once, with a human-readable message, if the sequence
        raises while stepping/reanalyzing/exposing, or if too few points
        remain for :func:`~focus.FocusSequence.fit_best_focus` to fit
        (e.g., after an early :func:`request_stop`). Mutually exclusive
        with ``sequenceFinished``/``singleExposureFinished``: exactly one
        of the three fires per run.
    singleExposureFinished : :class:`~PySide6.QtCore.Signal`
        Emitted once, with the resulting :class:`focus.StepResult`, when
        ``mode='single'`` completes successfully.
    """
    stepComplete = QtCore.Signal(object)
    sequenceFinished = QtCore.Signal(float, float)
    sequenceFailed = QtCore.Signal(str)
    singleExposureFinished = QtCore.Signal(object)

    def __init__(self, sequence, method='brightest', mode='step', exp_kwargs=None,
                 focus_value=None, parent=None):
        super().__init__(parent)
        if mode not in ('step', 'reanalyze', 'single'):
            raise ValueError(f"mode must be 'step', 'reanalyze', or 'single', got {mode!r}")
        if mode == 'single' and focus_value is None:
            raise ValueError("focus_value is required for mode='single'")
        self.sequence = sequence
        self.method = method
        self.mode = mode
        self.exp_kwargs = {} if exp_kwargs is None else exp_kwargs
        self.focus_value = focus_value
        self.stop_requested = False

    def request_stop(self):
        """
        Ask the sequence to stop after the exposure currently in progress
        finishes -- never mid-exposure (GUI_DESIGN.md §4.3). Has no effect
        once the run has already finished.
        """
        self.stop_requested = True

    def run(self):
        """
        `QThread`'s entry point: everything this method does runs on the
        new background thread, not the GUI thread, once :func:`start`
        is called. Overriding `run()` (rather than calling it directly)
        is what actually moves the work to that other thread -- Qt's
        machinery handles creating the OS-level thread and invoking this
        method on it. Signal ``.emit()`` calls made from here (e.g.
        ``self.stepComplete.emit(result)`` below) are still safe to
        connect to GUI-thread slots: Qt automatically queues the
        delivery back onto the GUI thread's event loop instead of
        calling the slot immediately, so widgets are never touched from
        the wrong thread.
        """
        self.stop_requested = False

        if self.mode == 'single':
            try:
                result = self.sequence.take_single_exposure(
                    self.focus_value, method=self.method, **self.exp_kwargs)
            except Exception as e:
                self.sequenceFailed.emit(f'Could not move to best focus: {e}')
                return
            self.singleExposureFinished.emit(result)
            return

        generator = (
            self.sequence.step(method=self.method) if self.mode == 'step'
            else self.sequence.reanalyze(method=self.method)
        )

        try:
            if self.mode == 'step' and self.sequence._exposure is not None:
                self.sequence._exposure.cfg.configure(**self.exp_kwargs)
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
