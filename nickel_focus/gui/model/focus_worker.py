"""
The :class:`PySide6.QtCore.QThread` engine that drives a
:class:`~nickel_focus.focus.FocusSequence` off the GUI thread.

See GUI_DESIGN.md §4.3: every hardware call in
:class:`~nickel_focus.scripts.focus.NickelFocus` is blocking, and a real
exposure takes seconds to tens of seconds, so running a sequence directly on the
GUI's event-loop thread would freeze the whole interface.
:class:`~nickel_focus.gui.model.focus_worker.FocusWorker` runs it here instead
and reports progress back via Qt signals.
"""
from nickel_focus.gui.qt import QtCore


class FocusWorker(QtCore.QThread):
    """
    Drive one :class:`~nickel_focus.focus.FocusSequence` run on a background
    thread.

    A new :class:`~nickel_focus.gui.model.focus_worker.FocusWorker` is expected
    to be created for each run (start a sequence, reanalyze one) rather than
    reused, matching :class:`PySide6.QtCore.QThread`'s own expectation that
    :meth:`PySide6.QtCore.QThread.start` is called at most once per instance.

    Parameters
    ----------
    focus_sequence : :class:`~nickel_focus.focus.FocusSequence`
        The sequence to drive.  Must already be constructed (e.g., a
        :class:`~nickel_focus.focus.GridFocusSequence` or
        :class:`~nickel_focus.focus.ArchiveFocusSequence`); this class only
        calls :meth:`~nickel_focus.focus.FocusSequence.step`,
        :meth:`~nickel_focus.focus.FocusSequence.reanalyze`, and
        :meth:`~nickel_focus.focus.FocusSequence.fit_best_focus` on it.
    method : :obj:`str`, :obj:`tuple`, optional
        The photometry method to use; see
        :func:`~nickel_focus.photometry.image_quality`.
    mode : :obj:`str`, optional
        ``'step'`` to advance the sequence
        (:meth:`~nickel_focus.focus.FocusSequence.step`, taking new exposures or
        replaying archived ones), ``'reanalyze'`` to re-run photometry on
        exposures already collected
        (:meth:`~nickel_focus.focus.FocusSequence.reanalyze`), without taking
        any new ones, or ``'single'`` to take one confirmation exposure at
        ``focus_value``
        (:meth:`~nickel_focus.focus.FocusSequence.take_single_exposure`) -- used
        for "Move to Best Focus" (GUI_DESIGN.md §5.4) and the standalone
        single-exposure workflow (§5.5).
    exp_kwargs : :obj:`dict`, optional
        Exposure settings (``record``, ``speed``, ``binning``, ``exptime``)
        applied via :meth:`~nickel_focus.focus.ExposureConfig.configure` before
        stepping, matching what the old CLI-only
        :meth:`~nickel_focus.focus.FocusSequence.execute` did.  Meaningful for
        ``mode='step'``/``'single'`` against a sequence with real (or fake)
        exposure hardware; ignored for ``mode='reanalyze'`` (no new exposures
        are taken) and harmless for archive/replay sequences (which have no
        exposure hardware to configure at all).
    focus_value : :obj:`int`, :obj:`float`, optional
        The focus value to move to before exposing.  Required for
        ``mode='single'``; ignored otherwise.

    Attributes
    ----------
    stepComplete : :class:`PySide6.QtCore.Signal`
        Emitted with one :class:`~nickel_focus.focus.StepResult` each time the
        driven generator yields.  Not emitted for ``mode='single'``, which has
        no generator to drive -- see
        :attr:`~nickel_focus.gui.model.focus_worker.FocusWorker.singleExposureFinished`.
    focusSequenceFinished : :class:`PySide6.QtCore.Signal`
        Emitted once, with ``(best_focus, best_fwhm)``, after the driven
        generator is exhausted (or stopped) and a quadratic fit to the results
        collected so far succeeds.  Not emitted for ``mode='single'``.
    focusSequenceFailed : :class:`PySide6.QtCore.Signal`
        Emitted once, with a human-readable message, if the sequence raises
        while stepping/reanalyzing/exposing, or if too few points remain for
        :meth:`~nickel_focus.focus.FocusSequence.fit_best_focus` to fit (e.g.,
        after an early
        :meth:`~nickel_focus.gui.model.focus_worker.FocusWorker.request_stop`).
        Mutually exclusive with
        :attr:`~nickel_focus.gui.model.focus_worker.FocusWorker.focusSequenceFinished`/
        :attr:`~nickel_focus.gui.model.focus_worker.FocusWorker.singleExposureFinished`:
        exactly one of the three fires per run.
    singleExposureFinished : :class:`PySide6.QtCore.Signal`
        Emitted once, with the resulting
        :class:`~nickel_focus.focus.StepResult`, when ``mode='single'``
        completes successfully.
    """
    stepComplete = QtCore.Signal(object)
    focusSequenceFinished = QtCore.Signal(float, float)
    focusSequenceFailed = QtCore.Signal(str)
    singleExposureFinished = QtCore.Signal(object)

    def __init__(
        self, focus_sequence, method='brightest', mode='step', exp_kwargs=None, focus_value=None,
        parent=None
    ):
        super().__init__(parent)
        if mode not in ('step', 'reanalyze', 'single'):
            raise ValueError(f"mode must be 'step', 'reanalyze', or 'single', got {mode!r}")
        if mode == 'single' and focus_value is None:
            raise ValueError("focus_value is required for mode='single'")
        self.focus_sequence = focus_sequence
        self.method = method
        self.mode = mode
        self.exp_kwargs = {} if exp_kwargs is None else exp_kwargs
        self.focus_value = focus_value
        self.stop_requested = False

    def request_stop(self):
        """
        Ask the sequence to stop after the exposure currently in progress
        finishes -- never mid-exposure (GUI_DESIGN.md §4.3). Has no effect once
        the run has already finished.
        """
        self.stop_requested = True

    def run(self):
        """
        :class:`PySide6.QtCore.QThread`'s entry point: everything this method
        does runs on the new background thread, not the GUI thread, once
        :meth:`PySide6.QtCore.QThread.start` is called.  Overriding
        :meth:`~nickel_focus.gui.model.focus_worker.FocusWorker.run` (rather
        than calling it directly) is what actually moves the work to that other
        thread -- Qt's machinery handles creating the OS-level thread and
        invoking this method on it.  Signal ``.emit()`` calls made from here
        (e.g.  ``self.stepComplete.emit(result)`` below) are still safe to
        connect to GUI-thread slots: Qt automatically queues the delivery back
        onto the GUI thread's event loop instead of calling the slot
        immediately, so widgets are never touched from the wrong thread.
        """
        self.stop_requested = False

        if self.mode == 'single':
            try:
                result = self.focus_sequence.take_single_exposure(
                    self.focus_value, method=self.method, **self.exp_kwargs
                )
            except Exception as e:
                self.focusSequenceFailed.emit(f'Could not move to best focus: {e}')
                return
            self.singleExposureFinished.emit(result)
            return

        generator = (
            self.focus_sequence.step(method=self.method) if self.mode == 'step'
            else self.focus_sequence.reanalyze(method=self.method)
        )

        try:
            if self.mode == 'step' and self.focus_sequence._exposure is not None:
                self.focus_sequence._exposure.cfg.configure(**self.exp_kwargs)
            for result in generator:
                self.stepComplete.emit(result)
                if self.stop_requested:
                    break
        except Exception as e:
            self.focusSequenceFailed.emit(f'Sequence failed: {e}')
            return

        try:
            best_focus, best_fwhm = self.focus_sequence.fit_best_focus(
                self.focus_sequence.observed_focus, self.focus_sequence.img_quality
            )
        except ValueError as e:
            self.focusSequenceFailed.emit(
                f'Stopped early -- not enough points for a focus fit: {e}'
            )
            return

        self.focusSequenceFinished.emit(best_focus, best_fwhm)
