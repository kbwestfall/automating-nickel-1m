"""
Wires the View (`gui.views.main_window.MainWindow` and its panels) to the
Model (`focus.FocusSequence` subclasses, via
`gui.model.sequence_worker.SequenceWorker`).

See GUI_DESIGN.md §4.3.
"""
import focus
from gui.model.sequence_worker import SequenceWorker
from gui.qt import QtCore


class Controller(QtCore.QObject):
    """
    Owns the currently active :class:`focus.FocusSequence` (if any) and
    the :class:`SequenceWorker` driving it, and mediates between it and a
    :class:`~gui.views.main_window.MainWindow`.

    Implements the "hardware exclusivity" state machine from §4.3: only
    one operation can be active at a time, and interactive source
    selection (`ImagePanel.sourceSelected`) is only enabled while nothing
    is running. As of this sub-phase, all three sequence types (Archive,
    Grid, Automated), reanalysis, "Move to Best Focus", and the
    standalone single-exposure workflow (§5.5) are all wired up.
    """

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.window = window
        self.sequence = None        # the current focus.FocusSequence, or None
        self.worker = None          # the SequenceWorker currently running, or None
        self.method = 'brightest'   # the photometry method currently in effect

        # A standalone single exposure (§5.5) not yet committed to
        # `sequence`, plus the throwaway `focus.FocusSequence` used only
        # to take/reanalyze it (its hardware handles, not its
        # bookkeeping, matter -- see `take_single_exposure`).
        self.pending_result = None
        self._standalone_sequence = None

        window.control_panel.startRequested.connect(self.start_sequence)
        window.control_panel.stopRequested.connect(self.stop)
        window.control_panel.methodChanged.connect(self.set_method)
        window.control_panel.moveToBestFocusRequested.connect(self.move_to_best_focus)
        window.control_panel.takeSingleExposureRequested.connect(self.take_single_exposure)
        window.control_panel.addToSequenceRequested.connect(self.add_pending_to_sequence)
        window.image_panel.sourceSelected.connect(self._on_source_selected)

        self._set_running(False)

    # -- actions the View can request ---------------------------------------

    def start_sequence(self):
        """Build a sequence from the control panel's configuration and run it."""
        if self.worker is not None:
            return  # hardware exclusivity: something is already running

        sequence_type = self.window.control_panel.get_sequence_type()
        config = self.window.control_panel.get_sequence_config()

        try:
            if sequence_type == 'archive':
                sequence = self._build_archive_sequence(config)
            elif sequence_type == 'grid':
                sequence = focus.GridFocusSequence(config['start'], config['step'],
                                                    nstep=config['nstep'])
            else:
                sequence = focus.AutomatedFocusSequence(config['start'], config['step'],
                                                         maxsteps=config['maxsteps'])
        except Exception as e:
            self.window.control_panel.show_failure(f'Could not start sequence: {e}')
            return

        if sequence_type != 'archive' and (sequence._focus is None or sequence._exposure is None):
            self.window.control_panel.show_failure(
                'Could not start sequence: no ktl connection is available for a live sequence.')
            return

        self.sequence = sequence
        self._clear_pending()
        self.window.image_panel.reset()
        self.window.curve_panel.reset()
        self.window.control_panel.reset()

        exp_kwargs = (None if sequence_type == 'archive'
                      else self.window.control_panel.get_exposure_config())
        self._start_worker(self.sequence, mode='step', exp_kwargs=exp_kwargs)

    def _build_archive_sequence(self, config):
        # GridFocusSequence is only used here for its focus-value
        # arithmetic (matches focus.py's own CLI archive-mode path); it
        # never touches hardware since ktl isn't connected for an archive
        # run.
        grid = focus.GridFocusSequence(config['start'], config['step'], nstep=config['nstep'])
        expected_files = [
            config['datadir'] / f"{config['prefix']}{config['obsnum'] + i}{config['suffix']}"
            for i in range(grid.nstep)
        ]
        missing = [f for f in expected_files if not f.is_file()]
        if missing:
            raise FileNotFoundError(
                'Expected to find the following files, but they are not available: '
                + ', '.join(str(f) for f in missing))
        return focus.ArchiveFocusSequence(list(grid.target_focus), expected_files)

    def reanalyze(self):
        """
        Re-run photometry on the already-collected exposures of whatever
        is currently loaded: the active `sequence`, or -- if none is
        loaded -- a pending standalone single exposure (§5.6's "no
        sequence loaded yet" case).
        """
        target = self.sequence if self.sequence is not None else self._standalone_sequence
        if self.worker is not None or target is None or not target.exposures:
            return
        self._start_worker(target, mode='reanalyze')

    def move_to_best_focus(self, focus_value):
        """
        Move to ``focus_value`` and take one confirmation exposure
        (:func:`focus.FocusSequence.take_single_exposure`) -- the action
        behind the "Move to Best Focus" button, which only emits
        `~gui.views.focus_control_panel.FocusControlPanel.moveToBestFocusRequested`
        once the user has confirmed a dialog and only when the finished
        sequence had a hardware connection to move.
        """
        if self.worker is not None or self.sequence is None:
            return  # hardware exclusivity: something is already running
        exp_kwargs = self.window.control_panel.get_exposure_config()
        self._start_worker(self.sequence, mode='single', exp_kwargs=exp_kwargs,
                            focus_value=focus_value)

    def take_single_exposure(self, focus_value):
        """
        Take one exposure at ``focus_value`` with no sequence bookkeeping
        (§5.5), e.g. to confirm the telescope landed on the intended
        field, or to mark a source via `ImagePanel.sourceSelected` before
        a real sequence starts. Uses a throwaway `focus.FocusSequence`
        for its hardware handles, independent of any loaded `sequence`.
        The result is held as `pending_result` until the user commits it
        (:func:`add_pending_to_sequence`) or discards it by starting a
        new sequence.
        """
        if self.worker is not None:
            return  # hardware exclusivity: something is already running
        standalone = focus.FocusSequence()
        if standalone._focus is None or standalone._exposure is None:
            self.window.control_panel.show_failure(
                'Could not take single exposure: no ktl connection is available.')
            return
        self._standalone_sequence = standalone
        exp_kwargs = self.window.control_panel.get_exposure_config()
        self._start_worker(standalone, mode='single', exp_kwargs=exp_kwargs,
                            focus_value=focus_value)

    def add_pending_to_sequence(self):
        """
        Commit `pending_result` -- a standalone single exposure (§5.5) --
        into the currently loaded `sequence`'s data, as though it had
        been collected as part of that sequence, and refit/redraw
        `~gui.views.focus_curve_panel.FocusCurvePanel`.
        """
        if self.pending_result is None or self.sequence is None:
            return
        result = self.pending_result
        seq = self.sequence
        seq.observed_focus.append(result.focus_value)
        seq.exposures.append(result.exposure)
        seq.img_quality.append(result.fwhm)
        seq.source_stamps.append(result.stamp)
        seq.centroids.append(result.centroid)
        seq.step_iter = len(seq.exposures)
        result.index = seq.step_iter - 1
        result.is_outlier = focus.FocusPlot.is_outlier(seq.centroids)
        self.window.image_panel.update_result(result)
        self.window.curve_panel.add_result(result)
        self._clear_pending()

    def stop(self):
        """Request that the running sequence stop between steps (§4.3)."""
        if self.worker is None:
            return
        self.worker.request_stop()
        self.window.control_panel.set_stopping()

    def set_method(self, method):
        """
        Update the photometry method used for future steps/reanalysis.
        ``method`` is ``'brightest'``, ``'weighted'``, or an ``(x, y)``
        coordinate tuple (from `ImagePanel.sourceSelected`).
        """
        self.method = method
        self.window.control_panel.set_method(method)

    # -- internals ------------------------------------------------------------

    def _start_worker(self, sequence, mode, exp_kwargs=None, focus_value=None):
        self.worker = SequenceWorker(sequence, method=self.method, mode=mode,
                                      exp_kwargs=exp_kwargs, focus_value=focus_value)
        self.worker.stepComplete.connect(self._on_step_complete)
        self.worker.sequenceFinished.connect(self._on_sequence_finished)
        self.worker.sequenceFailed.connect(self._on_sequence_failed)
        self.worker.singleExposureFinished.connect(self._on_single_exposure_finished)
        self.worker.finished.connect(self._on_worker_finished)
        self._set_running(True)
        self.worker.start()

    def _on_step_complete(self, result):
        reanalyzing_standalone = (
            self.worker is not None and self.worker.mode == 'reanalyze'
            and self._standalone_sequence is not None
            and self.worker.sequence is self._standalone_sequence
        )
        if self.worker is not None and self.worker.mode == 'reanalyze':
            self.window.image_panel.update_result(result)
            if reanalyzing_standalone:
                # Not on the curve panel -- a standalone exposure was
                # never added there either (§5.5: it isn't sequence data
                # until explicitly committed).
                self.pending_result = result
            else:
                self.window.curve_panel.update_result(result)
        else:
            self.window.image_panel.add_result(result)
            self.window.curve_panel.add_result(result)

        total = self.sequence.expected_steps if self.sequence is not None else None
        self.window.control_panel.update_step(result, total_expected=total)

    def _on_sequence_finished(self, best_focus, best_fwhm):
        can_move = (self.sequence is not None and self.sequence._focus is not None
                    and self.sequence._exposure is not None)
        self.window.control_panel.show_best_focus(best_focus, best_fwhm, can_move=can_move)

    def _on_sequence_failed(self, message):
        self.window.control_panel.show_failure(message)

    def _on_single_exposure_finished(self, result):
        self.window.image_panel.add_result(result)

        moving_to_best_focus = (self.worker is not None and self.sequence is not None
                                 and self.worker.sequence is self.sequence)
        if moving_to_best_focus:
            self.window.control_panel.show_confirmation(result)
            return

        # Standalone single-exposure workflow (§5.5): seed the throwaway
        # sequence's bookkeeping with this one exposure so `reanalyze()`
        # (via interactive source selection, §5.6) can run against it
        # even though no real sequence is loaded.
        seq = self._standalone_sequence
        seq.observed_focus.append(result.focus_value)
        seq.exposures.append(result.exposure)
        seq.img_quality.append(result.fwhm)
        seq.source_stamps.append(result.stamp)
        seq.centroids.append(result.centroid)
        seq.step_iter = 1
        self.pending_result = result
        self.window.control_panel.show_pending_exposure(result, can_add=self.sequence is not None)

    def _on_worker_finished(self):
        if self.worker is not None:
            self.worker.wait()
        self.worker = None
        self._set_running(False)

    def _on_source_selected(self, x, y):
        self.set_method((x, y))
        self.reanalyze()

    def _set_running(self, running):
        self.window.control_panel.set_running(running)
        self.window.image_panel.set_selection_enabled(not running)

    def _clear_pending(self):
        self.pending_result = None
        self._standalone_sequence = None
        self.window.control_panel.clear_pending_exposure()
