"""
Wires the View (`gui.views.main_window.MainWindow` and its panels) to the
Model (`focus.FocusSequence` subclasses, via
`gui.model.focus_worker.FocusWorker`; and
`slew.NickelTelescopePointing`, via `gui.model.slew_worker.SlewWorker`).

See GUI_DESIGN.md §4.3.
"""
from astropy.coordinates import SkyCoord

from nickel_focus import focus
from nickel_focus import slew
from nickel_focus.gui.model.focus_worker import FocusWorker
from nickel_focus.gui.model.slew_worker import SlewWorker
from nickel_focus.gui.qt import QtCore


class Controller(QtCore.QObject):
    """
    Owns the currently active :class:`focus.FocusSequence` (if any) and
    the :class:`FocusWorker` driving it, and mediates between it and a
    :class:`~gui.views.main_window.MainWindow`.

    Implements the "hardware exclusivity" state machine from §4.3: only
    one operation can be active at a time, and interactive source
    selection (`ImagePanel.sourceSelected`) is only enabled while nothing
    is running. All three sequence types (Grid, Automated, Archive/
    Replay), reanalysis, and the Single tab (which also serves as "move
    to best focus" -- see `gui.views.focus_control_panel.FocusControlPanel`)
    are wired up.

    Also owns the `slew.NickelTelescopePointing` handle used by the Slew
    tab: `move_to_target` (a slew, via `SlewWorker`) participates in the
    same hardware-exclusivity state as a focus sequence -- the telescope
    shouldn't move while it's mid-exposure, or vice versa -- while
    `find_nearest_object`/`find_nearest_pointing`/`find_nearest_focus`
    (a single fast ktl read plus a starlist search, not a move) do not.
    A `~PySide6.QtCore.QTimer` polls the telescope's current position
    for the tab's live display regardless of what else is running.
    """

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.window = window
        self.focus_sequence = None  # the current focus.FocusSequence, or None
        self.focus_worker = None    # the FocusWorker currently running, or None
        self.method = 'brightest'   # the photometry method currently in effect
        self.slew_worker = None     # the SlewWorker currently running, or None

        # The throwaway `focus.FocusSequence` used only to take/reanalyze
        # a standalone Single-tab exposure (its hardware handles, not its
        # bookkeeping, matter -- see `take_single_exposure`). Seeded with
        # that one exposure's data so `reanalyze()` (via interactive
        # source selection, §5.6) can run against it even though no real
        # sequence is loaded.
        self._standalone_focus_sequence = None

        # The telescope-pointing handle backing the Slew tab. Like
        # `focus.FocusSequence`'s `_focus`/`_exposure`, this is `None`
        # when no ktl connection is available (`NickelTelescopePointing`
        # raises `RuntimeError` in that case) -- every method below that
        # uses it checks for `None` first and reports a clear failure
        # instead.
        try:
            self.telescope = slew.NickelTelescopePointing()
        except RuntimeError:
            self.telescope = None

        # This is Qt's signal/slot mechanism, the backbone of how every
        # View talks to this Controller: each `connect()` call below
        # says "whenever this Signal is emitted (e.g.
        # `window.control_panel.startRequested.emit()`, from a button's
        # `clicked` handler), call this plain Python method (the
        # 'slot')." The View never calls `self.start_focus_sequence()`
        # directly -- it only knows how to `.emit()` its own signals, so
        # it stays usable even if what eventually happens in response
        # (or whether anything does) changes.
        window.control_panel.startRequested.connect(self.start_focus_sequence)
        window.control_panel.stopRequested.connect(self.stop)
        window.control_panel.takeSingleExposureRequested.connect(self.take_single_exposure)
        window.control_panel.moveToTargetRequested.connect(self.move_to_target)
        window.control_panel.findNearestObjectRequested.connect(self.find_nearest_object)
        window.control_panel.findNearestPointingRequested.connect(self.find_nearest_pointing)
        window.control_panel.findNearestFocusRequested.connect(self.find_nearest_focus)
        window.image_panel.sourceSelected.connect(self._on_source_selected)

        # Milliseconds, not seconds: a `QTimer` interval is always in ms.
        # Reading `telescope.current` is one quick ktl-cache read, not a
        # move -- safe (and useful) to keep polling no matter what else,
        # if anything, is running.
        self._position_timer = QtCore.QTimer(self)
        self._position_timer.timeout.connect(self._update_current_position)
        if self.telescope is None:
            self.window.control_panel.set_current_position('—', '—')
        else:
            self._update_current_position()
            self._position_timer.start(1000)

        self._set_running(False)

    # -- actions the View can request ---------------------------------------

    def start_focus_sequence(self):
        """Build a sequence from the control panel's configuration and run it."""
        if self.focus_worker is not None or self.slew_worker is not None:
            return  # hardware exclusivity: something is already running

        focus_sequence_type = self.window.control_panel.get_focus_sequence_type()
        config = self.window.control_panel.get_focus_sequence_config()

        try:
            if focus_sequence_type == 'archive':
                focus_sequence = self._build_archive_focus_sequence(config)
            elif focus_sequence_type == 'grid':
                focus_sequence = focus.GridFocusSequence(config['start'], config['step'],
                                                          nstep=config['nstep'])
            else:
                focus_sequence = focus.AutomatedFocusSequence(config['start'], config['step'],
                                                               maxsteps=config['maxsteps'])
        except Exception as e:
            self.window.control_panel.show_failure(f'Could not start sequence: {e}')
            return

        if focus_sequence_type != 'archive' and (
                focus_sequence._focus is None or focus_sequence._exposure is None):
            self.window.control_panel.show_failure(
                'Could not start sequence: no ktl connection is available for a live sequence.')
            return

        self.focus_sequence = focus_sequence
        self._standalone_focus_sequence = None
        self.window.image_panel.reset()
        self.window.curve_panel.reset()
        self.window.control_panel.reset()

        exp_kwargs = (None if focus_sequence_type == 'archive'
                      else self.window.control_panel.get_exposure_config())
        self._start_focus_worker(self.focus_sequence, mode='step', exp_kwargs=exp_kwargs)

    def _build_archive_focus_sequence(self, config):
        """Build a `focus.ArchiveFocusSequence` from the Replay tab's ``config``, or raise."""
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
        is currently loaded: the active `focus_sequence`, or -- if none
        is loaded -- a standalone Single-tab exposure (§5.6's "no
        sequence loaded yet" case).
        """
        target = (self.focus_sequence if self.focus_sequence is not None
                  else self._standalone_focus_sequence)
        if (self.focus_worker is not None or self.slew_worker is not None or target is None
                or not target.exposures):
            return
        if target is self.focus_sequence:
            # Clear the curve immediately rather than letting old points
            # be replaced one at a time as each re-measured result
            # streams in -- otherwise the plot would show a confusing
            # mix of old and new measurements while reanalysis runs.
            self.window.curve_panel.reset()
        self._start_focus_worker(target, mode='reanalyze')

    def take_single_exposure(self, focus_value):
        """
        Take one exposure at ``focus_value`` with no sequence bookkeeping
        -- the Single tab's action, which doubles as "move to best focus"
        when ``focus_value`` is left at its default (the most recent
        fitted best focus; see
        `~gui.views.focus_control_panel.FocusControlPanel.show_best_focus`).
        Uses a throwaway `focus.FocusSequence` for its hardware handles,
        independent of any loaded `focus_sequence`.
        """
        if self.focus_worker is not None or self.slew_worker is not None:
            return  # hardware exclusivity: something is already running
        standalone = focus.FocusSequence()
        if standalone._focus is None or standalone._exposure is None:
            self.window.control_panel.show_failure(
                'Could not take single exposure: no ktl connection is available.')
            return
        self._standalone_focus_sequence = standalone
        exp_kwargs = self.window.control_panel.get_exposure_config()
        self._start_focus_worker(standalone, mode='single', exp_kwargs=exp_kwargs,
                                  focus_value=focus_value)

    def stop(self):
        """Request that the running sequence stop between steps (§4.3)."""
        if self.focus_worker is None:
            return
        self.focus_worker.request_stop()
        self.window.control_panel.set_stopping()

    def set_method(self, method):
        """
        Update the photometry method used for future steps/reanalysis.
        ``method`` is ``'brightest'`` (the default) or an ``(x, y)``
        coordinate tuple (from `ImagePanel.sourceSelected`); the
        coordinates actually measured are reported per-exposure on the
        Log tab (`~gui.views.focus_control_panel.FocusControlPanel.update_step`/
        `~gui.views.focus_control_panel.FocusControlPanel.show_single_exposure_result`)
        rather than tracked as a separate "current method" display.
        """
        self.method = method

    def move_to_target(self, ra_text, dec_text):
        """
        Slew the telescope to the Slew tab's target RA/Dec -- on a
        background thread (`SlewWorker`), since
        `slew.NickelTelescopePointing.slew_to` can block for up to five
        minutes waiting for the telescope to arrive.

        Parameters
        ----------
        ra_text : :obj:`str`
            Target right ascension, as entered in the Slew tab's target
            RA field (a bare sexagesimal string, e.g. ``'05:30:00'``).
        dec_text : :obj:`str`
            Target declination, entered the same way (e.g.
            ``'+20:15:00'``).
        """
        if self.focus_worker is not None or self.slew_worker is not None:
            return  # hardware exclusivity: something is already running
        if self.telescope is None:
            self.window.control_panel.show_failure(
                'Could not move to target: no ktl connection is available.')
            return
        try:
            target = SkyCoord(ra=ra_text, dec=dec_text, unit=('hourangle', 'deg'))
        except ValueError as e:
            self.window.control_panel.show_failure(f'Could not parse target coordinates: {e}')
            return

        self.slew_worker = SlewWorker(self.telescope, target.ra, target.dec)
        self.slew_worker.slewFinished.connect(self._on_slew_finished)
        self.slew_worker.slewFailed.connect(self._on_slew_failed)
        self.slew_worker.finished.connect(self._on_slew_worker_finished)
        self._set_running(True)
        self.slew_worker.start()

    def find_nearest_object(self, file_text, search_text):
        """
        Handle the Slew tab's "Find nearest object" button: search
        ``file_text`` for the target nearest the telescope's current
        position, restricted to names containing ``search_text``
        (unrestricted if empty).

        Parameters
        ----------
        file_text : :obj:`str`
            Path to a starlist file entered in the Slew tab's file
            field. Unlike `find_nearest_pointing`/`find_nearest_focus`,
            there is no default here -- use one of those two for the
            packaged catalog.
        search_text : :obj:`str`
            Object-name search string entered in the Slew tab's search
            field, or an empty string to consider every target in the
            file.

        Raises
        ------
        ValueError
            Raised if ``file_text`` is empty.
        """
        if file_text == '':
            raise ValueError('A starlist file must be provided to find the nearest object.')
        search = None if search_text == '' else search_text
        self._find_nearest(search, file=file_text)

    def find_nearest_pointing(self):
        """
        Handle the Slew tab's "Find nearest pointing star" button:
        search the packaged default catalog for the nearest target whose
        name contains ``'Pointing'``, ignoring whatever is currently in
        the Slew tab's file/search fields (matching
        ``slew_to_nearest.py -s Pointing`` with no ``-f``).
        """
        self._find_nearest('Pointing')

    def find_nearest_focus(self):
        """
        Handle the Slew tab's "Find nearest focus star" button: like
        `find_nearest_pointing`, but for names containing ``'Focusing'``
        (matching ``slew_to_nearest.py -s Focusing`` with no ``-f``).
        """
        self._find_nearest('Focusing')

    # -- internals ------------------------------------------------------------

    def _find_nearest(self, obj_search_str, file=None):
        """
        Shared implementation for `find_nearest_object`/
        `find_nearest_pointing`/`find_nearest_focus`: locate the nearest
        target and populate/report it on the Slew tab, or report a clear
        failure.

        Parameters
        ----------
        obj_search_str
            Forwarded to `slew.find_nearest_target`.
        file
            Forwarded to `slew.find_nearest_target`; ``None`` (the
            default) searches the packaged default catalog.
        """
        if self.telescope is None:
            self.window.control_panel.show_failure(
                'Could not find nearest target: no ktl connection is available.')
            return
        try:
            name, ra, dec = slew.find_nearest_target(
                self.telescope.current, obj_search_str=obj_search_str, file=file)
        except (ValueError, FileNotFoundError) as e:
            self.window.control_panel.show_failure(f'Could not find nearest target: {e}')
            return
        ra_text = ra.to_string(unit='hourangle', sep=':', pad=True, precision=2)
        dec_text = dec.to_string(unit='deg', sep=':', pad=True, alwayssign=True, precision=2)
        self.window.control_panel.show_nearest_target(name, ra_text, dec_text)

    def _start_focus_worker(self, focus_sequence, mode, exp_kwargs=None, focus_value=None):
        """
        Create a `FocusWorker` for ``focus_sequence``, connect its
        signals to this Controller's ``_on_*`` handlers, and start it
        running on its own background thread (see `FocusWorker.run`).
        """
        self.focus_worker = FocusWorker(focus_sequence, method=self.method, mode=mode,
                                         exp_kwargs=exp_kwargs, focus_value=focus_value)
        self.focus_worker.stepComplete.connect(self._on_step_complete)
        self.focus_worker.focusSequenceFinished.connect(self._on_focus_sequence_finished)
        self.focus_worker.focusSequenceFailed.connect(self._on_focus_sequence_failed)
        self.focus_worker.singleExposureFinished.connect(self._on_single_exposure_finished)
        self.focus_worker.finished.connect(self._on_focus_worker_finished)
        self._set_running(True)
        self.focus_worker.start()

    def _on_step_complete(self, result):
        """
        Handle one `focus.StepResult` from the worker's
        `FocusWorker.stepComplete` signal: update the image/curve
        panels and the Log tab's step display.
        """
        is_reanalyze = self.focus_worker is not None and self.focus_worker.mode == 'reanalyze'
        reanalyzing_standalone = (
            is_reanalyze and self._standalone_focus_sequence is not None
            and self.focus_worker.focus_sequence is self._standalone_focus_sequence
        )
        if is_reanalyze:
            self.window.image_panel.update_result(result)
            if not reanalyzing_standalone:
                # The curve was already cleared in reanalyze(), so each
                # re-measured result is simply added back in as it
                # arrives, rather than replacing an old point in place --
                # that would otherwise show a confusing mix of old and
                # new measurements while reanalysis is still in progress.
                # A standalone exposure never appears on the curve panel
                # at all, so there's nothing to add for it.
                self.window.curve_panel.add_result(result)
        else:
            self.window.image_panel.add_result(result)
            self.window.curve_panel.add_result(result)

        total = self.focus_sequence.expected_steps if self.focus_sequence is not None else None
        self.window.control_panel.update_step(result, total_expected=total)

    def _on_focus_sequence_finished(self, best_focus, best_fwhm):
        """Handle `FocusWorker.focusSequenceFinished`: report the fitted result."""
        self.window.control_panel.show_best_focus(best_focus, best_fwhm)

    def _on_focus_sequence_failed(self, message):
        """Handle `FocusWorker.focusSequenceFailed`: report the failure message."""
        self.window.control_panel.show_failure(message)

    def _on_single_exposure_finished(self, result):
        """
        Handle `FocusWorker.singleExposureFinished`: display the
        exposure, seed `_standalone_focus_sequence`'s bookkeeping with it
        (so `reanalyze()` has something to work with), and report it.
        """
        self.window.image_panel.add_result(result)
        seq = self._standalone_focus_sequence
        seq.observed_focus.append(result.focus_value)
        seq.exposures.append(result.exposure)
        seq.img_quality.append(result.fwhm)
        seq.source_stamps.append(result.stamp)
        seq.centroids.append(result.centroid)
        seq.step_iter = 1
        self.window.control_panel.show_single_exposure_result(result)

    def _on_focus_worker_finished(self):
        """
        Handle `FocusWorker.finished` (a signal every `QThread` provides
        automatically, emitted once :func:`~FocusWorker.run` returns):
        release the worker and restore hardware exclusivity.
        """
        if self.focus_worker is not None:
            # QThread.wait() blocks (briefly) until the OS thread has
            # actually terminated. `finished` already means `run()` has
            # returned, so this is just tidying up bookkeeping, not
            # waiting for real work -- but it's the documented, correct
            # way to confirm a QThread is truly done before discarding
            # it.
            self.focus_worker.wait()
        self.focus_worker = None
        self._set_running(False)

    def _on_source_selected(self, x, y):
        """Handle `ImagePanel.sourceSelected`: mark that source, then reanalyze against it."""
        self.set_method((x, y))
        self.reanalyze()

    def _on_slew_finished(self):
        """Handle `SlewWorker.slewFinished`: report that the telescope reached its target."""
        self.window.control_panel.show_slew_result('Move to target complete.')

    def _on_slew_failed(self, message):
        """Handle `SlewWorker.slewFailed`: report the failure message."""
        self.window.control_panel.show_failure(message)

    def _on_slew_worker_finished(self):
        """
        Handle `SlewWorker.finished` (a signal every `QThread` provides
        automatically, emitted once :func:`~SlewWorker.run` returns):
        release the worker and restore hardware exclusivity. See
        `_on_focus_worker_finished` for why `wait()` is called here.
        """
        if self.slew_worker is not None:
            self.slew_worker.wait()
        self.slew_worker = None
        self._set_running(False)

    def _update_current_position(self):
        """
        Poll `telescope.current` and push the formatted result to the
        Slew tab's live "current position" display. A no-op if
        `telescope` is ``None`` (see the class docstring); the
        `~PySide6.QtCore.QTimer` driving this is never started in that
        case, but this also guards the one direct call made from
        `__init__` before that timer's first tick.
        """
        if self.telescope is None:
            return
        current = self.telescope.current
        ra_text = current.ra.to_string(unit='hourangle', sep=':', pad=True, precision=2)
        dec_text = current.dec.to_string(unit='deg', sep=':', pad=True, alwayssign=True,
                                          precision=2)
        self.window.control_panel.set_current_position(ra_text, dec_text)

    def _set_running(self, running):
        """Propagate the hardware-exclusivity state to the View (§4.3)."""
        self.window.control_panel.set_running(running)
        self.window.image_panel.set_selection_enabled(not running)
