"""
Wires the View (:class:`~nickel_focus.gui.views.main_window.MainWindow` and its
panels) to the Model (:class:`~nickel_focus.focus.FocusSequence` subclasses, via
:class:`~nickel_focus.gui.model.focus_worker.FocusWorker`; and
:class:`~nickel_focus.slew.NickelTelescopePointing`, via
:class:`~nickel_focus.gui.model.slew_worker.SlewWorker`).

See GUI_DESIGN.md §4.3.
"""
import logging

from astropy.coordinates import SkyCoord

from nickel_focus import focus
from nickel_focus import log
from nickel_focus import slew
from nickel_focus.gui.log_handler import QtLogHandler
from nickel_focus.gui.model.focus_worker import FocusWorker
from nickel_focus.gui.model.slew_worker import SlewWorker
from nickel_focus.gui.qt import QtCore
from nickel_focus.pkg.logger import GuiFormatter


class Controller(QtCore.QObject):
    """
    Owns the currently active :class:`~nickel_focus.focus.FocusSequence` (if
    any) and the :class:`~nickel_focus.gui.model.focus_worker.FocusWorker`
    driving it, and mediates between it and a
    :class:`~nickel_focus.gui.views.main_window.MainWindow`.

    Implements the "hardware exclusivity" state machine from §4.3: only one
    operation can be active at a time, and interactive source selection
    (:attr:`~nickel_focus.gui.views.image_panel.ImagePanel.sourceSelected`) is
    only enabled while nothing is running.  All three sequence types (Grid,
    Automated, Archive/ Replay), reanalysis, and the Single tab (which also
    serves as "move to best focus" -- see
    :class:`~nickel_focus.gui.views.focus_control_panel.FocusControlPanel`) are
    wired up.

    Also owns the :class:`~nickel_focus.slew.NickelTelescopePointing` handle
    used by the Slew tab:
    :meth:`~nickel_focus.gui.controller.Controller.move_to_target` (a slew, via
    :class:`~nickel_focus.gui.model.slew_worker.SlewWorker`) participates in the
    same hardware-exclusivity state as a focus sequence -- the telescope
    shouldn't move while it's mid-exposure, or vice versa -- while
    :meth:`~nickel_focus.gui.controller.Controller.find_nearest_object`/
    :meth:`~nickel_focus.gui.controller.Controller.find_nearest_pointing`/
    :meth:`~nickel_focus.gui.controller.Controller.find_nearest_focus` (a single
    fast ktl read plus a starlist search, not a move) do not.  A
    :class:`PySide6.QtCore.QTimer` polls the telescope's current position for
    the tab's live display regardless of what else is running.

    Without a live ``ktl`` connection, the Slew/Single/Grid/Auto tabs can't do
    anything but fail -- see
    :attr:`~nickel_focus.gui.controller.Controller.ktl_available` and
    :meth:`~nickel_focus.gui.views.focus_control_panel.FocusControlPanel.set_hardware_tabs_enabled`.

    Parameters
    ----------
    window : :class:`~nickel_focus.gui.views.main_window.MainWindow`
        The window to wire up.
    parent : :class:`PySide6.QtCore.QObject`, optional
        Passed through to the base :class:`PySide6.QtCore.QObject` constructor.
    force_enable_hardware_tabs : :obj:`bool`, optional
        Leave the Slew/Single/Grid/Auto tabs enabled even with no live ``ktl``
        connection, so a developer can still see and click through them while
        working on the GUI. Every action on those tabs still fails exactly as it
        would otherwise (each checks
        :attr:`~nickel_focus.gui.controller.Controller.telescope`/
        :attr:`~nickel_focus.focus.ktl` itself) -- this only skips graying them
        out.  Wired to the GUI script's suppressed ``--test`` flag; never set
        outside of that.
    """

    def __init__(self, window, parent=None, force_enable_hardware_tabs=False):
        super().__init__(parent)
        self.window = window

        # `log` is a module-level singleton, not per-Controller, so a new
        # Controller must remove any handler a previous one added before
        # adding its own -- otherwise a stale handler pointing at a
        # since-destroyed log_widget accumulates each time a Controller is
        # (re)constructed.
        log.remove_handlers_of_type(QtLogHandler)
        self._log_handler = QtLogHandler(self)
        self._log_handler.setFormatter(GuiFormatter())
        self._log_handler.setLevel(logging.INFO)
        self._log_handler.record_logged.connect(window.control_panel.append_log_line)
        log.addHandler(self._log_handler)

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

        # `focus.ktl` and `slew.ktl` are the same object in normal
        # operation (see `nickel_focus/pkg/ktl.py`), so checking one is
        # enough to know whether any of the four ktl-driven tabs can do
        # anything useful; `force_enable_hardware_tabs` overrides only
        # the View-level graying-out below, not this flag itself.
        self.ktl_available = focus.ktl is not None
        window.control_panel.set_hardware_tabs_enabled(
            self.ktl_available or force_enable_hardware_tabs)

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
        """
        Build a sequence from the control panel's configuration and run it.
        """
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
            self._fail(f'Could not start sequence: {e}')
            return

        if focus_sequence_type != 'archive' and (
                focus_sequence._focus is None or focus_sequence._exposure is None):
            self._fail(
                'Could not start sequence: no ktl connection is available for a live sequence.'
            )
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
        """
        Build a :class:`~nickel_focus.focus.ArchiveFocusSequence` from the
        Replay tab's ``config``, or raise.
        """
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
        Re-run photometry on the already-collected exposures of whatever is
        currently loaded: the active
        :attr:`~nickel_focus.gui.controller.Controller.focus_sequence`, or -- if
        none is loaded -- a standalone Single-tab exposure (§5.6's "no sequence
        loaded yet" case).
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
        Take one exposure at ``focus_value`` with no sequence bookkeeping -- the
        Single tab's action, which doubles as "move to best focus" when
        ``focus_value`` is left at its default (the most recent fitted best
        focus; see
        :meth:`~nickel_focus.gui.views.focus_control_panel.FocusControlPanel.show_best_focus`).
        Uses a throwaway :class:`~nickel_focus.focus.FocusSequence` for its
        hardware handles, independent of any loaded
        :attr:`~nickel_focus.gui.controller.Controller.focus_sequence`.
        """
        if self.focus_worker is not None or self.slew_worker is not None:
            return  # hardware exclusivity: something is already running
        standalone = focus.FocusSequence()
        if standalone._focus is None or standalone._exposure is None:
            self._fail('Could not take single exposure: no ktl connection is available.')
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
        ``method`` is ``'brightest'`` (the default) or an ``(x, y)`` coordinate
        tuple (from
        :attr:`~nickel_focus.gui.views.image_panel.ImagePanel.sourceSelected`);
        the coordinates actually measured are reported per-exposure on the Log
        tab
        (:meth:`~nickel_focus.gui.views.focus_control_panel.FocusControlPanel.update_step`/
        :meth:`~nickel_focus.gui.views.focus_control_panel.FocusControlPanel.show_single_exposure_result`)
        rather than tracked as a separate "current method" display.
        """
        self.method = method

    def move_to_target(self, ra_text, dec_text):
        """
        Slew the telescope to the Slew tab's target RA/Dec -- on a background
        thread (:class:`~nickel_focus.gui.model.slew_worker.SlewWorker`), since
        :meth:`~nickel_focus.slew.NickelTelescopePointing.slew_to` can block for
        up to five minutes waiting for the telescope to arrive.

        Parameters
        ----------
        ra_text : :obj:`str`
            Target right ascension, as entered in the Slew tab's target RA field
            (a bare sexagesimal string, e.g.  ``'05:30:00'``).
        dec_text : :obj:`str`
            Target declination, entered the same way (e.g.  ``'+20:15:00'``).
        """
        if self.focus_worker is not None or self.slew_worker is not None:
            return  # hardware exclusivity: something is already running
        if self.telescope is None:
            self._fail('Could not move to target: no ktl connection is available.')
            return
        try:
            target = SkyCoord(ra=ra_text, dec=dec_text, unit=('hourangle', 'deg'))
        except ValueError as e:
            self._fail(f'Could not parse target coordinates: {e}')
            return

        self.slew_worker = SlewWorker(self.telescope, target.ra, target.dec)
        self.slew_worker.slewFinished.connect(self._on_slew_finished)
        self.slew_worker.slewFailed.connect(self._fail)
        self.slew_worker.finished.connect(self._on_slew_worker_finished)
        self._set_running(True)
        self.slew_worker.start()

    def find_nearest_object(self, file_text, search_text):
        """
        Handle the Slew tab's "Find nearest object" button: search ``file_text``
        for the target nearest the telescope's current position, restricted to
        names containing ``search_text`` (unrestricted if empty).

        Parameters
        ----------
        file_text : :obj:`str`
            Path to a starlist file entered in the Slew tab's file field.
            Unlike
            :meth:`~nickel_focus.gui.controller.Controller.find_nearest_pointing`/
            :meth:`~nickel_focus.gui.controller.Controller.find_nearest_focus`,
            there is no default here -- use one of those two for the packaged
            catalog.
        search_text : :obj:`str`
            Object-name search string entered in the Slew tab's search field, or
            an empty string to consider every target in the file.

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
        Handle the Slew tab's "Find nearest pointing star" button: search the
        packaged default catalog for the nearest target whose name contains
        ``'Pointing'``, ignoring whatever is currently in the Slew tab's
        file/search fields (matching ``slew_to_nearest.py -s Pointing`` with no
        ``-f``).
        """
        self._find_nearest('Pointing')

    def find_nearest_focus(self):
        """
        Handle the Slew tab's "Find nearest focus star" button: like
        :meth:`~nickel_focus.gui.controller.Controller.find_nearest_pointing`,
        but for names containing ``'Focusing'`` (matching ``slew_to_nearest.py
        -s Focusing`` with no ``-f``).
        """
        self._find_nearest('Focusing')

    # -- internals ------------------------------------------------------------

    def _fail(self, message):
        """
        Log ``message`` as an error and report it as a failure on the Log tab.
        """
        log.error(message)
        self.window.control_panel.show_failure(message)

    def _find_nearest(self, obj_search_str, file=None):
        """
        Shared implementation for
        :meth:`~nickel_focus.gui.controller.Controller.find_nearest_object`/
        :meth:`~nickel_focus.gui.controller.Controller.find_nearest_pointing`/
        :meth:`~nickel_focus.gui.controller.Controller.find_nearest_focus`:
        locate the nearest target and populate/report it on the Slew tab, or
        report a clear failure.

        Parameters
        ----------
        obj_search_str
            Forwarded to :func:`~nickel_focus.slew.find_nearest_target`.
        file
            Forwarded to :func:`~nickel_focus.slew.find_nearest_target`;
            ``None`` (the default) searches the packaged default catalog.
        """
        if self.telescope is None:
            self._fail('Could not find nearest target: no ktl connection is available.')
            return
        try:
            name, ra, dec = slew.find_nearest_target(
                self.telescope.current, obj_search_str=obj_search_str, file=file)
        except (ValueError, FileNotFoundError) as e:
            self._fail(f'Could not find nearest target: {e}')
            return
        ra_text = ra.to_string(unit='hourangle', sep=':', pad=True, precision=2)
        dec_text = dec.to_string(unit='deg', sep=':', pad=True, alwayssign=True, precision=2)
        log.info(f'Nearest target: {name} (RA={ra_text}, Dec={dec_text})')
        self.window.control_panel.show_nearest_target(name, ra_text, dec_text)

    def _start_focus_worker(self, focus_sequence, mode, exp_kwargs=None, focus_value=None):
        """
        Create a :class:`~nickel_focus.gui.model.focus_worker.FocusWorker` for
        ``focus_sequence``, connect its signals to this Controller's ``_on_*``
        handlers, and start it running on its own background thread (see
        :meth:`~nickel_focus.gui.model.focus_worker.FocusWorker.run`).
        """
        self.focus_worker = FocusWorker(focus_sequence, method=self.method, mode=mode,
                                         exp_kwargs=exp_kwargs, focus_value=focus_value)
        self.focus_worker.stepComplete.connect(self._on_step_complete)
        self.focus_worker.focusSequenceFinished.connect(self._on_focus_sequence_finished)
        self.focus_worker.focusSequenceFailed.connect(self._fail)
        self.focus_worker.singleExposureFinished.connect(self._on_single_exposure_finished)
        self.focus_worker.finished.connect(self._on_focus_worker_finished)
        self._set_running(True)
        self.focus_worker.start()

    def _on_step_complete(self, result):
        """
        Handle one :class:`~nickel_focus.focus.StepResult` from the worker's
        :attr:`~nickel_focus.gui.model.focus_worker.FocusWorker.stepComplete`
        signal: update the image/curve panels and the Log tab's step display.
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
        step_text = f'Step {result.index + 1}'
        if total:
            step_text += f'/{total}'
        step_text += (f' — Focus {result.focus_value:.0f}, FWHM {result.fwhm:.2f}, '
                      f'Source ({result.centroid[0]:.1f}, {result.centroid[1]:.1f})')
        if result.is_outlier:
            step_text += '  [outlier]'
        log.info(step_text)
        self.window.control_panel.update_step(result, total_expected=total)

    def _on_focus_sequence_finished(self, best_focus, best_fwhm):
        """
        Handle
        :attr:`~nickel_focus.gui.model.focus_worker.FocusWorker.focusSequenceFinished`:
        report the fitted result.
        """
        log.info(f'Sequence finished: best focus {best_focus:.1f}, expected FWHM {best_fwhm:.2f}')
        self.window.control_panel.show_best_focus(best_focus, best_fwhm)

    def _on_single_exposure_finished(self, result):
        """
        Handle
        :attr:`~nickel_focus.gui.model.focus_worker.FocusWorker.singleExposureFinished`:
        display the exposure, seed
        :attr:`~nickel_focus.gui.controller.Controller._standalone_focus_sequence`'s
        bookkeeping with it (so
        :meth:`~nickel_focus.gui.controller.Controller.reanalyze` has something
        to work with), and report it.
        """
        self.window.image_panel.add_result(result)
        seq = self._standalone_focus_sequence
        seq.observed_focus.append(result.focus_value)
        seq.exposures.append(result.exposure)
        seq.img_quality.append(result.fwhm)
        seq.source_stamps.append(result.stamp)
        seq.centroids.append(result.centroid)
        seq.step_iter = 1
        log.info(f'Took exposure at focus {result.focus_value:.0f}: measured FWHM '
                 f'{result.fwhm:.2f}, source ({result.centroid[0]:.1f}, {result.centroid[1]:.1f})')
        self.window.control_panel.show_single_exposure_result(result)

    def _on_focus_worker_finished(self):
        """
        Handle :attr:`~nickel_focus.gui.model.focus_worker.FocusWorker.finished`
        (a signal every :class:`PySide6.QtCore.QThread` provides automatically,
        emitted once
        :meth:`~nickel_focus.gui.model.focus_worker.FocusWorker.run` returns):
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
        """
        Handle
        :attr:`~nickel_focus.gui.views.image_panel.ImagePanel.sourceSelected`:
        mark that source, then reanalyze against it.
        """
        self.set_method((x, y))
        self.reanalyze()

    def _on_slew_finished(self):
        """
        Handle
        :attr:`~nickel_focus.gui.model.slew_worker.SlewWorker.slewFinished`:
        report that the telescope reached its target.
        """
        log.info('Move to target complete.')
        self.window.control_panel.show_slew_result('Move to target complete.')

    def _on_slew_worker_finished(self):
        """
        Handle :attr:`~nickel_focus.gui.model.slew_worker.SlewWorker.finished`
        (a signal every :class:`PySide6.QtCore.QThread` provides automatically,
        emitted once :meth:`~nickel_focus.gui.model.slew_worker.SlewWorker.run`
        returns): release the worker and restore hardware exclusivity.  See
        :meth:`~nickel_focus.gui.controller.Controller._on_focus_worker_finished`
        for why :meth:`PySide6.QtCore.QThread.wait` is called here.
        """
        if self.slew_worker is not None:
            self.slew_worker.wait()
        self.slew_worker = None
        self._set_running(False)

    def _update_current_position(self):
        """
        Poll :attr:`~nickel_focus.slew.NickelTelescopePointing.current` and push
        the formatted result to the Slew tab's live "current position" display.
        A no-op if :attr:`~nickel_focus.gui.controller.Controller.telescope` is
        ``None`` (see the class docstring); the :class:`PySide6.QtCore.QTimer`
        driving this is never started in that case, but this also guards the one
        direct call made from
        :meth:`~nickel_focus.gui.controller.Controller.__init__` before that
        timer's first tick.
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
