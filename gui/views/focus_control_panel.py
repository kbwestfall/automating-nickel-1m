"""
Focus-sequence configuration and controls.

See GUI_DESIGN.md §5.4 and claude/GUI_IMPLEMENTATION.md's Phase 4 for the
tabbed redesign. This view only exposes configuration and emits signals
for requested actions; it never constructs a `focus.FocusSequence` or
talks to `gui.model.sequence_worker.SequenceWorker` itself -- that's the
Controller's job.
"""
from pathlib import Path

from gui.qt import QtCore, QtWidgets


def _int_spin(minimum, maximum, value):
    """A :class:`~PySide6.QtWidgets.QSpinBox` with no increment/decrement buttons."""
    spin = QtWidgets.QSpinBox()
    spin.setRange(minimum, maximum)
    spin.setValue(value)
    spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
    return spin


def _float_spin(minimum, maximum, value):
    """A :class:`~PySide6.QtWidgets.QDoubleSpinBox` with no increment/decrement buttons."""
    spin = QtWidgets.QDoubleSpinBox()
    spin.setRange(minimum, maximum)
    spin.setValue(value)
    spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
    return spin


def _exposure_field_rows():
    """
    Build exposure-time/speed/binning fields, used identically by the
    Single, Grid, and Auto tabs (all of which take real exposures).

    Returns
    -------
    tuple
        ``(rows, exptime_spin, speed_combo, binning_combo)``, where
        ``rows`` is a list of ``(label, widget)`` pairs suitable for
        :func:`_two_column_row`.
    """
    exptime_spin = _float_spin(0.1, 600., 5.)
    speed_combo = QtWidgets.QComboBox()
    speed_combo.addItems(['Slow', 'Fast'])
    binning_combo = QtWidgets.QComboBox()
    binning_combo.addItems(['1,1', '2,2', '4,4'])
    rows = [
        ('Exposure time (s):', exptime_spin),
        ('Speed:', speed_combo),
        ('Binning:', binning_combo),
    ]
    return rows, exptime_spin, speed_combo, binning_combo


def _two_column_row(left_rows, right_rows):
    """
    Lay ``left_rows``/``right_rows`` (each a list of ``(label, widget)``
    pairs) out as two side-by-side :class:`~PySide6.QtWidgets.QFormLayout`\\
    s, so a tab's fields use horizontal space instead of stacking
    everything in one tall column -- exposure parameters on the left,
    the focus value(s) that define the sequence on the right.

    Returns
    -------
    :class:`~PySide6.QtWidgets.QHBoxLayout`
    """
    left_form = QtWidgets.QFormLayout()
    for label, widget in left_rows:
        left_form.addRow(label, widget)
    right_form = QtWidgets.QFormLayout()
    for label, widget in right_rows:
        right_form.addRow(label, widget)

    row = QtWidgets.QHBoxLayout()
    row.addLayout(left_form)
    row.addLayout(right_form)
    return row


class FocusControlPanel(QtWidgets.QWidget):
    """
    Sequence configuration/action tabs and Start/Stop -- plus live
    status/history on the Log tab and short usage reminders on the Help
    tab.

    One tab per action -- Single, Grid, Auto, Replay -- each showing only
    the fields relevant to it, with exposure parameters and the focus
    value(s) that define the sequence laid out in two side-by-side
    columns to keep each tab short (see :func:`_two_column_row`).
    Start/Stop are shared below the tabs rather than duplicated on each
    one, since they mean the same thing regardless of which tab is
    active: "Start" reads whichever tab is currently selected via
    :func:`get_sequence_type`/:func:`get_sequence_config` (or, for
    Single, emits :attr:`takeSingleExposureRequested` directly) and
    "Stop" always targets whatever's currently running.

    The Single tab doubles as "move to best focus": its focus field
    defaults to the most recent fitted best focus (via
    :func:`show_best_focus`), so acquiring a single exposure there with
    the default value *is* moving to best focus; changing the value
    first tests any other focus instead. There is no separate
    confirmation step -- reviewing/editing the value before pressing
    Start is the confirmation.

    There is no photometry-method selector: every measurement uses the
    brightest detected source by default, and the coordinates actually
    used (whether the automatic brightest source or one picked by
    clicking the image and pressing 'm' -- see `ImagePanel`) are
    reported on the Log tab via :func:`update_step`/
    :func:`show_single_exposure_result` rather than shown as a persistent
    "current method" state.

    Attributes
    ----------
    startRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when "Start" is clicked while Grid/Auto/Replay is active;
        the Controller should read :func:`get_sequence_type` and
        :func:`get_sequence_config` (and :func:`get_exposure_config` for
        Grid/Auto) to build and run the sequence.
    stopRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when "Stop" is clicked.
    takeSingleExposureRequested : :class:`~PySide6.QtCore.Signal`
        Emitted with a focus value when "Start" is clicked while Single
        is active.
    """
    startRequested = QtCore.Signal()
    stopRequested = QtCore.Signal()
    takeSingleExposureRequested = QtCore.Signal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._running = False

        self.single_tab = self._build_single_tab()
        self.grid_tab = self._build_grid_tab()
        self.auto_tab = self._build_auto_tab()
        self.replay_tab = self._build_replay_tab()
        self.log_tab = self._build_log_tab()
        self.help_tab = self._build_help_tab()

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.single_tab, 'Single')
        self.tabs.addTab(self.grid_tab, 'Grid')
        self.tabs.addTab(self.auto_tab, 'Auto')
        self.tabs.addTab(self.replay_tab, 'Replay')
        self.tabs.addTab(self.log_tab, 'Log')
        self.tabs.addTab(self.help_tab, 'Help')
        self.tabs.currentChanged.connect(self._on_tab_changed)

        # Every input widget across the four actionable tabs, locked
        # uniformly while anything is running (hardware exclusivity:
        # only one operation at a time, so there's no reason to allow
        # editing a different tab's configuration meanwhile). The tab
        # bar itself is deliberately never disabled, so the Log tab
        # stays reachable while a sequence runs.
        self._config_widgets = [
            self.single_focus_spin, self.single_exptime_spin, self.single_speed_combo,
            self.single_binning_combo,
            self.grid_start_spin, self.grid_step_spin, self.grid_nstep_spin,
            self.grid_exptime_spin, self.grid_speed_combo, self.grid_binning_combo,
            self.auto_start_spin, self.auto_step_spin, self.auto_maxsteps_spin,
            self.auto_exptime_spin, self.auto_speed_combo, self.auto_binning_combo,
            self.replay_datadir_edit, self.replay_browse_button, self.replay_prefix_edit,
            self.replay_suffix_edit, self.replay_obsnum_spin, self.replay_start_spin,
            self.replay_step_spin, self.replay_nstep_spin,
        ]

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.tabs)
        layout.addLayout(self._build_start_stop_row())

        self._update_start_button()

    # -- widget construction ----------------------------------------------

    def _build_single_tab(self):
        self.single_focus_spin = _int_spin(165, 500, 340)
        exposure_rows, self.single_exptime_spin, self.single_speed_combo, \
            self.single_binning_combo = _exposure_field_rows()
        focus_rows = [('Focus value:', self.single_focus_spin)]

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(_two_column_row(exposure_rows, focus_rows))
        layout.addStretch(1)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        return widget

    def _build_grid_tab(self):
        self.grid_start_spin = _int_spin(165, 500, 340)
        self.grid_step_spin = _int_spin(1, 100, 5)
        self.grid_nstep_spin = _int_spin(3, 100, 5)
        exposure_rows, self.grid_exptime_spin, self.grid_speed_combo, \
            self.grid_binning_combo = _exposure_field_rows()
        focus_rows = [
            ('Start focus:', self.grid_start_spin),
            ('Step size:', self.grid_step_spin),
            ('Number of steps:', self.grid_nstep_spin),
        ]

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(_two_column_row(exposure_rows, focus_rows))
        layout.addStretch(1)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        return widget

    def _build_auto_tab(self):
        self.auto_start_spin = _int_spin(165, 500, 340)
        self.auto_step_spin = _int_spin(1, 100, 5)
        self.auto_maxsteps_spin = _int_spin(2, 100, 12)
        exposure_rows, self.auto_exptime_spin, self.auto_speed_combo, \
            self.auto_binning_combo = _exposure_field_rows()
        focus_rows = [
            ('Start focus:', self.auto_start_spin),
            ('Step size:', self.auto_step_spin),
            ('Max steps:', self.auto_maxsteps_spin),
        ]

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(_two_column_row(exposure_rows, focus_rows))
        layout.addStretch(1)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        return widget

    def _build_replay_tab(self):
        self.replay_datadir_edit = QtWidgets.QLineEdit('.')
        self.replay_browse_button = QtWidgets.QPushButton('Browse…')
        self.replay_browse_button.clicked.connect(self._on_browse_clicked)
        self.replay_prefix_edit = QtWidgets.QLineEdit('n')
        self.replay_suffix_edit = QtWidgets.QLineEdit('.fits')
        self.replay_obsnum_spin = _int_spin(0, 999999, 0)
        self.replay_start_spin = _int_spin(165, 500, 340)
        self.replay_step_spin = _int_spin(1, 100, 5)
        self.replay_nstep_spin = _int_spin(3, 100, 5)

        datadir_row = QtWidgets.QHBoxLayout()
        datadir_row.addWidget(self.replay_datadir_edit, 1)
        datadir_row.addWidget(self.replay_browse_button)
        datadir_form = QtWidgets.QFormLayout()
        datadir_form.addRow('Data directory:', datadir_row)

        filename_rows = [
            ('Prefix:', self.replay_prefix_edit),
            ('Suffix:', self.replay_suffix_edit),
            ('Obsnum:', self.replay_obsnum_spin),
        ]
        focus_rows = [
            ('Start focus:', self.replay_start_spin),
            ('Step size:', self.replay_step_spin),
            ('Number of steps:', self.replay_nstep_spin),
        ]

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(datadir_form)
        layout.addLayout(_two_column_row(filename_rows, focus_rows))
        layout.addStretch(1)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        return widget

    def _build_log_tab(self):
        self.status_label = QtWidgets.QLabel('')
        self.step_label = QtWidgets.QLabel('Step: —')
        self.log_widget = QtWidgets.QPlainTextEdit()
        self.log_widget.setReadOnly(True)
        self.log_widget.setMaximumBlockCount(500)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.status_label)
        layout.addWidget(self.step_label)
        layout.addWidget(self.log_widget)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        return widget

    def _build_help_tab(self):
        text = (
            '<b>Single</b> — take one exposure at the given focus. '
            'Defaults to the most recent fitted best focus, so acquiring '
            'it as-is moves to best focus; change the value first to '
            'test any other focus.<br><br>'
            '<b>Grid</b> — step through an evenly spaced focus grid.<br>'
            '<b>Auto</b> — adaptively search for the best focus.<br>'
            '<b>Replay</b> — reprocess an existing set of exposures from disk.<br>'
            '<b>Log</b> — live status and history.<br><br>'
            '<b>Source selection</b> — the brightest star is measured '
            "automatically by default (see the Log for its coordinates); "
            "hover over a different star in the image and press 'm' to "
            'measure that one instead.<br><br>'
            '<b>Start / Stop</b> — run or stop the action configured in '
            'the current tab.'
        )
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)
        label.setTextFormat(QtCore.Qt.TextFormat.RichText)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addStretch(1)

        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        return widget

    def _build_start_stop_row(self):
        self.start_button = QtWidgets.QPushButton('Start')
        self.stop_button = QtWidgets.QPushButton('Stop')
        self.stop_button.setEnabled(False)
        self.start_button.clicked.connect(self._on_start_clicked)
        self.stop_button.clicked.connect(self.stopRequested.emit)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.start_button)
        row.addWidget(self.stop_button)
        return row

    # -- public API used by the Controller ---------------------------------

    def get_sequence_type(self):
        """
        The sequence type for the currently active tab, as ``'grid'``,
        ``'automated'``, or ``'archive'`` -- or ``None`` if Single, Log,
        or Help is active (in which case "Start" is disabled anyway).
        """
        tab = self.tabs.currentWidget()
        if tab is self.grid_tab:
            return 'grid'
        if tab is self.auto_tab:
            return 'automated'
        if tab is self.replay_tab:
            return 'archive'
        return None

    def get_sequence_config(self):
        """
        Return the currently active tab's configuration as a
        :obj:`dict`. Only includes keys relevant to that tab -- e.g. a
        Grid config has no ``datadir``, an Auto config has ``maxsteps``
        instead of ``nstep`` -- matching `focus.py`'s own CLI arguments.
        Empty for Single/Log/Help.
        """
        tab = self.tabs.currentWidget()
        if tab is self.grid_tab:
            return {
                'start': self.grid_start_spin.value(),
                'step': self.grid_step_spin.value(),
                'nstep': self.grid_nstep_spin.value(),
            }
        if tab is self.auto_tab:
            return {
                'start': self.auto_start_spin.value(),
                'step': self.auto_step_spin.value(),
                'maxsteps': self.auto_maxsteps_spin.value(),
            }
        if tab is self.replay_tab:
            return {
                'datadir': Path(self.replay_datadir_edit.text()),
                'prefix': self.replay_prefix_edit.text(),
                'suffix': self.replay_suffix_edit.text(),
                'obsnum': self.replay_obsnum_spin.value(),
                'start': self.replay_start_spin.value(),
                'step': self.replay_step_spin.value(),
                'nstep': self.replay_nstep_spin.value(),
            }
        return {}

    def get_exposure_config(self):
        """
        Return the currently active tab's exposure settings as a
        :obj:`dict` with keys ``exptime`` (:obj:`float`), ``speed``,
        ``binning`` (:obj:`str`) -- passed to `ExposureConfig.configure`
        before a live exposure. Empty for Replay/Log/Help, which take no
        new exposures.
        """
        tab = self.tabs.currentWidget()
        if tab is self.single_tab:
            exptime, speed, binning = (self.single_exptime_spin, self.single_speed_combo,
                                        self.single_binning_combo)
        elif tab is self.grid_tab:
            exptime, speed, binning = (self.grid_exptime_spin, self.grid_speed_combo,
                                        self.grid_binning_combo)
        elif tab is self.auto_tab:
            exptime, speed, binning = (self.auto_exptime_spin, self.auto_speed_combo,
                                        self.auto_binning_combo)
        else:
            return {}
        return {'exptime': exptime.value(), 'speed': speed.currentText(),
                'binning': binning.currentText()}

    def set_running(self, running):
        """Toggle widget states for whether an action is currently running."""
        self._running = running
        self.stop_button.setEnabled(running)
        for widget in self._config_widgets:
            widget.setEnabled(not running)

        if running:
            self.start_button.setEnabled(False)
        else:
            self._update_start_button()
            # Only clear a stale "Stopping..." message -- not a
            # confirmation/failure the worker's finishing signal may have
            # just set moments before this.
            if self.status_label.text().startswith('Stopping'):
                self.status_label.setText('')

    def set_stopping(self):
        """Reflect a Stop request that hasn't taken effect yet (§4.3)."""
        self.stop_button.setEnabled(False)
        self.status_label.setText('Stopping — waiting for current exposure to finish...')

    def update_step(self, result, total_expected=None):
        """Update the live step/status display for one new `focus.StepResult`."""
        text = f'Step {result.index + 1}'
        if total_expected:
            text += f'/{total_expected}'
        text += (f' — Focus {result.focus_value:.0f}, FWHM {result.fwhm:.2f}, '
                  f'Source ({result.centroid[0]:.1f}, {result.centroid[1]:.1f})')
        if result.is_outlier:
            text += '  [outlier]'
        self.step_label.setText(text)
        self.log_widget.appendPlainText(text)

    def show_best_focus(self, best_focus, best_fwhm):
        """
        Report a completed sequence's fitted result, and set it as the
        Single tab's default focus value -- rounded, since a focus value
        is a whole-unit telescope position even though the fit itself can
        land on a fractional one.
        """
        text = f'Sequence finished: best focus {best_focus:.1f}, expected FWHM {best_fwhm:.2f}'
        self.status_label.setText(text)
        self.log_widget.appendPlainText(text)
        self.single_focus_spin.setValue(round(best_focus))

    def show_single_exposure_result(self, result):
        """Report one exposure taken from the Single tab."""
        text = (f'Took exposure at focus {result.focus_value:.0f}: measured FWHM '
                f'{result.fwhm:.2f}, source ({result.centroid[0]:.1f}, {result.centroid[1]:.1f})')
        self.status_label.setText(text)
        self.log_widget.appendPlainText(text)

    def show_failure(self, message):
        """Display a sequence failure (e.g., from `SequenceWorker.sequenceFailed`)."""
        self.status_label.setText(message)
        self.log_widget.appendPlainText(f'ERROR: {message}')

    def reset(self):
        """Clear live status and log for a new run."""
        self.step_label.setText('Step: —')
        self.status_label.setText('')
        self.log_widget.clear()

    # -- internals ----------------------------------------------------------

    def _on_tab_changed(self, _index):
        if self._running:
            return  # set_running() alone governs Start/Stop while running
        self._update_start_button()

    def _update_start_button(self):
        tab = self.tabs.currentWidget()
        self.start_button.setEnabled(
            tab in (self.single_tab, self.grid_tab, self.auto_tab, self.replay_tab))

    def _on_start_clicked(self):
        if self.tabs.currentWidget() is self.single_tab:
            self.takeSingleExposureRequested.emit(self.single_focus_spin.value())
        else:
            self.startRequested.emit()

    def _on_browse_clicked(self):
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Select archive directory', self.replay_datadir_edit.text())
        if directory:
            self.replay_datadir_edit.setText(directory)
