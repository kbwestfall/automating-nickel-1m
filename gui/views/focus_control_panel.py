"""
Focus-sequence configuration and controls.

See GUI_DESIGN.md §5.4. This view only exposes configuration and emits
signals for requested actions; it never constructs a `focus.FocusSequence`
or talks to `gui.model.sequence_worker.SequenceWorker` itself -- that's
the Controller's job (sub-phase 7).
"""
from pathlib import Path

from gui.qt import QtCore, QtWidgets

_NO_MOVE_TOOLTIP = 'Only available after a live (Grid/Automated) sequence finishes'
_NO_PENDING_TOOLTIP = 'Take a single exposure first, with a sequence already loaded'


class FocusControlPanel(QtWidgets.QWidget):
    """
    Sequence configuration, Start/Stop, live status, and results.

    All three sequence types -- Archive/Replay, Grid, and Automated -- are
    selectable; which configuration fields are enabled follows the
    selected type (:func:`_on_sequence_type_changed`): Archive needs the
    data directory/prefix/suffix/obsnum fields and exposure settings are
    irrelevant (no new exposures are taken); Grid/Automated need exposure
    settings and have no archive fields to fill in; Grid uses "number of
    steps," Automated uses "max steps" instead.

    Attributes
    ----------
    startRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when "Start" is clicked; the Controller should read
        :func:`get_sequence_type` and :func:`get_sequence_config` (and,
        for a live sequence type, :func:`get_exposure_config`) to build
        and run the sequence.
    stopRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when "Stop" is clicked.
    methodChanged : :class:`~PySide6.QtCore.Signal`
        Emitted with ``'brightest'`` or ``'weighted'`` when the user picks
        one from the method combo box.
    moveToBestFocusRequested : :class:`~PySide6.QtCore.Signal`
        Emitted with the target focus value once the user confirms the
        "Move to best focus" dialog. Only enabled after a sequence
        finishes with a hardware connection to move (see
        :func:`show_best_focus`'s ``can_move`` argument) -- archive/replay
        sequences have no hardware to move.
    takeSingleExposureRequested : :class:`~PySide6.QtCore.Signal`
        Emitted with a focus value when "Take Single Exposure" is
        clicked (§5.5) -- a standalone exposure with no sequence
        bookkeeping, e.g. to confirm the field or mark a source before a
        real sequence starts.
    addToSequenceRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when "Add to Existing Sequence" is clicked, to commit the
        pending single exposure (see :func:`show_pending_exposure`) into
        the currently loaded sequence's data. Starting a new sequence
        instead (§5.5) needs no separate signal -- it's just
        ``startRequested`` -- the pending exposure is simply discarded.
    """
    startRequested = QtCore.Signal()
    stopRequested = QtCore.Signal()
    methodChanged = QtCore.Signal(str)
    moveToBestFocusRequested = QtCore.Signal(float)
    takeSingleExposureRequested = QtCore.Signal(float)
    addToSequenceRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._best_focus = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self._build_sequence_type_group())
        layout.addWidget(self._build_sequence_config_group())
        layout.addWidget(self._build_exposure_settings_group())
        layout.addWidget(self._build_method_group())
        layout.addLayout(self._build_start_stop_row())
        layout.addWidget(self._build_single_exposure_group())
        layout.addWidget(self._build_status_group())
        layout.addWidget(self._build_result_group())
        layout.addStretch(1)

        for radio in (self.archive_radio, self.grid_radio, self.automated_radio):
            radio.toggled.connect(lambda checked: self._on_sequence_type_changed())
        self._on_sequence_type_changed()

    # -- widget construction ----------------------------------------------

    def _build_sequence_type_group(self):
        group = QtWidgets.QGroupBox('Sequence Type')
        self.archive_radio = QtWidgets.QRadioButton('Archive / Replay')
        self.grid_radio = QtWidgets.QRadioButton('Grid')
        self.automated_radio = QtWidgets.QRadioButton('Automated')
        self.archive_radio.setChecked(True)

        self.sequence_type_buttons = QtWidgets.QButtonGroup(self)
        for radio in (self.archive_radio, self.grid_radio, self.automated_radio):
            self.sequence_type_buttons.addButton(radio)

        row = QtWidgets.QHBoxLayout(group)
        row.addWidget(self.archive_radio)
        row.addWidget(self.grid_radio)
        row.addWidget(self.automated_radio)
        return group

    def _build_sequence_config_group(self):
        group = QtWidgets.QGroupBox('Sequence Configuration')

        self.datadir_edit = QtWidgets.QLineEdit('.')
        self.browse_button = QtWidgets.QPushButton('Browse…')
        self.browse_button.clicked.connect(self._on_browse_clicked)
        self.prefix_edit = QtWidgets.QLineEdit('n')
        self.suffix_edit = QtWidgets.QLineEdit('.fits')
        self.obsnum_spin = QtWidgets.QSpinBox()
        self.obsnum_spin.setRange(0, 999999)
        self.obsnum_spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        # Focus values and step sizes are both expressed in the same
        # physical telescope-position units (whole units), even though a
        # best-fit focus (§ move-to-best-focus) can land on a fractional
        # value -- that gets rounded before it's ever offered here, so
        # these entry fields are always integers.
        self.start_spin = QtWidgets.QSpinBox()
        self.start_spin.setRange(165, 500)
        self.start_spin.setValue(340)
        self.start_spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.step_spin = QtWidgets.QSpinBox()
        self.step_spin.setRange(1, 100)
        self.step_spin.setValue(5)
        self.step_spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.nstep_spin = QtWidgets.QSpinBox()
        self.nstep_spin.setRange(3, 100)
        self.nstep_spin.setValue(5)
        self.nstep_spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.maxsteps_spin = QtWidgets.QSpinBox()
        self.maxsteps_spin.setRange(2, 100)
        self.maxsteps_spin.setValue(12)
        self.maxsteps_spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)

        form = QtWidgets.QFormLayout(group)
        datadir_row = QtWidgets.QHBoxLayout()
        datadir_row.addWidget(self.datadir_edit, 1)
        datadir_row.addWidget(self.browse_button)
        form.addRow('Data directory:', datadir_row)
        form.addRow('Prefix:', self.prefix_edit)
        form.addRow('Suffix:', self.suffix_edit)
        form.addRow('Obsnum:', self.obsnum_spin)
        form.addRow('Start focus:', self.start_spin)
        form.addRow('Step size:', self.step_spin)
        form.addRow('Number of steps (Archive/Grid):', self.nstep_spin)
        form.addRow('Max steps (Automated):', self.maxsteps_spin)
        return group

    def _build_exposure_settings_group(self):
        group = QtWidgets.QGroupBox('Exposure Settings')
        self.exptime_spin = QtWidgets.QDoubleSpinBox()
        self.exptime_spin.setRange(0.1, 600.)
        self.exptime_spin.setValue(5.)
        self.exptime_spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.speed_combo = QtWidgets.QComboBox()
        self.speed_combo.addItems(['Slow', 'Fast'])
        self.binning_combo = QtWidgets.QComboBox()
        self.binning_combo.addItems(['1,1', '2,2', '4,4'])

        form = QtWidgets.QFormLayout(group)
        form.addRow('Exposure time (s):', self.exptime_spin)
        form.addRow('Speed:', self.speed_combo)
        form.addRow('Binning:', self.binning_combo)
        return group

    def _build_method_group(self):
        group = QtWidgets.QGroupBox('Photometry Method')
        self.method_combo = QtWidgets.QComboBox()
        self.method_combo.addItems(['Brightest', 'Weighted'])
        self.method_combo.currentTextChanged.connect(
            lambda text: self.methodChanged.emit(text.lower()))
        self.method_label = QtWidgets.QLabel('Current method: Brightest')

        col = QtWidgets.QVBoxLayout(group)
        col.addWidget(self.method_combo)
        col.addWidget(self.method_label)
        return group

    def _build_start_stop_row(self):
        self.start_button = QtWidgets.QPushButton('Start')
        self.stop_button = QtWidgets.QPushButton('Stop')
        self.stop_button.setEnabled(False)
        self.start_button.clicked.connect(self.startRequested.emit)
        self.stop_button.clicked.connect(self.stopRequested.emit)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.start_button)
        row.addWidget(self.stop_button)
        return row

    def _build_single_exposure_group(self):
        group = QtWidgets.QGroupBox('Single Exposure')
        self.single_focus_spin = QtWidgets.QSpinBox()
        self.single_focus_spin.setRange(165, 500)
        self.single_focus_spin.setValue(340)
        self.single_focus_spin.setButtonSymbols(QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.take_single_exposure_button = QtWidgets.QPushButton('Take Single Exposure')
        self.take_single_exposure_button.clicked.connect(
            lambda: self.takeSingleExposureRequested.emit(self.single_focus_spin.value()))

        self.pending_label = QtWidgets.QLabel('No pending exposure')
        self.add_to_sequence_button = QtWidgets.QPushButton('Add to Existing Sequence')
        self.add_to_sequence_button.setEnabled(False)
        self.add_to_sequence_button.setToolTip(_NO_PENDING_TOOLTIP)
        self.add_to_sequence_button.clicked.connect(self.addToSequenceRequested.emit)

        form = QtWidgets.QFormLayout()
        form.addRow('Focus value:', self.single_focus_spin)
        col = QtWidgets.QVBoxLayout(group)
        col.addLayout(form)
        col.addWidget(self.take_single_exposure_button)
        col.addWidget(self.pending_label)
        col.addWidget(self.add_to_sequence_button)
        return group

    def _build_status_group(self):
        group = QtWidgets.QGroupBox('Status')
        self.status_label = QtWidgets.QLabel('')
        self.step_label = QtWidgets.QLabel('Step: —')
        self.log_widget = QtWidgets.QPlainTextEdit()
        self.log_widget.setReadOnly(True)
        self.log_widget.setMaximumBlockCount(500)

        col = QtWidgets.QVBoxLayout(group)
        col.addWidget(self.status_label)
        col.addWidget(self.step_label)
        col.addWidget(self.log_widget)
        return group

    def _build_result_group(self):
        group = QtWidgets.QGroupBox('Result')
        self.result_label = QtWidgets.QLabel('Best focus: —')
        self.move_to_best_focus_button = QtWidgets.QPushButton('Move to Best Focus')
        self.move_to_best_focus_button.setEnabled(False)
        self.move_to_best_focus_button.setToolTip(_NO_MOVE_TOOLTIP)
        self.move_to_best_focus_button.clicked.connect(self._on_move_to_best_focus_clicked)

        col = QtWidgets.QVBoxLayout(group)
        col.addWidget(self.result_label)
        col.addWidget(self.move_to_best_focus_button)
        return group

    # -- public API used by the Controller ---------------------------------

    def get_sequence_config(self):
        """
        Return the current sequence configuration as a :obj:`dict`. Not
        every key is relevant to every sequence type -- the Controller
        picks out what it needs based on :func:`get_sequence_type`.

        Keys: ``datadir`` (:class:`pathlib.Path`), ``prefix``, ``suffix``
        (:obj:`str`, Archive only), ``obsnum`` (:obj:`int`, Archive
        only), ``start``, ``step`` (:obj:`int`, focus values, all
        types), ``nstep`` (:obj:`int`, Archive/Grid), and ``maxsteps``
        (:obj:`int`, Automated only) -- matching `focus.py`'s own CLI
        arguments.
        """
        return {
            'datadir': Path(self.datadir_edit.text()),
            'prefix': self.prefix_edit.text(),
            'suffix': self.suffix_edit.text(),
            'obsnum': self.obsnum_spin.value(),
            'start': self.start_spin.value(),
            'step': self.step_spin.value(),
            'nstep': self.nstep_spin.value(),
            'maxsteps': self.maxsteps_spin.value(),
        }

    def get_sequence_type(self):
        """The selected sequence type, as ``'archive'``, ``'grid'``, or ``'automated'``."""
        if self.grid_radio.isChecked():
            return 'grid'
        if self.automated_radio.isChecked():
            return 'automated'
        return 'archive'

    def get_exposure_config(self):
        """
        Return the current exposure settings as a :obj:`dict` with keys
        ``exptime`` (:obj:`float`), ``speed``, ``binning`` (:obj:`str`) --
        passed to `ExposureConfig.configure` before a live sequence
        steps. Not applicable in archive/replay mode.
        """
        return {
            'exptime': self.exptime_spin.value(),
            'speed': self.speed_combo.currentText(),
            'binning': self.binning_combo.currentText(),
        }

    def get_selected_method(self):
        """The method combo's current selection, as ``'brightest'`` or ``'weighted'``."""
        return self.method_combo.currentText().lower()

    def set_method(self, method):
        """
        Update the read-only "current method" display. ``method`` is
        ``'brightest'``, ``'weighted'``, or an ``(x, y)`` coordinate tuple
        (set after an image click is selected; see `ImagePanel.sourceSelected`).
        """
        if isinstance(method, tuple):
            text = f'Selected source ({method[0]:.1f}, {method[1]:.1f})'
        else:
            text = str(method).capitalize()
        self.method_label.setText(f'Current method: {text}')

    def set_running(self, running):
        """Toggle widget states for whether a sequence is currently running."""
        self.start_button.setEnabled(not running)
        self.stop_button.setEnabled(running)
        for radio in (self.archive_radio, self.grid_radio, self.automated_radio):
            radio.setEnabled(not running)
        for widget in (self.start_spin, self.step_spin, self.method_combo,
                       self.single_focus_spin, self.take_single_exposure_button):
            widget.setEnabled(not running)

        if running:
            # Force every type-specific field off; which ones matter
            # depends on the sequence type, and none of them should be
            # editable mid-run regardless.
            for widget in (self.datadir_edit, self.browse_button, self.prefix_edit,
                           self.suffix_edit, self.obsnum_spin, self.nstep_spin,
                           self.maxsteps_spin, self.exptime_spin, self.speed_combo,
                           self.binning_combo, self.move_to_best_focus_button,
                           self.add_to_sequence_button):
                widget.setEnabled(False)
        else:
            # Restore the correct per-type enabled state, not just "all on".
            self._on_sequence_type_changed()
            # Only clear a stale "Stopping..." message -- not a
            # confirmation/failure the worker's finishing signal may have
            # just set moments before this (e.g. show_confirmation() from
            # a completed "Move to Best Focus").
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
        text += f' — Focus {result.focus_value:.0f}, FWHM {result.fwhm:.2f}'
        if result.is_outlier:
            text += '  [outlier]'
        self.step_label.setText(text)
        self.log_widget.appendPlainText(text)

    def show_best_focus(self, best_focus, best_fwhm, can_move=False):
        """
        Display a completed sequence's fitted result. ``can_move`` enables
        "Move to Best Focus" -- the Controller sets it based on whether
        the sequence that produced this result has a hardware connection
        to move (archive/replay sequences never do).
        """
        self._best_focus = best_focus
        self.result_label.setText(f'Best focus: {best_focus:.1f}   Expected FWHM: {best_fwhm:.2f}')
        self.log_widget.appendPlainText(
            f'Sequence finished: best focus {best_focus:.1f}, expected FWHM {best_fwhm:.2f}')
        self.move_to_best_focus_button.setEnabled(can_move)
        self.move_to_best_focus_button.setToolTip('' if can_move else _NO_MOVE_TOOLTIP)

    def show_confirmation(self, result):
        """Report the confirmation exposure taken by "Move to Best Focus"."""
        text = f'Moved to focus {result.focus_value:.0f}: measured FWHM {result.fwhm:.2f}'
        self.status_label.setText(text)
        self.log_widget.appendPlainText(text)

    def show_failure(self, message):
        """Display a sequence failure (e.g., from `SequenceWorker.sequenceFailed`)."""
        self.status_label.setText(message)
        self.log_widget.appendPlainText(f'ERROR: {message}')

    def show_pending_exposure(self, result, can_add):
        """
        Report a standalone single exposure (§5.5) awaiting a decision --
        "Add to Existing Sequence" or discard by starting a new one.
        ``can_add`` reflects whether a sequence is actually loaded to add
        it to.
        """
        self.pending_label.setText(
            f'Pending: Focus {result.focus_value:.0f}, FWHM {result.fwhm:.2f}')
        self.log_widget.appendPlainText(
            f'Took single exposure at focus {result.focus_value:.0f}: FWHM {result.fwhm:.2f}')
        self.add_to_sequence_button.setEnabled(can_add)
        self.add_to_sequence_button.setToolTip('' if can_add else _NO_PENDING_TOOLTIP)

    def clear_pending_exposure(self):
        """Clear the pending single-exposure display, e.g. once committed or discarded."""
        self.pending_label.setText('No pending exposure')
        self.add_to_sequence_button.setEnabled(False)
        self.add_to_sequence_button.setToolTip(_NO_PENDING_TOOLTIP)

    def reset(self):
        """Clear live status, log, and result display for a new run."""
        self._best_focus = None
        self.step_label.setText('Step: —')
        self.result_label.setText('Best focus: —')
        self.status_label.setText('')
        self.log_widget.clear()
        self.move_to_best_focus_button.setEnabled(False)
        self.move_to_best_focus_button.setToolTip(_NO_MOVE_TOOLTIP)
        self.clear_pending_exposure()

    # -- internals ----------------------------------------------------------

    def _on_sequence_type_changed(self):
        """Enable/disable configuration fields to match the selected sequence type."""
        is_archive = self.archive_radio.isChecked()
        is_automated = self.automated_radio.isChecked()

        for widget in (self.datadir_edit, self.browse_button, self.prefix_edit,
                       self.suffix_edit, self.obsnum_spin):
            widget.setEnabled(is_archive)
        self.nstep_spin.setEnabled(not is_automated)
        self.maxsteps_spin.setEnabled(is_automated)

        is_live = not is_archive
        for widget in (self.exptime_spin, self.speed_combo, self.binning_combo):
            widget.setEnabled(is_live)

    def _on_browse_clicked(self):
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Select archive directory', self.datadir_edit.text())
        if directory:
            self.datadir_edit.setText(directory)

    def _on_move_to_best_focus_clicked(self):
        if self._best_focus is None:
            return
        # The fit can land on a fractional focus value, but a focus
        # target is a whole-unit telescope position -- round once here,
        # so the confirmation dialog shows the observer exactly the
        # value that will actually be commanded.
        target = round(self._best_focus)
        reply = QtWidgets.QMessageBox.question(
            self, 'Move to best focus',
            f'Move the telescope focus to {target}?',
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No,
        )
        if reply == QtWidgets.QMessageBox.StandardButton.Yes:
            self.moveToBestFocusRequested.emit(target)
