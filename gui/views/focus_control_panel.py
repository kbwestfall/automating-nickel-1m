"""
Focus-sequence configuration and controls.

See GUI_DESIGN.md §5.4. This view only exposes configuration and emits
signals for requested actions; it never constructs a `focus.FocusSequence`
or talks to `gui.model.sequence_worker.SequenceWorker` itself -- that's
the Controller's job (sub-phase 7).
"""
from pathlib import Path

from gui.qt import QtCore, QtWidgets

_PHASE_3_TOOLTIP = 'Requires a live ktl connection (Phase 3)'
_NOT_APPLICABLE_TOOLTIP = 'Not applicable in archive/replay mode'


class FocusControlPanel(QtWidgets.QWidget):
    """
    Sequence configuration, Start/Stop, live status, and results.

    Only "Archive / Replay" is enabled in this phase (`GUI_DESIGN.md` §9
    phased plan, sub-phase 6): "Grid" and "Automated" are present but
    disabled, so this widget's shape doesn't need to change once Phase 3
    wires up live sequences. Likewise, exposure settings (exposure time,
    speed, binning) are shown but disabled, since archive/replay mode
    never takes new exposures.

    Attributes
    ----------
    startRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when "Start" is clicked; the Controller should read
        :func:`get_archive_config` to build the sequence.
    stopRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when "Stop" is clicked.
    methodChanged : :class:`~PySide6.QtCore.Signal`
        Emitted with ``'brightest'`` or ``'weighted'`` when the user picks
        one from the method combo box.
    moveToBestFocusRequested : :class:`~PySide6.QtCore.Signal`
        Emitted with the target focus value once the user confirms the
        "Move to best focus" dialog. The button that triggers this is
        disabled unconditionally in this phase -- there's no hardware to
        move yet -- but the confirmation logic is already wired up so
        Phase 3 only needs to enable it.
    """
    startRequested = QtCore.Signal()
    stopRequested = QtCore.Signal()
    methodChanged = QtCore.Signal(str)
    moveToBestFocusRequested = QtCore.Signal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._best_focus = None

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self._build_sequence_type_group())
        layout.addWidget(self._build_archive_config_group())
        layout.addWidget(self._build_exposure_settings_group())
        layout.addWidget(self._build_method_group())
        layout.addLayout(self._build_start_stop_row())
        layout.addWidget(self._build_status_group())
        layout.addWidget(self._build_result_group())
        layout.addStretch(1)

    # -- widget construction ----------------------------------------------

    def _build_sequence_type_group(self):
        group = QtWidgets.QGroupBox('Sequence Type')
        self.archive_radio = QtWidgets.QRadioButton('Archive / Replay')
        self.grid_radio = QtWidgets.QRadioButton('Grid')
        self.automated_radio = QtWidgets.QRadioButton('Automated')
        self.archive_radio.setChecked(True)
        for radio in (self.grid_radio, self.automated_radio):
            radio.setEnabled(False)
            radio.setToolTip(_PHASE_3_TOOLTIP)

        self.sequence_type_buttons = QtWidgets.QButtonGroup(self)
        for radio in (self.archive_radio, self.grid_radio, self.automated_radio):
            self.sequence_type_buttons.addButton(radio)

        row = QtWidgets.QHBoxLayout(group)
        row.addWidget(self.archive_radio)
        row.addWidget(self.grid_radio)
        row.addWidget(self.automated_radio)
        return group

    def _build_archive_config_group(self):
        group = QtWidgets.QGroupBox('Archive Configuration')

        self.datadir_edit = QtWidgets.QLineEdit('.')
        self.browse_button = QtWidgets.QPushButton('Browse…')
        self.browse_button.clicked.connect(self._on_browse_clicked)
        self.prefix_edit = QtWidgets.QLineEdit('n')
        self.suffix_edit = QtWidgets.QLineEdit('.fits')
        self.obsnum_spin = QtWidgets.QSpinBox()
        self.obsnum_spin.setRange(0, 999999)
        self.start_spin = QtWidgets.QDoubleSpinBox()
        self.start_spin.setRange(165., 500.)
        self.start_spin.setValue(340.)
        self.step_spin = QtWidgets.QDoubleSpinBox()
        self.step_spin.setRange(0.1, 100.)
        self.step_spin.setValue(5.)
        self.nstep_spin = QtWidgets.QSpinBox()
        self.nstep_spin.setRange(3, 100)
        self.nstep_spin.setValue(5)

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
        form.addRow('Number of steps:', self.nstep_spin)
        return group

    def _build_exposure_settings_group(self):
        group = QtWidgets.QGroupBox('Exposure Settings')
        self.exptime_spin = QtWidgets.QDoubleSpinBox()
        self.exptime_spin.setRange(0.1, 600.)
        self.exptime_spin.setValue(5.)
        self.speed_combo = QtWidgets.QComboBox()
        self.speed_combo.addItems(['Slow', 'Fast'])
        self.binning_combo = QtWidgets.QComboBox()
        self.binning_combo.addItems(['1,1', '2,2', '4,4'])
        for widget in (self.exptime_spin, self.speed_combo, self.binning_combo):
            widget.setEnabled(False)
            widget.setToolTip(_NOT_APPLICABLE_TOOLTIP)

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
        self.move_to_best_focus_button.setToolTip(_PHASE_3_TOOLTIP)
        self.move_to_best_focus_button.clicked.connect(self._on_move_to_best_focus_clicked)

        col = QtWidgets.QVBoxLayout(group)
        col.addWidget(self.result_label)
        col.addWidget(self.move_to_best_focus_button)
        return group

    # -- public API used by the Controller ---------------------------------

    def get_archive_config(self):
        """
        Return the current archive-mode configuration as a :obj:`dict`
        with keys ``datadir`` (:class:`pathlib.Path`), ``prefix``,
        ``suffix`` (:obj:`str`), ``obsnum`` (:obj:`int`), and ``start``,
        ``step`` (:obj:`float`), ``nstep`` (:obj:`int`) -- matching
        `focus.py`'s own `--datadir`/`--prefix`/`--suffix`/`--obsnum`
        and positional ``focus``/`--nstep` CLI arguments.
        """
        return {
            'datadir': Path(self.datadir_edit.text()),
            'prefix': self.prefix_edit.text(),
            'suffix': self.suffix_edit.text(),
            'obsnum': self.obsnum_spin.value(),
            'start': self.start_spin.value(),
            'step': self.step_spin.value(),
            'nstep': self.nstep_spin.value(),
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
        self.archive_radio.setEnabled(not running)
        for widget in (self.datadir_edit, self.browse_button, self.prefix_edit,
                       self.suffix_edit, self.obsnum_spin, self.start_spin,
                       self.step_spin, self.nstep_spin, self.method_combo):
            widget.setEnabled(not running)
        if not running:
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
        text += f' — Focus {result.focus_value:.1f}, FWHM {result.fwhm:.2f}'
        if result.is_outlier:
            text += '  [outlier]'
        self.step_label.setText(text)
        self.log_widget.appendPlainText(text)

    def show_best_focus(self, best_focus, best_fwhm):
        """Display a completed sequence's fitted result."""
        self._best_focus = best_focus
        self.result_label.setText(f'Best focus: {best_focus:.1f}   Expected FWHM: {best_fwhm:.2f}')
        self.log_widget.appendPlainText(
            f'Sequence finished: best focus {best_focus:.1f}, expected FWHM {best_fwhm:.2f}')

    def show_failure(self, message):
        """Display a sequence failure (e.g., from `SequenceWorker.sequenceFailed`)."""
        self.status_label.setText(message)
        self.log_widget.appendPlainText(f'ERROR: {message}')

    def reset(self):
        """Clear live status, log, and result display for a new run."""
        self._best_focus = None
        self.step_label.setText('Step: —')
        self.result_label.setText('Best focus: —')
        self.status_label.setText('')
        self.log_widget.clear()

    # -- internals ----------------------------------------------------------

    def _on_browse_clicked(self):
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Select archive directory', self.datadir_edit.text())
        if directory:
            self.datadir_edit.setText(directory)

    def _on_move_to_best_focus_clicked(self):
        if self._best_focus is None:
            return
        reply = QtWidgets.QMessageBox.question(
            self, 'Move to best focus',
            f'Move the telescope focus to {self._best_focus:.1f}?',
            QtWidgets.QMessageBox.StandardButton.Yes | QtWidgets.QMessageBox.StandardButton.No,
            QtWidgets.QMessageBox.StandardButton.No,
        )
        if reply == QtWidgets.QMessageBox.StandardButton.Yes:
            self.moveToBestFocusRequested.emit(self._best_focus)
