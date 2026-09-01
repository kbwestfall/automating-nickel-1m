"""
Focus-sequence configuration and controls.

See GUI_DESIGN.md §5.4 and claude/GUI_IMPLEMENTATION.md's Phase 4 for the
tabbed redesign. This view only exposes configuration and emits signals
for requested actions; it never constructs a `focus.FocusSequence` or
talks to `gui.model.focus_worker.FocusWorker` itself -- that's the
Controller's job.
"""
from pathlib import Path

from nickel_focus.gui.qt import QtCore, QtWidgets


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
        :func:`_group_box`.
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
    # QFormLayout is a layout specialized for exactly this
    # "label: field" arrangement -- addRow() places each label to the
    # left of its widget and keeps every label in the form the same
    # width, so the fields line up in a column even though the labels
    # themselves are different lengths.
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


def _tighten_margins(widget):
    """
    Shrink ``widget``'s layout margins and inter-row spacing from the
    style's defaults (11px margins, 6px spacing -- sized to comfortably
    fit a titled group box's title, and to keep rows from ever looking
    crowded even when a form has just one or two of them) down to a
    tighter 6px/2px -- used throughout this tab's group boxes and tab
    pages to cut down on unused white space, most of which was that
    default padding repeated on every box and between every row rather
    than anything the content itself needed.

    Parameters
    ----------
    widget : :class:`~PySide6.QtWidgets.QWidget`
        The widget whose layout's margins/spacing should be tightened;
        must already have a layout set.
    """
    layout = widget.layout()
    layout.setContentsMargins(6, 6, 6, 6)
    layout.setSpacing(2)


def _group_box(rows):
    """
    A titleless :class:`~PySide6.QtWidgets.QGroupBox` containing ``rows``
    laid out as a :class:`~PySide6.QtWidgets.QFormLayout` -- the border
    alone groups the fields, the same untitled convention used by the
    Slew tab's sub-panels (see `FocusControlPanel._build_slew_tab`'s
    docstring for why no title is needed there).

    Parameters
    ----------
    rows : list
        ``(label, field)`` pairs, in the order they should appear;
        ``field`` may be a widget or a layout (a
        :class:`~PySide6.QtWidgets.QFormLayout` row accepts either).

    Returns
    -------
    :class:`~PySide6.QtWidgets.QGroupBox`
    """
    form = QtWidgets.QFormLayout()
    for label, field in rows:
        form.addRow(label, field)
    group = QtWidgets.QGroupBox()
    group.setLayout(form)
    _tighten_margins(group)
    return group


def _button_row(*buttons):
    """A :class:`~PySide6.QtWidgets.QHBoxLayout` containing each of ``buttons``, in order."""
    row = QtWidgets.QHBoxLayout()
    for button in buttons:
        row.addWidget(button)
    return row


def _add_centered_label(layout, text):
    """
    Add a :class:`~PySide6.QtWidgets.QLabel` to ``layout``, horizontally
    centered over whatever ``layout`` holds below it, rather than
    stretched to its full width or left-aligned.

    Parameters
    ----------
    layout : :class:`~PySide6.QtWidgets.QBoxLayout`
        The layout to add the label to.
    text : :obj:`str`
        The label's text.
    """
    layout.addWidget(QtWidgets.QLabel(text), alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)


def _tab_widget(layout):
    """
    Wrap ``layout`` in a scrollable page, ready to hand to
    :func:`~PySide6.QtWidgets.QTabWidget.addTab` -- every
    ``_build_*_tab`` method below ends by returning this.

    Scrollable (rather than a bare :class:`~PySide6.QtWidgets.QWidget`)
    so a tab whose content doesn't fit the available height scrolls
    internally, keeping the tab bar itself fixed in place instead of
    scrolling out of view along with everything else -- `MainWindow` no
    longer needs to size its own control-panel area around the tallest
    page (Help) for exactly this reason. The frame is turned off since
    the surrounding `QTabWidget` already visually frames each page;
    without this, a page would show two nested borders.
    """
    widget = QtWidgets.QWidget()
    widget.setLayout(layout)
    _tighten_margins(widget)

    scroll = QtWidgets.QScrollArea()
    scroll.setWidget(widget)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
    return scroll


class FocusControlPanel(QtWidgets.QWidget):
    """
    Sequence configuration/action tabs -- plus live status/history on
    the Log tab and short usage reminders on the Help tab.

    One tab per action -- Slew, Single, Grid, Auto, Replay -- each
    showing only the fields relevant to it. Single/Grid/Auto/Replay each
    lay their fields out as two side-by-side, titleless
    :class:`~PySide6.QtWidgets.QGroupBox` sub-panels (see
    :func:`_group_box`) -- exposure parameters (or, for Replay,
    on-disk file input) on the left, the focus value(s) that define the
    sequence on the right -- with its own acquisition button(s) below
    both, spanning the full tab width rather than belonging to either
    sub-panel: Single has "Acquire"; Grid and Auto have "Acquire" and
    "Interrupt"; Replay has "Load". Log and Help have no buttons at all.
    Because a `QTabWidget` only lets the user interact with the
    currently visible page, each button unambiguously means "do this
    tab's action" -- there's no need to track which tab is active to
    decide what a click means, only to decide what
    :func:`get_focus_sequence_type`/:func:`get_focus_sequence_config`
    should return once a click has already happened.

    The Slew tab (first in the list) is a small control panel around
    `slew.NickelTelescopePointing`/`slew.find_nearest_target`, mirroring
    `scripts/slew_to_nearest.py`: a live "current position" display and
    a manually-editable RA/Dec target with its own "Move to Target"
    button on the left; a starlist file/browse field, an object-name
    search string, and the three "Find nearest..." buttons (which
    populate the target fields rather than moving the telescope
    themselves) on the right. Unlike the other tabs, it has no
    `get_*_config` counterpart -- each button emits its own request
    signal with exactly the text fields it needs, since there is no
    single "the currently active tab's config" to read (this tab is
    never in competition with Single/Grid/Auto/Replay for what
    `startRequested` means).

    The Single tab doubles as "move to best focus": its focus field
    defaults to the most recent fitted best focus (via
    :func:`show_best_focus`), so acquiring a single exposure there with
    the default value *is* moving to best focus; changing the value
    first tests any other focus instead. There is no separate
    confirmation step -- reviewing/editing the value before pressing
    Acquire is the confirmation.

    There is no photometry-method selector: every measurement uses the
    brightest detected source by default, and the coordinates actually
    used (whether the automatic brightest source or one picked by
    clicking the image and pressing 'm' -- see `ImagePanel`) are
    reported on the Log tab via :func:`update_step`/
    :func:`show_single_exposure_result` rather than shown as a persistent
    "current method" state.

    The Options tab holds app-wide preferences -- currently just
    "Remember settings between sessions," which opts into (and, when
    unchecked, immediately erases) the settings persistence handled by
    `~gui.views.main_window.MainWindow`. It's meant to grow as more
    such preferences become useful, rather than being specific to any
    one sequence action the way the other tabs are.

    Attributes
    ----------
    startRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when Grid/Auto's "Acquire" or Replay's "Load" is
        clicked; the Controller should read
        :func:`get_focus_sequence_type` and
        :func:`get_focus_sequence_config` (and :func:`get_exposure_config`
        for Grid/Auto) to build and run the sequence.
    stopRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when Grid/Auto's "Interrupt" is clicked.
    takeSingleExposureRequested : :class:`~PySide6.QtCore.Signal`
        Emitted with a focus value when Single's "Acquire" is clicked.
    moveToTargetRequested : :class:`~PySide6.QtCore.Signal`
        Emitted with ``(ra_text, dec_text)`` when the Slew tab's "Move to
        Target" is clicked; the Controller should parse these and call
        `slew.NickelTelescopePointing.slew_to`.
    findNearestObjectRequested : :class:`~PySide6.QtCore.Signal`
        Emitted with ``(file_text, search_text)`` when the Slew tab's
        "Find nearest object" is clicked -- either may be an empty
        string, meaning "use the default" (see
        `slew.find_nearest_target`).
    findNearestPointingRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when the Slew tab's "Find nearest pointing star" is
        clicked; the Controller should search the packaged default
        catalog with ``obj_search_str='Pointing'``, ignoring whatever is
        currently in the file/search fields.
    findNearestFocusRequested : :class:`~PySide6.QtCore.Signal`
        Emitted when the Slew tab's "Find nearest focus star" is
        clicked; like ``findNearestPointingRequested``, but with
        ``obj_search_str='Focusing'``.
    """
    startRequested = QtCore.Signal()
    stopRequested = QtCore.Signal()
    takeSingleExposureRequested = QtCore.Signal(float)
    moveToTargetRequested = QtCore.Signal(str, str)
    findNearestObjectRequested = QtCore.Signal(str, str)
    findNearestPointingRequested = QtCore.Signal()
    findNearestFocusRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        self.slew_tab = self._build_slew_tab()
        self.single_tab = self._build_single_tab()
        self.grid_tab = self._build_grid_tab()
        self.auto_tab = self._build_auto_tab()
        self.replay_tab = self._build_replay_tab()
        self.log_tab = self._build_log_tab()
        self.options_tab = self._build_options_tab()
        self.help_tab = self._build_help_tab()

        # QTabWidget shows exactly one of its pages (each an ordinary
        # widget, built above) at a time, switched via the row of tabs
        # it draws automatically from the label given to addTab().
        # Pages that aren't currently showing still exist as real
        # widgets in memory -- their state (a spin box's value, etc.)
        # isn't lost when you switch away -- they're just not painted or
        # interactive until their tab is selected again.
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.slew_tab, 'Slew')
        self.tabs.addTab(self.single_tab, 'Single')
        self.tabs.addTab(self.grid_tab, 'Grid')
        self.tabs.addTab(self.auto_tab, 'Auto')
        self.tabs.addTab(self.replay_tab, 'Replay')
        self.tabs.addTab(self.log_tab, 'Log')
        self.tabs.addTab(self.options_tab, 'Options')
        self.tabs.addTab(self.help_tab, 'Help')

        # Every input widget across the four actionable tabs, locked
        # uniformly while anything is running (hardware exclusivity:
        # only one operation at a time, so there's no reason to allow
        # editing a different tab's configuration meanwhile). The tab
        # bar itself is deliberately never disabled, so the Log tab
        # stays reachable while a sequence runs.
        self._config_widgets = [
            self.slew_target_ra_edit, self.slew_target_dec_edit, self.slew_file_edit,
            self.slew_browse_button, self.slew_search_edit,
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
        # The "go" buttons on the five actionable tabs, disabled while
        # anything is running.
        self._acquire_buttons = [
            self.slew_move_button, self.slew_find_object_button,
            self.slew_find_pointing_button, self.slew_find_focus_button,
            self.single_acquire_button, self.grid_acquire_button,
            self.auto_acquire_button, self.replay_load_button,
        ]
        # The "stop" buttons on Grid/Auto, enabled only while something
        # is running -- either one works, since `stopRequested` always
        # targets whatever's actually running, regardless of which tab
        # is currently visible.
        self._interrupt_buttons = [self.grid_interrupt_button, self.auto_interrupt_button]

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.tabs)
        layout.setContentsMargins(0, 0, 0, 0)

    # -- widget construction ----------------------------------------------

    def _build_acquire_interrupt_row(self):
        """
        An "Acquire"/"Interrupt" button pair wired to
        :attr:`startRequested`/:attr:`stopRequested`, used identically
        by the Grid and Auto tabs (the only two tabs with a running,
        interruptible sequence). Interrupt starts disabled; see
        :func:`set_running`/:func:`set_stopping` for how the two
        buttons' enabled states are actually managed once the panel is
        in use.

        Returns
        -------
        tuple
            ``(button_row, acquire_button, interrupt_button)``.
        """
        acquire_button = QtWidgets.QPushButton('Acquire')
        # Connecting a button's `clicked` signal straight to another
        # signal's `.emit` -- rather than to a Python method that then
        # calls `.emit()` itself -- is a common Qt shorthand: `.emit` is
        # itself just a callable, so anything that can be a slot
        # (including it) can be connected directly.
        acquire_button.clicked.connect(self.startRequested.emit)
        interrupt_button = QtWidgets.QPushButton('Interrupt')
        interrupt_button.setEnabled(False)
        interrupt_button.clicked.connect(self.stopRequested.emit)
        return _button_row(acquire_button, interrupt_button), acquire_button, interrupt_button

    def _build_slew_tab(self):
        """
        Build the Slew tab as two untitled `~PySide6.QtWidgets.QGroupBox`
        sub-panels: current/target RA-Dec and "Move to Target" on the
        left, starlist file/browse, a target-name search pattern, and a
        single row of three "Find nearest ..." buttons ("Object",
        "Pointing *", "Focusing *" -- ``*`` standing in for "star") on
        the right. Unlike every other tab here (deliberately flat since
        the tab title alone says what its fields are for -- see
        claude/GUI_IMPLEMENTATION.md's Phase 4), Slew covers two
        distinct sub-workflows that its one tab title can't itself
        disambiguate, which is what makes the extra grouping worth its
        cost here specifically; each panel's own content (a coordinate
        pair vs. a search field and buttons) already identifies its
        purpose, so neither needs its own title on top of that border.
        """
        # Short "RA:"/"Dec:" field labels, not "Current RA:"/"Target RA:"
        # etc. -- those are wide enough (with two side-by-side columns
        # per row) to noticeably widen this tab beyond every other one;
        # a short section header above each row disambiguates instead,
        # at a fraction of the width cost (GUI_DESIGN.md's §9-phase-5
        # scroll-area floor otherwise has to grow to match).
        self.slew_current_ra_label = QtWidgets.QLabel('—')
        self.slew_current_dec_label = QtWidgets.QLabel('—')
        self.slew_target_ra_edit = QtWidgets.QLineEdit()
        self.slew_target_dec_edit = QtWidgets.QLineEdit()
        # Wide enough for a full sexagesimal value (e.g. '+00:00:00.0')
        # to show without scrolling -- a bare default QLineEdit is only
        # wide enough for a few characters. `sizeHint()` on a QLineEdit
        # already containing that text turned out to reserve much more
        # than the text itself needs (a fixed, style-dependent margin,
        # not something proportional to the text); the font metrics'
        # bare advance for the text plus a small fixed pad is tighter
        # and still leaves room to spare (confirmed against
        # `cursorRect()` with the text entered).
        coord_width = (
            self.slew_target_ra_edit.fontMetrics().horizontalAdvance('+00:00:00.0') + 10
        )
        self.slew_target_ra_edit.setMinimumWidth(coord_width)
        self.slew_target_dec_edit.setMinimumWidth(coord_width)
        current_row = _two_column_row(
            [('RA:', self.slew_current_ra_label)],
            [('Dec:', self.slew_current_dec_label)],
        )
        target_row = _two_column_row(
            [('RA:', self.slew_target_ra_edit)],
            [('Dec:', self.slew_target_dec_edit)],
        )
        self.slew_move_button = QtWidgets.QPushButton('Move to Target')
        self.slew_move_button.clicked.connect(
            lambda: self.moveToTargetRequested.emit(
                self.slew_target_ra_edit.text(), self.slew_target_dec_edit.text()))

        position_layout = QtWidgets.QVBoxLayout()
        _add_centered_label(position_layout, 'Current coordinates')
        position_layout.addLayout(current_row)
        _add_centered_label(position_layout, 'Target coordinates')
        position_layout.addLayout(target_row)
        position_layout.addLayout(_button_row(self.slew_move_button))
        position_layout.addStretch(1)
        position_group = QtWidgets.QGroupBox()
        position_group.setLayout(position_layout)
        _tighten_margins(position_group)

        self.slew_file_edit = QtWidgets.QLineEdit('')
        self.slew_browse_button = QtWidgets.QPushButton('Browse…')
        self.slew_browse_button.clicked.connect(self._on_slew_browse_clicked)
        file_row = QtWidgets.QHBoxLayout()
        file_row.addWidget(self.slew_file_edit, 1)
        file_row.addWidget(self.slew_browse_button)

        self.slew_search_edit = QtWidgets.QLineEdit('')

        # Short button labels ("Object"/"Pointing */"Focusing *", the
        # asterisk standing in for "star") so the three fit in one row
        # instead of stacking -- the "Find nearest ..." label above
        # supplies the words those labels leave out.
        self.slew_find_object_button = QtWidgets.QPushButton('Object')
        self.slew_find_object_button.clicked.connect(
            lambda: self.findNearestObjectRequested.emit(
                self.slew_file_edit.text(), self.slew_search_edit.text()))
        self.slew_find_pointing_button = QtWidgets.QPushButton('Pointing *')
        self.slew_find_pointing_button.clicked.connect(self.findNearestPointingRequested.emit)
        self.slew_find_focus_button = QtWidgets.QPushButton('Focusing *')
        self.slew_find_focus_button.clicked.connect(self.findNearestFocusRequested.emit)

        find_row = _button_row(
            self.slew_find_object_button, self.slew_find_pointing_button,
            self.slew_find_focus_button)
        # Tighter than the default gap between widgets in a layout --
        # frees up width for the Position panel next to it without
        # actually shrinking any button.
        find_row.setSpacing(2)

        # Narrow enough that neither dictates the Search panel's width
        # on its own -- that's the button row above, now the widest
        # thing here. Still wide enough to show a short search pattern
        # or filename without immediately scrolling.
        text_width = QtWidgets.QLineEdit('M' * 8).sizeHint().width()
        self.slew_file_edit.setMinimumWidth(text_width)
        self.slew_search_edit.setMinimumWidth(text_width)

        search_layout = QtWidgets.QVBoxLayout()
        _add_centered_label(search_layout, 'Starlist')
        search_layout.addLayout(file_row)
        _add_centered_label(search_layout, 'Target pattern')
        search_layout.addWidget(self.slew_search_edit)
        _add_centered_label(search_layout, 'Find nearest ...')
        search_layout.addLayout(find_row)
        search_layout.addStretch(1)
        search_group = QtWidgets.QGroupBox()
        search_group.setLayout(search_layout)
        _tighten_margins(search_group)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(position_group, 1)
        layout.addWidget(search_group, 1)
        return _tab_widget(layout)

    def _build_single_tab(self):
        """
        Build the Single tab as two side-by-side sub-panels -- exposure
        settings, focus value -- with "Acquire" below both, spanning the
        full tab width rather than belonging to either one.
        """
        self.single_focus_spin = _int_spin(165, 500, 340)
        exposure_rows, self.single_exptime_spin, self.single_speed_combo, \
            self.single_binning_combo = _exposure_field_rows()
        focus_rows = [('Focus value:', self.single_focus_spin)]

        self.single_acquire_button = QtWidgets.QPushButton('Acquire')
        self.single_acquire_button.clicked.connect(
            lambda: self.takeSingleExposureRequested.emit(self.single_focus_spin.value()))

        columns = QtWidgets.QHBoxLayout()
        columns.addWidget(_group_box(exposure_rows), 1)
        columns.addWidget(_group_box(focus_rows), 1)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(columns)
        layout.addLayout(_button_row(self.single_acquire_button))
        layout.addStretch(1)
        return _tab_widget(layout)

    def _build_grid_tab(self):
        """
        Build the Grid tab as two side-by-side sub-panels -- exposure
        settings, start focus/step size/number of steps -- with
        Acquire/Interrupt below both, spanning the full tab width.
        """
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
        button_row, self.grid_acquire_button, self.grid_interrupt_button = \
            self._build_acquire_interrupt_row()

        columns = QtWidgets.QHBoxLayout()
        columns.addWidget(_group_box(exposure_rows), 1)
        columns.addWidget(_group_box(focus_rows), 1)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(columns)
        layout.addLayout(button_row)
        layout.addStretch(1)
        return _tab_widget(layout)

    def _build_auto_tab(self):
        """
        Build the Auto tab as two side-by-side sub-panels -- exposure
        settings, start focus/step size/max steps -- with
        Acquire/Interrupt below both, spanning the full tab width.
        """
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
        button_row, self.auto_acquire_button, self.auto_interrupt_button = \
            self._build_acquire_interrupt_row()

        columns = QtWidgets.QHBoxLayout()
        columns.addWidget(_group_box(exposure_rows), 1)
        columns.addWidget(_group_box(focus_rows), 1)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(columns)
        layout.addLayout(button_row)
        layout.addStretch(1)
        return _tab_widget(layout)

    def _build_replay_tab(self):
        """
        Build the Replay tab as two side-by-side sub-panels -- on-disk
        file input (data directory, filename fields), the focus grid
        describing the sequence to replay -- with "Load" below both,
        spanning the full tab width.
        """
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

        # "Data directory" centered above the text box/browse button,
        # rather than a QFormLayout label beside them -- matching the
        # Slew tab's "Starlist" field (see `_build_slew_tab`).
        file_layout = QtWidgets.QVBoxLayout()
        _add_centered_label(file_layout, 'Data directory')
        file_layout.addLayout(datadir_row)
        filename_form = QtWidgets.QFormLayout()
        filename_form.addRow('Prefix:', self.replay_prefix_edit)
        filename_form.addRow('Suffix:', self.replay_suffix_edit)
        filename_form.addRow('Obsnum:', self.replay_obsnum_spin)
        file_layout.addLayout(filename_form)
        file_group = QtWidgets.QGroupBox()
        file_group.setLayout(file_layout)
        _tighten_margins(file_group)

        focus_rows = [
            ('Start focus:', self.replay_start_spin),
            ('Step size:', self.replay_step_spin),
            ('Number of steps:', self.replay_nstep_spin),
        ]

        self.replay_load_button = QtWidgets.QPushButton('Load')
        self.replay_load_button.clicked.connect(self.startRequested.emit)

        columns = QtWidgets.QHBoxLayout()
        columns.addWidget(file_group, 1)
        columns.addWidget(_group_box(focus_rows), 1)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(columns)
        layout.addLayout(_button_row(self.replay_load_button))
        layout.addStretch(1)
        return _tab_widget(layout)

    def _build_log_tab(self):
        """Build the Log tab: the status line, step line, and scrolling history."""
        # Word wrap is essential here, not cosmetic: a long message (e.g.
        # a failure listing several missing file paths) in an unwrapped
        # QLabel reports its *entire single-line* width as the label's
        # minimum size, which balloons this whole panel's width to match.
        # Wrapping lets it grow in height instead, which the surrounding
        # QScrollArea (see MainWindow) handles gracefully.
        self.status_label = QtWidgets.QLabel('')
        self.status_label.setWordWrap(True)
        self.step_label = QtWidgets.QLabel('Step: —')
        self.step_label.setWordWrap(True)
        self.log_widget = QtWidgets.QPlainTextEdit()
        self.log_widget.setReadOnly(True)
        # A QPlainTextEdit grows without bound by default, appending
        # one line at a time all session -- capping it here is what
        # keeps a very long run from slowly consuming ever more memory.
        self.log_widget.setMaximumBlockCount(500)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.status_label)
        layout.addWidget(self.step_label)
        layout.addWidget(self.log_widget)
        return _tab_widget(layout)

    def _build_options_tab(self):
        """Build the Options tab: currently just the "remember settings" checkbox."""
        # Meant to grow as more app-wide (as opposed to per-sequence)
        # preferences become useful -- currently just the one checkbox.
        self.remember_settings_checkbox = QtWidgets.QCheckBox(
            'Remember settings between sessions')

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.remember_settings_checkbox)
        layout.addStretch(1)
        return _tab_widget(layout)

    def _build_help_tab(self):
        """Build the Help tab: a single block of short, rich-text usage reminders."""
        text = (
            '<b>Slew</b> — move the telescope to a target. Enter RA/Dec '
            'directly, or use "Find nearest..." to search a starlist '
            '(the packaged pointing/focus catalog by default) and fill '
            'them in, then press "Move to Target."<br><br>'
            '<b>Single</b> — take one exposure at the given focus. '
            'Defaults to the most recent fitted best focus, so acquiring '
            'it as-is moves to best focus; change the value first to '
            'test any other focus.<br><br>'
            '<b>Grid</b> — step through an evenly spaced focus grid.<br>'
            '<b>Auto</b> — adaptively search for the best focus.<br>'
            '<b>Replay</b> — reprocess an existing set of exposures from disk.<br>'
            '<b>Log</b> — live status and history.<br>'
            '<b>Options</b> — app-wide preferences.<br><br>'
            '<b>Source selection</b> — the brightest star is measured '
            "automatically by default (see the Log for its coordinates); "
            "hover over a different star in the image and press 'm' to "
            'measure that one instead.<br><br>'
            '<b>Acquire / Load</b> — run the action configured in that '
            'tab.<br>'
            '<b>Interrupt</b> — stop a running Grid/Auto sequence between '
            'steps (does not abort the current exposure).'
        )
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)
        # A QLabel treats its text as plain text by default -- this is
        # what makes it interpret the <b>/<br> tags above as HTML markup
        # instead of displaying them as literal characters.
        label.setTextFormat(QtCore.Qt.TextFormat.RichText)
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(label)
        layout.addStretch(1)
        return _tab_widget(layout)

    # -- public API used by the Controller/MainWindow -----------------------

    def preferred_height_excluding_help(self):
        """
        The tallest *minimum* height among every tab page except Help,
        plus the tab bar -- i.e. how tall this panel should start out so
        every *interactive* tab displays without scrolling by default.

        Not a hard floor: every tab page is individually scrollable (see
        `_tab_widget`), so this panel can be shrunk well below this value
        -- any tab that no longer fits just scrolls internally, and the
        tab bar itself (outside all of that scrolling) stays put either
        way. Help is excluded because it's a large block of static
        reference text that's taller than every other tab combined with
        the tab bar; sizing the *initial* height around it too would
        waste space under the six shorter, interactive tabs every time
        the window opens, for the sake of a tab that's fine to need a
        scrollbar immediately.

        Each candidate tab's *content widget* `minimumSizeHint` is used
        here, not the tab page's own (a `~PySide6.QtWidgets.QScrollArea`,
        whose `sizeHint` reflects a generically "comfortable" size that
        can run well past what its content actually needs -- e.g. the
        Log tab's `QPlainTextEdit`, which is happy to be much shorter
        than its `sizeHint` and was otherwise the tallest tab here,
        ahead of Replay).

        Returns
        -------
        :obj:`int`
            Preferred height, in pixels, that fits every tab page except
            Help without scrolling.
        """
        other_tabs = [self.slew_tab, self.single_tab, self.grid_tab, self.auto_tab,
                      self.replay_tab, self.log_tab, self.options_tab]
        tallest_other = max(tab.widget().minimumSizeHint().height() for tab in other_tabs)
        tab_bar_height = self.tabs.tabBar().sizeHint().height()
        # `QTabWidget` reserves a couple more pixels around the page
        # itself (its pane frame) beyond just the tab bar -- without
        # this, the tallest tab ends up a hair short of fitting and gets
        # an unnecessary, barely-there scrollbar (confirmed against the
        # actual rendered tab: a couple of pixels was all it took).
        return tallest_other + tab_bar_height + 4

    def get_focus_sequence_type(self):
        """
        The sequence type for the currently active tab, as ``'grid'``,
        ``'automated'``, or ``'archive'`` -- or ``None`` if Single, Log,
        or Help is active (in which case no button can have triggered
        `startRequested` in the first place).
        """
        tab = self.tabs.currentWidget()
        if tab is self.grid_tab:
            return 'grid'
        if tab is self.auto_tab:
            return 'automated'
        if tab is self.replay_tab:
            return 'archive'
        return None

    def get_focus_sequence_config(self):
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

    def get_settings_state(self):
        """
        Return every persistable field's current value as a flat
        :obj:`dict` keyed by ``'<tab>/<field>'`` (e.g. ``'grid/start'``),
        suitable for a settings store like
        :class:`~PySide6.QtCore.QSettings`. The Single tab's focus value
        is deliberately excluded -- it tracks the most recent fit (see
        :func:`show_best_focus`), not a persisted default.
        """
        state = {}
        for key, widget in self._settings_fields().items():
            if isinstance(widget, QtWidgets.QComboBox):
                state[key] = widget.currentText()
            elif isinstance(widget, QtWidgets.QLineEdit):
                state[key] = widget.text()
            else:
                state[key] = widget.value()
        return state

    def set_settings_state(self, state):
        """
        Apply a settings :obj:`dict` from :func:`get_settings_state` (or
        any subset of it -- e.g. an older saved state missing a field
        added later). Keys not recognized, or missing entirely, are
        left at their current value.
        """
        for key, widget in self._settings_fields().items():
            if key not in state:
                continue
            value = state[key]
            if isinstance(widget, QtWidgets.QComboBox):
                widget.setCurrentText(value)
            elif isinstance(widget, QtWidgets.QLineEdit):
                widget.setText(value)
            else:
                widget.setValue(value)

    def _settings_fields(self):
        """
        Map of persistable setting key to the widget holding it. The
        single source of truth for both :func:`get_settings_state` and
        :func:`set_settings_state`, so the two can't drift apart into
        covering different fields.
        """
        return {
            'slew/file': self.slew_file_edit,
            'slew/search': self.slew_search_edit,
            'single/exptime': self.single_exptime_spin,
            'single/speed': self.single_speed_combo,
            'single/binning': self.single_binning_combo,
            'grid/start': self.grid_start_spin,
            'grid/step': self.grid_step_spin,
            'grid/nstep': self.grid_nstep_spin,
            'grid/exptime': self.grid_exptime_spin,
            'grid/speed': self.grid_speed_combo,
            'grid/binning': self.grid_binning_combo,
            'auto/start': self.auto_start_spin,
            'auto/step': self.auto_step_spin,
            'auto/maxsteps': self.auto_maxsteps_spin,
            'auto/exptime': self.auto_exptime_spin,
            'auto/speed': self.auto_speed_combo,
            'auto/binning': self.auto_binning_combo,
            'replay/datadir': self.replay_datadir_edit,
            'replay/prefix': self.replay_prefix_edit,
            'replay/suffix': self.replay_suffix_edit,
            'replay/obsnum': self.replay_obsnum_spin,
            'replay/start': self.replay_start_spin,
            'replay/step': self.replay_step_spin,
            'replay/nstep': self.replay_nstep_spin,
        }

    def set_hardware_tabs_enabled(self, enabled):
        """
        Enable or disable the four tabs that need a live ``ktl``
        connection to do anything (Slew, Single, Grid, Auto) -- called
        once from `Controller.__init__` (see there for how ``enabled``
        is decided). This is a one-time, ``ktl``-availability toggle,
        independent of :func:`set_running`'s per-action running state.

        Disabling a tab via
        `~PySide6.QtWidgets.QTabWidget.setTabEnabled` grays it out and
        disables every widget on it, so it stays visible but can't be
        clicked into or interacted with; the Replay tab -- the only one
        that works with no ``ktl`` connection at all -- is left alone,
        and becomes the tab shown when the window first opens.

        Parameters
        ----------
        enabled : :obj:`bool`
            Whether the four ``ktl``-driven tabs should be enabled.
        """
        for tab in (self.slew_tab, self.single_tab, self.grid_tab, self.auto_tab):
            self.tabs.setTabEnabled(self.tabs.indexOf(tab), enabled)
        if not enabled:
            self.tabs.setCurrentWidget(self.replay_tab)

    def set_running(self, running):
        """Toggle widget states for whether an action is currently running."""
        for widget in self._config_widgets:
            widget.setEnabled(not running)
        for button in self._acquire_buttons:
            button.setEnabled(not running)
        for button in self._interrupt_buttons:
            button.setEnabled(running)

        if not running and self.status_label.text().startswith('Stopping'):
            # Only clear a stale "Stopping..." message -- not a
            # confirmation/failure the worker's finishing signal may have
            # just set moments before this.
            self.status_label.setText('')

    def set_stopping(self):
        """Reflect an Interrupt request that hasn't taken effect yet (§4.3)."""
        for button in self._interrupt_buttons:
            button.setEnabled(False)
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

    def set_current_position(self, ra_text, dec_text):
        """
        Update the Slew tab's live "current position" display.

        Parameters
        ----------
        ra_text : :obj:`str`
            The telescope's current right ascension, already formatted
            for display (e.g. ``'05:30:00'``) -- formatting an `Angle`
            is the Controller's job, not this passive view's.
        dec_text : :obj:`str`
            The telescope's current declination, formatted the same way.
        """
        self.slew_current_ra_label.setText(ra_text)
        self.slew_current_dec_label.setText(dec_text)

    def show_nearest_target(self, name, ra_text, dec_text):
        """
        Populate the Slew tab's target RA/Dec fields with a found
        nearest target, and report it (mirrors
        `scripts/slew_to_nearest.py`'s printed message).

        Parameters
        ----------
        name : :obj:`str`
            The nearest target's name.
        ra_text : :obj:`str`
            The nearest target's right ascension, already formatted for
            display.
        dec_text : :obj:`str`
            The nearest target's declination, already formatted for
            display.
        """
        self.slew_target_ra_edit.setText(ra_text)
        self.slew_target_dec_edit.setText(dec_text)
        self.status_label.setText(f'Nearest target: {name} (RA={ra_text}, Dec={dec_text})')

    def show_slew_result(self, message):
        """Report a completed "Move to Target" slew (e.g., from `SlewWorker.slewFinished`)."""
        self.status_label.setText(message)

    def show_best_focus(self, best_focus, best_fwhm):
        """
        Report a completed sequence's fitted result, and set it as the
        Single tab's default focus value -- rounded, since a focus value
        is a whole-unit telescope position even though the fit itself can
        land on a fractional one.
        """
        self.status_label.setText(
            f'Sequence finished: best focus {best_focus:.1f}, expected FWHM {best_fwhm:.2f}')
        self.single_focus_spin.setValue(round(best_focus))

    def show_single_exposure_result(self, result):
        """Report one exposure taken from the Single tab."""
        self.status_label.setText(
            f'Took exposure at focus {result.focus_value:.0f}: measured FWHM '
            f'{result.fwhm:.2f}, source ({result.centroid[0]:.1f}, {result.centroid[1]:.1f})')

    def show_failure(self, message):
        """Display a sequence failure (e.g., from `FocusWorker.focusSequenceFailed`)."""
        self.status_label.setText(message)

    def append_log_line(self, text):
        """Append one line of text to the Log tab's scrolling history."""
        self.log_widget.appendPlainText(text)

    def reset(self):
        """Clear live status and log for a new run."""
        self.step_label.setText('Step: —')
        self.status_label.setText('')
        self.log_widget.clear()

    # -- internals ----------------------------------------------------------

    def _on_browse_clicked(self):
        """
        Slot for `replay_browse_button`'s ``clicked``: open the native
        OS directory picker and, unless the user cancels (signaled by
        an empty string, not `None`), fill it into `replay_datadir_edit`.
        """
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Select archive directory', self.replay_datadir_edit.text())
        if directory != '':
            self.replay_datadir_edit.setText(directory)

    def _on_slew_browse_clicked(self):
        """
        Slot for `slew_browse_button`'s ``clicked``: open the native OS
        file picker and, unless the user cancels (signaled by an empty
        string, not `None`), fill the selected path into
        `slew_file_edit`.
        """
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select starlist file', self.slew_file_edit.text())
        if path != '':
            self.slew_file_edit.setText(path)
