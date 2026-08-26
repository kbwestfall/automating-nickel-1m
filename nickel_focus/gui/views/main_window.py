"""
Top-level window laying out the three panels.

See GUI_DESIGN.md §5.1.
"""
from nickel_focus.gui.qt import QtCore, QtWidgets
from nickel_focus.gui.views.image_panel import ImagePanel
from nickel_focus.gui.views.focus_curve_panel import FocusCurvePanel
from nickel_focus.gui.views.focus_control_panel import FocusControlPanel


class MainWindow(QtWidgets.QMainWindow):
    """
    Lays out `ImagePanel` (left), `FocusCurvePanel` (top right), and
    `FocusControlPanel` (bottom right) per §5.1's sketch. Purely
    structural -- wiring the panels together is `gui.controller.Controller`'s
    job, not this class's.

    Also owns settings persistence (GUI_DESIGN.md §9 phase 5): each
    panel's last-used configuration is restored on construction and
    saved on close, via `~PySide6.QtCore.QSettings`. Deliberately not
    saved on every change -- these are convenience defaults for next
    time, not data that needs crash-safety.

    Persistence is opt-in, via the Options tab's "Remember settings
    between sessions" checkbox
    (`~gui.views.focus_control_panel.FocusControlPanel.remember_settings_checkbox`),
    and that opt-in is never itself written to disk as a separate flag --
    doing so would mean *every* user silently gets a settings file
    written the moment they close the window, whether they ever checked
    the box or not. Instead, whether anything was ever saved *is* the
    opt-in signal: on load, the checkbox is checked if (and only if) a
    previous session actually saved something; on close, nothing is
    written at all unless the checkbox is checked, and any
    previously-saved configuration is actively erased the moment it's
    unchecked (rather than merely left unsaved and lingering unused, or
    the checkbox itself lingering on disk with nothing behind it).
    """

    #: Organization/application name identifying this app's `QSettings`
    #: store (an macOS plist, a Windows registry key, or a Linux ini
    #: file, depending on platform).
    _SETTINGS_ORG = 'LickObservatory'
    _SETTINGS_APP = 'NickelFocusGUI'

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Nickel Focus GUI')

        # Each of these is a QWidget subclass -- constructing one here
        # doesn't show anything on screen by itself; a widget only
        # becomes visible once it's placed into a visible parent's
        # layout (below) and that parent is eventually shown.
        self.image_panel = ImagePanel()
        self.curve_panel = FocusCurvePanel()
        self.control_panel = FocusControlPanel()

        # FocusControlPanel stacks enough group boxes that its natural
        # minimum height can exceed a smaller screen's available height;
        # a QScrollArea's own minimum size hint doesn't inherit that, so
        # wrapping it here is what actually lets the window shrink to fit
        # (see _size_to_screen) instead of being held open by the layout.
        # A QScrollArea is itself a widget: it displays exactly one
        # child widget (set via setWidget) inside a scrollable viewport,
        # adding scrollbars only when that child doesn't fit.
        control_scroll = QtWidgets.QScrollArea()
        control_scroll.setWidget(self.control_panel)
        control_scroll.setWidgetResizable(True)
        # `setWidgetResizable(True)` has a real bug: if the splitter ever
        # squeezes the viewport narrower than the panel's own minimum
        # width, the panel doesn't get a horizontal scrollbar as you'd
        # expect -- it snaps to some unrelated, much larger width (seen
        # in practice: a ~340px-minimum panel jumping to ~640px) and gets
        # stuck there. An explicit minimum width on the scroll area itself
        # (a hard floor, unlike a sizeHint) keeps the splitter from ever
        # asking for less than the panel can actually provide, so that
        # code path never triggers. The small buffer accounts for the
        # vertical scrollbar's own width once one appears.
        control_scroll.setMinimumWidth(self.control_panel.minimumSizeHint().width() + 20)

        # QSplitter arranges its child widgets in a row or column with a
        # user-draggable handle between each pair, unlike a plain layout
        # (QVBoxLayout/QHBoxLayout elsewhere in this codebase), whose
        # relative sizes the user can't adjust interactively. Splitters
        # nest like any other widget: `right` (curve plot over the
        # scrollable control panel) becomes one side of `central`
        # (image panel beside that whole right-hand column).
        right = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        right.addWidget(self.curve_panel)
        right.addWidget(control_scroll)

        central = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        central.addWidget(self.image_panel)
        central.addWidget(right)
        # Give both sides of the horizontal splitter equal claim on any
        # extra space when the window is resized (an equal, nonzero
        # stretch factor on each), rather than one side absorbing all of
        # it while the other stays fixed.
        central.setStretchFactor(0, 1)
        central.setStretchFactor(1, 1)

        # A QMainWindow manages one designated central widget (menus,
        # toolbars, and dockable panels -- none used here -- go around
        # it); this is what actually makes `central` (and everything
        # nested inside it) appear when the window is shown.
        self.setCentralWidget(central)
        self._size_to_screen(1200, 600)
        self._load_settings()

    def closeEvent(self, event):
        """
        Save settings before the window actually closes. This overrides
        a `QWidget` *event handler* rather than connecting to a signal
        (contrast `gui.controller.Controller`, which is all signal/slot
        connections) -- Qt calls this method itself whenever the window
        is about to close, however that was triggered (the close button,
        `QWidget.close()`, or quitting the application), so overriding
        it is the standard way to run code at that specific moment.
        """
        self._save_settings()
        super().closeEvent(event)

    def _settings(self):
        """The `QSettings` store this window persists to/restores from."""
        return QtCore.QSettings(self._SETTINGS_ORG, self._SETTINGS_APP)

    def _save_settings(self):
        """
        If "Remember settings between sessions" is unchecked, touch
        nothing -- unless a previous session actually saved something
        (in which case it's erased, honoring an opt-*out*). A fresh
        install/session that never opts in must never cause a settings
        file to be written at all.
        """
        settings = self._settings()
        if not self.control_panel.remember_settings_checkbox.isChecked():
            if 'control_panel' in settings.childGroups():
                settings.remove('control_panel')
                settings.remove('image_panel')
            return
        for key, value in self.control_panel.get_settings_state().items():
            settings.setValue(f'control_panel/{key}', value)
        for key, value in self.image_panel.get_settings_state().items():
            settings.setValue(f'image_panel/{key}', value)

    def _load_settings(self):
        """
        Check the "Remember settings between sessions" box if (and only
        if) a previous session actually saved something -- that presence
        *is* the opt-in state, rather than a separately-persisted flag
        (see the class docstring) -- and only then restore it, falling
        back to whatever the widgets were already constructed with for
        any key never saved before (a field added since the last save).
        Passing that current value as `QSettings.value`'s ``defaultValue``
        also tells it what type to coerce the stored value to, so no
        separate type map is needed here. This is read-only: merely
        checking/reading must never itself create a settings file.
        """
        settings = self._settings()
        remember = 'control_panel' in settings.childGroups()
        self.control_panel.remember_settings_checkbox.setChecked(remember)
        if not remember:
            return

        control_defaults = self.control_panel.get_settings_state()
        control_state = {key: settings.value(f'control_panel/{key}', default)
                          for key, default in control_defaults.items()}
        self.control_panel.set_settings_state(control_state)

        image_defaults = self.image_panel.get_settings_state()
        image_state = {key: settings.value(f'image_panel/{key}', default)
                        for key, default in image_defaults.items()}
        self.image_panel.set_settings_state(image_state)

    def _size_to_screen(self, preferred_width, preferred_height):
        """
        Size the window to ``(preferred_width, preferred_height)``, or
        smaller if the screen's available area (excluding docks/taskbars)
        can't fit that -- so the resize handles are never pushed off
        screen -- then center it.
        """
        screen = self.screen() or QtWidgets.QApplication.primaryScreen()
        if screen is None:
            self.resize(preferred_width, preferred_height)
            return
        available = screen.availableGeometry()
        width = min(preferred_width, int(available.width() * 0.9))
        height = min(preferred_height, int(available.height() * 0.9))
        self.resize(width, height)
        self.move(available.center() - self.rect().center())
