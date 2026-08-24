"""
Interactive display of the exposures in a focus sequence.

See GUI_DESIGN.md §5.2 and §5.6.
"""
import bisect

from matplotlib import patches
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from astropy.visualization import (
    ImageNormalize, ZScaleInterval, MinMaxInterval, LinearStretch, SqrtStretch, LogStretch,
    AsinhStretch,
)

from gui.qt import QtCore, QtWidgets


class ImagePanel(QtWidgets.QWidget):
    """
    Shows the exposures collected by a focus sequence, one at a time.

    Baseline interactions (GUI_DESIGN.md §5.2): pan via scroll bars, wheel
    zoom, a "Recenter" button that returns to a fit-to-view zoom, and a
    stretch selector. The measured source is boxed on every frame -- red
    normally, yellow (with an "Outlier centroid" label) if
    :attr:`focus.StepResult.is_outlier` is set. A drop-down lists every
    exposure added so far, always kept sorted by focus value, and tracks
    the most recently *added* exposure until the user manually picks a
    different one (§5.2's "auto-follow latest" decision) -- note that
    "most recently added" is not the same as "highest focus value," since
    an adaptive (`AutomatedFocusSequence`) run doesn't acquire focus
    values in order.

    Hovering over a displayed source and pressing 'm' (§5.6) emits
    :attr:`sourceSelected` with the data coordinates under the cursor,
    but only while :func:`set_selection_enabled` has been called with
    ``True`` -- gating that on whether a sequence is currently running is
    the Controller's job, not this panel's. A key press rather than a
    plain click is deliberate: on macOS, clicking an unfocused window
    both refocuses it *and* is delivered to whatever's under the cursor,
    so a stray click meant only to bring the window forward was enough
    to kick off a reanalysis.

    Does not implement the design doc's optional "center on measured
    source" recenter variant; "Recenter" always fits the whole frame.
    """

    #: Emitted with the clicked ``(x, y)`` data coordinates, only while
    #: selection is enabled (see :func:`set_selection_enabled`).
    sourceSelected = QtCore.Signal(float, float)

    #: Available stretch options, by name, as callables returning an
    #: `astropy.visualization.ImageNormalize` for a given data array.
    #: The non-linear options (Sqrt/Log/Asinh) pair `ZScaleInterval`'s
    #: already-good limits with a different stretch curve, to bring out
    #: faint structure (Sqrt/Asinh) or compress a large dynamic range
    #: (Log), without needing a separate interval choice for each. Kept
    #: to short, single-word names -- a longer combo entry (e.g. the
    #: originally-tried "ZScale + Sqrt") widens the whole combo box
    #: enough to push the window past a smaller screen's available width
    #: (see claude/GUI_IMPLEMENTATION.md's Phase 5 sub-phase 3 log).
    STRETCHES = {
        'ZScale': lambda data: ImageNormalize(data, interval=ZScaleInterval(),
                                               stretch=LinearStretch()),
        'Min/Max': lambda data: ImageNormalize(data, interval=MinMaxInterval(),
                                                stretch=LinearStretch()),
        'Sqrt': lambda data: ImageNormalize(data, interval=ZScaleInterval(),
                                             stretch=SqrtStretch()),
        'Log': lambda data: ImageNormalize(data, interval=ZScaleInterval(),
                                            stretch=LogStretch()),
        'Asinh': lambda data: ImageNormalize(data, interval=ZScaleInterval(),
                                              stretch=AsinhStretch()),
    }

    def __init__(self, parent=None):
        super().__init__(parent)

        self._results = []          # StepResults, kept sorted by focus_value
        self._latest_result = None  # most recently *added* result (acquisition order)
        self._current = None        # currently displayed result
        self._data_shape = None     # (ny, nx) of the currently displayed frame
        self._zoom = 1.0
        self._stretch_name = next(iter(self.STRETCHES))
        self._selection_enabled = False

        self.figure = Figure(figsize=(5, 5))
        # No ticks/labels/title are ever drawn (the exposure/focus value
        # is already in the drop-down above), so the image can fill the
        # whole figure instead of leaving matplotlib's default margins.
        self.figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
        self.canvas = FigureCanvasQTAgg(self.figure)
        # Needed to actually receive the 'm' key press below.
        self.canvas.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        self.exposure_combo = QtWidgets.QComboBox()
        self.stretch_combo = QtWidgets.QComboBox()
        self.stretch_combo.addItems(list(self.STRETCHES))
        self.recenter_button = QtWidgets.QPushButton('Recenter')

        self.h_scroll = QtWidgets.QScrollBar(QtCore.Qt.Orientation.Horizontal)
        self.v_scroll = QtWidgets.QScrollBar(QtCore.Qt.Orientation.Vertical)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel('Exposure:'))
        # The trailing 1 is a *stretch factor*: when this row is given
        # more width than its widgets need, QHBoxLayout hands the extra
        # space to whichever widgets have a nonzero factor, in
        # proportion to it. Every other widget here defaults to 0 (no
        # extra space), so the combo box is the only one that grows.
        top.addWidget(self.exposure_combo, 1)
        top.addWidget(QtWidgets.QLabel('Stretch:'))
        top.addWidget(self.stretch_combo)
        top.addWidget(self.recenter_button)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.canvas, 0, 0)
        grid.addWidget(self.v_scroll, 0, 1)
        grid.addWidget(self.h_scroll, 1, 0)

        # Passing `self` to the constructor both creates this layout and
        # immediately calls `self.setLayout(...)` with it -- equivalent
        # to (and just a shorthand for) creating it plain and calling
        # `setLayout` afterward, as the panels built from bare `QWidget`s
        # elsewhere in this codebase do. `top`/`grid` above are added as
        # *nested* sub-layouts, not widgets -- layouts, like widgets,
        # can contain other layouts to build up more complex
        # arrangements than a single row/column/grid alone could.
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(top)
        layout.addLayout(grid)

        # Ordinary Qt signal/slot connections (see gui/controller.py for
        # more on what that means).
        self.exposure_combo.currentIndexChanged.connect(self._on_combo_changed)
        self.stretch_combo.currentTextChanged.connect(self._on_stretch_changed)
        self.recenter_button.clicked.connect(self.recenter)
        self.h_scroll.valueChanged.connect(self._on_scroll)
        self.v_scroll.valueChanged.connect(self._on_scroll)
        # `mpl_connect` is matplotlib's *own* callback registry, not a
        # Qt signal -- the embedded figure's mouse/keyboard events don't
        # become Qt signals, so matplotlib-specific events (like a plot
        # scroll or a key press while the canvas has focus) need
        # matplotlib's own connection mechanism instead of `.connect()`.
        self.canvas.mpl_connect('scroll_event', self._on_wheel_zoom)
        self.canvas.mpl_connect('key_press_event', self._on_key_press)

    # -- public API used by the Controller -------------------------------

    def set_selection_enabled(self, enabled):
        """
        Enable or disable click-to-select-source (§5.6). The Controller
        is responsible for calling this appropriately (enabled before a
        sequence/single exposure starts or after one finishes, disabled
        while one is running).
        """
        self._selection_enabled = bool(enabled)

    def add_result(self, result):
        """
        Add a newly-collected :class:`focus.StepResult` to the exposure
        drop-down, inserted to keep it sorted by focus value, and display
        it if this panel is currently tracking the most recently added
        exposure (see the class docstring).
        """
        tracking_latest = self._current is None or self._current is self._latest_result

        focus_values = [r.focus_value for r in self._results]
        insert_at = bisect.bisect_right(focus_values, result.focus_value)
        self._results.insert(insert_at, result)
        self._latest_result = result

        label = f'{result.exposure.stem} — Focus {result.focus_value:.0f}'
        self.exposure_combo.blockSignals(True)
        self.exposure_combo.insertItem(insert_at, label)
        if tracking_latest:
            self.exposure_combo.setCurrentIndex(insert_at)
        self.exposure_combo.blockSignals(False)

        if tracking_latest:
            self.show_result(result)

    def update_result(self, result):
        """
        Replace an existing entry -- matched by ``result.exposure`` -- with
        an updated measurement, without changing its position in the
        drop-down (its focus value, and thus sort position, doesn't
        change during a reanalysis). Falls back to :func:`add_result` if
        no entry for that exposure exists yet. Used when
        :func:`focus.FocusSequence.reanalyze` re-measures exposures
        already on display (§5.6).
        """
        existing_index = next(
            (i for i, r in enumerate(self._results) if r.exposure == result.exposure), None)
        if existing_index is None:
            self.add_result(result)
            return

        was_current = self._current is self._results[existing_index]
        was_latest = self._latest_result is self._results[existing_index]

        self._results[existing_index] = result
        if was_latest:
            self._latest_result = result

        label = f'{result.exposure.stem} — Focus {result.focus_value:.0f}'
        self.exposure_combo.blockSignals(True)
        self.exposure_combo.setItemText(existing_index, label)
        self.exposure_combo.blockSignals(False)

        if was_current:
            self.show_result(result)

    def show_result(self, result):
        """Display ``result``, resetting to a fit-to-view zoom."""
        self._current = result
        self._zoom = 1.0
        self._render(reset_view=True)

    def _show_result_preserving_view(self, result):
        """
        Display ``result`` keeping the current zoom and view center fixed
        -- used when the user picks a different exposure from the
        drop-down, so browsing between already-collected frames doesn't
        keep resetting back to a fit-to-view zoom.
        """
        if self._current is None or self._data_shape is None:
            self.show_result(result)
            return

        old_view_w, old_view_h = self._view_extent()
        center_x = self.h_scroll.value() + old_view_w / 2
        center_y = self.v_scroll.value() + old_view_h / 2

        self._current = result
        self._data_shape = result.frame.shape
        self._update_scrollbar_ranges()

        new_view_w, new_view_h = self._view_extent()
        new_x0 = round(center_x - new_view_w / 2)
        new_y0 = round(center_y - new_view_h / 2)
        self._set_scroll_position(new_x0, new_y0)
        self._render(reset_view=False)

    def reset(self):
        """Clear the panel: no exposures, nothing displayed."""
        self._results = []
        self._latest_result = None
        self._current = None
        self._data_shape = None
        self.exposure_combo.blockSignals(True)
        self.exposure_combo.clear()
        self.exposure_combo.blockSignals(False)
        self.ax.clear()
        self.canvas.draw_idle()

    def recenter(self):
        """Reset to a fit-to-view zoom (GUI_DESIGN.md §5.2)."""
        self.set_zoom(1.0)

    def get_settings_state(self):
        """
        Return the persistable display preference (currently just the
        stretch choice) as a flat :obj:`dict`, suitable for a settings
        store like :class:`~PySide6.QtCore.QSettings`.
        """
        return {'stretch': self._stretch_name}

    def set_settings_state(self, state):
        """Apply a settings :obj:`dict` from :func:`get_settings_state`."""
        name = state.get('stretch')
        if name in self.STRETCHES:
            self.stretch_combo.setCurrentText(name)

    def set_zoom(self, zoom):
        """Set the display magnification (``1.0`` fits the whole frame)."""
        self._zoom = max(1.0, float(zoom))
        if self._current is not None:
            self._update_scrollbar_ranges()
            self._apply_view_limits()
            self.canvas.draw_idle()

    # -- internals --------------------------------------------------------

    def _on_combo_changed(self, index):
        """Slot for `exposure_combo`'s ``currentIndexChanged``: display the newly picked entry."""
        if 0 <= index < len(self._results):
            self._show_result_preserving_view(self._results[index])

    def _on_stretch_changed(self, name):
        """Slot for `stretch_combo`'s ``currentTextChanged``: re-render with the new stretch."""
        self._stretch_name = name
        if self._current is not None:
            self._render(reset_view=False)

    def _on_scroll(self, _value):
        """
        Slot for both scrollbars' ``valueChanged``. The new value itself
        (``_value``) isn't needed -- :func:`_apply_view_limits` reads
        both bars' current positions directly -- but a slot connected to
        a signal that carries an argument must still accept it
        positionally; the leading underscore is just this codebase's way
        of flagging "intentionally unused."
        """
        if self._current is not None:
            self._apply_view_limits()
            self.canvas.draw_idle()

    def _on_wheel_zoom(self, event):
        """
        matplotlib ``scroll_event`` callback (mouse wheel or trackpad
        scroll over the canvas): zoom in/out, anchored on the cursor
        (see :func:`_zoom_at`).
        """
        if event.inaxes != self.ax or self._current is None:
            return
        if event.button == 'up':
            factor = 1.2
        elif event.button == 'down':
            factor = 1 / 1.2
        else:
            return
        self._zoom_at(event.xdata, event.ydata, self._zoom * factor)

    def _zoom_at(self, xdata, ydata, new_zoom):
        """
        Like :func:`set_zoom`, but keeps the data point ``(xdata, ydata)``
        fixed under the cursor rather than anchored at the view's
        top-left corner -- otherwise every scroll/trackpad zoom event
        visibly drifts the image out from under the cursor.
        """
        new_zoom = max(1.0, float(new_zoom))
        if new_zoom == self._zoom:
            return
        old_view_w, old_view_h = self._view_extent()
        frac_x = (xdata - self.h_scroll.value()) / old_view_w
        frac_y = (ydata - self.v_scroll.value()) / old_view_h

        self._zoom = new_zoom
        self._update_scrollbar_ranges()
        new_view_w, new_view_h = self._view_extent()
        new_x0 = round(xdata - frac_x * new_view_w)
        new_y0 = round(ydata - frac_y * new_view_h)
        self._set_scroll_position(new_x0, new_y0)
        self._apply_view_limits()
        self.canvas.draw_idle()

    def _on_key_press(self, event):
        """
        matplotlib ``key_press_event`` callback: pressing 'm' over the
        canvas, while selection is enabled, emits :attr:`sourceSelected`
        with the data coordinates under the cursor (see the class
        docstring for why a key press rather than a click).
        """
        if event.key != 'm' or not self._selection_enabled or event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        self.sourceSelected.emit(float(event.xdata), float(event.ydata))

    def _render(self, reset_view):
        """
        Redraw the currently displayed result from scratch: image,
        stretch, and the measured-source box (or nothing, if
        `_current` is `None`). ``reset_view`` controls whether the
        zoom/scroll state is recomputed for a new frame shape (a
        genuinely new image) or left alone (e.g. just changing the
        stretch on the same image).
        """
        result = self._current
        self.ax.clear()
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        if result is None:
            self.canvas.draw_idle()
            return

        data = result.frame
        normalize = self.STRETCHES[self._stretch_name](data)
        self.ax.imshow(data, cmap='gray', origin='lower', norm=normalize)

        stamp_size = int(result.fwhm * 10)
        half = stamp_size / 2
        color = 'yellow' if result.is_outlier else 'red'
        rect = patches.Rectangle(
            (result.centroid[0] - half, result.centroid[1] - half), stamp_size, stamp_size,
            linewidth=2, edgecolor=color, facecolor='none')
        self.ax.add_patch(rect)
        if result.is_outlier:
            self.ax.text(result.centroid[0], result.centroid[1] + half + 5, 'Outlier centroid',
                          color='yellow', fontsize=10, ha='center')

        if reset_view:
            self._data_shape = data.shape
            self._update_scrollbar_ranges()
        self._apply_view_limits()
        self.canvas.draw_idle()

    def _view_extent(self):
        """Return ``(view_w, view_h)`` in data pixels at the current zoom."""
        ny, nx = self._data_shape
        return nx / self._zoom, ny / self._zoom

    def _set_scroll_position(self, x0, y0):
        """
        Set both scrollbars' values at once, clamped to their valid
        range, without triggering :func:`_on_scroll` for either -- used
        whenever this panel repositions the view itself (recentering on
        a data point, panning by right-clicking, etc.) rather than in
        response to the user actually dragging a scrollbar. Signals are
        blocked because :func:`_on_scroll` would otherwise fire once per
        bar and redraw twice for what's really one logical view change;
        the caller is expected to refresh the display itself afterward.
        """
        for bar, value in ((self.h_scroll, x0), (self.v_scroll, y0)):
            bar.blockSignals(True)
            bar.setValue(max(0, min(value, bar.maximum())))
            bar.blockSignals(False)

    def _update_scrollbar_ranges(self):
        """
        Recompute each scrollbar's valid range/page size for the
        current frame shape and zoom level (e.g. after loading a new
        frame or changing zoom) -- signals are blocked for the same
        reason as :func:`_set_scroll_position`.
        """
        ny, nx = self._data_shape
        view_w, view_h = self._view_extent()
        for bar, extent, full in ((self.h_scroll, view_w, nx), (self.v_scroll, view_h, ny)):
            bar.blockSignals(True)
            bar.setRange(0, max(0, round(full - extent)))
            bar.setPageStep(max(1, round(extent)))
            bar.setValue(min(bar.value(), bar.maximum()))
            bar.blockSignals(False)

    def _apply_view_limits(self):
        """
        Set the matplotlib axes' visible data range from the current
        scrollbar positions and zoom level. Doesn't redraw the canvas
        itself -- callers are responsible for that (`draw_idle`), since
        this is often one of several changes made before a single
        redraw.
        """
        view_w, view_h = self._view_extent()
        x0 = self.h_scroll.value()
        y0 = self.v_scroll.value()
        self.ax.set_xlim(x0, x0 + view_w)
        self.ax.set_ylim(y0, y0 + view_h)
