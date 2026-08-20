"""
Interactive display of the exposures in a focus sequence.

See GUI_DESIGN.md §5.2 and §5.6.
"""
import bisect

from matplotlib import patches
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from astropy.visualization import ImageNormalize, ZScaleInterval, MinMaxInterval, LinearStretch

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

    Clicking on a displayed source (§5.6) emits :attr:`sourceSelected`
    with the clicked data coordinates, but only while
    :func:`set_selection_enabled` has been called with ``True`` --
    gating that on whether a sequence is currently running is the
    Controller's job, not this panel's.

    Does not implement the design doc's optional "center on measured
    source" recenter variant; "Recenter" always fits the whole frame.
    """

    #: Emitted with the clicked ``(x, y)`` data coordinates, only while
    #: selection is enabled (see :func:`set_selection_enabled`).
    sourceSelected = QtCore.Signal(float, float)

    #: Available stretch options, by name, as callables returning an
    #: `astropy.visualization.ImageNormalize` for a given data array.
    STRETCHES = {
        'ZScale': lambda data: ImageNormalize(data, interval=ZScaleInterval(),
                                               stretch=LinearStretch()),
        'Min/Max': lambda data: ImageNormalize(data, interval=MinMaxInterval(),
                                                stretch=LinearStretch()),
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
        self.canvas = FigureCanvasQTAgg(self.figure)
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
        top.addWidget(self.exposure_combo, 1)
        top.addWidget(QtWidgets.QLabel('Stretch:'))
        top.addWidget(self.stretch_combo)
        top.addWidget(self.recenter_button)

        grid = QtWidgets.QGridLayout()
        grid.addWidget(self.canvas, 0, 0)
        grid.addWidget(self.v_scroll, 0, 1)
        grid.addWidget(self.h_scroll, 1, 0)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(top)
        layout.addLayout(grid)

        self.exposure_combo.currentIndexChanged.connect(self._on_combo_changed)
        self.stretch_combo.currentTextChanged.connect(self._on_stretch_changed)
        self.recenter_button.clicked.connect(self.recenter)
        self.h_scroll.valueChanged.connect(self._on_scroll)
        self.v_scroll.valueChanged.connect(self._on_scroll)
        self.canvas.mpl_connect('scroll_event', self._on_wheel_zoom)
        self.canvas.mpl_connect('button_press_event', self._on_click)

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

        label = f'{result.exposure.stem} — Focus {result.focus_value:.1f}'
        self.exposure_combo.blockSignals(True)
        self.exposure_combo.insertItem(insert_at, label)
        if tracking_latest:
            self.exposure_combo.setCurrentIndex(insert_at)
        self.exposure_combo.blockSignals(False)

        if tracking_latest:
            self.show_result(result)

    def show_result(self, result):
        """Display ``result``, resetting to a fit-to-view zoom."""
        self._current = result
        self._zoom = 1.0
        self._render(reset_view=True)

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

    def set_zoom(self, zoom):
        """Set the display magnification (``1.0`` fits the whole frame)."""
        self._zoom = max(1.0, float(zoom))
        if self._current is not None:
            self._update_scrollbar_ranges()
            self._apply_view_limits()
            self.canvas.draw_idle()

    # -- internals --------------------------------------------------------

    def _on_combo_changed(self, index):
        if 0 <= index < len(self._results):
            self.show_result(self._results[index])

    def _on_stretch_changed(self, name):
        self._stretch_name = name
        if self._current is not None:
            self._render(reset_view=False)

    def _on_scroll(self, _value):
        if self._current is not None:
            self._apply_view_limits()
            self.canvas.draw_idle()

    def _on_wheel_zoom(self, event):
        if event.inaxes != self.ax or self._current is None:
            return
        if event.button == 'up':
            self.set_zoom(self._zoom * 1.2)
        elif event.button == 'down':
            self.set_zoom(self._zoom / 1.2)

    def _on_click(self, event):
        if not self._selection_enabled or event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        self.sourceSelected.emit(float(event.xdata), float(event.ydata))

    def _render(self, reset_view):
        result = self._current
        self.ax.clear()
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        if result is None:
            self.canvas.draw_idle()
            return

        data = result.frame
        normalize = self.STRETCHES[self._stretch_name](data)
        self.ax.imshow(data, cmap='viridis', origin='lower', norm=normalize)
        self.ax.set_title(f'{result.exposure.stem}  Focus: {result.focus_value:.1f}')

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

    def _update_scrollbar_ranges(self):
        ny, nx = self._data_shape
        view_w, view_h = self._view_extent()
        for bar, extent, full in ((self.h_scroll, view_w, nx), (self.v_scroll, view_h, ny)):
            bar.blockSignals(True)
            bar.setRange(0, max(0, round(full - extent)))
            bar.setPageStep(max(1, round(extent)))
            bar.setValue(min(bar.value(), bar.maximum()))
            bar.blockSignals(False)

    def _apply_view_limits(self):
        view_w, view_h = self._view_extent()
        x0 = self.h_scroll.value()
        y0 = self.v_scroll.value()
        self.ax.set_xlim(x0, x0 + view_w)
        self.ax.set_ylim(y0, y0 + view_h)
