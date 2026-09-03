"""
Live FWHM-vs-focus plot for a focus sequence.

See GUI_DESIGN.md §5.3.
"""
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from nickel_focus import quadratic
from nickel_focus.gui.qt import QtWidgets


class FocusCurvePanel(QtWidgets.QWidget):
    """
    Scatter of FWHM versus focus value for a sequence's results so far, with the
    best-fit quadratic and its vertex overlaid once there are enough points to
    fit one (matching :meth:`~nickel_focus.focus.FocusSequence.fit_best_focus`'s
    own 3-point minimum). Outlier points
    (:attr:`~nickel_focus.focus.StepResult.is_outlier`) are drawn as a distinct
    color/marker from normal ones.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._results = []

        # 'constrained' layout recomputes margins on every draw so the
        # axis labels/title always have enough room, regardless of how
        # small the splitter/window makes this panel -- unlike a fixed
        # margin, which gets clipped once the panel shrinks below it.
        self.figure = Figure(figsize=(5, 4), layout='constrained')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot(111)
        self._reset_axis()

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.canvas)

    def add_result(self, result):
        """Add one more :class:`~nickel_focus.focus.StepResult` and redraw."""
        self._results.append(result)
        self._render()

    def reset(self):
        """Clear the plot: no points, no fitted curve."""
        self._results = []
        self._reset_axis()
        self.canvas.draw_idle()

    def _reset_axis(self):
        """
        Clear the plot and redraw its static decorations (labels, title, grid)
        -- everything
        :meth:`~nickel_focus.gui.views.focus_curve_panel.FocusCurvePanel._render`
        would otherwise have to re-apply after every
        :meth:`matplotlib.axes.Axes.clear`, factored out since both
        :meth:`~nickel_focus.gui.views.focus_curve_panel.FocusCurvePanel._render`
        and
        :meth:`~nickel_focus.gui.views.focus_curve_panel.FocusCurvePanel.reset`
        need exactly this.
        """
        self.ax.clear()
        self.ax.set_xlabel('Focus Value')
        self.ax.set_ylabel('FWHM (pixels)')
        self.ax.set_title('Focus Curve')
        self.ax.grid(True, alpha=0.3)

    def _render(self):
        """
        Redraw the whole plot from
        :attr:`~nickel_focus.gui.views.focus_curve_panel.FocusCurvePanel._results`:
        points, plus the fitted curve/vertex if enough exist.
        """
        self._reset_axis()
        if not self._results:
            self.canvas.draw_idle()
            return

        normal = [r for r in self._results if not r.is_outlier]
        outliers = [r for r in self._results if r.is_outlier]

        if normal:
            self.ax.scatter([r.focus_value for r in normal], [r.fwhm for r in normal],
                             color='tab:blue', marker='o', label='Measured')
        if outliers:
            self.ax.scatter([r.focus_value for r in outliers], [r.fwhm for r in outliers],
                             color='gold', marker='x', label='Outlier')

        focus_values = [r.focus_value for r in self._results]
        fwhm_values = [r.fwhm for r in self._results]
        if len(self._results) >= 3:
            try:
                a, b, c = quadratic.fit_quadratic(focus_values, fwhm_values)
                x_vertex, y_vertex = quadratic.vertex(a, b, c)
            except (ValueError, ZeroDivisionError):
                pass
            else:
                x_smooth = np.linspace(min(focus_values), max(focus_values), 50)
                y_smooth = a*x_smooth**2 + b*x_smooth + c
                self.ax.plot(x_smooth, y_smooth, 'r-')
                self.ax.scatter([x_vertex], [y_vertex], color='green', zorder=3,
                                 label=f'Best focus: {x_vertex:.1f}')

        self.ax.legend(loc='best')
        self.canvas.draw_idle()
