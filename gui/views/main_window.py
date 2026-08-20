"""
Top-level window laying out the three panels.

See GUI_DESIGN.md §5.1.
"""
from gui.qt import QtCore, QtWidgets
from gui.views.image_panel import ImagePanel
from gui.views.focus_curve_panel import FocusCurvePanel
from gui.views.focus_control_panel import FocusControlPanel


class MainWindow(QtWidgets.QMainWindow):
    """
    Lays out `ImagePanel` (left), `FocusCurvePanel` (top right), and
    `FocusControlPanel` (bottom right) per §5.1's sketch. Purely
    structural -- wiring the panels together is `gui.controller.Controller`'s
    job, not this class's.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Nickel Focus GUI')

        self.image_panel = ImagePanel()
        self.curve_panel = FocusCurvePanel()
        self.control_panel = FocusControlPanel()

        right = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        right.addWidget(self.curve_panel)
        right.addWidget(self.control_panel)

        central = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        central.addWidget(self.image_panel)
        central.addWidget(right)
        central.setStretchFactor(0, 1)
        central.setStretchFactor(1, 1)

        self.setCentralWidget(central)
        self.resize(1200, 800)
