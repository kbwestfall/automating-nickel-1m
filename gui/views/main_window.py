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

        # FocusControlPanel stacks enough group boxes that its natural
        # minimum height can exceed a smaller screen's available height;
        # a QScrollArea's own minimum size hint doesn't inherit that, so
        # wrapping it here is what actually lets the window shrink to fit
        # (see _size_to_screen) instead of being held open by the layout.
        control_scroll = QtWidgets.QScrollArea()
        control_scroll.setWidget(self.control_panel)
        control_scroll.setWidgetResizable(True)

        right = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        right.addWidget(self.curve_panel)
        right.addWidget(control_scroll)

        central = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        central.addWidget(self.image_panel)
        central.addWidget(right)
        central.setStretchFactor(0, 1)
        central.setStretchFactor(1, 1)

        self.setCentralWidget(central)
        self._size_to_screen(1200, 600)

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
