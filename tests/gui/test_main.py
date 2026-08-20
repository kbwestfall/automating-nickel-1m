"""Smoke tests for :mod:`gui.main`'s scaffolding."""
import gui.main
from gui.qt import QtWidgets


def test_build_window_opens_without_error(qapp):
    window = gui.main.build_window()
    window.show()
    qapp.processEvents()

    assert window.windowTitle() == 'Nickel Focus GUI', 'window title not set as expected'

    window.close()


def test_build_window_fits_within_the_screen(qapp):
    # A fixed 1200x800 window can be larger than the available screen
    # (excluding docks/taskbars), pushing its resize handles off screen
    # and leaving the user stuck -- the window should never open larger
    # than what's actually available.
    window = gui.main.build_window()
    window.show()
    qapp.processEvents()

    available = window.screen().availableGeometry()
    assert window.width() <= available.width(), 'window should not be wider than the screen'
    assert window.height() <= available.height(), 'window should not be taller than the screen'

    window.close()


def test_control_scroll_area_has_a_floor_at_the_panels_minimum_width(qapp):
    # QScrollArea.setWidgetResizable(True) has a real bug: if its
    # viewport ever gets squeezed narrower than the scrolled widget's own
    # minimum width, the widget doesn't get a horizontal scrollbar as
    # you'd expect -- it snaps to some unrelated, much larger width (seen
    # in practice: a ~340px-minimum panel jumping to ~640px) and gets
    # stuck there, buttons and all. The mitigation is giving the scroll
    # area itself an explicit minimum width -- a hard floor, unlike a
    # sizeHint -- that the splitter can't ask it to go below, so that
    # code path never triggers. This confirms that floor is actually in
    # place and covers the panel's current minimum.
    window = gui.main.build_window()
    window.show()
    qapp.processEvents()

    control_scroll = window.centralWidget().widget(1).widget(1)
    assert isinstance(control_scroll, QtWidgets.QScrollArea), \
        'setup: expected the control panel to be wrapped in a QScrollArea'
    assert control_scroll.minimumWidth() >= window.control_panel.minimumSizeHint().width(), \
        "the scroll area's minimum width must never be less than the panel's own minimum"

    window.close()
