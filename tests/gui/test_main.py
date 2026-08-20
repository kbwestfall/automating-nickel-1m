"""Smoke tests for :mod:`gui.main`'s scaffolding."""
import gui.main


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
