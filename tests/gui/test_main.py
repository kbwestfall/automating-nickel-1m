"""Smoke tests for :mod:`gui.main`'s scaffolding."""
import gui.main


def test_build_window_opens_without_error(qapp):
    window = gui.main.build_window()
    window.show()
    qapp.processEvents()

    assert window.windowTitle() == 'Nickel Focus GUI', 'window title not set as expected'

    window.close()
