#! @KPYTHON@

from nickel_focus.scripts import scriptbase

class NickelFocusGUI(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Focusing GUI for the 1m Nickel Telescope.', width=width
        )
        return parser

    @classmethod
    def main(cls, args):
        """
        Build the window and its `Controller`, show the window, and run the
        Qt event loop until the user quits.
        """

        from nickel_focus.gui import launcher
        from nickel_focus.gui.controller import Controller

        app = launcher.build_app()
        window = launcher.build_window()
        # `controller` looks unused (hence the noqa), but it must stay alive
        # for as long as the window does: PySide6 only keeps a signal
        # connection working while *something* in Python still holds a
        # reference to the connected object. Once `main()`'s local variable
        # `controller` would normally go out of scope, nothing in Python
        # would reference this Controller instance any more (Qt's C++ side
        # doesn't parent it to the window automatically), and it -- along
        # with every signal connection it made in `Controller.__init__` --
        # could be garbage-collected out from under the running window.
        controller = Controller(window)  # noqa: F841 (kept alive for the life of main())
        window.show()
        # QApplication.exec() starts Qt's event loop: it blocks here,
        # dispatching mouse clicks, key presses, timers, and cross-thread
        # signal deliveries (e.g. from SequenceWorker) to the right widgets
        # one at a time, until the application quits (e.g. the last window
        # closes). Nothing in this GUI runs unless this loop is running --
        # that's true of every Qt application, not just this one.
        app.exec()
