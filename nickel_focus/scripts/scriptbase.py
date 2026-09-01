"""
Implements base classes for use with scripts.
"""
import argparse
import datetime
from pathlib import Path

from IPython import embed

from nickel_focus import log

class ScriptBase:
    """
    Provides a base class for all scripts.
    """
    @classmethod
    def entry_point(cls):
        """
        Defines the main script entry point.
        """
        cls.main(cls.parse_args())

    @classmethod
    def name(cls):
        """
        Provide the name of the script.  By default, this is the name of the
        module with "nickel" prepended.
        """
        return f"nickel_{cls.__module__.split('.')[-1]}"

    @classmethod
    def parse_args(cls, options=None):
        """
        Parse the command-line arguments.
        """
        parser = cls.get_parser()
        cls._fill_parser_cwd(parser)
        return parser.parse_args() if options is None else parser.parse_args(options)

    @staticmethod
    def _fill_parser_cwd(parser):
        """
        Replace the default of any action that is exactly ``'current working
        directory'`` with the value of ``Path.cwd()``.

        The ``parser`` is edited *in place*.

        Args:
            parser (argparse.ArgumentParser):
                The argument parsing object to edit.
        """
        for action in parser._actions:
            if action.default == 'current working directory':
                action.default = str(Path.cwd())

    # Base classes should override this
    @classmethod
    def main(cls, args):
        """
        Execute the script.
        """
        pass

    @classmethod
    def get_parser(
        cls, description=None, width=None, formatter=argparse.ArgumentDefaultsHelpFormatter,
        include_log_options=True, default_log_file=False
    ):
        """
        Construct the command-line argument parser.

        Derived classes should override this.  Ideally they should use this
        base-class method to instantiate the ArgumentParser object and then fill
        in the relevant parser arguments

        .. warning::

            *Any* argument that defaults to the
            string ``'current working directory'`` will be replaced by the
            result of ``Path.cwd()`` when the script is executed.  This means
            help dialogs will include this replacement, and parsing of the
            command line will use ``Path.cwd()`` as the default.

        Parameters
        ----------
        description : :obj:`str`, optional
            A short description of the purpose of the script.
        width : :obj:`int`, optional
            Restrict the width of the formatted help output to be no longer than
            this number of characters, if possible given the help formatter.  If
            None, the width is the same as the terminal width.
        formatter : `argparse.HelpFormatter`_
            Class used to format the help output.
        include_log_options : :obj:`bool`, optional
            Include options that define the logging level(s) and log file.
        default_log_file : :obj:`bool`, optional
            If true, script will use the default log file name if none is
            provided.  Ignored if ``include_log_options`` is False.

        Returns
        -------
        `argparse.ArgumentParser`_
            Command-line interpreter.
        """
        parser = argparse.ArgumentParser(
            description=description, formatter_class=lambda prog: formatter(prog, width=width)
        )
        if not include_log_options:
            return parser
        # Add the logging options
        parser.add_argument(
            '-v', '--verbosity', type=int, default=1,
            help='Verbosity level, which must be 0, 1, or 2.  Level 0 includes warning and error '
                 'messages, level 1 adds informational messages, and level 2 adds debugging '
                 'messages and the calling sequence.'
        )
        parser.add_argument(
            '--log_file', type=str, default='default' if default_log_file else None,
            help='Name for the log file.  If set to "default", a default name is used.  If None, '
                 'a log file is not produced.'
        )
        parser.add_argument(
            '--log_level', type=int, default=None,
            help='Verbosity level for the log file.  If a log file is produce and this is None, '
                 'the file log will match the console stream log.'
        )
        return parser

    @classmethod
    def init_log(cls, args):
        """
        Initialize the logger provided the command-line arguments.
        """
        level = log.convert_verbosity_to_logging_level(args.verbosity)
        log_file_level = None if args.log_level is None else \
            log.convert_verbosity_to_logging_level(args.log_level)
        if args.log_file == 'default':
            _log_file = cls.default_log_file()
        elif args.log_file in ['None', None]:
            _log_file = None
        else:
            _log_file = args.log_file
        log.init(level=level,
                 log_file=_log_file,
                 log_file_level=log_file_level)

    @classmethod
    def default_log_file(cls):
        """
        Set the default name for the log file.
        """
        # Create a UT timestamp (to the minute) for the log filename
        timestamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%d-%H%M")
        return f'{cls.name()}_{timestamp}.log'
