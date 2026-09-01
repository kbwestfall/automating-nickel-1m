"""
Script that slews the Nickel telescope to the nearest known target in a
starlist.
"""

from nickel_focus.scripts import scriptbase


class NickelSlewToNearest(scriptbase.ScriptBase):
    """
    Slew the telescope to the nearest pointing/focus star, or the
    nearest target in a user-supplied starlist.
    """

    @classmethod
    def get_parser(cls, width=None):
        """
        Construct the command-line argument parser.

        Parameters
        ----------
        width
            Restrict the width of the formatted help output to be no
            longer than this number of characters, if possible; see
            `scriptbase.ScriptBase.get_parser`.

        Returns
        -------
        argparse.ArgumentParser
            Command-line interpreter.
        """
        parser = super().get_parser(
            description='Slew the telescope to the nearest target', width=width
        )

        parser.add_argument(
            '-f', '--starlist', default=None, type=str,
            help=(
                'Name of a starlist file with the possible targets.  If not provided, defaults '
                'to the list of pointing and focusing stars.'
            )
        )
        parser.add_argument(
            '-s', '--search', default=None,
            help=(
                'Search string used to down-select the target list.  For example, to find the '
                'nearest pointing star in the default star list, use "-s Pointing".'
            )
        )
        parser.add_argument(
            '--dry_run', default=False, action='store_true',
            help=(
                'Run in dry run mode; i.e., the script will print the nearest object, but the '
                'telescope will not move.'
            )
        )
        return parser

    @classmethod
    def main(cls, args):
        """
        Execute the script: find the nearest target and, unless running
        in dry-run mode, slew the telescope to it.

        Parameters
        ----------
        args
            Parsed command-line arguments, as returned by
            `get_parser`/`~nickel_focus.scripts.scriptbase.ScriptBase.parse_args`.
        """
        from nickel_focus import slew

        cls.init_log(args)

        telescope = slew.NickelTelescopePointing()

        target_name, target_ra, target_dec = slew.find_nearest_target(
            telescope.current, obj_search_str=args.search, file=args.starlist
        )
        print(f'Nearest object is: {target_name} at RA={target_ra}, Dec={target_dec}')

        if not args.dry_run:
            telescope.slew_to(target_ra, target_dec)
