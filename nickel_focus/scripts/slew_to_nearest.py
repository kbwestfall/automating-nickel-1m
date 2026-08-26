
from nickel_focus.scripts import scriptbase


class NickelSlewToNearest(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
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

        from nickel_focus import log
        from nickel_focus import slew

        cls.init_log(args)

        telescope = slew.NickelTelescopePointing()

        target_name, target_ra, target_dec = slew.find_nearest_star(
            telescope.current, obj_search_str=args.search, file=args.starlist
        )
        log.info(f'Nearest object is: {target_name} at RA={target_ra}, Dec={target_dec}')

        if not args.dry_run:
            telescope.slew_to(target_ra, target_dec)
