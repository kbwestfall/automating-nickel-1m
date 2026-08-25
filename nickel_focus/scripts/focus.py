#! @KPYTHON@

from nickel_focus.scripts import scriptbase

class NickelFocus(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(
            description='Automated focus measurement for the 1m Nickel Telescope.', width=width
        )
        parser.add_argument(
            'focus', nargs='+', type=float, help=(
                'Focus starting value and step size.  You can also provide the last focus value.  '
                'Use --n to set the number of steps instead of the last focus value.  If neither '
                'are provided (last focus or number of focus steps), the code performs an '
                'automated number of steps to fill out the focus curve.'
            )
        )
        parser.add_argument(
            '-n', '--nstep', default=None, type=int, help=(
                'The number of focus steps to perform.  Ignored if the ending focus value is '
                'provided (see --focus).'
            )
        )
#        parser.add_argument(
#            '--method', type=str, nargs='+', default='brightest', help=(
#                'The method used for calculating the image quality.  Must be "brightest" to use '
#                'the brightest detected source in the field, "weighted" to use a weighted mean of '
#                'all detected sources, or provide two pixel coordinates (x,y or column,row) to '
#                'use a the detected source closest to the provided coordinates.'
#            )
#        )
        parser.add_argument(
            '--maxsteps', type=int, default=12, help=(
                'If using the automated focus sequence, this is the maximum of steps that are '
                'allowed.'
            )
        )
        parser.add_argument(
            '--obsnum', type=int, default=None, help=(
                'Re-analyze a focus sequence starting with the provided focus values and this '
                'observation number.  The number of available images must match the focus '
                'sequence requested.  Requires --datadir, --prefix, and --suffix to locate the '
                'exposures; does not require a ktl connection.'
            )
        )
        parser.add_argument(
            '--datadir', default='.', type=str, help=(
                'Directory holding the archived exposures to re-analyze.  Only used with --obsnum.'
            )
        )
        parser.add_argument(
            '--prefix', default='n', type=str, help=(
                'Filename prefix of the archived exposures, e.g. "n" for files named like '
                '"n2165.fits".  Only used with --obsnum.'
            )
        )
        parser.add_argument(
            '--suffix', default='.fits', type=str, help=(
                'Filename suffix of the archived exposures, e.g. ".fits".  Only used with '
                '--obsnum.'
            )
        )
        parser.add_argument(
            '-t', '--exptime', default=5, type=float,
            help='Exposure time in seconds for each exposure.'
        )
        parser.add_argument(
            '-b', '--binning', default=None, choices=['1,1', '2,2', '4,4'],
            help='BinningExposure speed.  Must be "1,1", "2,2", or "4.4".'
        )
        parser.add_argument(
            '-s', '--speed', default=None, choices=['Slow', 'Fast'],
            help='Exposure speed.  Must be Slow or Fast.  If None'
        )
        parser.add_argument(
            '-o', '--ofile', default=None, type=str,help=(
                'Output file for the measured focus data.  This can be used to exclude and refit '
                'the best focus.'
            )
        )
        parser.add_argument(
            '--refit', action='store_true',
            help='Refit the focus curve.  The output file (see --ofile) must be provided.'
        )
        parser.add_argument(
            '--omit', type=int, nargs='+', default=None,
            help='List of observation numbers to omit from the curve fitting.'
        )
        parser.add_argument(
            '--config', action='store_true',
            help='Only print the current exposure configuration and focus'
        )
        parser.add_argument(
            '--verbose', action='store_true',
            help='Enable verbose output for debugging'
        )
        parser.add_argument(
            '--no-plot', action='store_true',
            help='Disable the live focus-sequence plot.'
        )
        return parser

    @classmethod
    def main(cls, args):

#        import warnings
        from pathlib import Path

        from astropy.table import Table
        from IPython import embed
        from matplotlib import pyplot
        import numpy as np

        from nickel_focus import focus

        # Set the read speed
        _speed = args.speed
#        if _speed == 'Fast':
#            warnings.warn('Fast does not work!  Setting to slow.')
#            _speed == 'Slow'
#        if _speed == 'Slow':
#            _speed = '0.05MHz'
#        elif _speed == 'Fast':
#            _speed = '1.0MHz'

        if args.refit:
            raise NotImplementedError('Not ready to refit.')
            if args.ofile is None:
                raise ValueError(
                    'To refit, must provide output file name from a previous focus sequence.'
                )
            _ofile = Path(args.ofile).absolute()
            if not _ofile.is_file():
                raise FileNotFoundError(f'{_ofile.name} not found!  Correct the output file name.')
            tbl = Table.read(_ofile, format='ecsv')
            if args.omit is not None:
                # Find the entries in the table and remove them
                pass
            best_focus, best_img_quality = FocusSequence.fit_best_focus(tbl['FOCUS'], tbl['SIGMA'])
            # And then finish this out
            return
    
        # Check if output file exists

        if len(args.focus) == 3 or args.nstep is not None:
            # Perform a grid or archive sequence
            end = None if len(args.focus) == 2 else args.focus[2]
            seq = focus.GridFocusSequence(args.focus[0], args.focus[1], end=end, nstep=args.nstep)

            if args.obsnum is not None:
                # Construct the expected exposure filenames directly, instead of
                # asking ktl for the recording directory/prefix/suffix, so that
                # archived sequences can be re-analyzed without a ktl connection.
                datadir = Path(args.datadir).absolute()
                expected_files = np.array([
                    datadir / f'{args.prefix}{args.obsnum + i}{args.suffix}'
                    for i in range(seq.nstep)
                ])
                indx = np.array([f.is_file() for f in expected_files])
                if not np.all(indx):
                    missing = expected_files[np.logical_not(indx)]
                    raise FileNotFoundError(
                        'Expected to find the following files, but they are not available: '
                        f'{", ".join(str(f) for f in missing)}'
                    )
                seq = focus.ArchiveFocusSequence(seq.target_focus, expected_files)
        else:
            seq = focus.AutomatedFocusSequence(
                args.focus[0], args.focus[1], maxsteps=args.maxsteps
            )

        best_focus, best_img_quality = seq.execute(
            goto=False, method='brightest', record=True, speed=_speed, exptime=args.exptime,
            binning=args.binning, plot=not args.no_plot
        )
        print(f'Best focus: {best_focus:.1f}')
        print(f'Expected sigma: {best_img_quality:.1f} pixels')

        if args.ofile is not None:
            # Outlier status is recomputed the same way step() does it live:
            # relative to every centroid observed up to and including that
            # one, not the full, final set.
            outlier = [
                focus.FocusPlot.is_outlier(seq.centroids[:i+1]) for i in range(len(seq.centroids))
            ]
            tbl = Table({
                'EXPOSURE': [str(f) for f in seq.exposures],
                'FOCUS': seq.observed_focus,
                'SIGMA': seq.img_quality,
                'OUTLIER': outlier,
            })
            tbl.write(args.ofile, format='ecsv', overwrite=True)
            print(f'Wrote focus data to {args.ofile}')

        if seq.plot is not None:
            pyplot.ioff()
            pyplot.show(block=True)
