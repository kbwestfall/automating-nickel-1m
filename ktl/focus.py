#! @KPYTHON@

import warnings
from IPython import embed

from pathlib import Path
import logging
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import quadratic

from astropy.table import Table

from photometry import image_quality
from photometry import Grid

import ktl


def setup_logging(log_level='INFO', log_file=None):
    """
    Setup logging configuration
    
    Parameters:
    log_level: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    log_file: Optional filename to save logs to file
    """
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s'
    )
    
    # Setup root logger
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, log_level.upper()))
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (optional)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)  
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


class Focus:
    def __init__(self):
        self.pocstop = ktl.cache('nickelpoco', 'POCSTOP')
        self.secpa = ktl.cache('nickelpoco', 'POCSECPA')
        self.secpd = ktl.cache('nickelpoco', 'POCSECPD')
        self.seclk = ktl.cache('nickelpoco', 'POCSECLK')
        self.expstate = ktl.cache('nscicam', 'EXPSTATE')

    @property
    def current(self):
        """
        Return the current focus position
        """
        return self.secpa.read()

    def set_to(self, focus_value):
        """
        Set the focus to the provided value.

        Movement must be enabled and there must not be an ongoing exposure.

        Parameters
        ----------
        focus_value : :obj:`int`
            The requested focus value.  Must be between 165 and 500, inclusive.
            If the requested position is already within 0.1 of the current
            value, no change is made.

        Raises
        ------
        ValueError
            Raised if the focus value is outside the allowed range, if movement
            is disabled, or if the exposure state is anything except that the
            camera is ready for another exposure to begin.
        """

        # Check that the requested focus value is valid
        if focus_value < 165 or focus_value > 500:
            raise ValueError(f'Focus value {focus_value} is out of range (165-500).')
        
        # Make sure movement is enabled.  Do NOT enable movement via this
        # script!
        if not self.pocstop.waitFor('== allowed', timeout=0.5):
            raise ValueError('Telescope movement is disabled!')

        # Check that an exposure isn't currently happening
        if not self.expstate.waitFor('== Ready', timeout=0.5):
            raise ValueError('Camera exposure state not ready. Cannot change focus.')

        print(f'Current POCSECPD: {self.secpd.read()}')
        print(f'Current POCSECPA: {self.secpa.read()}')

        if abs(float(self.secpd.read()) - focus_value) < .1:
            print(f'POCSECPD already set to {focus_value}. No change needed.')
            return
        
        print(f'Current POCSECLK: {self.seclk.read()}')
        self.seclk.write('off')
        self.seclk.read()
        print('Set POCSECLK to off')

        print(f'POCSECPA: {self.secpa.read()}')
        self.secpd.write(focus_value)
        print(f'POCSECPD: {self.secpd.read()}')

        if not self.seclk.waitFor('== on', timeout=30):
            # TODO: Explicitly set the lock to on?
            raise ValueError("POCSECLK did not turn on. Focus change failed.")
            
        print(f"Successfully changed focus to {focus_value}")


class ExposurePath:
    def __init__(self):
        self.exprec = ktl.cache('nscicam', 'RECORD')
        self.expresult = ktl.cache('nscicam', 'EXPRESULT')
        self.scratchdir = ktl.cache('nscicam', 'SCRATCHDIR')
        self.recorddir = ktl.cache('nscicam', 'RECORDDIR')
        self.prefix = ktl.cache('nscicam', 'FITSPREFIX')
        self.obsnum = ktl.cache('nscicam', 'OBSNUM')
        self.suffix = ktl.cache('nscicam', 'FITSSUFFIX')

    @property
    def previous(self):
        return self.expresult.read()

    @property
    def next(self):
        return self.for_obsnum(self.obsnum.read())

    def for_obsnum(self, obsnum, assume_recorded=False):
        # TODO: Is this the correct way to check the keyword has a given value?
        record = self.exprec.read() == 'Yes'
        if record or assume_recorded:
            path = Path(self.recorddir.read()).absolute()
        else:
            path = Path(self.scratchdir.read()).absolute()
        return str(path / f'{self.prefix.read()}{obsnum}{self.suffix.read()}')
    

class ExposureConfig:
    def __init__(self):
        self.exprec = ktl.cache('nscicam', 'RECORD')
        self.inttime = ktl.cache('nscicam', 'EXPOSURE')
        self.expspeed  = ktl.cache('nscicam', 'AMPCONF')
        self.expbin  = ktl.cache('nscicam', 'BINNING')

    def configure(self, record=None, speed=None, binning=None, exptime=None):
        if record is not None:
            self.exprec.write(record)
        if speed is not None:
            self.expspeed.write(speed)
        if binning is not None:
            self.expbin.write(binning)
        if exptime is not None:
            self.inttime.write(exptime)

    @property
    def exptime(self):
        return self.inttime.read()

    def __repr__(self):
        return (
            'Exposure settings:\n'
            f'    Record: {self.exprec.read()}\n'
            f'    Time: {self.inttime.read()}\n'
            f'    Speed: {self.expspeed.read()}\n'
            f'    Binning: {self.expbin.read()}\n'
        )
    
class Exposure:

    def __init__(self):
    
        self.path = ExposurePath()
        self.cfg = ExposureConfig()

        # SCICAM exposure keywords
        self.expstate = ktl.cache('nscicam', 'EXPSTATE')
        self.expstate.callback(self._expstate_callback)
        self.expstate.monitor()
        self.expstate_value = None

        self.expstart = ktl.cache('nscicam', 'EXPOSE')
        self.expstart.callback(self._filepath_callback)
        self.expstart.monitor()
        self.filepath = None

    def _expstate_callback(self, keyword):
        self.expstate_value = keyword.read()
#        print(f'EXPSTATE updated: {self.expstate_value}')

    def _filepath_callback(self, keyword):
        self.filepath = self.path.next
#        print(f'FILEPATH updated: {self.filepath}')

    def expose(self, record=None, speed=None, binning=None, exptime=None):

        self.cfg.configure(record=record, speed=speed, binning=binning, exptime=exptime)

        # Check that an exposure isn't currently happening
        if not self.expstate.waitFor('== Ready', timeout=15):
            raise ValueError('Camera exposure state not ready. Cannot take exposure.')

        # Start the exposure
        self.expstart.write('StartX')

        # Wait for it to start
        if not self.expstate.waitFor('== Start', timeout=30):
            raise ValueError('Exposure start (EXPSTATE == Start) not detected within timeout')

        # Then wait for it to be ready again
        if self.expstate.waitFor('== Ready', timeout=round(float(self.cfg.exptime) + 90.)):
            print('Exposure completed successfully')
        else:
            raise ValueError('Exposure EXPSTATE=Ready not detected within timeout')
        

class FocusSequence:
    """
    Perform a focus sequence.

    This is the base class for other derived classes that specify how the focus
    sequence proceeds.
    """
    def __init__(self):
        self._focus = Focus()
        self._exposure = Exposure()
        self.method = None
        self.reset()

    def reset(self):
        self.observed_focus = []
        self.exposures = []
        self.img_quality = []
        self.step_iter = 0
        self.source_stamps = []

    def execute(self, verbose=True, goto=True, method='brightest', **exp_kwargs):

        self._exposure.configure(**exp_kwargs)
        self.reset()

        while self.continue_sequence():
            self.observed_focus += [self.step_focus()]
            self.exposures += [self.take_exposure()]
            data, bkg, src_data, img_quality, source_stamp \
                = image_quality(self.exposures[-1], method=method)
#            self.show(data-bkg, source_stamp, src_data, self._focus.current, img_quality)
            self.source_stamps += [source_stamp]
            self.img_quality += [img_quality]
            self.step_iter += 1

        best_focus, best_img_quality = self.fit_best_focus(self.observed_focus, self.img_quality)
        if goto:
            self._focus.set_to(best_focus)
            self.take_exposure()
            sigma_at_best_focus = self.measure_fwhm(self.exposures[-1])

        return best_focus, best_img_quality

    def continue_sequence(self):
        raise NotImplementedError(
            'Method used to check if the sequence should continue is not implemented by the '
            'base class.'
        )

    def step_focus(self):
        raise NotImplementedError(
            'Method used to increment the focus setting is not implemented by the base class!'
        )

    def take_exposure(self):
        self._exposure.expose()
        return self._exposure.path.previous

    @staticmethod
    def fit_best_focus(focus, img_quality):
        if len(focus) < 3:
            raise ValueError('Insufficient number of focus values observed.')
        a, b, c = quadratic.fit_quadratic(focus, img_quality)
        return quadratic.vertex(a, b, c)


class GridFocusSequence(FocusSequence):
    def __init__(self, start, step, end=None, nstep=None):

        super().__init__()

        if end is None and nstep is None:
            raise ValueError('Must provide either the ending value or the number of steps.')
        self.nstep = nstep if end is None else int((end-start)/step+1)
        _end = start + (self.nstep-1)*step
        self.target_focus = np.round(np.linspace(start, _end, self.nstep))
        embed()
        exit()

    def continue_sequence(self):
        return self.step_iter < self.nstep
    
    def step_focus(self):
        self._focus.set_to(self.focus_values[self.step_iter])
        return self._focus.current


class AutomatedFocusSequence(FocusSequence):
    def __init__(self, start, step, maxsteps=12):

        if maxsteps is None or maxsteps < 2:
            raise ValueError('maxsteps cannot be None and must be at least 2.')

        super().__init__()

        self.start = start
        self.step = step
        self.direction = None
        self.last = None
        self.maxsteps = maxsteps

    def reset(self):
        super().reset()
        self.direction = None
        self.last = None

    def continue_sequence(self):
        return (
            self.step_iter < maxsteps
            and (self.step_iter < 2 
                or (self.last is not None and self.step_iter < self.last)
            )
        )
    
    def step_focus(self):
        if self.step_iter == 0:
            next_focus = self.start
        elif self.step_iter == 1:
            next_focus = self.start + self.step
        elif self.step_iter == 2 and self.img_quality[0] > self.img_quality[1]:
            self.direction = 1
            next_focus = self.focus_values[1] + self.step
        elif self.step_iter == 2 and self.img_quality[0] < self.img_quality[1]:
            self.direction = -1
            next_focus = self.focus_values[0] - self.step
        elif self.step_iter > 2 and self.img_quality[-1] > self.img_quality[-2]:
            self.last = self.step_iter + 2
            if self.last > self.maxsteps:
                warnings.warn(
                    f'Number of steps to fulfill sequence ({self.last}) is more than the '
                    f'maximum number of steps requested ({self.maxsteps}).')
            next_focus = self.focus_values[-1] + self.direction * self.step
        else:
            next_focus = self.focus_values[-1] + self.direction * self.step
        self._focus.set_to(next_focus)
        return self._focus.current


class ArchiveFocusSequence(FocusSequence):
    def __init__(self, observed_focus, exposures):
        super().__init__()

        self._observed_focus = observed_focus
        self._exposures = exposures

    def continue_sequence(self):
        return self.step_iter < self.nstep
    
    def step_focus(self):
        return self._observed_focus[self.step_iter]

    def take_exposure(self):
        return self._exposures[self.step_iter]


def detect_outliers(curve, plot, threshold=2.0):
    
    centroid_x = np.array([img.focus_star['Centroid'][0] for img in curve], dtype=float)
    centroid_y = np.array([img.focus_star['Centroid'][1] for img in curve], dtype=float)

    median_x = np.median(centroid_x)
    median_y = np.median(centroid_y)

    distances = np.sqrt((centroid_x - median_x)**2 + (centroid_y - median_y)**2)
    mean_distance = np.mean(distances)
    std_distance = np.std(distances)
    outlier_threshold = mean_distance + threshold * std_distance

    outliers = []
    for i, img in enumerate(curve):
        if distances[i] > outlier_threshold:
            outliers.append(img)

            cutout_size = int(3 * img.focus_star['FWHM'])

            half_size = cutout_size // 2
    
            # Calculate cutout boundaries
            x_center = int(round(img.focus_star['Centroid'][0]))
            y_center = int(round(img.focus_star['Centroid'][1]))
            # WHEN MOVING TO KTLPRACTICE, MAKE SURE TO USE CORRECT DICT KEY

            x_min = max(0, x_center - half_size)
            x_max = min(1024, x_center + half_size + 1)
            y_min = max(0, y_center - half_size)
            y_max = min(1024, y_center + half_size + 1)
            # USING HARDCODED 1024 FOR NOW
    
            # highlight cutout region
            rect = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, linewidth=2, edgecolor='yellow', facecolor='none')
            plot.ax_left.add_patch(rect)

    plot.fig.suptitle(f'Median Centroid ({median_x:.2f}×{median_y:.2f})', fontsize=16, fontweight='bold')

    return outliers, float(median_x), float(median_y)


def refit_curve(omit):
    try:
        data = Table.read("focus_data.ecsv")
    except Exception as e:
        raise Exception(f"Error reading focus data: {e}")

    images = []
    for line in data:
        images.append({
            'ObsNum': line['ObsNum'],
            'Focus': line['Focus'],
            'FWHM': line['FWHM'],
            'Centroid': (line['CentroidX'], line['CentroidY'])
        })

    if omit:
        curve = [img for img in images if int(img['ObsNum']) not in omit]
        print(f'Refitting curve with the following observations omitted: {omit}')
        if len(curve) <= 2:
            print('Not enough observations left to refit the curve. Exiting.')
            return None
    else:
        print('Refitting curve with all observations included.')

    x_values = [float(img['Focus']) for img in curve]
    y_values = [float(img['FWHM']) for img in curve]

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)

    print(f"Refitted curve: Optimal focus: {x_vertex}, FWHM: {y_vertex}")

    plt.figure(figsize=(12, 12))


    x_min, x_max = min(x_values), max(x_values)
    x_smooth = np.linspace(x_min, x_max, 50)
    y_smooth = a * x_smooth**2 + b * x_smooth + c

    plt.scatter(x_values, y_values, label='Measured FWHM', color='blue')
    plt.plot(x_smooth, y_smooth, 'r-', label='Fitted Quadratic')
    plt.scatter([x_vertex], [y_vertex], color='green', label='Optimal focus', zorder=3)
    plt.axvline(x=x_vertex, color='green', linestyle='--')

    plt.xlabel('Focus Value', fontsize=12)
    plt.ylabel('FWHM (pixels)', fontsize=12)
    plt.title('Focus Curve Analysis', fontsize=14, fontweight='bold')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)

    return curve

class TempFocusObject:
    def __init__(self, focus_star_dict):
        self.focus_star = focus_star_dict

def reevaluate(focus_coords, keyword, verbose):

    try:
        data = Table.read("focus_data.ecsv")
    except Exception as e:
        raise Exception(f"Error reading focus data: {e}")

    images = []
    for line in data:
        images.append({
            'ObsNum': line['ObsNum'],
            'Focus': line['Focus'],
            'FWHM': line['FWHM'],
            'Centroid': (line['CentroidX'], line['CentroidY'])
        })

    grid = Grid()
    plt.ion()
    plt.show(block=False)
    plt.pause(0.1)

    curve = set()

    for obs_index, obs in enumerate(images):
        obs_num = obs['ObsNum']
        filename = f"{keyword.dir_key.read()}/nickel/{keyword.file_key.read()}{obs_num}.{keyword.suffix_key.read()}"
        focus_value = obs['Focus']
        grid.index = obs_index

        focus_star = photometry(filename, obs_num, focus_value, grid, focus_coords, verbose=verbose)

        if focus_star is None:
            print(f"Image {filename} has no FWHM. Skipping.")
            continue

        temp_obj = TempFocusObject(focus_star)
        curve.add(temp_obj)

    curve = sorted(curve, key=lambda img: img.focus_star['Focus'])

    outliers, median_x, median_y = detect_outliers(curve, grid, threshold=2.0)
    if outliers:
        print(f"\n Median centroid: {centroid_median_x:.2f}, {centroid_median_y:.2f}")
        print(f"Outlier observations:")
        for outlier in outliers:
            x, y = outlier.focus_star['Centroid']
            print(f"   d{outlier.focus_star['ObsNum']}: ({x:.2f}, {y:.2f}) - Focus: {outlier.focus_star['Focus']}, FWHM: {outlier.focus_star['FWHM']:.3f}")

    print(curve)
    x_values = []
    y_values = []
    obs_values = []
    print("Reevaluated curve with the following focus values:")
    for img in curve:
        print(f"OBS: {img.focus_star['ObsNum']} Focus: {img.focus_star['Focus']}, FWHM: {img.focus_star['FWHM']}")
        x_values.append(img.focus_star['Focus'])
        y_values.append(img.focus_star['FWHM'])
        obs_values.append(img.focus_star['ObsNum'])

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)
    print(f"\nOptimal focus: {x_vertex}, FWHM: {y_vertex}")

    grid.set_right_axis(x_values, y_values, a, b, c, x_vertex, y_vertex, obs_values, outliers)

    return curve

import argparse

def main():

#    log_filename = f'focus_finding.log'
#    logger = setup_logging(log_level='INFO', log_file=log_filename)
#    logger.info("Starting focus finding process")

    parser = argparse.ArgumentParser(description="Automate the focus finding process.")

    parser.add_argument('focus', nargs='+', type=float,
        help='Focus starting value and step size.  You can also provide the last focus value.  '
             'Use --n to set the number of steps instead of the last focus value.  If neither '
             'are provided (last focus or number of focus steps), the code performs an automated '
             'number of steps to fill out the focus curve.'
    )
    parser.add_argument('-n', '--nstep', default=None, type=int,
        help='The number of focus steps to perform.  Ignored if the ending focus value is '
             'provided (see --focus).'
    )
    parser.add_argument('--method', type=str, nargs='+', default='brightest',
        help='The method used for calculating the image quality.  Must be "brightest" to use the '
             'brightest detected source in the field, "weighted" to use a weighted mean of all '
             'detected sources, or provide two pixel coordinates (x,y or column,row) to use a '
             'the detected source closest to the provided coordinates.'
    )
    parser.add_argument('--maxsteps', type=int, default=12,
        help='If using the automated focus sequence, this is the maximum of steps that are '
             'allowed.'
    )
    parser.add_argument('--obsnum', type=int, default=None,
        help='Re-analyze a focus sequence starting with the provided focus values and this '
             'observation number.  The number of available images must match the focus sequence '
             'requested.'
    )
    parser.add_argument('-t', '--exptime', default=5, type=float,
                        help='Exposure time in seconds for each exposure.')
    parser.add_argument('-b', '--binning', default=None, choices=['1,1', '2,2', '4,4'],
                        help='BinningExposure speed.  Must be "1,1", "2,2", or "4.4".')
    parser.add_argument('-s', '--speed', default=None, choices=['Slow', 'Fast'],
                        help='Exposure speed.  Must be Slow or Fast.  If None')
    parser.add_argument('-o', '--ofile', default=None, type=str,
        help='Output file for the measured focus data.  This can be used to exclude and refit the '
             'best focus.'
    )
    parser.add_argument('--refit', action='store_true',
        help='Refit the focus curve.  The output file (see --ofile) must be provided.'
    )
    parser.add_argument('--omit', type=int, nargs='+', default=None,
        help='List of observation numbers to omit from the curve fitting.'
    )
    parser.add_argument('--verbose', action='store_true',
        help='Enable verbose output for debugging'
    )
    args = parser.parse_args()

    # Set the read speed
    _speed = args.speed
    if _speed == 'Fast':
        warnings.warn('Fast does not work!  Setting to slow.')
        _speed == 'Slow'
    if _speed == 'Slow':
        _speed = '0.05MHz'
    elif _speed == 'Fast':
        _speed = '1.0MHz'

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
        seq = GridFocusSequence(args.focus[0], args.focus[1], end=end, nstep=args.nstep)

        if args.obsnum is not None:
            # Use GridFocusSequence to set the expected focus values
            expected_files = np.array([
                seq._exposure.path.from_obsnum(args.obsnum + i, assumed_recorded=True)
                for i in range(seq.nsteps)
            ])
            indx = np.array([Path(f).absolute().is_file() for f in expected_files])
            if not np.all(indx):
                raise FileNotFoundError('Expected to find the following files, but they are not '
                                        f'available: {", ".join(expected_files[indx].tolist())}')
            seq = ArchiveFocusSequence(seq.target_focus, expected_files)
    else:
        seq = AutomatedFocusSequence(args.focus[0], args.focus[1], maxsteps=args.maxsteps)

    embed()
    exit()

    best_focus, best_img_quality = seq.execute(goto=False, method=args.method, record=True,
                                               speed=_speed, exptime=args.exptime,
                                               binning=args.binning)
    
    # - Print the best focus and the img quality
    # - Write the output file if provided

            


#    logger.info(f"Arguments: {vars(args)}")

    keyword = Keyword('Yes', args.length_exposure, args.exposure_speed)
    # Record: 0 No 1 Yes, Exposure: seconds, Speed: 0 Slow 1 Medium 2 Fast

    event = Event(keyword, args.verbose)
    keyword.verbose = args.verbose

    filepath = f"{keyword.dir_key.read()}/{keyword.file_key.read()}{keyword.obs_key.read()}.{keyword.suffix_key.read()}"

    if args.reevaluate is True:

        curve = reevaluate(args.focus_coords, keyword, args.verbose)

        data = Table()
        data['ObsNum'] = [np.array(img.focus_star['ObsNum'], dtype=int) for img in curve]
        data['Focus'] = [np.array(img.focus_star['Focus'], dtype=float) for img in curve]
        data['FWHM'] = [np.array(img.focus_star['FWHM'], dtype=float) for img in curve]
        data['CentroidX'] = [np.array(img.focus_star['Centroid'][0], dtype=float) for img in curve]
        data['CentroidY'] = [np.array(img.focus_star['Centroid'][1], dtype=float) for img in curve]
        data.description = 'Focus curve data with observations, focus values, FWHM, and centroids'
        data.meta['CentroidX'] = f'median: {np.median(data["CentroidX"]):.2f}, std: {np.std(data["CentroidX"]):.2f}'
        data.meta['CentroidY'] = f'median: {np.median(data["CentroidY"]):.2f}, std: {np.std(data["CentroidY"]):.2f}'
        data.write("focus_data.ecsv", overwrite=True)

    elif args.refit is True:
        curve = refit_curve(args.omit)

    if keyword.event_key.read() != 'ControllerReady' and keyword.event_key.read() != 'ExposeSequenceDone':
        raise Exception("Controller not ready. Cannot take exposure.")

    #check if pocstop is enabled
    if keyword.stop_key.read() == 1:
        print("POCSTOP is 'disabled'. Waiting for 'enabled' to allow motion")
    if not keyword.stop_key.waitFor('== allowed', timeout=30):
        raise Exception("POCSTOP is 'disabled'. Set to 'enabled' to allow motion")
    
    if args.refit is False and args.reevaluate is False:

        if args.focus_end:
            curve = manual_focus_finder(args.focus_start, args.focus_end, args.step_size, args.focus_coords, keyword)
        else:
            curve = auto_focus_finder(args.focus_start, args.step_size, args.focus_coords, keyword)

        data = Table()
        data['ObsNum'] = [np.array(img.focus_star['ObsNum'], dtype=int) for img in curve]
        data['Focus'] = [np.array(img.focus_star['Focus'], dtype=float) for img in curve]
        data['FWHM'] = [np.array(img.focus_star['FWHM'], dtype=float) for img in curve]
        data['CentroidX'] = [np.array(img.focus_star['Centroid'][0], dtype=float) for img in curve]
        data['CentroidY'] = [np.array(img.focus_star['Centroid'][1], dtype=float) for img in curve]
        data.description = 'Focus curve data with observations, focus values, FWHM, and centroids'
        data.meta['CentroidX'] = f'median: {np.median(data["CentroidX"]):.2f}, std: {np.std(data["CentroidX"]):.2f}'
        data.meta['CentroidY'] = f'median: {np.median(data["CentroidY"]):.2f}, std: {np.std(data["CentroidY"]):.2f}'
        data.write("focus_data.ecsv", overwrite=True)
   

    plt.ioff()
    plt.show(block=True)


if __name__ == "__main__":
    main()


# /data/nickel
# 8:30

