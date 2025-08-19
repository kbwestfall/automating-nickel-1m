#! @KPYTHON@

from IPython import embed

from pathlib import Path
import logging
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import quadratic

from astropy.io import fits
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
    def __init__(self, scratch=False):
        self.expresult = ktl.cache('nscicam', 'EXPRESULT')
        self.recorddir = ktl.cache('nscicam', 'SCRATCHDIR' if scratch else 'RECORDDIR')
        self.prefix = ktl.cache('nscicam', 'FITSPREFIX')
        self.obsnum = ktl.cache('nscicam', 'OBSNUM')
        self.suffix = ktl.cache('nscicam', 'FITSSUFFIX')

    @property
    def previous(self):
        return self.expresult.read()

    @property
    def next(self):
        return self.for_obsnum(self.obsnum.read())

    def for_obsnum(self, obsnum):
        path = Path(self.recorddir.read()).absolute()
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
    
        self._exppath = ExposurePath()
        self._expcfg = ExposureConfig()

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
        print(f'EXPSTATE updated: {self.expstate_value}')

    def _filepath_callback(self, keyword):
        self.filepath = self._exppath.next()
        print(f'FILEPATH updated: {self.filepath}')

    def expose(self, record=None, speed=None, binning=None, exptime=None):

        self._expcfg.configure(record=record, speed=speed, binning=binning, exptime=exptime)

        # Check that an exposure isn't currently happening
        if not self.expstate.waitFor('== Ready', timeout=15):
            error_msg = "Camera exposure state not ready. Cannot take exposure."
            self.logger.error(error_msg)
            raise ValueError(error_msg)
        
        self.expstart.write('StartX')

        if not self.expstate.waitFor('== Start', timeout=30):
            raise ValueError('Exposure start (EXPSTATE == Start) not detected within timeout')

        if self.expstate.waitFor('== Ready', timeout=round(self._expcfg.exptime + 90.)):
            print('Exposure completed successfully')
        else:
            raise ValueError('Exposure EXPSTATE=Ready not detected within timeout')
        

#class FocusData:
    #def __init__(self, image):



class FocusSequence:
    """
    Perform a focus sequence.

    Parameters
    ----------
    focus_start : :obj:`int`, :obj`float`
        Starting focus value; e.g., 330
    focus_step : :obj:`int`, :obj:`float`
        Change in focus between different images; e.g., 5
    focus_end : :obj:`int`, :obj:`float`, optional
        The focus at which to stop, inclusive.  If None and ``nstep`` is None,
        use the automated procedure to find the best focus.  If both are
        provided, ``focus_end`` takes precedence.
    nstep : :obj:`int`, optional
        The number of focus steps to take.  If None and ``focus_end`` is None, 
        use the automated procedure to find the best focus.  If both are
        provided, ``focus_end`` takes precedence.
    obsnum : :obj:`int`, optional
        Observation number to start with for a previously observed sequence.  If
        None, new images are taken.
    verbose : :obj:`bool`, optional
        Verbose output the screen.
    """
    def __init__(self):

        self.logger = logging.getLogger(__name__)
        self.verbose = verbose

        self.data = None

        # Object used to change the focus
        self._focus = Focus()
        self._exposure = Exposure()

    def execute(self, speed='0.05MHz', binning='2,2', exptime=5, verbose=True):

        self._exposure.configure(record=True, speed=speed, binning=binning, exptime=exptime)

        while not self.done():
            self.step_focus()
            self._exposure.expose()
            self.measure_fwhm()

        self.fit_best_focus()

    def measure_fwhm(self):

        self.change_focus(focus_value)
        self.focus_value = focus_value
        self.exposure()

        filepath = self.keyword.filepath
        self.logger.debug(f"Exposure saved at: {filepath}")
        print(f"Exposure being saved at: {filepath}")

        obs_num = self.keyword.obs_key.read()

        hdu = fits.open(filepath)

        self.focus_star = photometry(filepath, obs_num, focus_value, grid, focus_coords, verbose=self.verbose)
        self.fwhm = self.focus_star['FWHM']
        self.logger.debug(f"Photometry complete - FWHM: {self.fwhm}")
        print(f" {filepath} FWHM: {self.fwhm} \n")



def curve_finder(image1, image2, seen, focus_coords, keyword, grid, direction=None):
    seen.add(image1)
    seen.add(image2)

    if image1.fwhm is None or image2.fwhm is None:
        raise Exception("A source could not be detected in one of the images. Cannot continue with focus finding.")
    
    if image1.fwhm > image2.fwhm:
        if direction is None:
            direction = 'right'
        if direction == 'left':
            return curve_helper(image1, image2, seen, focus_coords, keyword, grid)
        focus3 = image2.focus_value + (image2.focus_value - image1.focus_value)
        image3 = Event(keyword)
        image3.sequence(focus3, focus_coords, grid)
        grid.index += 1
        return curve_finder(image2, image3, seen, focus_coords, keyword, grid, direction)
    elif image1.fwhm < image2.fwhm:
        if direction is None:
            direction = 'left'
        if direction == 'right':
            return curve_helper(image1, image2, seen, focus_coords, keyword, grid)
        focus3 = image1.focus_value - (image2.focus_value - image1.focus_value)
        image3 = Event(keyword)
        image3.sequence(focus3, focus_coords, grid)
        grid.index += 1
        return curve_finder(image3, image1, seen, focus_coords, keyword, grid, direction)

def curve_helper(image1, image2, seen, focus_coords, keyword, grid, iterations=2):
    if iterations > 0:
        if image1.fwhm is None or image2.fwhm is None:
            raise Exception("A source could not be detected in one of the images. Cannot continue with focus finding.")

        if image1.fwhm > image2.fwhm:
            focus3 = image1.focus_value - (image2.focus_value - image1.focus_value)
            image3 = Event(keyword)
            image3.sequence(focus3, focus_coords, grid)
            grid.index += 1
            seen.add(image3)
            return curve_helper(image3, image1, seen, focus_coords, keyword, grid, iterations-1)
        elif image1.fwhm < image2.fwhm:
            focus3 = image2.focus_value + (image2.focus_value - image1.focus_value)
            image3 = Event(keyword)
            image3.sequence(focus3, focus_coords, grid)
            grid.index += 1
            seen.add(image3)
            return curve_helper(image2, image3, seen, focus_coords, keyword, grid, iterations-1)
    else:
        return seen

def auto_focus_finder(initial_focus, step_size, focus_coords, keyword):

    grid = Grid()
    plt.ion()
    plt.show(block=False)
    plt.pause(0.1)

    seen = set()
    
    image1 = Event(keyword)
    image1.sequence(initial_focus, focus_coords, grid)
    grid.index += 1

    image2 = Event(keyword)
    image2.sequence(initial_focus + step_size, focus_coords, grid)
    grid.index += 1

    curve = curve_finder(image1, image2, seen, focus_coords, keyword, grid)
    if len(curve) < 3:
        raise Exception("Not enough images. Cannot find optimal focus.")
    curve = sorted(curve, key=lambda img: img.focus_value)

    outliers, median_x, median_y = detect_outliers(curve, grid, threshold=2.0)
    if outliers:
        print(f"\n Median centroid: {centroid_median_x:.2f}, {centroid_median_y:.2f}")
        print(f"Outlier observations:")
        for outlier in outliers:
            x, y = outlier.focus_star['Centroid']
            print(f"   d{outlier.focus_star['ObsNum']}: ({x:.2f}, {y:.2f}) - Focus: {outlier.focus_star['Focus']}, FWHM: {outlier.focus_star['FWHM']:.3f}")

    x_values = []
    y_values = []
    obs_values = []
    print("Curve found with the following focus values:")
    for img in curve:
        print(f"OBS: {img.focus_star['ObsNum']} Focus: {img.focus_value}, FWHM: {img.fwhm}")
        x_values.append(img.focus_value)
        y_values.append(img.fwhm)
        obs_values.append(img.focus_star['ObsNum'])

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)
    print(f"\nOptimal focus: {x_vertex}, FWHM: {y_vertex}")

    grid.set_right_axis(x_values, y_values, a, b, c, x_vertex, y_vertex, obs_values, outliers)

    return curve

def manual_focus_finder(initial_focus, end_focus, step_size, focus_coords, keyword):
    if (end_focus - initial_focus) // step_size < 2:
        raise Exception("Not enough images to find a focus curve. Increase the range or decrease the step size.")

    grid = Grid()
    plt.ion()
    plt.show(block=False)
    plt.pause(0.1)

    curve = set()
    
    focus = initial_focus
    while focus <= end_focus:
        image = Event(keyword, keyword.verbose)
        image.sequence(focus, focus_coords, grid)
        grid.index += 1
        focus += step_size
        if image.fwhm is None:
            print(f"Image at focus {focus} has no FWHM. Skipping.")
            continue
        curve.add(image)
        plt.draw()

    outliers, median_x, median_y = detect_outliers(curve, grid, threshold=2.0)
    if outliers:
        print(f"\n Median centroid: {median_x:.2f}, {median_y:.2f}")
        print(f"Outlier observations:")
        for outlier in outliers:
            x, y = outlier.focus_star['Centroid']
            print(f"   d{outlier.focus_star['ObsNum']}: ({x:.2f}, {y:.2f}) - Focus: {outlier.focus_star['Focus']}, FWHM: {outlier.focus_star['FWHM']:.3f}")
    

    x_values = []
    y_values = []
    obs_values = []
    print("Curve found with the following focus values:")
    for img in curve:
        print(f"OBS: {img.focus_star['ObsNum']} Focus: {img.focus_value}, FWHM: {img.fwhm}")
        x_values.append(img.focus_value)
        y_values.append(img.fwhm)
        obs_values.append(img.focus_star['ObsNum'])

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)
    print(f"\nOptimal focus: {x_vertex}, FWHM: {y_vertex}")

    grid.set_right_axis(x_values, y_values, a, b, c, x_vertex, y_vertex, obs_values, outliers)

    return curve

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

    focus = Focus()

    embed()
    exit()

    fseq = FocusSequence(340, 5)
    #fseq.set_focus(360)
    #fseq.set_focus(351)
    fseq.take_exposure(record=False)

    embed()
    exit()

    # timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f'focus_finding.log'
    logger = setup_logging(log_level='INFO', log_file=log_filename)
    logger.info("Starting focus finding process")

    parser = argparse.ArgumentParser(description="Automate the focus finding process.")
    parser.add_argument('-fs', '--focus_start', default=350, type=int, help='Start Focus value')
    parser.add_argument('-fe', '--focus_end', default=None, type=int, help='End Focus value')
    parser.add_argument('-s', '--step_size', default=5, type=int, help='Step size for focus increments')
    parser.add_argument('-el', '--length_exposure', default=1, type=float, help='Exposure length in seconds')
    parser.add_argument('-es', '--exposure_speed', default='Fast', choices=['Slow', 'Medium', 'Fast'], help='Exposure speed: Slow, Medium, Fast')
    parser.add_argument('--focus_coords', type=float, nargs='*', default=None, help='Focus star coordinates (x, y) for manual focus finding')
    parser.add_argument('--refit', action='store_true', help='Refit the focus curve with omitted outliers')
    parser.add_argument('--omit', type=int, nargs='*', default=None, help='List of observation numbers to omit from the curve fitting')
    parser.add_argument('--reevaluate', action='store_true', help='Reevaluate the focus curve with the last recorded focus data')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output for debugging')
    args = parser.parse_args()

    logger.info(f"Arguments: {vars(args)}")

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

