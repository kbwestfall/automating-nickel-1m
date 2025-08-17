#! @KPYTHON@

from IPython import embed

import logging
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import quadratic

from astropy.io import fits
from astropy.table import Table

from photometry import photometry, Grid

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


class Keyword:
    def __init__(self, record, exposure, speed):

        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Initializing Keyword with record={record}, exposure={exposure}, speed={speed}")

        self.verbose = False

        self.record = record
        self.exposure = exposure
        self.speed = speed
        self.start = 'Yes'

        self.secpa_key = ktl.cache('nickelpoco', 'POCSECPA')
        self.secpd_key = ktl.cache('nickelpoco', 'POCSECPD')
        self.seclk_key = ktl.cache('nickelpoco', 'POCSECLK')
        self.obs_key = ktl.cache('nickucam', 'OBSNUM')
        self.dir_key = ktl.cache('nickucam', 'OUTDIR')
        self.file_key = ktl.cache('nickucam', 'OUTFILE')
        self.suffix_key = ktl.cache('nickucam', 'RECORD_SUFFIX')
        self.record_key = ktl.cache('nickucam', 'RECORD')
        self.exposure_key = ktl.cache('nickucam', 'EXPOSURE')
        self.speed_key = ktl.cache('nickucam', 'READSPEED')
        self.start_key = ktl.cache('nickucam', 'START')

        self.stop_key = ktl.cache('nickelpoco', 'POCSTOP')

        self.event_key = ktl.cache('nickucam', 'EVENT')
        self.event_key.callback(self.event_callback)
        self.event_key.monitor()

        self.start_key.callback(self.filepath_callback)
        self.start_key.monitor()

    def event_callback(self, keyword):
        self.event_value = keyword.read()
        self.logger.debug(f'EVENT updated: {self.event_value}')

    def filepath_callback(self, keyword):
        self.filepath = f"{self.dir_key.read()}/nickel/{self.file_key.read()}{self.obs_key.read()}.{self.suffix_key.read()}"
        self.logger.debug(f'FILEPATH updated: {self.filepath}')
        # print(f'update FILEPATH: {self.filepath}')

#    def general_callback(self, keyword):
#        self.general = keyword.read()
    
    # def wait_for(self, keyword, expected_value, timeout=15):
    #     keyword.callback(self.general_callback)
    #     keyword.monitor(start=True, prime=True)
    #     start_time = time.time()
    #     while time.time() - start_time < timeout:
    #         if self.general in expected_value:
    #             return True
    #     return False

    # def wait_until(self, keyword, expected_value, timeout=15):
    #     self.logger.debug(f"Waiting for {keyword} to be in {expected_value} (timeout: {timeout}s)")
    #     start_time = time.time()
    #     while time.time() - start_time < timeout:
    #         if keyword.read() in expected_value:
    #             return True
    #     return False

        

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
    def __init__(self, focus_start, focus_step, focus_end=None, nstep=None, obsnum=None,
                 verbose=True):
        self.start = focus_start
        self.step = focus_step
        self.logger = logging.getLogger(__name__)
#        self.keyword = keyword
        self.verbose = verbose

        # Keys
        self._secpa = None
        self._secpd = None
        self._seclk = None

        self._expstate = None

    @property
    def secpa(self):
        if self._secpa is None:
            self._secpa = ktl.cache('nickelpoco', 'POCSECPA')
        return self._secpa

    @property
    def secpd(self):
        if self._secpd is None:
            self._secpd = ktl.cache('nickelpoco', 'POCSECPD')
        return self._secpd

    @property
    def seclk(self):
        if self._seclk is None:
            self._seclk = ktl.cache('nickelpoco', 'POCSECLK')
        return self._seclk
    
    @property
    def expstate(self):
        if self._expstate is None:
            self._expstate = ktl.cache('nscicam', 'EXPSTATE')
            self._expstate.callback(self.expstate_callback)
            self._expstate.monitor()
        return self._expstate

    def expstate_callback(self, keyword):
        self.expstate_value = keyword.read()
        self.logger.debug(f'EXPSTATE updated: {self.expstate_value}')
    
    def set_focus(self, focus_value):

        self.logger.info(f"Changing focus to {focus_value}")
        
        self.logger.debug(f'Current POCSECPD: {self.secpd.read()}')
        self.logger.debug(f'Current POCSECPA: {self.secpa.read()}')

        if abs(float(self.secpd.read()) - focus_value) < .1:
            print(f'POCSECPD already set to {focus_value}. No change needed.')
            return
        
        if focus_value < 165 or focus_value > 500:
            error_msg = f"Focus value {focus_value} is out of range (165-500)."
            self.logger.error(error_msg)
            raise ValueError(error_msg)
        
        if not self.expstate.waitFor('== Ready', timeout=15):
            error_msg = "Camera exposure state not ready. Cannot change focus."
            self.logger.error(error_msg)
            raise Exception(error_msg)

        self.logger.debug(f'Current POCSECLK: {self.seclk.read()}')
        self.seclk.write('off')
        self.seclk.read()
        self.logger.debug("Set POCSECLK to 'off'")

        self.seclk.write('on')
        self.seclk.read()
        self.logger.debug("Set POCSECLK to 'on'")
        embed()
        exit()

        print(f'POCSECPA: {self.keyword.secpa_key.read()}')
        self.keyword.secpd_key.write(focus_value)
        self.logger.debug(f'Set POCSECPD to {focus_value}')
        print(f'POCSECPD: {self.keyword.secpd_key.read()}')

        if not self.keyword.seclk_key.waitFor('== on', timeout=30):
            error_msg = "POCSECLK did not turn on. Focus change failed."
            self.logger.error(error_msg)
            raise Exception(error_msg)
            
        self.logger.info(f"Successfully changed focus to {focus_value}")
    ### CHANGE FOCUS ###

    ### TAKE EXPOSURE ###
    def exposure(self):
        self.logger.info(f"Starting exposure: {self.keyword.exposure}s at {self.keyword.speed} speed")

        if not self.keyword.event_key.waitFor('== ControllerReady', timeout=15):
            error_msg = "Controller not ready. Cannot take exposure."
            self.logger.error(error_msg)
            raise Exception(error_msg)

        self.keyword.record_key.write(self.keyword.record)
        self.logger.debug(f'Set RECORD: {self.keyword.record_key.read()}')

        self.keyword.exposure_key.write(self.keyword.exposure)
        self.logger.debug(f'Set EXPOSURE: {self.keyword.exposure_key.read()}')

        self.keyword.speed_key.write(self.keyword.speed)
        self.logger.debug(f'Set SPEED: {self.keyword.speed_key.read()}')

        self.keyword.start_key.write(self.keyword.start)
        self.logger.debug(f'SET START: {self.keyword.start_key.read()}')

        if not self.keyword.event_key.waitFor('== ExposureBegin', timeout=30):
            error_msg = "ExposureBegin not detected within timeout"
            self.logger.error(error_msg)
            raise Exception(error_msg)  

        self.logger.debug("Waiting for ExposureEnd...")
        if not self.keyword.event_key.waitFor('== ExposureEnd', timeout=round(self.keyword.exposure * 1.2)):
            error_msg = "ExposureEnd not detected within timeout"
            self.logger.error(error_msg)
            raise Exception(error_msg)

        if self.keyword.event_key.waitFor('== ControllerReady', timeout=30):
            self.logger.debug("Exposure completed successfully")
        else:
            error_msg = "ControllerReady not detected within timeout"
            self.logger.warning(error_msg)
            raise Exception(error_msg)
    ### TAKE EXPOSURE ###

    def sequence(self, focus_value, focus_coords, grid):

        # if self.keyword.record == 'Yes':
        #     self.keyword.dir_key.write("/data")
        #     self.keyword.file_key.write(f"{count}focus_")
        #     self.keyword.obs_key.write(str(focus_value))
        #     self.keyword.suffix_key.write("fits")


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

    fseq = FocusSequence(340, 5)
    fseq.set_focus(360)

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

