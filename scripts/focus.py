#! @KPYTHON@

import argparse
import warnings
from IPython import embed

from pathlib import Path
from dataclasses import dataclass
import logging

import numpy as np
from matplotlib import pyplot, patches

from astropy.table import Table
from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch

from photometry import image_quality
import quadratic

try:
    import ktl
except ModuleNotFoundError:
    warnings.warn('Unable to import ktl package.  Functionality will be limited.')
    ktl = None


#def setup_logging(log_level='INFO', log_file=None):
#    """
#    Setup logging configuration
#    
#    Parameters:
#    log_level: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
#    log_file: Optional filename to save logs to file
#    """
#    # Create formatter
#    formatter = logging.Formatter(
#        '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s'
#    )
#    
#    # Setup root logger
#    logger = logging.getLogger()
#    logger.setLevel(getattr(logging, log_level.upper()))
#    
#    # Clear any existing handlers
#    logger.handlers.clear()
#    
#    # Console handler
#    console_handler = logging.StreamHandler(sys.stdout)
#    console_handler.setLevel(getattr(logging, log_level.upper()))
#    console_handler.setFormatter(formatter)
#    logger.addHandler(console_handler)
#    
#    # File handler (optional)
#    if log_file:
#        file_handler = logging.FileHandler(log_file)
#        file_handler.setLevel(logging.DEBUG)  
#        file_handler.setFormatter(formatter)
#        logger.addHandler(file_handler)
#    
#    return logger


class Focus:
    """
    KTL object used to manipulate the telescope focus.

    Attributes
    ----------
    pocstop : ktl.Keyword

    """
    def __init__(self):
        if ktl is None:
            raise RuntimeError(
                'The ktl package is not available.  Controlling the telescope focus requires a '
                'kpython environment with ktl installed.'
            )
        # Controls overall movement of mechanisms
        self.pocstop = ktl.cache('nickelpoco', 'POCSTOP')
        # Actual secondary position
        self.secpa = ktl.cache('nickelpoco', 'POCSECPA')
        # Desired secondary position
        self.secpd = ktl.cache('nickelpoco', 'POCSECPD')
        # Controls the secondary locking mechanism
        self.seclk = ktl.cache('nickelpoco', 'POCSECLK')
        # Provides the state of the exposures
        self.expstate = ktl.cache('nscicam', 'EXPSTATE')

    @property
    def current(self):
        """
        The current focus position (floating-point)
        """
        return float(self.secpa.read())

    def set_to(self, focus_value):
        """
        Set the focus to the provided value.

        Movement must be enabled and there must not be an ongoing exposure.

        Parameters
        ----------
        focus_value : :obj:`int`, :obj:`float`
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
            raise ValueError(f'Requested focus ({focus_value}) is out of range (165-500).')
        
        # Make sure movement is enabled.  Do NOT enable movement via this
        # script!
        if not self.pocstop.waitFor('== allowed', timeout=0.5):
            raise ValueError('Telescope movement is disabled!')

        # Check that an exposure isn't currently happening
        if not self.expstate.waitFor('== Ready', timeout=0.5):
            raise ValueError('Camera exposure state not ready. Cannot change focus.')

        # Check if the focus is already at the requested position.
        if abs(float(self.secpd.read()) - focus_value) < .1:
            print(f'POCSECPD already set to {focus_value}. No change needed.')
            return

        # Change the focus       
        print('Unlocking secondary')
        self.seclk.write('off')
        self.seclk.read()

        print(f'Actual position: {self.secpa.read()}')
        self.secpd.write(focus_value)
        print(f'Desired position: {self.secpd.read()}')

        if not self.seclk.waitFor('== on', timeout=30):
            # TODO: Explicitly set the lock to on?
#            self.seclk.write('on')
#            self.seclk.read()
            raise ValueError("POCSECLK did not turn on. Focus change failed.")
            
        print(f"Successfully changed focus to {focus_value}")


class ExposurePath:
    """
    KTL object used to read and set the paths for the image exposures.
    """
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
        return Path(self.expresult.read())

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
        return path / f'{self.prefix.read()}{obsnum}{self.suffix.read()}'
    

class ExposureConfig:
    def __init__(self):
        self.exprec = ktl.cache('nscicam', 'RECORD')
        self.inttime = ktl.cache('nscicam', 'EXPOSURE')
        # NOTE: AMPCONFLIST gives the options
        self.expspeed  = ktl.cache('nscicam', 'AMPCONF')
        self.expspeed_options = ktl.cache('nscicam', 'AMPCONFLIST').read().split(',')
        # TODO: Is there a way to programmatically get the binning options
        self.expbin  = ktl.cache('nscicam', 'BINNING')
        self.expbin_options = ['1,1', '2,2', '4,4']

    def configure(self, record=None, speed=None, binning=None, exptime=None):
        if record is not None:
            self.exprec.write(record)
        if speed is not None:
            if speed not in self.expspeed_options:
                raise ValueError(
                    f'{speed} is not a valid read speed.  Must be one of {self.expspeed_options}.'
                )
            self.expspeed.write(speed)
        if binning is not None:
            if binning not in self.expbin_options:
                raise ValueError(
                    f'{binning} is not a valid binning.  Must be one of {self.expbin_options}.'
                )
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
        if ktl is None:
            raise RuntimeError(
                'The ktl package is not available.  Taking new exposures requires a kpython '
                'environment with ktl installed.'
            )
        self.path = ExposurePath()
        self.cfg = ExposureConfig()

        # SCICAM exposure keywords
        self.expstate = ktl.cache('nscicam', 'EXPSTATE')
        self.expstate.monitor()

        self.expstart = ktl.cache('nscicam', 'EXPOSE')
        self.expstart.monitor()

    def expose(self, record=None, speed=None, binning=None, exptime=None):

        # Check that an exposure isn't currently happening
        if not self.expstate.waitFor('== Ready', timeout=15):
            raise ValueError('Camera exposure state not ready. Cannot take exposure.')

        self.cfg.configure(record=record, speed=speed, binning=binning, exptime=exptime)

        # Start the exposure
        # TODO: This needs to be confirmed once Will settles on the naming
        self.expstart.write('StartX')

        # Wait for it to start
        if not self.expstate.waitFor('== Start', timeout=30):
            raise ValueError('Exposure start (EXPSTATE == Start) not detected within timeout')

        # Then wait for it to be ready again
        if self.expstate.waitFor('== Ready', timeout=round(float(self.cfg.exptime) + 90.)):
            print('Exposure completed successfully')
        else:
            raise ValueError('Exposure EXPSTATE=Ready not detected within timeout')
        

class FocusPlot:
    """
    Live view of a focus sequence in progress.

    Shows the most recent exposure (with the measured source boxed, and
    flagged if its centroid looks like an outlier relative to the rest of
    the sequence collected so far), a grid of the per-exposure source
    stamps, and an evolving plot of FWHM versus focus value.

    Parameters
    ----------
    nstamps : :obj:`int`
        Number of stamps to reserve space for (an upper bound on the number
        of exposures expected in the sequence).
    ncols : :obj:`int`, optional
        Number of columns in the stamp grid.
    """
    def __init__(self, nstamps, ncols=4):
        self.ncols = ncols
        self.nrows = int(np.ceil(nstamps / ncols))

        self.fig = pyplot.figure(figsize=(6 + 3*self.ncols, max(3*self.nrows, 6)),
                                  constrained_layout=True)
        subfigs = self.fig.subfigures(1, 3, width_ratios=[1.2, self.ncols, 1.2])

        self.frame_ax = subfigs[0].subplots(1, 1)
        self.frame_ax.set_aspect('equal')

        self.stamp_axes = subfigs[1].subplots(self.nrows, self.ncols, squeeze=False)
        for ax in self.stamp_axes.flat:
            ax.axis('off')

        self.curve_ax = subfigs[2].subplots(1, 1)
        self._reset_curve_axis()

        pyplot.ion()
        pyplot.show(block=False)
        pyplot.pause(0.1)

    def _reset_curve_axis(self):
        self.curve_ax.clear()
        self.curve_ax.set_xlabel('Focus Value')
        self.curve_ax.set_ylabel('FWHM (pixels)')
        self.curve_ax.set_title('Focus Curve')
        self.curve_ax.grid(True, alpha=0.3)

    def update_frame(self, data, coords, stamp_size, label, focus_value, is_outlier):
        """Show the most recent full frame, with the measured source boxed."""
        self.frame_ax.clear()
        norm = ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
        self.frame_ax.imshow(data, cmap='viridis', origin='lower', norm=norm)
        self.frame_ax.set_xticks([])
        self.frame_ax.set_yticks([])
        self.frame_ax.set_title(f'{label}  Focus: {focus_value:.1f}')

        half = stamp_size / 2
        color = 'yellow' if is_outlier else 'red'
        rect = patches.Rectangle((coords[0]-half, coords[1]-half), stamp_size, stamp_size,
                                  linewidth=2, edgecolor=color, facecolor='none')
        self.frame_ax.add_patch(rect)
        if is_outlier:
            self.frame_ax.text(coords[0], coords[1]+half+5, 'Outlier centroid', color='yellow',
                                fontsize=10, ha='center')
        self._draw()

    def add_stamp(self, index, stamp, focus_value, fwhm, is_outlier=False):
        """Show the source stamp for the index-th exposure taken so far."""
        row, col = divmod(index, self.ncols)
        if row >= self.nrows:
            warnings.warn(f'No plot space left for stamp {index}; skipping.')
            return
        ax = self.stamp_axes[row, col]
        ax.clear()
        ax.axis('off')
        norm = ImageNormalize(stamp, interval=ZScaleInterval(), stretch=LinearStretch())
        ax.imshow(stamp, origin='lower', cmap='viridis', norm=norm)
        title_color = 'darkorange' if is_outlier else 'black'
        ax.set_title(f'Focus {focus_value:.1f}, FWHM {fwhm:.2f}', fontsize=9, color=title_color)
        self._draw()

    def update_curve(self, focus_values, fwhm_values):
        """
        Redraw the focus-vs-FWHM curve, including the best-fit quadratic
        once there are enough points to fit one.
        """
        self._reset_curve_axis()
        self.curve_ax.scatter(focus_values, fwhm_values, color='tab:blue')

        if len(focus_values) >= 3:
            a, b, c = quadratic.fit_quadratic(focus_values, fwhm_values)
            x_vertex, y_vertex = quadratic.vertex(a, b, c)
            x_smooth = np.linspace(min(focus_values), max(focus_values), 50)
            y_smooth = a*x_smooth**2 + b*x_smooth + c
            self.curve_ax.plot(x_smooth, y_smooth, 'r-')
            self.curve_ax.scatter([x_vertex], [y_vertex], color='green', zorder=3,
                                   label=f'Best focus: {x_vertex:.1f}')
            self.curve_ax.legend(loc='best')

        self._draw()

    def _draw(self):
        self.fig.canvas.draw_idle()
        pyplot.pause(0.001)

    @staticmethod
    def is_outlier(centroids, threshold=2.0):
        """
        Check whether the most recently added centroid in ``centroids`` is
        an outlier relative to the full set collected so far, using the same
        median-distance statistic as the original focus-finding code.
        """
        _centroids = np.atleast_2d(centroids).astype(float)
        if len(_centroids) < 3:
            return False
        median = np.median(_centroids, axis=0)
        distances = np.sqrt(np.sum((_centroids - median)**2, axis=1))
        outlier_threshold = np.mean(distances) + threshold * np.std(distances)
        return distances[-1] > outlier_threshold


@dataclass
class StepResult:
    """
    The measurement produced by one step of a focus sequence.

    Attributes
    ----------
    index : :obj:`int`
        0-based position of this step within the sequence.
    focus_value : :obj:`float`
        The focus value used for this exposure.
    exposure : :class:`pathlib.Path`
        Path to the exposure's FITS file.
    fwhm : :obj:`float`
        Measured image quality (FWHM) of the selected source.
    frame : :class:`numpy.ndarray`
        Background-subtracted image data for the full exposure.
    stamp : :class:`numpy.ndarray`
        Cutout of ``frame`` around the selected source.
    centroid : :obj:`tuple`
        The (x,y) centroid of the selected source.
    is_outlier : :obj:`bool`
        Whether this step's centroid is an outlier relative to the rest of
        the sequence collected so far; see :func:`FocusPlot.is_outlier`.
    """
    index: int
    focus_value: float
    exposure: Path
    fwhm: float
    frame: np.ndarray
    stamp: np.ndarray
    centroid: tuple
    is_outlier: bool


class FocusSequence:
    """
    Perform a focus sequence.

    This is the base class for other derived classes that specify how the focus
    sequence proceeds.
    """
    def __init__(self):
        # `_focus` and `_exposure` control the telescope and camera hardware.
        # They require a live ktl connection, which is not needed by
        # sequences (e.g., ArchiveFocusSequence) that only reprocess
        # exposures that already exist on disk.
        self._focus = Focus() if ktl is not None else None
        self._exposure = Exposure() if ktl is not None else None
        self.method = None
        self.reset()

    def reset(self):
        self.observed_focus = []
        self.exposures = []
        self.img_quality = []
        self.step_iter = 0
        self.source_stamps = []
        self.centroids = []

    @property
    def expected_steps(self):
        raise NotImplementedError(
            'Subclasses must report the (maximum) number of steps in the sequence, used to size '
            'the live plot.'
        )

    def step(self, method='brightest'):
        """
        Advance the focus sequence one exposure at a time.

        This is the engine shared by :func:`execute` (used by the CLI
        script) and any other driver of the sequence (e.g., a GUI worker
        thread): each iteration sets the focus, takes or retrieves an
        exposure, measures its image quality, and yields the result. It
        does not decide when the sequence should stop (see
        :func:`continue_sequence`) or what to do once it has; callers are
        responsible for calling :func:`reset` beforehand and for reacting
        to each :class:`StepResult` as it's yielded (e.g., to update a
        plot).

        Parameters
        ----------
        method : :obj:`str`, :obj:`tuple`, optional
            The photometry method passed to
            :func:`photometry.image_quality`.

        Yields
        ------
        StepResult
            The result of the most recently completed step.
        """
        self.method = method
        while self.continue_sequence():
            self.observed_focus += [self.step_focus()]
            self.exposures += [self.take_exposure()]
            data, bkg, src_data, img_quality, source_stamp, coords \
                = image_quality(self.exposures[-1], method=method)
            self.source_stamps += [source_stamp]
            self.img_quality += [img_quality]
            self.centroids += [coords]

            result = StepResult(
                index=self.step_iter,
                focus_value=self.observed_focus[-1],
                exposure=self.exposures[-1],
                fwhm=img_quality,
                frame=data - bkg,
                stamp=source_stamp,
                centroid=coords,
                is_outlier=FocusPlot.is_outlier(self.centroids),
            )
            self.step_iter += 1
            yield result

    def reanalyze(self, method='brightest'):
        """
        Re-run photometry on every exposure this sequence has already
        collected, using a new ``method``, without taking any new
        exposures.

        This replaces the stored ``img_quality``/``source_stamps``/
        ``centroids`` in place; ``observed_focus`` and ``exposures`` (and
        thus which exposures exist and what focus values they were taken
        at) are left untouched. It's meant for interactively changing
        which source is used for the FWHM measurement (e.g., after a user
        clicks a different star in the displayed image) and seeing the
        effect on every exposure already taken, live or archived, without
        re-observing. If nothing has been collected yet (``exposures`` is
        empty), this yields nothing.

        Parameters
        ----------
        method : :obj:`str`, :obj:`tuple`, optional
            The photometry method passed to
            :func:`photometry.image_quality`.

        Yields
        ------
        StepResult
            The updated result for each already-collected exposure, in
            the order they were originally taken.
        """
        self.method = method
        focus_values = list(self.observed_focus)
        exposures = list(self.exposures)

        self.img_quality = []
        self.source_stamps = []
        self.centroids = []

        for i, (focus_value, exposure) in enumerate(zip(focus_values, exposures)):
            data, bkg, src_data, img_quality, source_stamp, coords \
                = image_quality(exposure, method=method)
            self.source_stamps += [source_stamp]
            self.img_quality += [img_quality]
            self.centroids += [coords]

            yield StepResult(
                index=i,
                focus_value=focus_value,
                exposure=exposure,
                fwhm=img_quality,
                frame=data - bkg,
                stamp=source_stamp,
                centroid=coords,
                is_outlier=FocusPlot.is_outlier(self.centroids),
            )

    def execute(self, verbose=True, goto=True, method='brightest', plot=True, **exp_kwargs):

        if self._exposure is not None:
            self._exposure.cfg.configure(**exp_kwargs)
        self.reset()

        self.plot = FocusPlot(self.expected_steps) if plot else None

        for result in self.step(method=method):
            if self.plot is not None:
                self.plot.update_frame(result.frame, result.centroid, int(result.fwhm*10),
                                        result.exposure.stem, result.focus_value,
                                        result.is_outlier)
                self.plot.add_stamp(result.index, result.stamp, result.focus_value,
                                     result.fwhm, is_outlier=result.is_outlier)
                self.plot.update_curve(self.observed_focus, self.img_quality)

        best_focus, best_img_quality = self.fit_best_focus(self.observed_focus, self.img_quality)

        if goto:
            if self._focus is None:
                raise ValueError(
                    'No telescope focus control is available (ktl not connected); cannot move '
                    'to the best-fit focus.'
                )
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

    @property
    def expected_steps(self):
        return self.nstep

    def continue_sequence(self):
        return self.step_iter < self.nstep

    def step_focus(self):
        self._focus.set_to(self.target_focus[self.step_iter])
        return self._focus.current


class AutomatedFocusSequence(FocusSequence):
    def __init__(self, start, step, maxsteps=12):

        if maxsteps is None or maxsteps < 2:
            raise ValueError('maxsteps cannot be None and must be at least 2.')

        super().__init__()

        self.start = start
        # NOTE: Named `step_size`, not `step`, because `step` is the name
        # of the generator method inherited from FocusSequence -- an
        # instance attribute of the same name would shadow it entirely
        # (seq.step(...) would try to call a float).
        self.step_size = step
        self.direction = None
        self.last = None
        self.maxsteps = maxsteps

    def reset(self):
        super().reset()
        self.direction = None
        self.last = None

    @property
    def expected_steps(self):
        return self.maxsteps

    def continue_sequence(self):
        return (
            self.step_iter < self.maxsteps
            and (self.step_iter < 2
                 or self.last is None
                 or (self.last is not None and self.step_iter < self.last)
            )
        )

    def step_focus(self):
        if self.step_iter == 0:
            next_focus = self.start
        elif self.step_iter == 1:
            next_focus = self.start + self.step_size
        elif self.step_iter == 2 and self.img_quality[0] > self.img_quality[1]:
            self.direction = 1
            next_focus = self.observed_focus[1] + self.step_size
        elif self.step_iter == 2 and self.img_quality[0] < self.img_quality[1]:
            self.direction = -1
            next_focus = self.observed_focus[0] - self.step_size
        elif self.last is None and self.step_iter > 2 and self.img_quality[-1] > self.img_quality[-2]:
            self.last = self.step_iter + 2
            if self.last > self.maxsteps:
                warnings.warn(
                    f'Number of steps to fulfill sequence ({self.last}) is more than the '
                    f'maximum number of steps requested ({self.maxsteps}).')
            next_focus = self.observed_focus[-1] + self.direction * self.step_size
        else:
            next_focus = self.observed_focus[-1] + self.direction * self.step_size
        self._focus.set_to(next_focus)
        return self._focus.current


class ArchiveFocusSequence(FocusSequence):
    def __init__(self, observed_focus, exposures):
        super().__init__()

        self._observed_focus = observed_focus
        self._exposures = exposures
        self.nstep = len(self._observed_focus)

    @property
    def expected_steps(self):
        return self.nstep

    def continue_sequence(self):
        return self.step_iter < self.nstep

    def step_focus(self):
        return self._observed_focus[self.step_iter]

    def take_exposure(self):
        return self._exposures[self.step_iter]


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
             'requested.  Requires --datadir, --prefix, and --suffix to locate the exposures; '
             'does not require a ktl connection.'
    )
    parser.add_argument('--datadir', default='.', type=str,
        help='Directory holding the archived exposures to re-analyze.  Only used with --obsnum.'
    )
    parser.add_argument('--prefix', default='n', type=str,
        help='Filename prefix of the archived exposures, e.g. "n" for files named like '
             '"n2165.fits".  Only used with --obsnum.'
    )
    parser.add_argument('--suffix', default='.fits', type=str,
        help='Filename suffix of the archived exposures, e.g. ".fits".  Only used with --obsnum.'
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
    parser.add_argument('--config', action='store_true',
        help='Only print the current exposure configuration and focus'
    )
    parser.add_argument('--verbose', action='store_true',
        help='Enable verbose output for debugging'
    )
    parser.add_argument('--no-plot', action='store_true',
        help='Disable the live focus-sequence plot.'
    )
    args = parser.parse_args()

    # Set the read speed
    _speed = args.speed
#    if _speed == 'Fast':
#        warnings.warn('Fast does not work!  Setting to slow.')
#        _speed == 'Slow'
#    if _speed == 'Slow':
#        _speed = '0.05MHz'
#    elif _speed == 'Fast':
#        _speed = '1.0MHz'

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
                raise FileNotFoundError('Expected to find the following files, but they are not '
                                        f'available: {", ".join(str(f) for f in missing)}')
            seq = ArchiveFocusSequence(seq.target_focus, expected_files)
    else:
        seq = AutomatedFocusSequence(args.focus[0], args.focus[1], maxsteps=args.maxsteps)

    best_focus, best_img_quality = seq.execute(goto=False, method=args.method, record=True,
                                               speed=_speed, exptime=args.exptime,
                                               binning=args.binning, plot=not args.no_plot)
    print(f'Best focus: {best_focus:.1f}')
    print(f'Expected sigma: {best_img_quality:.1f} pixels')

    if seq.plot is not None:
        pyplot.ioff()
        pyplot.show(block=True)

    # TODO:
    # - Write the output file if provided


if __name__ == "__main__":
    main()

