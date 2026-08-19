import warnings
from IPython import embed

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.ndimage import binary_dilation

from astropy.io import fits
from astropy import table
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils import segmentation
from astropy.visualization import (ImageNormalize, ZScaleInterval, LinearStretch)

# detect threshold
# use threshold to detect sources
# grow sources with binary dilation
# estimate background on source subtracted image
# inject sources back into image
# iterate while threshold converges

# use moments to determine shape of source(s)
# decide whether shape is reasonable
    # ratio between x and y aren't too different and x and y aren't too big/small
# return sources and their attributes

# find source closest to provided coordinates else brightest source
# return fwhm of that source

class Grid:

    def __init__(self, num_rows=4, num_cols=3):
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.index = 0

        self.fig = plt.figure(figsize=(20, 10), constrained_layout=True)
    
        # Create subfigures: left, center, right
        subfigs = self.fig.subfigures(1, 3, width_ratios=[1, 2, 1])

        # Left panel
        self.ax_left = subfigs[0].subplots(1, 1)
        self.ax_left.set_aspect('equal')
        self.ax_left.set_title('Source Location')

        # Center grid
        self.center_axes = subfigs[1].subplots(num_rows, num_cols)
        if num_rows == 1 and num_cols == 1:
            self.center_axes = [self.center_axes]
        elif num_rows == 1:
            pass
        elif num_cols == 1:
            self.center_axes = [[ax] for ax in self.center_axes]

        # Format center subplots
        for i in range(num_rows):
            for j in range(num_cols):
                ax = self.get_center_axis(i, j)
                ax.set_aspect('equal')
                

        # Right panel
        self.ax_right = subfigs[2].subplots(1, 1)
        self.ax_right.set_aspect('equal')
        self.ax_right.set_title('Focus Curve')


    def set_left_axis(self, data):
        self.ax_left.clear()
        norm = ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
        self.ax_left.imshow(data, cmap='viridis', origin='lower', norm=norm)
        self.ax_left.set_xticks([])
        self.ax_left.set_yticks([])
        self.ax_left.set_title('Full Image')

    def set_center_axis(self, data, obs_num, focus_value, rect):
        ax = self.get_center_axis_by_index(self.index)
        if ax is None:
            print(f"Index {self.index} is out of bounds for the center axes grid.")
            return
        
        ax.clear()
        ax.imshow(data, origin='lower', cmap='viridis')
        ax.set_title(f'd{obs_num} Focus: {focus_value}')
        ax.set_xticks([])
        ax.set_yticks([])

        self.ax_left.add_patch(rect)

        plt.draw()
        plt.pause(1)

        # Add red rectangle around the cutout region
        # from matplotlib.patches import Rectangle
        # rect = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
        #                 linewidth=2, edgecolor='red', facecolor='none')
        # ax.add_patch(rect)

        # cutout_center_x = focus_star['Centroid'][0] - x_min
        # cutout_center_y = focus_star['Centroid'][1] - y_min
        # ax.plot(cutout_center_x, cutout_center_y, 'r+', markersize=15, markeredgewidth=3)

    def get_center_axis(self, row, col):
        """Get center axis by row/col coordinates"""
        if self.num_rows == 1 and self.num_cols == 1:
            return self.center_axes[0]
        elif self.num_rows == 1:
            return self.center_axes[col]
        elif self.num_cols == 1:
            return self.center_axes[row][0]
        else:
            return self.center_axes[row][col]
    
    def get_center_axis_by_index(self, linear_index):
        """Get center axis by linear index"""
        row = linear_index // self.num_cols
        col = linear_index % self.num_cols
        
        if row >= self.num_rows:
            return None  # Index out of bounds
        
        return self.get_center_axis(row, col)

    def set_right_axis(self, x_values, y_values, a, b, c, x_vertex, y_vertex, obs_values, outliers):
        self.ax_right.clear()

        x_min, x_max = min(x_values), max(x_values)
        x_smooth = np.linspace(x_min, x_max, 50)
        y_smooth = a * x_smooth**2 + b * x_smooth + c

        self.ax_right.scatter(x_values, y_values, color='blue')
        for i, obs in enumerate(obs_values):
            self.ax_right.text(x_values[i] + 1, y_values[i], f'{obs}', color='blue', fontsize=12)
        if (outliers):
            for outlier in outliers:
                self.ax_right.scatter(outlier['focus_value'], outlier['fwhm'], color='yellow', label=f'Outlier Obs {outlier["obs"]}')
                self.ax_right.text(outlier['focus_value'] + 1, outlier['fwhm'], f'{outlier["obs"]}', color='yellow', fontsize=12)
        self.ax_right.plot(x_smooth, y_smooth, 'r-')
        self.ax_right.scatter([x_vertex], [y_vertex], color='green', label='Optimal focus', zorder=3)
        self.ax_right.axvline(x=x_vertex, color='green', linestyle='--')

        self.ax_right.set_xlabel('Focus Value', fontsize=12)
        self.ax_right.set_ylabel('FWHM (pixels)', fontsize=12)
        self.ax_right.set_title('Focus Curve Analysis', fontsize=14, fontweight='bold')
        self.ax_right.legend(loc='upper right')
        self.ax_right.grid(True, alpha=0.3)

        # Calculate aspect ratio to make plot square
        x_range = max(x_values) - min(x_values)
        y_range = max(y_values) - min(y_values)
        x_padding = x_range * 0.1
        y_padding = y_range * 0.1
        self.ax_right.set_xlim(min(x_values) - x_padding, max(x_values) + x_padding)
        self.ax_right.set_ylim(min(y_values) - y_padding, max(y_values) + y_padding)
        aspect_ratio = (x_range + 2*x_padding) / (y_range + 2*y_padding)
        self.ax_right.set_aspect(aspect_ratio, adjustable='box')


def find_sources(data, max_iterations=5, grow=7, atol=0.1, rtol=0.01, verbose=False):
    """
    Find sources in an image.

    Iteratively uses `photutils.segmentation.detect_threshold` to determine the
    image detection threshold and `photutils.segmentation.detect_sources` to
    identify pixels with sources.  Source pixels masks are grown and the
    non-source pixels are used to measure and subtract the background.
    Iterations stop when subsequent measurements of the threshold are within the
    provided relative tolerance.

    Parameters
    ----------
    data : numpy.ndarray
        Image data
    max_iteration : :obj:`int`, optional
        Maximum number of iterations
    grow : :obj:`int`, optional
        Number of pixels to grow the source mask.
    atol : :obj:`float`, optional
        Absolute tolerance used to test for convergence of the detection
        threshold.  See `numpy.isclose`.
    rtol : :obj:`float`, optional
        Relative tolerance used to test for convergence of the detection
        threshold.  See `numpy.isclose`.
    verbose : :obj:`bool`, optional
        Print progress

    Returns
    -------
    background : :obj:`float`
        The estimated background in the image
    source_mask : `photutils.segmentation.core.SegmentationImage`
        Segmentation image.
    """
    previous_threshold = None               # Previous threshold value
    bkg = 0.0                               # Background
    # Sigma-clipping object
    sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
    # Structure used to grow the mask
    structure = np.ones((grow, grow), dtype=bool)

    for iteration in range(max_iterations):

        # Subtract the background
        _data = data - bkg
        # Get the threshold image     
        threshold = segmentation.detect_threshold(_data, nsigma=5.0, sigma_clip=sigma_clip)
        # Detect sources above the threshold
        sources = segmentation.detect_sources(_data, threshold, npixels=10)
        # Grow the mask
        grown_source_mask = binary_dilation(sources.data_ma, structure=structure)
        # Get the background and add it to the total
        bkg += sigma_clipped_stats(_data, sigma=3.0, mask=grown_source_mask)[1]
        if verbose:
            print(f'Updated background: {bkg:.1f}')

        # Calculate the median of the threshold image
        med_threshold = np.median(threshold)
        if previous_threshold is None or \
                not np.isclose(med_threshold, previous_threshold, atol=atol, rtol=rtol):
            if verbose:
                print(f'Updated threshold: {med_threshold:.1f}')
            # This is the first iteration
            previous_threshold = med_threshold
            continue

        # Converged
        if verbose:
            print(f'Source detection converged after {iteration+1} iterations')
        break

    return bkg, sources


def evaluate_shape(data, source_mask, verbose=False):

    source_data = data * source_mask
    M0 = np.sum(source_data)

    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])

    meshgrid_x, meshgrid_y = np.meshgrid(x, y)

    M1_x = np.sum(meshgrid_x * source_data) / M0
    M1_y = np.sum(meshgrid_y * source_data) / M0

    M2_x = np.sum((meshgrid_x)** 2 * source_data) / M0
    M2_y = np.sum((meshgrid_y)** 2 * source_data) / M0

    sigma_x = np.sqrt(M2_x - M1_x**2)
    sigma_y = np.sqrt(M2_y - M1_y**2)

    if verbose:
        print(f"M0: {M0}; M1: {M1_x:.2f}, {M1_y:.2f}; M2: {M2_x:.2f}, {M2_y:.2f}")
        print(f"Sigma_x: {sigma_x:.10e}, Sigma_y: {sigma_y:.10e}")

    if sigma_x <= 1e-8 or sigma_y <= 1e-8:
        if verbose:
            print("Warning: sigma_x or sigma_y is too small, cannot compute FWHM accurately.\n")

        shape = {
            'M0': 0,
            'Centroid': (M1_x, M1_y),
            'FWHM': None
        }
        return shape
    elif min(sigma_x, sigma_y) / max(sigma_x, sigma_y) < 0.5:
        if verbose:
            print("Warning: sigma_x and sigma_y are too different, cannot compute FWHM accurately.\n")
        shape = {
            'M0': 0,
            'Centroid': (M1_x, M1_y),
            'FWHM': None
        }
        return shape

    FWHM_x = 2 * np.sqrt(2 * np.log(2)) * sigma_x
    FWHM_y = 2 * np.sqrt(2 * np.log(2)) * sigma_y
    average_FWHM = (FWHM_x + FWHM_y) / 2

    shape = {
        'M0': M0,
        'Centroid': (float(M1_x), float(M1_y)),
        'FWHM': average_FWHM
    }

    return shape

def cutout(data, focus_star, obs_num, focus_value, plot, verbose=False):
    cutout_size = int(3 * focus_star['FWHM'])
    if cutout_size % 2 == 0:
        cutout_size += 1

    half_size = cutout_size // 2
    
    # Calculate cutout boundaries
    x_center = int(round(focus_star['Centroid'][0]))
    y_center = int(round(focus_star['Centroid'][1]))

    x_min = max(0, x_center - half_size)
    x_max = min(data.shape[1], x_center + half_size + 1)
    y_min = max(0, y_center - half_size)
    y_max = min(data.shape[0], y_center + half_size + 1)
    
    # Extract cutout from data
    cutout_data = data[y_min:y_max, x_min:x_max]

     # highlight cutout region
    rect = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, linewidth=2, edgecolor='red', facecolor='none')

    if verbose:
        print(f"Cutout size: {cutout_size}x{cutout_size} pixels")
        print(f"Cutout bounds: x=[{x_min}:{x_max}], y=[{y_min}:{y_max}]")
        print(f"Actual cutout shape: {cutout_data.shape}")

#    fit = fit_fwhm(cutout_data)
        
    plot.set_center_axis(cutout_data, obs_num, focus_value, rect)

#    return fit


def moment2d(x, y, z):
    """
    Calculate moments of a 2D dataset.

    Parameters
    ----------
    x : array-like
        First coordinate of the data
    y : array-like
        Second coordinate of the data
    z : array-like
        Value of the data at each provide x and y coordinate.

    Returns:
    tot : :obj:`float`
        Sum of ``z``
    cx : :obj:`float`
        Weighted mean of ``x``
    cy : :obj:`float`
        Weighted mean of ``y``
    sx : :obj:`float`
        Weighted standard deviation of ``x``
    sy : :obj:`float`
        Weighted standard deviation of ``y``
    """
    _x = np.asarray(x)
    _y = np.asarray(y)
    _z = np.asarray(z)

    tot = np.sum(_z)
    if np.absolute(tot) < 1e-6:
        raise ValueError('Sum of the data is too close to 0.')

    cx = np.sum(_x*_z)/tot
    cy = np.sum(_y*_z)/tot

    sx = np.sqrt(np.sum(_x**2*_z)/tot - cx**2)
    sy = np.sqrt(np.sum(_y**2*_z)/tot - cy**2)

    return tot, cx, cy, sx, sy


def empty_source_table(length):
    return table.Table([
        table.Column(name='ID', dtype=int, length=length, description='Source ID'),
        table.Column(name='CNTS', dtype=float, length=length, description='Total counts'),
        table.Column(name='CENX', dtype=float, length=length, description='X centroid'),
        table.Column(name='CENY', dtype=float, length=length, description='Y centroid'),
        table.Column(name='SIGX', dtype=float, length=length, description='X sigma'),
        table.Column(name='SIGY', dtype=float, length=length, description='Y sigma'),
    ])


def evaluate_sources(data, sources, verbose=False):
    """
    Provided a source segmentation image, measure the first 3 moments of all
    sources.

    Parameters
    ----------
    data : `numpy.ndarray`
        Background subtracted, raw image data
    sources : `photutils.segmentation.core.SegmentationImage`
        Source segmentation image.
    verbose : :obj:`bool`, optional
        Print progress messages

    Returns
    -------
    `astropy.table.Table`
        Table with the source measurements.  Columns are:
            - ID: Source ID number
            - CENX : center along X (column)
            - CENY : center along Y (row)
            - SIGX : dispersion along X (column)
            - SIGY : dispersion along Y (row)
    """
    if verbose:
        print(f"Number of sources detected: {sources.nlabels}\n")
    if sources.nlabels == 0:
        raise ValueError('No sources found.')

    # Construct the coordinate images
    img_x, img_y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))

    # For each source, get the total flux, and the 1st and 2nd moments along each axis.
    src_data = empty_source_table(sources.nlabels)
    for i, source in enumerate(sources.segments):
        if verbose:
            print(f"Evaluating source with label {source.label}")
        src_data['ID'][i] = source.label
        indx = sources.data == source.label
        if np.sum(indx) == 0:
            warnings.warn(f'No data assocated with source {source.label}')
        src_data['CNTS'][i], src_data['CENX'][i], src_data['CENY'][i], \
            src_data['SIGX'][i], src_data['SIGY'][i] \
                = moment2d(img_x[indx], img_y[indx], data[indx])

    # Check if sources are real        
    #   - Bad sigma measurements
    small_sources = (src_data['SIGX'] < 1e-8) |  (src_data['SIGY'] < 1e-8)
    if np.any(small_sources):
        if verbose:
            warnings.warn(f'Removing {np.sum(small_sources)} sources with errantly small widths.')
        src_data = src_data[np.logical_not(small_sources)]

    if len(src_data) == 0:
        raise ValueError('No good sources found')

    axis_ratio = src_data['SIGX'] / src_data['SIGY']
    ellip_sources = (axis_ratio < 0.5) | (axis_ratio > 2)
    if np.any(ellip_sources):
        if verbose:
            warnings.warn(f'Removing {np.sum(small_sources)} sources with large ellipticity.')
        src_data = src_data[np.logical_not(ellip_sources)]

    if len(src_data) == 0:
        raise ValueError('No good sources found')

    return src_data


def image_quality(fits_file, method='brightest', verbose=False):
    """
    Evaluate the image quality of the provided data.

    Parameters
    ----------
    fits_file : :obj:`str`, `Path`
        File with raw image data
    method : :obj:`str`, :obj:`tuple`, optional
        Method used to measure the image quality.  Must be:

            - ``brightest``: Return image quality measurement for the brightest
              source in the field

            - ``weighted``: Return the flux-weighted mean of the image-quality
              measurements for all sources.

            - :obj:`tuple`: Provide a tuple with coordinates and use the source
              closest to the provides coordinates.
    verbose : :obj:`bool`, optional
        Print status messages

    Returns
    -------

    """
    with fits.open(fits_file) as hdu:
        data = hdu[0].data.astype(float)
    bkg, sources = find_sources(data)
    src_data = evaluate_sources(data-bkg, sources, verbose=verbose)

    if method == 'brightest':
        target_source = np.argmax(src_data['CNTS'])
        img_quality = (src_data['SIGX'][target_source] + src_data['SIGY'][target_source])/2
        coords = (src_data['CENX'][target_source], src_data['CENY'][target_source])
        stamp = extract_stamp(data-bkg, coords, int(img_quality*10))
    elif method == 'weighted':
        img_quality = np.sum(src_data['CNTS'] * (src_data['SIGX'] + src_data['SIGY'])/2) \
                        / np.sum(src_data['CNTS'])
        target_source = np.argmax(src_data['CNTS'])
        coords = (src_data['CENX'][target_source], src_data['CENY'][target_source])
        stamp = extract_stamp(data-bkg, coords, int(img_quality*10))
    elif not isinstance(method, tuple):
        raise ValueError('image_quality method must be brightest, weighted, or a tuple of '
                         'coordinates')
    else:
        try:
            dist = (src_data['CENX'] - method[0])**2 + (src_data['CENY'] - method[1])**2
        except Exception as e:
            raise ValueError('Could not use tuple provided to method keyword to find nearest '
                             'source to use for image quality measurement.  Original excception '
                             f'message: {e}.')
        target_source = np.argmin(dist)
        img_quality = (src_data['SIGX'][target_source] + src_data['SIGY'][target_source])/2
        coords = (src_data['CENX'][target_source], src_data['CENY'][target_source])
        stamp = extract_stamp(data-bkg, coords, int(img_quality*10))

    return data, bkg, src_data, img_quality, stamp


def extract_stamp(data, coords, size):
    sx = int(np.floor(coords[0] - size / 2))
    ex = int(np.ceil(coords[0] + size / 2)) + 1
    sy = int(np.floor(coords[1] - size / 2))
    ey = int(np.ceil(coords[1] + size / 2)) + 1
    return data[sy:ey,sx:ex]


