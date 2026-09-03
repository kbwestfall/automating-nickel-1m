from astropy.io import fits
from astropy import table
from astropy.stats import SigmaClip, sigma_clipped_stats
from IPython import embed
import numpy as np
from photutils import segmentation
from scipy.ndimage import binary_dilation

from nickel_focus import log

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


def find_sources(data, max_iterations=5, grow=7, atol=0.1, rtol=0.01):
    """
    Find sources in an image.

    Iteratively uses :func:`photutils.segmentation.detect_threshold` to
    determine the image detection threshold and
    :func:`photutils.segmentation.detect_sources` to identify pixels with
    sources.  Source pixels masks are grown and the non-source pixels are used
    to measure and subtract the background.  Iterations stop when subsequent
    measurements of the threshold are within the provided relative tolerance.

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
        threshold.  See :func:`numpy.isclose`.
    rtol : :obj:`float`, optional
        Relative tolerance used to test for convergence of the detection
        threshold.  See :func:`numpy.isclose`.

    Returns
    -------
    background : :obj:`float`
        The estimated background in the image
    source_mask : :class:`photutils.segmentation.core.SegmentationImage`
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
        threshold = segmentation.detect_threshold(_data, n_sigma=5.0, sigma_clip=sigma_clip)
        # Detect sources above the threshold
        sources = segmentation.detect_sources(_data, threshold, n_pixels=10)
        # Grow the mask
        grown_source_mask = binary_dilation(sources.data_masked, structure=structure)
        # Get the background and add it to the total
        bkg += sigma_clipped_stats(_data, sigma=3.0, mask=grown_source_mask)[1]
        log.debug(f'Updated background: {bkg:.1f}')

        # Calculate the median of the threshold image
        med_threshold = np.median(threshold)
        if previous_threshold is None or \
                not np.isclose(med_threshold, previous_threshold, atol=atol, rtol=rtol):
            log.debug(f'Updated threshold: {med_threshold:.1f}')
            # This is the first iteration
            previous_threshold = med_threshold
            continue

        # Converged
        log.debug(f'Source detection converged after {iteration+1} iterations')
        break

    return bkg, sources


def evaluate_shape(data, source_mask):

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

    log.debug(f"M0: {M0}; M1: {M1_x:.2f}, {M1_y:.2f}; M2: {M2_x:.2f}, {M2_y:.2f}")
    log.debug(f"Sigma_x: {sigma_x:.10e}, Sigma_y: {sigma_y:.10e}")

    if sigma_x <= 1e-8 or sigma_y <= 1e-8:
        log.warning('sigma_x or sigma_y is too small, cannot compute FWHM accurately.')

        shape = {
            'M0': 0,
            'Centroid': (M1_x, M1_y),
            'FWHM': None
        }
        return shape
    elif min(sigma_x, sigma_y) / max(sigma_x, sigma_y) < 0.5:
        log.warning('sigma_x and sigma_y are too different, cannot compute FWHM accurately.')
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

    Returns: tot : :obj:`float`
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


def evaluate_sources(data, sources):
    """
    Provided a source segmentation image, measure the first 3 moments of all
    sources.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        Background subtracted, raw image data
    sources : :class:`photutils.segmentation.core.SegmentationImage`
        Source segmentation image.

    Returns
    -------
    :class:`astropy.table.Table`
        Table with the source measurements.  Columns are:
            - ID: Source ID number
            - CENX : center along X (column)
            - CENY : center along Y (row)
            - SIGX : dispersion along X (column)
            - SIGY : dispersion along Y (row)
    """
    log.info(f"Number of sources detected: {sources.n_labels}")
    if sources.n_labels == 0:
        raise ValueError('No sources found.')

    # Construct the coordinate images
    img_x, img_y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))

    # For each source, get the total flux, and the 1st and 2nd moments along each axis.
    src_data = empty_source_table(sources.n_labels)
    for i, source in enumerate(sources.segments):
        log.info(f"Evaluating source with label {source.label}")
        src_data['ID'][i] = source.label
        indx = sources.data == source.label
        if np.sum(indx) == 0:
            log.warning(f'No data assocated with source {source.label}')
        src_data['CNTS'][i], src_data['CENX'][i], src_data['CENY'][i], \
            src_data['SIGX'][i], src_data['SIGY'][i] \
                = moment2d(img_x[indx], img_y[indx], data[indx])

    # Check if sources are real
    #   - Bad sigma measurements
    small_sources = (src_data['SIGX'] < 1e-8) |  (src_data['SIGY'] < 1e-8)
    if np.any(small_sources):
        log.warning(f'Removing {np.sum(small_sources)} sources with errantly small widths.')
        src_data = src_data[np.logical_not(small_sources)]

    if len(src_data) == 0:
        raise ValueError('No good sources found')

    axis_ratio = src_data['SIGX'] / src_data['SIGY']
    ellip_sources = (axis_ratio < 0.5) | (axis_ratio > 2)
    if np.any(ellip_sources):
        log.warning(f'Removing {np.sum(small_sources)} sources with large ellipticity.')
        src_data = src_data[np.logical_not(ellip_sources)]

    if len(src_data) == 0:
        raise ValueError('No good sources found')

    return src_data


def image_quality(fits_file, method='brightest'):
    """
    Evaluate the image quality of the provided data.

    Parameters
    ----------
    fits_file : :obj:`str`, :class:`pathlib.Path`
        File with raw image data
    method : :obj:`str`, :obj:`tuple`, optional
        Method used to measure the image quality.  Must be:

            - ``brightest``: Return image quality measurement for the brightest
              source in the field

            - ``weighted``: Return the flux-weighted mean of the image-quality
              measurements for all sources.

            - :obj:`tuple`: Provide a tuple with coordinates and use the source
              closest to the provides coordinates.

    Returns
    -------
    data : :class:`numpy.ndarray`
        Raw image data.
    bkg : :obj:`float`
        Estimated background.
    src_data : :class:`astropy.table.Table`
        Table with the source measurements; see
        :func:`~nickel_focus.photometry.evaluate_sources`.
    img_quality : :obj:`float`
        The image-quality measurement (mean of the source's x and y sigma).
    stamp : :class:`numpy.ndarray`
        Cutout of the background-subtracted data around the selected source.
    coords : :obj:`tuple`
        The (x,y) centroid of the selected source.
    """
    with fits.open(fits_file) as hdu:
        data = hdu[0].data.astype(float)
    bkg, sources = find_sources(data)
    src_data = evaluate_sources(data-bkg, sources)

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

    return data, bkg, src_data, img_quality, stamp, coords


def extract_stamp(data, coords, size):
    sx = int(np.floor(coords[0] - size / 2))
    ex = int(np.ceil(coords[0] + size / 2)) + 1
    sy = int(np.floor(coords[1] - size / 2))
    ey = int(np.ceil(coords[1] + size / 2)) + 1
    return data[sy:ey,sx:ex]


