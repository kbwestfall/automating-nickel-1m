"""
Shared pytest configuration and fixtures for the Nickel automation test
suite.

Adds ``scripts/`` to :obj:`sys.path` so that ``focus``, ``photometry``, and
``quadratic`` import the same way they import each other today (a flat
namespace, matching how they run at the telescope), and forces a
non-interactive matplotlib backend before anything imports :mod:`focus`
(which imports :mod:`matplotlib.pyplot`), so tests never require a
display.

Per the project's testing strategy, everything here is expected to run
without the ``ktl`` package installed: fixtures only ever produce data for
:class:`focus.ArchiveFocusSequence`, which never touches
:class:`focus.Focus`/:class:`focus.Exposure` hardware.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / 'scripts'
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

import numpy as np
import pytest
from astropy.io import fits


def gaussian_frame(fwhm, amplitude=8000., background=200., noise_sigma=5.,
                    shape=(101, 101), center=None, seed=0):
    """
    Generate a small synthetic image with a single 2D Gaussian source.

    A Gaussian (rather than, e.g., a Moffat profile) is used because
    :func:`photometry.evaluate_shape` converts its moment-based sigma to a
    FWHM using the Gaussian relationship
    (``FWHM = 2*sqrt(2*ln(2))*sigma``), so a Gaussian source is what lets
    the requested ``fwhm`` be recovered (up to segmentation-mask
    truncation effects) rather than introducing an unrelated,
    profile-dependent bias.

    Parameters
    ----------
    fwhm : :obj:`float`
        Full width at half maximum of the source, in pixels.
    amplitude : :obj:`float`, optional
        Peak source flux, in counts.
    background : :obj:`float`, optional
        Constant background level, in counts.
    noise_sigma : :obj:`float`, optional
        Standard deviation of the Gaussian noise added to every pixel.
    shape : :obj:`tuple`, optional
        ``(ny, nx)`` shape of the returned array.
    center : :obj:`tuple`, optional
        ``(x, y)`` center of the source.  Defaults to the center of the
        frame.
    seed : :obj:`int`, optional
        Seed for the noise random number generator, so tests are
        deterministic.

    Returns
    -------
    :class:`numpy.ndarray`
        The synthetic image, as ``float32``.
    """
    ny, nx = shape
    if center is None:
        center = (nx / 2, ny / 2)
    y, x = np.mgrid[0:ny, 0:nx]

    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    source = amplitude * np.exp(-((x - center[0])**2 + (y - center[1])**2) / (2 * sigma**2))

    rng = np.random.default_rng(seed)
    data = background + source + rng.normal(scale=noise_sigma, size=shape)
    return data.astype(np.float32)


@pytest.fixture
def focus_sweep(tmp_path):
    """
    Write a small synthetic focus sequence to disk.

    The sequence follows a known quadratic,
    ``fwhm(focus) = fwhm_min + curvature*(focus - best_focus)**2``, so
    tests can check that fitting the sequence recovers ``best_focus`` and
    ``fwhm_min``. Filenames follow the ``{prefix}{obsnum}{suffix}``
    convention :func:`focus.main` expects for ``--obsnum``/``--datadir``/
    ``--prefix``/``--suffix`` archive-mode replay.

    Returns
    -------
    :obj:`dict`
        Keys: ``focus_values`` (:obj:`list` of :obj:`float`), ``files``
        (:obj:`list` of :class:`pathlib.Path`, one per focus value),
        ``best_focus``, ``fwhm_min`` (:obj:`float`, the true values used to
        generate the sweep), ``datadir`` (:class:`pathlib.Path`),
        ``prefix``, ``suffix`` (:obj:`str`), and ``obsnum`` (:obj:`int`,
        the observation number of the first file).
    """
    best_focus = 350.
    fwhm_min = 3.0
    curvature = 0.01
    focus_values = [340., 345., 350., 355., 360.]
    obsnum = 2000
    prefix, suffix = 'n', '.fits'

    files = []
    for i, focus_value in enumerate(focus_values):
        fwhm = fwhm_min + curvature * (focus_value - best_focus)**2
        data = gaussian_frame(fwhm, seed=i)
        path = tmp_path / f'{prefix}{obsnum + i}{suffix}'
        fits.writeto(path, data, overwrite=True)
        files.append(path)

    return {
        'focus_values': focus_values,
        'files': files,
        'best_focus': best_focus,
        'fwhm_min': fwhm_min,
        'datadir': tmp_path,
        'prefix': prefix,
        'suffix': suffix,
        'obsnum': obsnum,
    }
