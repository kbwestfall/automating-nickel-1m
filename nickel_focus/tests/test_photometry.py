"""Tests for :mod:`photometry`."""
import logging

import numpy as np

from nickel_focus import photometry


def _source_frame():
    """A small synthetic image: one bright square source on a flat background."""
    data = np.full((51, 51), 200., dtype=float)
    data[20:31, 20:31] = 8000.
    return data


def test_find_sources_debug_logging(caplog):
    data = _source_frame()

    with caplog.at_level(logging.DEBUG, logger='nickel_focus'):
        photometry.find_sources(data)

    messages = [r.getMessage() for r in caplog.records if r.name == 'nickel_focus']
    assert any('Updated background' in m for m in messages), \
        'find_sources should log its background updates at debug level'
    assert any('converged' in m for m in messages), \
        'find_sources should log convergence progress at debug level'


def test_find_sources_silent_at_info_level(caplog):
    data = _source_frame()

    with caplog.at_level(logging.INFO, logger='nickel_focus'):
        photometry.find_sources(data)

    nickel_focus_records = [r for r in caplog.records if r.name == 'nickel_focus']
    assert nickel_focus_records == [], \
        'find_sources should produce no records at all once debug logging is off'


def test_evaluate_shape_warns_on_degenerate_sigma(caplog):
    # A single hot pixel has zero moment spread in both axes.
    data = np.zeros((10, 10))
    data[5, 5] = 100.
    source_mask = data > 0

    with caplog.at_level(logging.WARNING, logger='nickel_focus'):
        shape = photometry.evaluate_shape(data, source_mask)

    assert shape['FWHM'] is None, 'a degenerate source should report no FWHM'
    assert any(
        r.levelname == 'WARNING' and 'too small' in r.getMessage()
        for r in caplog.records if r.name == 'nickel_focus'
    ), 'a degenerate-sigma source should log a warning'


def test_evaluate_shape_warns_on_dissimilar_sigma(caplog):
    # A source spread widely in x but confined to two adjacent rows in y
    # gives a much larger sigma_x than sigma_y.
    data = np.zeros((12, 12))
    data[5, 0:10] = 100.
    data[6, 0:10] = 100.
    source_mask = data > 0

    with caplog.at_level(logging.WARNING, logger='nickel_focus'):
        shape = photometry.evaluate_shape(data, source_mask)

    assert shape['FWHM'] is None, 'a too-elliptical source should report no FWHM'
    assert any(
        r.levelname == 'WARNING' and 'too different' in r.getMessage()
        for r in caplog.records if r.name == 'nickel_focus'
    ), 'a too-elliptical source should log a warning'
