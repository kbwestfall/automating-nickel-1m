"""Tests for :mod:`quadratic`."""
import numpy as np
import pytest

import quadratic


def test_fit_quadratic_recovers_known_coefficients():
    a_true, b_true, c_true = 0.02, -1.5, 25.
    x = np.array([10., 20., 30., 40., 50.])
    y = a_true * x**2 + b_true * x + c_true

    a, b, c = quadratic.fit_quadratic(x, y)

    assert a == pytest.approx(a_true, abs=1e-8), 'wrong quadratic coefficient'
    assert b == pytest.approx(b_true, abs=1e-6), 'wrong linear coefficient'
    assert c == pytest.approx(c_true, abs=1e-4), 'wrong constant coefficient'


def test_fit_quadratic_requires_at_least_three_points():
    with pytest.raises(ValueError):
        quadratic.fit_quadratic([1., 2.], [1., 2.])


def test_fit_quadratic_requires_matching_lengths():
    with pytest.raises(ValueError):
        quadratic.fit_quadratic([1., 2., 3.], [1., 2.])


def test_vertex():
    # y = 2*(x-5)**2 + 3 = 2x^2 - 20x + 53
    x_vertex, y_vertex = quadratic.vertex(2., -20., 53.)
    assert x_vertex == pytest.approx(5.), 'wrong vertex x-coordinate'
    assert y_vertex == pytest.approx(3.), 'wrong vertex y-coordinate'
