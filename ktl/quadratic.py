import numpy as np


def fit_quadratic(x, y):
    _x = np.asarray(x)
    _y = np.asarray(y)
    if _x.size < 3:
        raise ValueError("At least three data points are required to fit a quadratic function.")
    if _y.size != _x.size:
        raise ValueError('Mismatch between the number of x and y values.')
    return np.polyfit(_x, _y, 2)


def vertex(a, b, c):
    return -b / (2 * a), c - (b**2 / (4 * a))

