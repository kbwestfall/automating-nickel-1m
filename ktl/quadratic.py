import numpy as np


def fit_quadratic(x_values, y_values):
    x = np.array(x_values)
    y = np.array(y_values)

    if len(x) < 3:
        raise ValueError("At least three data points are required to fit a quadratic function.")
    
    coefficients = np.polyfit(x, y, 2)
    a, b, c = coefficients

    return a, b, c

def vertex(a, b, c):
    x_vertex = -b / (2 * a)
    y_vertex = c - (b**2 / (4 * a))
    return x_vertex, y_vertex

