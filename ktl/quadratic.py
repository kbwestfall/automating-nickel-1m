import numpy as np
from photometry import photometry
from simulation import simulation
from focus import Image

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

# x_values = [335, 340, 345, 350, 355, 360]
# y_values = [16, 12.4, 9.6, 7.6, 6.4, 6.0]

# a, b, c = fit_quadratic(x_values, y_values)
# print(f"coefficients: a={a}, b={b}, c={c}")

# x_vertex, y_vertex = vertex(a, b, c)
# print(f"Vertex: ({x_vertex}, {y_vertex})")

def main():
    focus_values = np.arange(335, 390, 5)
    fwhm_values = []

    for focus in focus_values:
        image = Image(focus)
        fwhm_values.append(image.fwhm)

    a, b, c = fit_quadratic(focus_values, fwhm_values)
    x_vertex, y_vertex = vertex(a, b, c)

    print(f"Optimal focus: {x_vertex}, FWHM: {y_vertex}")

if __name__ == "__main__":
    main()
