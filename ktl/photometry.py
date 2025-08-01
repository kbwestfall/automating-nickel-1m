from astropy.io import fits
from photutils.psf import fit_fwhm
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.segmentation import detect_threshold, detect_sources
from scipy.ndimage import binary_dilation
from astropy.visualization import (ImageNormalize, ZScaleInterval, LinearStretch)


import quadratic

import numpy as np
import matplotlib.pyplot as plt
import argparse

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

        fig = plt.figure(figsize=(20, 10), constrained_layout=True)
    
        # Create subfigures: left, center, right
        subfigs = fig.subfigures(1, 3, width_ratios=[1, 2, 1])
        
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

        # Add title
        fig.suptitle(f'Focus Sequence Analysis ({num_rows}×{num_cols})', 
                    fontsize=16, fontweight='bold')

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


def find_sources(data, max_iterations=5, verbose=False):

    previous_threshold = None

    for iteration in range(max_iterations):
        
        sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
        threshold = detect_threshold(data, nsigma=5.0, sigma_clip=sigma_clip)
        
        sources = detect_sources(data, threshold, npixels=10)

        source_mask = sources.data_ma
        structure = np.ones((7, 7), dtype=bool)
        grown_source_mask = binary_dilation(source_mask, structure=structure)


        mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=grown_source_mask)

        data = data - median

        if previous_threshold is not None:
            threshold_median = np.median(threshold)
            previous_threshold_median = np.median(previous_threshold)
            threshold_change = abs(threshold_median - previous_threshold_median) / previous_threshold_median
            if verbose:
                print(f"Threshold change: {threshold_change:.6f}")
            
            if threshold_change < 0.01:
                print(f"Converged after {iteration + 1} iterations\n")
                return data, sources
        
        previous_threshold = threshold


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

    if verbose:
        print(f"Cutout size: {cutout_size}x{cutout_size} pixels")
        print(f"Cutout bounds: x=[{x_min}:{x_max}], y=[{y_min}:{y_max}]")
        print(f"Actual cutout shape: {cutout_data.shape}")

    fit = fit_fwhm(cutout_data)
        
    plot.set_center_axis(cutout_data, obs_num, focus_value)

    return fit


def evaluate_sources(data, sources, focus_coords, verbose=False):
    print(f"Number of sources detected: {sources.nlabels}\n")
    results = []

    for source in sources.segments:

        source_mask = sources.data == source.label
        if verbose:
            print(f"Evaluating source with label {source.label}")
        shape = evaluate_shape(data, source_mask, verbose=verbose)

        source_attr = {
            'Label': source.label,
            'M0': shape['M0'],
            'Centroid': shape['Centroid'],
            'FWHM': shape['FWHM']
        }

        results.append(source_attr)

    min_dist = float('inf')
    focus_star = None

    if focus_coords is not None:
        for source in results:
            centroid_x, centroid_y = source['Centroid']
            dist = np.sqrt((centroid_x - focus_coords[0])**2 + (centroid_y - focus_coords[1])**2)

            if dist < min_dist:
                min_dist = dist
                focus_star = source
    else:
        focus_star = max(results, key=lambda x: x['M0'])


    return focus_star


def photometry(fits_file, obs_num, focus_value, plot, focus_coords, verbose=False):

    if fits_file:
        hdu = fits.open(fits_file)
        data = hdu[0].data
        plot.set_left_axis(data)
    else:
        raise ValueError("No FITS file provided. Please specify a file.")

    data, sources = find_sources(data)
    focus_star = evaluate_sources(data, sources, focus_coords, verbose=verbose)

    if focus_star['FWHM'] is not None:
        focus_star['Focus'] = focus_value
        focus_star['ObsNum'] = obs_num

        fit = cutout(data, focus_star, obs_num, focus_value, plot, verbose=verbose)

        print(f"Using source with label {focus_star['Label']} at ({focus_star['Centroid'][0]}, {focus_star['Centroid'][1]}) with FWHM: {focus_star['FWHM']} and fit_fwhm: {fit}\n")
        return focus_star
    else:
        return None





