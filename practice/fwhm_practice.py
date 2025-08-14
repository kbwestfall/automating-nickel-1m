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

        fig = plt.figure(figsize=(20, 10), constrained_layout=True)
    
        # Create subfigures: left, center, right
        subfigs = fig.subfigures(1, 3, width_ratios=[1, 2, 1])
        
        # Left panel
        self.ax_left = subfigs[0].subplots(1, 1)
        self.ax_left.set_aspect('equal')
        # self.ax_left.text(0.5, 0.5, 'Full Image\nwith Cutout Box', ha='center', va='center', fontsize=14)
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
                # ax.text(0.5, 0.5, f'Cutout\n{i*num_cols + j + 1}', ha='center', va='center', fontsize=10)

        # Right panel
        self.ax_right = subfigs[2].subplots(1, 1)
        # self.ax_right.set_aspect('equal', adjustable='box')
        # self.ax_right.text(0.5, 0.5, 'FWHM vs\nFocus Value', ha='center', va='center', fontsize=14)
        self.ax_right.set_title('Focus Curve')

        # Add title
        fig.suptitle(f'Focus Sequence Analysis ({num_rows}×{num_cols})', 
                    fontsize=16, fontweight='bold')

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


def find_sources(data, max_iterations=5, verbose=False):

    previous_threshold = None

    for iteration in range(max_iterations):
        
        sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
        threshold = detect_threshold(data, nsigma=5.0, sigma_clip=sigma_clip)
        
        sources = detect_sources(data, threshold, npixels=10)

        source_mask = sources.data_ma
        structure = np.ones((7, 7), dtype=bool)
        grown_source_mask = binary_dilation(source_mask, structure=structure)

        if verbose:
            print(f"Original data pixel count: {np.sum(data)}")

        # fig, axes = plt.subplots(2, 2, figsize=(14, 14))

        # axes[0][0].imshow(source_mask, origin='lower', cmap='gray')
        # axes[0][0].set_title('Original Source Mask')

        # axes[0][1].imshow(grown_source_mask, origin='lower', cmap='gray')
        # axes[0][1].set_title('Grown Source Mask (3x3)')

        # axes[1][0].imshow(data, origin='lower', cmap='gray')
        # axes[1][0].set_title('Original Data')

        mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=grown_source_mask)

        data = data - median

        if verbose:
            print(f"Background subtracted data pixel count: {np.sum(data)}")

        # axes[1][1].imshow(data, origin='lower', cmap='gray')
        # axes[1][1].set_title('Background Subtracted Data')

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

    # plt.figure(figsize=(8, 8))
    # plt.imshow(source_data, origin='lower', cmap='gray')
    # plt.show()

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
        'Centroid': (M1_x, M1_y),
        'FWHM': average_FWHM
    }

    return shape

def cutout(data, focus_star, plot, index, verbose=False):
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

    if not verbose:
        # fig, axes = plt.subplots(1, 2, figsize=(12, 8))
        ax = plot.get_center_axis_by_index(index)

        if ax is None:
            print(f"Index {index} is out of bounds for the center axes grid.")
            return fit

        ax.clear()
        ax.imshow(cutout_data, origin='lower', cmap='viridis')
        
        # Add red rectangle around the cutout region
        # from matplotlib.patches import Rectangle
        # rect = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
        #                 linewidth=2, edgecolor='red', facecolor='none')
        # ax.add_patch(rect)


        # Mark center in cutout coordinates
        cutout_center_x = focus_star['Centroid'][0] - x_min
        cutout_center_y = focus_star['Centroid'][1] - y_min
        ax.plot(cutout_center_x, cutout_center_y, 'r+', markersize=15, markeredgewidth=3)

        ax.set_title(f'Cutout ({cutout_data.shape[1]}×{cutout_data.shape[0]} pixels)')
        ax.set_xticks([])
        ax.set_yticks([])

    return fit


def evaluate_sources(data, sources, focus_x, focus_y, verbose=False):
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

    if focus_x is not None and focus_y is not None:
        for source in results:
            centroid_x, centroid_y = source['Centroid']
            dist = np.sqrt((centroid_x - focus_x)**2 + (centroid_y - focus_y)**2)

            if dist < min_dist:
                min_dist = dist
                focus_star = source
    else:
        focus_star = max(results, key=lambda x: x['M0'])


    return focus_star


def photometry(fits_file, plot, index, focus_x=None, focus_y=None, verbose=False):

    if fits_file:
        hdu = fits.open(fits_file)
        data = hdu[0].data
        norm = ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
        plot.ax_left.imshow(data, cmap='viridis', origin='lower', norm=norm)
        plot.ax_left.set_xticks([])
        plot.ax_left.set_yticks([])
        plot.ax_left.set_title('Full Image')
    else:
        raise ValueError("No FITS file provided. Please specify a file.")

    data, sources = find_sources(data)
    focus_star = evaluate_sources(data, sources, focus_x=focus_x, focus_y=focus_y, verbose=verbose)

    fit = cutout(data, focus_star, plot, index, verbose=verbose)

    if focus_star is not None:
        print(f"Using source with label {focus_star['Label']} at ({focus_star['Centroid'][0]}, {focus_star['Centroid'][1]}) with FWHM: {focus_star['FWHM']} and fit_fwhm: {fit}\n")
        return focus_star['FWHM']
            

# fits_file = "data/nickel/focus9_6.fits"
# # fits_file = "../obs_images/raw/d1110.fits"

# photometry(fits_file, verbose=True)





# continue with detecting bad sources
# start with fine grain focus
# go through sequence and find best focus
    # do that with the brightest source and specify multiple sources and see if they all converge to similar values
    # also do that with fit_fwhm and what it converges to
# we also want to work on visual representations: box showing which star is being used, plot of FWHM vs focus value

# after that we can work on coarse grain focus
# figure out a way to measure sources when way out of focus

# field fine
# 1109    355
# 1110    357       # if brightest source is used d1110 chooses different source at (69, 986) with flux 70,219 vs usual source (486, 590) with flux 69,904
# 1111    359
# 1112    361
# 1113    363
# 1114    365
# 1115    367
# 1116    369
# 1117    371
# 1118    373
# 1119    375



def focus_finder(obs_start, obs_end, initial_focus, step_size, focus_x=None, focus_y=None, verbose=False):

    plot = Grid(num_rows=4, num_cols=3)

    images = []
    for obs_index, obs in enumerate(range(obs_start, obs_end + 1)):
        filename = f"../obs_images/raw/d{obs}.fits"
        print(f"Processing {filename}...")

        fwhm = photometry(filename, plot=plot, index=obs_index, focus_x=focus_x, focus_y=focus_y, verbose=verbose)

        attr = {
            'obs': obs,
            'fwhm': fwhm,
            'focus_value': initial_focus + (obs - obs_start) * step_size
        }
        images.append(attr)

    curve = sorted(images, key=lambda img: img['focus_value'])
    x_values = []
    y_values = []
    print("Curve found with the following focus values:")
    for img in curve:
        print(f"OBS: {img['obs']}, Focus: {img['focus_value']}, FWHM: {img['fwhm']}")
        x_values.append(img['focus_value'])
        y_values.append(img['fwhm'])

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)
    print(f"\nOptimal focus: {x_vertex}, FWHM: {y_vertex}")

    x_min, x_max = min(x_values), max(x_values)
    x_smooth = np.linspace(x_min, x_max, 50)
    y_smooth = a * x_smooth**2 + b * x_smooth + c

    plot.ax_right.scatter(x_values, y_values, label='Measured FWHM', color='blue')
    plot.ax_right.plot(x_smooth, y_smooth, 'r-', label='Fitted Quadratic')
    plot.ax_right.scatter([x_vertex], [y_vertex], color='green', label='Optimal focus', zorder=3)
    plot.ax_right.axvline(x=x_vertex, color='green', linestyle='--')

    plot.ax_right.set_xlabel('Focus Value', fontsize=12)
    plot.ax_right.set_ylabel('FWHM (pixels)', fontsize=12)
    plot.ax_right.set_title('Focus Curve Analysis', fontsize=14, fontweight='bold')
    plot.ax_right.legend(loc='best')
    plot.ax_right.grid(True, alpha=0.3)

    x_range = max(x_values) - min(x_values)
    y_range = max(y_values) - min(y_values)
    x_padding = x_range * 0.1
    y_padding = y_range * 0.1
    plot.ax_right.set_xlim(min(x_values) - x_padding, max(x_values) + x_padding)
    plot.ax_right.set_ylim(min(y_values) - y_padding, max(y_values) + y_padding)
    aspect_ratio = (x_range + 2*x_padding) / (y_range + 2*y_padding)
    plot.ax_right.set_aspect(aspect_ratio, adjustable='box')

    plt.show()

    return x_values, y_values


x_values, y_values = focus_finder(1109, 1119, 355, 2, verbose=False)

filename = 'focus_data.txt'
with open(filename, 'w') as f:
        # Write x_values on first line
        f.write(' '.join(map(str, x_values)) + '\n')
        # Write y_values on second line
        f.write(' '.join(map(str, y_values)) + '\n')



# Fine Field
# (485, 585): Optimal focus: 364.14860580366167, FWHM: 4.903481438876952
# (70, 980): Optimal focus: 364.0175987649038, FWHM: 4.86989269634887
# (873, 800): Optimal focus: 363.82007058960636, FWHM: 4.8494302870685715
# (734, 322): Optimal focus: 363.2149298264376, FWHM: 4.590148570307974
# (200, 345): Optimal focus: 363.7223486423801, FWHM: 4.760352676480579
    
# small focus stars
# (530, 396): Optimal focus: 362.60179806239313, FWHM: 4.147551719544026
# (402, 565): Optimal focus: 357.7769910685564, FWHM: 4.208899930379744
# (701, 490): d1110 is none 7:1 ratio, but d1111 is nearly 1:1


# Fine vband
# 1077    355     12.0
# 1078    357
# 1079    359     10.1
# 1080    361      8.4
# 1081    363      8.5
# 1082    365      8.7
# 1083    367
# 1084    369
# 1085    371
# 1086    373


# fine rband
# 1033    368      9.4
# 1034    366      9.9
# 1035    364      8.1
# 1036    362      7.6
# 1037    360      8.8
# 1038    358     17.4
# 1039    356      9.2
# 1040    354     10.3
# 1041    352     11.3
# 1042    350     12.9




# output a file with the focus values and fwhm to be used for refitting 
    # give the user the option to ommit bad data points
# iterative fitting ommiting outlier points (e.g. 3 sigma from the fit) (don't do this yet)

# restructure the image output
    # one image on the left showing the full image with the cutout box around the source
        # if the centroid is outside the cutout region, show the outlier centroid with its id
    # grid of images in the center showing the cutouts of each source at different focus values
    # one image on the right showing the plot of the fwhm vs focus value

    # make it so that as the program runs and points are added, the plot is updated


