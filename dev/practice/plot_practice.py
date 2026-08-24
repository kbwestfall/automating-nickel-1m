from matplotlib import pyplot as plt
from matplotlib import gridspec

import numpy as np
from astropy.io import fits
import quadratic
from astropy.visualization import (ImageNormalize, ZScaleInterval, LinearStretch)


# fig = plt.figure(figsize=(15, 10))
# num_exposures = 9
# plt.subplot2grid((num_exposures+2, num_exposures), (0, 0), rowspan=num_exposures)


# def identify_axes(ax_dict, fontsize=48):
#     """
#     Helper to identify the Axes in the examples below.

#     Draws the label in a large font in the center of the Axes.

#     Parameters
#     ----------
#     ax_dict : dict[str, Axes]
#         Mapping between the title / label and the Axes.
#     fontsize : int, optional
#         How big the label should be.
#     """
#     kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
#     for k, ax in ax_dict.items():
#         ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)

# mosaic = "ABCDK;AEFGK;AHIJK"

# fig = plt.figure(figsize=(15, 10))
# ax_dict = fig.subplot_mosaic(mosaic)
# identify_axes(ax_dict)


def create_focus_layout_constrained(num_rows=3, num_cols=3):
    """
    Create layout using constrained layout for better automatic sizing
    """
    
    fig = plt.figure(figsize=(20, 10), constrained_layout=True)
    
    # Create subfigures: left, center, right
    subfigs = fig.subfigures(1, 3, width_ratios=[1, 2, 1])
    
    # Left panel
    ax_left = subfigs[0].subplots(1, 1)
    ax_left.set_aspect('equal', adjustable='box')
    # ax_left.text(0.5, 0.5, 'Full Image\nwith Cutout Box', ha='center', va='center', fontsize=14)
    ax_left.set_title('Source Location')
    
    # Center grid
    center_axes = subfigs[1].subplots(num_rows, num_cols)
    if num_rows == 1 and num_cols == 1:
        center_axes = [center_axes]
    elif num_rows == 1:
        pass
    elif num_cols == 1:
        center_axes = [[ax] for ax in center_axes]
    
    # Format center subplots
    for i in range(num_rows):
        for j in range(num_cols):
            if num_rows == 1 and num_cols == 1:
                ax = center_axes[0]
            elif num_rows == 1:
                ax = center_axes[j]
            elif num_cols == 1:
                ax = center_axes[i][0]
            else:
                ax = center_axes[i][j]
            
            ax.set_aspect('equal')
            ax.text(0.5, 0.5, f'Cutout\n{i*num_cols + j + 1}', ha='center', va='center', fontsize=10)
            
    
    # Right panel
    ax_right = subfigs[2].subplots(1, 1)
    ax_right.set_aspect('equal', adjustable='box')
    # ax_right.text(0.5, 0.5, 'FWHM vs\nFocus Value', ha='center', va='center', fontsize=14)
    ax_right.set_title('Focus Curve')
    
    # Add title
    fig.suptitle(f'Focus Sequence Analysis ({num_rows}×{num_cols})', 
                fontsize=16, fontweight='bold')
    
    return fig, ax_left, center_axes, ax_right

fig, ax_left, center_axes, ax_right = create_focus_layout_constrained(num_rows=4, num_cols=3)


fits_file = '../obs_images/raw/d1091.fits'
hdu = fits.open(fits_file)
data = hdu[0].data


norm = ImageNormalize(data, interval=ZScaleInterval(), stretch=LinearStretch())
ax_left.imshow(data, cmap='viridis', origin='lower', norm=norm)
ax_left.set_xticks([])
ax_left.set_yticks([])
ax_left.set_title('Full Image')


filename = 'focus_data.txt'
with open(filename, 'r') as f:
    lines = f.readlines()
    
    x_values = np.array([float(x) for x in lines[0].strip().split()])
    y_values = np.array([float(y) for y in lines[1].strip().split()])

a, b, c = quadratic.fit_quadratic(x_values, y_values)
x_vertex, y_vertex = quadratic.vertex(a, b, c)
print(f"\nOptimal focus: {x_vertex}, FWHM: {y_vertex}")

x_min, x_max = min(x_values), max(x_values)
x_smooth = np.linspace(x_min, x_max, 50)
y_smooth = a * x_smooth**2 + b * x_smooth + c

ax_right.scatter(x_values, y_values, label='Measured FWHM', color='blue')
ax_right.plot(x_smooth, y_smooth, 'r-')
ax_right.scatter([x_vertex], [y_vertex], color='green', label='Optimal focus', zorder=3)
ax_right.axvline(x=x_vertex, color='green', linestyle='--')

ax_right.set_xlabel('Focus Value', fontsize=12)
ax_right.set_ylabel('FWHM (pixels)', fontsize=12)
ax_right.set_title('Focus Curve Analysis', fontsize=14, fontweight='bold')
ax_right.legend(loc='best')
ax_right.grid(True, alpha=0.3)

# Calculate aspect ratio to make plot square but with meaningful scaling
x_range = max(x_values) - min(x_values)
y_range = max(y_values) - min(y_values)

# Add some padding to the ranges
x_padding = x_range * 0.1
y_padding = y_range * 0.1

ax_right.set_xlim(min(x_values) - x_padding, max(x_values) + x_padding)
ax_right.set_ylim(min(y_values) - y_padding, max(y_values) + y_padding)

# Set aspect ratio to make it square
aspect_ratio = (x_range + 2*x_padding) / (y_range + 2*y_padding)
ax_right.set_aspect(aspect_ratio, adjustable='box')
    


plt.show()