import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.psf import fit_fwhm
from photutils.detection import DAOStarFinder
from photutils.centroids import centroid_1dg, centroid_2dg, centroid_com, centroid_quadratic
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint


import argparse

def main():
    parser = argparse.ArgumentParser(description='Photometry of a FITS file.')
    parser.add_argument('-f', '--file', type=str, help='Path to the FITS file')
    parser.add_argument('--fwhm', type=float, default=10, help='Estimated FWHW)')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Enable verbose output')
    args = parser.parse_args()

    if args.verbose:
        plt.ion()

    photometry(args.file, args.fwhm, args.verbose)

    if args.verbose:
        plt.ioff()
        input("Press Enter to exit...")

def photometry(fits_file, est_fwhm=10, verbose=False):

    hdu = fits.open(fits_file)
    data = hdu[0].data
    if verbose:
        print(hdu.info(), "\n")

    if verbose:
        max_index = np.argmax(data)
        y_max, x_max = np.unravel_index(max_index, data.shape)
        print(f"Brightest pixel at position: x={x_max}, y={y_max}, value={data[y_max, x_max]}")
        
    sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
    threshold = detect_threshold(data, nsigma=5.0, sigma_clip=sigma_clip)
    segment_img = detect_sources(data, threshold, npixels=10)
    footprint = circular_footprint(radius=int(est_fwhm))
    mask = segment_img.make_source_mask(footprint=footprint)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)

    if verbose:
        print(f"Mean: {mean}, Median: {median}, Std: {std}")

    # background = data.copy()
    # background[mask] = median

    # if verbose:
    #     plt.figure(figsize=(15, 8))
    #     plt.subplot(1, 3, 1)
    #     plt.imshow(data, origin='lower', cmap='viridis')
    #     plt.title('Original Image')

    #     plt.subplot(1, 3, 2)
    #     plt.imshow(background, origin='lower', cmap='viridis')
    #     plt.title('Noise Image')

    #     plt.subplot(1, 3, 3)
    #     plt.imshow(data - median, origin='lower', cmap='viridis')
    #     plt.title('Removed Noise Image')

    #     plt.tight_layout()
    #     plt.draw()
    #     # plt.savefig('psf-practice/noise.png')

    data = data - median

    sources = DAOStarFinder(100, est_fwhm).find_stars(data)
    sources.sort('flux', reverse=True)
    focus_star = sources[0]
    if verbose: 
        print(focus_star, "\n")

    region_size = est_fwhm * 3
    y_min = int(max(0, focus_star['ycentroid'] - region_size//2))
    y_max = int(min(data.shape[0], focus_star['ycentroid'] + region_size//2))
    x_min = int(max(0, focus_star['xcentroid'] - region_size//2))
    x_max = int(min(data.shape[1], focus_star['xcentroid'] + region_size//2))
    data = data[y_min:y_max, x_min:x_max]

    xcent, ycent = centroid_2dg(data)

    if verbose:
        plt.figure(figsize=(8, 8))
        plt.imshow(data, origin='lower', cmap='grey')
        plt.plot(xcent, ycent, 'r+', markersize=10)
        plt.annotate('Centroid', (xcent, ycent), xytext=(5, 5), textcoords='offset points', color='white')
        plt.draw()  
        # plt.savefig('psf-practice/cutout.png')


    fit_x = data.shape[0]
    fit_y = data.shape[1]
    if fit_x%2 == 0:
        fit_x -= 1
    if fit_y%2 == 0:
        fit_y -= 1
    fit_shape = (fit_x, fit_y)

    if verbose:
        est_fwhm = fit_fwhm(data, fit_shape=fit_shape)
        # est_fwhm.sort()
        # est_fwhm = est_fwhm[-1]
        print(f"Estimated FWHM: {est_fwhm}")

    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])

    meshgrid_x, meshgrid_y = np.meshgrid(x, y)

    M1_x = np.sum(meshgrid_x * data) / np.sum(data)
    M1_y = np.sum(meshgrid_y * data) / np.sum(data)

    M2_x = np.sum((meshgrid_x)** 2 * data) / np.sum(data)
    M2_y = np.sum((meshgrid_y)** 2 * data) / np.sum(data)

    sigma_x = np.sqrt(M2_x - M1_x**2)
    sigma_y = np.sqrt(M2_y - M1_y**2)

    FWHM_x = 2 * np.sqrt(2 * np.log(2)) * sigma_x
    FWHM_y = 2 * np.sqrt(2 * np.log(2)) * sigma_y
    average_FWHM = (FWHM_x + FWHM_y) / 2

    if verbose:
        print(f"FWHM_x: {FWHM_x}, FWHM_y: {FWHM_y}, Average FWHM: {average_FWHM}")

    return average_FWHM


# true_fwhm = []
# obs_fwhm = []
# ratio_fwhm = []

# for i in range(6, 13, 1):
#     true = float(i)
#     true_fwhm.append(true)
#     file_name = str(float(i))
#     file_name = file_name.replace('.', '_')
#     obs = photometry(f'psf-practice/images/focus_{file_name}.fits')
#     obs_fwhm.append(obs)
#     ratio_fwhm.append(obs/true)

# plt.figure(figsize=(8, 8))
# plt.plot(true_fwhm, ratio_fwhm, 'ro', label=f'True FWHM: {true_fwhm}')
# plt.xlabel('True FWHM')
# plt.ylabel('Observed FWHM')
# plt.title('True vs Observed FWHM')
# plt.savefig(f'psf-practice/images/true_vs_obs.png')

# photometry('psf-practice/images/focus_24.fits')

# photometry('psf-practice/images/n1043.fits')

if __name__ == "__main__":
    main()