import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import simple_norm
from photutils.datasets import (make_model_image, make_model_params,
                                make_noise_image, make_random_models_table)
from astropy.modeling.models import Moffat2D



def gamma_to_fwhm(gamma, alpha):
    return 2 * gamma * (2 ** (1 / alpha) - 1) ** 0.5

def fwhm_to_gamma(fwhm, alpha):
    return (fwhm / 2) * (1 / ((2 ** (1 / alpha) - 1) ** 0.5))

def simulation(filepath, model=Moffat2D(), n_sources=1, shape=(100,100), alpha=3.5, fwhm=6.0):
    gamma = fwhm_to_gamma(fwhm, alpha)
    param_ranges = {'amplitude': [100, 500],
                    'x_0': [50, 50],
                    'y_0': [50, 50],
                    'gamma': [gamma, gamma],
                    'alpha': [alpha, alpha]}
    params = make_random_models_table(n_sources, param_ranges, seed=123)

    model_shape = (int(fwhm*3), int(fwhm*3))
    clean_data = make_model_image(shape, model, params, model_shape=model_shape)

    background_stddev = 10
    background_mean = 0
    noise = make_noise_image(shape, distribution='gaussian', mean=background_mean, stddev=background_stddev, seed=123)
    noise_data = clean_data + noise

    hdu = fits.PrimaryHDU(noise_data)
    file_name = str(fwhm)
    file_name = file_name.replace('.', '_')
    hdu.writeto(f'{filepath}.fits', overwrite=True)
    plt.figure(figsize=(8, 8))
    plt.imshow(noise_data, origin='lower', cmap='viridis')
    plt.savefig(f'{filepath}.png')


# for fwhm in [6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]:
#     simulation(fwhm=fwhm)

# simulation(fwhm=6.0)

# plt.imshow(noise_data, origin='lower')
# plt.tight_layout()
# plt.savefig('psf-practice/simulated_image.png')

# from astropy.io import fits
# hdu = fits.PrimaryHDU(noise_data)
# hdu.writeto('psf-practice/images/simulated_image.fits', overwrite=True)