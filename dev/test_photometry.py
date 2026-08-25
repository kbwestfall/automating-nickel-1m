
from pathlib import Path

from IPython import embed
from astropy.io import fits

import photometry

def test_find_sources():

    root = Path('/Users/westfall/Work/Lick/Nickel/science_camera/obs/2025-08-14')
    tst_img = root / 'n5' / 'n5053.fits'

    with fits.open(tst_img) as hdu:
        data = hdu[0].data.astype(float)

    bkg, src_img = photometry.find_sources(data, verbose=True)

def test_image_quality():

    root = Path('/Users/westfall/Work/Lick/Nickel/science_camera/obs/2025-08-14')
    tst_img = root / 'n5' / 'n5053.fits'

    data, bkg, sources, img_quality \
        = photometry.image_quality(tst_img, method=(516, 520), verbose=True)

    data, bkg, sources, img_quality \
        = photometry.image_quality(tst_img, method='brightest', verbose=True)

    data, bkg, sources, img_quality \
        = photometry.image_quality(tst_img, method='weighted', verbose=True)

    embed()
    exit()

test_image_quality()


