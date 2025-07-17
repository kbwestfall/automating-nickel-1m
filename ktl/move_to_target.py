#! @KPYTHON@
import ktl

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time 
from astropy import units as u
import datetime
import numpy as np

lick = EarthLocation.of_site('Lick Observatory')
print(lick)

t = Time(datetime.datetime.now(datetime.UTC), location=lick)
print(f'Current time: {t}')

lst = t.sidereal_time('mean')

ra_key = ktl.cache('nickelpoco', 'POCORAA')
dec_key = ktl.cache('nickelpoco', 'POCODECA')
lick_ra = ra_key.read()
lick_dec = dec_key.read()
print(f'Lick RA: {lick_ra}')
print(f'Lick DEC: {lick_dec}')

def find_nearest_focus_star(telescope_ra, telescope_dec, file):
    telescope_coord = SkyCoord(ra=telescope_ra, dec=telescope_dec, unit=('hourangle', 'deg'))

    min_separation = float('inf')
    nearest_star = None

    focus_stars = np.loadtxt(file, dtype=str, usecols=(0, 1, 2), comments='Pointing')

    focus_coords = SkyCoord(ra=focus_stars[:, 1], dec=focus_stars[:, 2], unit=('hourangle', 'deg'))

    separations = telescope_coord.separation(focus_coords)

    min_index = np.argmin(separations)
    nearest_star = focus_stars[min_index]

    return nearest_star[0], nearest_star[1], nearest_star[2]

    # with open(file, 'r') as f:
    #     for line in f:
    #         line = line.strip()
    #         if line and line.startswith('Focusing'):
    #             parts = line.split()
    #             star_name = parts[0]
    #             star_ra = parts[1]
    #             star_dec = parts[2]
    #             star_coord = SkyCoord(ra=star_ra, dec=star_dec, unit=('hourangle', 'deg'))

    #             separation = telescope_coord.separation(star_coord).deg
    #             if separation < min_separation:
    #                 min_separation = separation
    #                 nearest_star = {
    #                     'name': star_name,
    #                     'ra': star_ra,
    #                     'dec': star_dec,
    #                 }

    # return nearest_star['name'], nearest_star['ra'], nearest_star['dec']


# Focusing00   00:30:48.780  +28:06:42.42   2000.0   f
# target_ra = '00:30:48.780'
# target_dec = '+28:06:42.42'

target_name, target_ra, target_dec = find_nearest_focus_star(lick_ra, lick_dec, 'point_focus.txt')
print(f'Target Name: {target_name}')
print(f'Target RA: {target_ra}')
print(f'Target DEC: {target_dec}')

target_coords = SkyCoord(ra=target_ra, dec=target_dec, unit=('hourangle', 'deg'))

ha_key = ktl.cache('nickelpoco', 'POCOHAA')
print(f'Lick HA: {ha_key.read()}')

ha = lst - target_coords.ra
print(f'Target HA: {ha.to_string(unit=u.hourangle, sep=":")}')


