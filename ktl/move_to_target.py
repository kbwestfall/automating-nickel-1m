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

telescope_coord = SkyCoord(ra=lick_ra, dec=lick_dec, unit=('hourangle', 'deg'))

def find_nearest_focus_star(telescope, file):
    

    min_separation = float('inf')
    nearest_star = None

    focus_stars = np.loadtxt(file, dtype=str, usecols=(0, 1, 2), comments='Pointing')

    focus_coords = SkyCoord(ra=focus_stars[:, 1], dec=focus_stars[:, 2], unit=('hourangle', 'deg'))

    separations = telescope.separation(focus_coords)

    min_index = np.argmin(separations)
    nearest_star = focus_stars[min_index]

    return nearest_star[0], nearest_star[1], nearest_star[2]

target_name, target_ra, target_dec = find_nearest_focus_star(telescope_coord, 'point_focus.txt')
print(f'Target Name: {target_name}')
print(f'Target RA: {target_ra}')
print(f'Target DEC: {target_dec}')

target_coords = SkyCoord(ra=target_ra, dec=target_dec, unit=('hourangle', 'deg'))
print(f'Target Coordinates: {target_coords}')
print(f'Target Coordinates (RA, DEC): {target_coords.ra}, {target_coords.dec}')



lick_ha = lst - telescope.ra
ha_key = ktl.cache('nickelpoco', 'POCOHAA')
print(f'Lick Read HA: {ha_key.read()}')
print(f'Lick HA: {lick_ha.to_string(unit=u.hourangle, sep=":")}')

ha = lst - target_coords.ra
print(f'Target HA: {ha.to_string(unit=u.hourangle, sep=":")}')



track_key = ktl.cache('nickelpoco', 'POCTRCK')
ra_desired = ktl.cache('nickelpoco', 'POCRAD')
dec_desired = ktl.cache('nickelpoco', 'POCDECD')
target_key = ktl.cache('nickelpoco', 'POCOT')

if track_key.read() == 'off':
    raise Exception("Tracking is not enabled. Please enable tracking before moving to target.")

if target_key.read() != '0':
    print("Telescope not ready to move to target. Waiting for ready (POCOT to be 0)")
if not target_key.waitFor('== 0', timeout=30):
    raise Exception("Telescope is not ready to move to target. Please wait until the telescope is ready.")

ra_desired.write(target_coords.ra.to('hourangle').value)
dec_desired.write(target_coords.dec.to('deg').value)


