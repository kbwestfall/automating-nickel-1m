#! @KPYTHON@
import ktl

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time 
from astropy import units as u
import datetime
import numpy as np

import warnings
warnings.filterwarnings("ignore", message="Input line .* contained no data")

def wait_until(self, keyword, expected_value, timeout=15):
        start_time = time.time()
        while time.time() - start_time < timeout:
            if keyword.read() in expected_value:
                return True
        return False

def find_nearest_star(telescope, star_type, file):

    if star_type == 'focus':
        print("Finding nearest focus star...")
        star_type = 'Pointing'
    elif star_type == 'pointing':
        print("Finding nearest pointing star...")
        star_type = 'Focusing'
    else:
        raise ValueError("Invalid star type. Choose either 'focus' or 'pointing'.")
    

    min_separation = float('inf')
    nearest_star = None

    focus_stars = np.loadtxt(file, dtype=str, usecols=(0, 1, 2), comments=star_type)

    focus_coords = SkyCoord(ra=focus_stars[:, 1], dec=focus_stars[:, 2], unit=('hourangle', 'deg'))

    separations = telescope.separation(focus_coords)

    min_index = np.argmin(separations)
    nearest_star = focus_stars[min_index]

    target_name = nearest_star[0]
    target_ra = nearest_star[1]
    target_dec = nearest_star[2]

    print(f'Target Name: {target_name}')
    print(f'Target RA: {target_ra}')
    print(f'Target DEC: {target_dec}')

    return target_name, target_ra, target_dec


def move_to_target(target_coords):
    stop_key = ktl.cache('nickelpoco', 'POCSTOP')
    target_key = ktl.cache('nickelpoco', 'POCOT')
    track_key = ktl.cache('nickelpoco', 'POCTRCK')
    ra_desired = ktl.cache('nickelpoco', 'POCRAD')
    dec_desired = ktl.cache('nickelpoco', 'POCDECD')

    if stop_key.read() != 0:
        print("POCSTOP is 'disabled'. Waiting for 'allowed' to allow motion")
    if not wait_until(stop_key, 'allowed', timeout=30):
        raise Exception("POCSTOP is 'disabled'. Set to 'allowed' to allow motion")

    if track_key.read() == 'off':
        raise Exception("Tracking is not enabled. Please enable tracking before moving to target.")
        # 0->off, 1->on

    if target_key.read() != '0':
        print("Telescope not ready to move to target. Waiting for ready (POCOT to be 0)")
        # 1->have not reached target, 0-> on target and stable, -1-> failed to reach target within time limit and not on target.
    if not wait_until(target_key, '0', timeout=30):
        raise Exception("Telescope is not ready to move to target. Please wait until the telescope is ready.")

    ra_desired.write(target_coords.ra.to('hourangle').value)
    dec_desired.write(target_coords.dec.to('deg').value)

    if wait_until(target_key, '0', timeout=60):
        print("Telescope is now on target.")
    else:
        if target_key.read() == '-1':
            raise Exception("Failed to reach target within time limit. Telescope is not on target.")
        else:
            raise Exception("Telescope is not on target. Please check the status.")

# POCOBJNA = Pointing16

import argparse

def main():
    parser = argparse.ArgumentParser(description='Move telescope to target')
    parser.add_argument('--dry_run', action='store_true', help='Run in dry run mode (no actual movement)')
    parser.add_argument('--star_type', choices=['focus', 'pointing'], default='focus', help='Type of star to use for target selection (focus or pointing)')
    args = parser.parse_args()

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

    target_name, target_ra, target_dec = find_nearest_star(telescope_coord, args.star_type, 'point_focus.txt')

    target_coords = SkyCoord(ra=target_ra, dec=target_dec, unit=('hourangle', 'deg'))

    if not args.dry_run:
        move_to_target(target_coords)

if __name__ == '__main__':
    main()

# note that if first time acquiring target and not idling on pocot=0, then it will timeout waiting for pocot=0