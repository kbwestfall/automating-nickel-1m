import datetime
from importlib import resources
from pathlib import Path
import warnings

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time 
import numpy as np

try:
    import ktl
except ModuleNotFoundError:
    warnings.warn('Unable to import ktl package.  Functionality will be limited.')
    ktl = None

from nickel_focus import log


class TelescopePointing:
    """
    KTL object used to manipulate the telescope pointing.
    """
    def __init__(self):
        if ktl is None:
            raise RuntimeError(
                'The ktl package is not available.  Controlling the telescope pointing requires a '
                'kpython environment with ktl installed.'
            )

        # Controls overall movement of mechanisms
        self.pocstop = ktl.cache('nickelpoco', 'POCSTOP')
        # Indicates whether or not the telescope is tracking
        self.track = ktl.cache('nickelpoco', 'POCTRCK')
        # Indicates whether or not the telescope is ready to move to a target
        self.target = ktl.cache('nickelpoco', 'POCOT')
        # Provides the actual RA of the telescope pointing
        self.raa = ktl.cache('nickelpoco', 'POCORAA')
        # Sets the desired RA of the telescope pointing
        self.rad = ktl.cache('nickelpoco', 'POCORAD')
        # Provides the actual Dec of the telescope pointing
        self.deca = ktl.cache('nickelpoco', 'POCODECA') 
        # Sets the desired Dec of the telescope pointing
        self.decd = ktl.cache('nickelpoco', 'POCODECD') 

    def slew_to(self, ra, dec):

        if not self.pocstop.waitFor('== allowed', timeout=0.5):
            raise ValueError('Telescope movement is disabled!')

        if not self.track.waitFor('== off', timeout=0.5):
            raise ValueError('Tracking is disabled!')

        if not self.target.waitFor('== 0', timeout=0.5):
            raise ValueError('Telescope is not ready to move to a new target!')

        log.info(f'Moving to RA={ra}, Dec={dec}')
        target_coo = SkyCoord(ra=ra, dec=dec, unit=('hourangle', 'deg'))
        self.rad.write(target_coo.ra.to('hourangle').value)
        self.decd.write(target_coo.dec.to('deg').value)

        # TODO: Need a better way to determine if the move fails!  Timeout
        # currently set at 5 min.
        if not self.target.waitFor('== 0', timeout=300):
            # TODO: Is there a way to appropriately handle this?
            raise ValueError('Telescope failed to make it to target within 5 min.')


def find_nearest_target(telescope_coo, obj_search_str=None, file=None):

    _file = (
        resources.files('nickel_focus') / 'data' / 'point_focus.txt'
        if file is None else Path(file).absolute()
    )
    coo_tbl = np.loadtxt(_file, dtype=str)

    if obj_search_str is not None:
        selected = [obj_search_str in obj_name for obj_name in coo_tbl[:,0]]
        if np.sum(selected) == 0:
            raise ValueError(f'{_file} has no object names containing {obj_search_str}')
        coo_tbl = coo_tbl[selected]

    obj_coo = SkyCoord(ra=coo_tbl[:,1], dec=coo_tbl[:,2], unit=('hourangle', 'deg'))
    offsets = telescope_coo.separation(obj_coo)

    nearest_obj = coo_tbl[np.argmin(offsets)]
    log.info(f'Nearest object is: {nearest_obj[0]} at RA={nearest_obj[1]}, Dec={nearest_obj[2]}')
    return nearest_obj[0], nearest_obj[1], nearest_obj[2]


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

# note: if first time acquiring target and not idling on pocot=0, it will timeout waiting for pocot=0 unless the user moves to a target and is successful within the timeout period.