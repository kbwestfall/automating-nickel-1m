
from importlib import resources
from pathlib import Path
import warnings

from astropy.coordinates import SkyCoord
import numpy as np

try:
    import ktl
except ModuleNotFoundError:
    warnings.warn('Unable to import ktl package.  Functionality will be limited.')
    ktl = None

from nickel_focus import log
from nickel_focus import starlist


class NickelTelescopePointing:
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

    @property
    def current(self):
        return SkyCoord(ra=self.raa.read(), dec=self.deca.read(), unit=('hourangle', 'deg'))

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
    tbl = starlist.parse_starlist(_file)

    if obj_search_str is not None:
        selected = [obj_search_str in obj_name for obj_name in tbl['name']]
        if np.sum(selected) == 0:
            raise ValueError(f'{_file} has no object names containing {obj_search_str}')
        tbl = tbl[selected]

    obj_coo = SkyCoord(ra=tbl['ra'], dec=tbl['dec'], unit=('hourangle', 'deg'))
    offsets = telescope_coo.separation(obj_coo)
    nearest_obj = tbl[np.argmin(offsets)]
    return nearest_obj['name'][0], nearest_obj['ra'][0], nearest_obj['dec'][0]
