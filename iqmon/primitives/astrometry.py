from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
from astropy import units as u
from astropy import stats
from astropy.wcs import WCS
from astropy.time import Time
from astropy.table import Table, Column


from keckdrpframework.primitives.base_primitive import BasePrimitive

from .utils import pre_condition, post_condition



##-----------------------------------------------------------------------------
## Primitive: SolveAstrometry
##-----------------------------------------------------------------------------
class SolveAstrometry(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'Skip image is not set',
                                not self.action.args.skip),
                  pre_condition(self, 'WCS not already solved',
                                self.action.args.wcs is None),
                 ]
        force_solve = self.cfg['Astrometry'].getboolean('force_solve', False)
        return np.all(checks) if force_solve is False else True

    def _post_condition(self):
        """
        Check for conditions necessary to verify that the process ran
        correctly.
        """
        checks = []
        return np.all(checks)

    def _perform(self):
        """
        Returns an Argument() with the parameters that depend on this
        operation.
        """
        self.log.info(f"Running {self.__class__.__name__} action")

        nx, ny = self.action.args.ccddata.shape
        estimated_pixel_scale = self.cfg['Telescope'].getfloat('pixel_scale', 1)
        search_radius = self.cfg['Astrometry'].getfloat('search_radius', 1)
        solve_field = self.cfg['Astrometry'].get('solve_field', '/usr/local/bin/solve-field')
        wcs_output_file = Path('~/tmp.wcs').expanduser()
        axy_output_file = Path('~/tmp.axy').expanduser()
        solvetimeout = self.cfg['Astrometry'].getint('solve_timeout', 120)
        cmd = [f'{solve_field}', '-p', '-O', '-N', 'none', '-B', 'none',
               '-U', 'none', '-S', 'none', '-M', 'none', '-R', 'none',
               '--axy', f'{axy_output_file}', '-W', f'{wcs_output_file}',
               '-z', '4',
               '-L', f'{0.9*estimated_pixel_scale}',
               '-H', f'{1.1*estimated_pixel_scale}',
               '-u', 'arcsecperpix',
               '-t', f"{self.cfg['Astrometry'].getfloat('tweak_order', 2)}",
               '-3', f'{self.action.args.header_pointing.ra.deg}',
               '-4', f'{self.action.args.header_pointing.dec.deg}',
               '-5', f'{search_radius}',
               '-l', f'{solvetimeout}',
               ]
        if self.cfg['Astrometry'].get('astrometry_cfg_file', None) is not None:
            cmd.extend(['-b', self.cfg['Astrometry'].get('astrometry_cfg_file')])

        if self.action.args.fitsfilepath.suffix == '.fz':
            self.log.info(f'  Uncompressing {self.action.args.fitsfilepath.name}')
            tmpfile = Path('~/tmp.fts').expanduser().absolute()
            if tmpfile.exists():
                self.log.debug(f'  Removing old {tmpfile}')
                tmpfile.unlink()
            funpackcmd = ['funpack', '-O', f'{tmpfile}',
                          f'{self.action.args.fitsfilepath}']
            self.log.debug(f'  funpack command: {funpackcmd}')
            funpack_result = subprocess.run(funpackcmd,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
            cmd.append(f'{tmpfile}')
        else:
            cmd.append(f'{self.action.args.fitsfilepath}')

        self.log.debug(f"  Solve astrometry command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, timeout=solvetimeout,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        except subprocess.TimeoutExpired as e:
            self.log.warning('The solve-field process timed out')
            return self.action.args

        self.log.debug(f"  Returncode = {result.returncode}")
        for line in result.stdout.decode().split('\n'):
            self.log.debug(line)
        if result.returncode != 0:
            self.log.warning(f'Astrometry solve failed:')
            for line in result.stdout.decode().split('\n'):
                self.log.warning(line)
            for line in result.stderr.decode().split('\n'):
                self.log.warning(line)
            return self.action.args

        if wcs_output_file.exists() is True:
            self.log.debug(f"  Found {wcs_output_file}")
        else:
            self.log.warning(f"Could not find {wcs_output_file}")
            return self.action.args
        if axy_output_file.exists() is True:
            self.log.debug(f"  Found {axy_output_file}. Deleteing.")
            axy_output_file.unlink()
        # Open wcs output
        self.log.debug(f"  Creating astropy.wcs.WCS object")
        output_wcs = WCS(f'{wcs_output_file}')
        self.log.debug(f"Deleteing {wcs_output_file}")
        wcs_output_file.unlink()
        if output_wcs.is_celestial is False:
            self.log.info(f"  Could not parse resulting WCS as celestial")
            return self.action.args

        self.log.info(f"  Solve complete")
        self.action.args.meta['wcs_header'] = output_wcs.to_header_string()
        self.action.args.wcs = output_wcs

        # Determine Pointing Error
        r, d = self.action.args.wcs.all_pix2world([nx/2.], [ny/2.], 1)
        self.action.args.wcs_pointing = c.SkyCoord(r[0], d[0], frame='fk5',
                                          equinox='J2000',
                                          unit=(u.deg, u.deg),
                                          obstime=self.action.args.kd.obstime())
        self.action.args.meta['perr'] = self.action.args.wcs_pointing.separation(
                                        self.action.args.header_pointing).to(u.arcmin).value
        self.log.info(f'Pointing error = {self.action.args.meta.get("perr"):.1f}')

        return self.action.args

