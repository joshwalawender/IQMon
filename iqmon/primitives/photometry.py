from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
from astropy import units as u
from astropy import stats
from astropy.time import Time
from astropy.table import Table, Column
import photutils
import sep


from keckdrpframework.primitives.base_primitive import BasePrimitive

from .utils import pre_condition, post_condition, mode


##-----------------------------------------------------------------------------
## Primitive: MakeSourceMask
##-----------------------------------------------------------------------------
class MakeSourceMask(BasePrimitive):
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
                  pre_condition(self, 'Extraction requested',
                                self.cfg['Extract'].getboolean('do_extraction', False) is True),
                 ]
        return np.all(checks)

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

        snr = self.cfg['Extract'].getfloat('source_mask_snr', 5)
        self.log.debug(f"  Using snr = {snr}")
        npixels = self.cfg['Extract'].getfloat('source_mask_npixels', 5)
        self.log.debug(f"  Using npixels = {npixels}")
        source_mask = photutils.make_source_mask(self.action.args.ccddata,
                                snr, npixels)
        self.action.args.source_mask = source_mask

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: CreateBackground
##-----------------------------------------------------------------------------
class CreateBackground(BasePrimitive):
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
                  pre_condition(self, 'Extraction requested',
                                self.cfg['Extract'].getboolean('do_extraction', False) is True),
                 ]
        return np.all(checks)

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

        box_size = self.cfg['Extract'].getint('background_box_size', 128)
        self.log.debug(f"  Using box size = {box_size} pixels")

        bkg = photutils.Background2D(self.action.args.ccddata,
                                     box_size=box_size,
                                     mask=self.action.args.source_mask,
                                     sigma_clip=stats.SigmaClip())
        self.action.args.background = bkg

        return self.action.args



##-----------------------------------------------------------------------------
## Primitive: ExtractStars
##-----------------------------------------------------------------------------
class ExtractStars(BasePrimitive):
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
                  pre_condition(self, 'Background has been generated',
                                self.action.args.background is not None),
                  pre_condition(self, 'Extraction requested',
                                self.cfg['Extract'].getboolean('do_extraction', False) is True),
                 ]
        return np.all(checks)

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

        exptime = self.action.args.meta.get('exptime')
        pixel_scale = self.cfg['Telescope'].getfloat('pixel_scale', 1)
        thresh = self.cfg['Extract'].getint('extract_threshold', 9)
        minarea = self.cfg['Extract'].getint('extract_minarea', 7)
        mina = self.cfg['Extract'].getfloat('fwhm_mina', 1)
        minb = self.cfg['Extract'].getfloat('fwhm_minb', 1)
        faint_limit_pct = self.cfg['Extract'].getfloat('faint_limit_percentile', 0)
        bright_limit_pct = self.cfg['Extract'].getfloat('bright_limit_percentile', 100)
        radius_limit = self.cfg['Extract'].getfloat('radius_limit_pix', 4000)

        bsub = self.action.args.ccddata.data - self.action.args.background.background
        seperr = self.action.args.ccddata.uncertainty.array
        sepmask = self.action.args.ccddata.mask

        # Define quick utility function
        def run_sep(bsub, seperr, sepmask, thresh, minarea):
            try:
                objects = sep.extract(bsub, err=seperr, mask=sepmask,
                                      thresh=float(thresh), minarea=minarea)
                return objects
            except Exception as e:
                if str(e)[:27] == 'internal pixel buffer full:':
                    return None
                else:
                    raise SEPError(str(e))

        objects = None
        while objects is None:
            try:
                self.log.info(f'Invoking SEP with threshold: {thresh}')
                objects = run_sep(bsub, seperr, sepmask, thresh, minarea)
                thresh += 9
            except SEPError as e:
                self.log.error('Source extractor failed:')
                self.log.error(e)
                return self.action.args

        t = Table(objects)
        t['flux'] /= exptime

        # Add radius (of star from center of image) to table
        ny, nx = bsub.shape
        r = np.sqrt((t['x']-nx/2.)**2 + (t['y']-ny/2.)**2)
        t.add_column(Column(data=r.data, name='r', dtype=np.float))

        # Add FWHM to table
        coef = 2*np.sqrt(2*np.log(2))
        fwhm = np.sqrt((coef*t['a'])**2 + (coef*t['b'])**2)
        t.add_column(Column(data=fwhm.data, name='FWHM', dtype=np.float))

        # Add ellipticities to table
        ellipticities = t['a']/t['b']
        t.add_column(Column(data=ellipticities.data, name='ellipticity', dtype=np.float))

        # Filter out stars based on bright and faint limits
        faint_limit = np.percentile(t['flux'], faint_limit_pct)
        bright_limit = np.percentile(t['flux'], bright_limit_pct)
        self.log.info(f'  Faintest {faint_limit_pct:.1f}% flux {faint_limit:f}')
        self.log.info(f'  Brightest {bright_limit_pct:.1f}% flux {bright_limit:f}')
        filtered = (t['a'] < mina) | (t['b'] < minb) | (t['flag'] > 0) | (t['flux'] > bright_limit) | (t['flux'] < faint_limit) | (t['r'] > radius_limit)
        self.log.debug(f'  Removing {np.sum(filtered):d}/{len(filtered):d}'\
                       f' extractions from FWHM calculation')
        self.log.debug(f"    {np.sum( (t['a'] < mina) )} removed for fwhm_mina limit")
        self.log.debug(f"    {np.sum( (t['b'] < minb) )} removed for fwhm_minb limit")
        self.log.debug(f"    {np.sum( (t['flag'] > 0) )} removed for source extractor flags")
        self.log.debug(f"    {np.sum( (t['flux'] < faint_limit) )} removed for faint limit")
        self.log.debug(f"    {np.sum( (t['flux'] > bright_limit) )} removed for bright limit")

        self.action.args.meta['n_objects'] = len(t[~filtered])
        self.log.info(f'  Found {self.action.args.meta.get("n_objects"):d} stars')

        self.action.args.objects = t[~filtered]
        self.action.args.objects.sort('flux')
        self.action.args.objects.reverse()

        if self.action.args.meta.get("n_objects") == 0:
            self.log.warning('No stars found')
            return self.action.args
        else:
            FWHM_pix = np.median(t['FWHM'][~filtered])
            FWHM_mode_bin = pixel_scale*0.25
            FWHM_pix_mode = mode(t['FWHM'][~filtered]/FWHM_mode_bin)*FWHM_mode_bin
            self.log.info(f'  Median FWHM = {FWHM_pix:.1f} pix ({FWHM_pix*pixel_scale:.2f} arcsec)')
            self.log.info(f'  Mode FWHM = {FWHM_pix_mode:.1f} pix ({FWHM_pix_mode*pixel_scale:.2f} arcsec)')
            ellipticity = np.median(t['ellipticity'][~filtered])
            ellipticity_mode_bin = 0.05
            ellipticity_mode = mode(t['ellipticity'][~filtered]/ellipticity_mode_bin)*ellipticity_mode_bin
            self.log.info(f'  Median ellipticity = {ellipticity:.2f}')
            self.log.info(f'  Mode ellipticity = {ellipticity_mode:.2f}')
            self.action.args.meta['fwhm'] = FWHM_pix_mode
            self.action.args.meta['ellipticity'] = ellipticity_mode

        ## Do photutils photometry measurement
        positions = [(det['x'], det['y']) for det in self.action.args.objects]
        ap_radius = self.cfg['Photometry'].getfloat('aperture_radius', 2)*FWHM_pix
        star_apertures = photutils.CircularAperture(positions, ap_radius)
        sky_apertures = photutils.CircularAnnulus(positions,
                                                  r_in=int(np.ceil(1.5*ap_radius)),
                                                  r_out=int(np.ceil(2.0*ap_radius)))
        phot_table = photutils.aperture_photometry(
                               self.action.args.ccddata,
                               [star_apertures, sky_apertures])
        phot_table['sky'] = phot_table['aperture_sum_1'] / sky_apertures.area()
        med_sky = np.median(phot_table['sky'])
        self.log.info(f'  Median Sky = {med_sky.value:.0f} e-/pix')
        self.action.args.meta['sky_background'] = med_sky.value
        self.action.args.objects.add_column(phot_table['sky'])
        bkg_sum = phot_table['aperture_sum_1'] / sky_apertures.area() * star_apertures.area()
        final_sum = (phot_table['aperture_sum_0'] - bkg_sum)
        final_uncert = (bkg_sum + final_sum)**0.5 * u.electron**0.5
        phot_table['apflux'] = final_sum/exptime
        self.action.args.objects.add_column(phot_table['apflux'])
        phot_table['apuncert'] = final_uncert/exptime
        self.action.args.objects.add_column(phot_table['apuncert'])
        phot_table['snr'] = final_sum/final_uncert
        self.action.args.objects.add_column(phot_table['snr'])

        where_bad = (final_sum <= 0)
        self.log.info(f'  {np.sum(where_bad)} stars rejected for flux < 0')
        self.action.args.objects = self.action.args.objects[~where_bad]
        self.action.args.meta['n_fluxes'] = len(self.action.args.objects)
        self.log.info(f'  Fluxes for {self.action.args.meta["n_fluxes"]:d} stars')

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: GetCalibrationStars
##-----------------------------------------------------------------------------
class GetCalibrationStars(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = []
        return np.all(checks)

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

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: AssociateCalibratorStars
##-----------------------------------------------------------------------------
class AssociateCalibratorStars(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = []
        return np.all(checks)

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

        return self.action.args

