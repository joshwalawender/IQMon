from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
from astropy import units as u
from astropy import stats
from astropy.time import Time
from astropy.nddata import CCDData
import ccdproc


from keckdrpframework.primitives.base_primitive import BasePrimitive

from .utils import pre_condition, post_condition, find_master


##-----------------------------------------------------------------------------
## Primitive: SubtractBias
##-----------------------------------------------------------------------------
class SubtractBias(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument
        master_dir = self.cfg['Calibrations'].get('DirectoryForMasters', None)
        self.master_bias_file = find_master(master_dir, 'Bias', action.args.meta)

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'Skip image is not set',
                                not self.action.args.skip),
                  pre_condition(self, 'master bias is available',
                                self.master_bias_file is not None),
                  pre_condition(self, 'Image type is not BIAS',
                                self.action.args.imtype != 'BIAS'),
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
        self.log.info(f"  Found master bias file: {self.master_bias_file.name}")
        master_bias_ccddata = CCDData.read(self.master_bias_file, unit="adu")
        self.log.info(f"  Subtracting bias")
        self.action.args.ccddata = ccdproc.subtract_bias(self.action.args.ccddata,
                                                         master_bias_ccddata)

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: SubtractDark
##-----------------------------------------------------------------------------
class SubtractDark(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument
        master_dir = self.cfg['Calibrations'].get('DirectoryForMasters', None)
        self.master_dark_file = find_master(master_dir, 'Dark', action.args.meta)

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'Skip image is not set',
                                not self.action.args.skip),
                  pre_condition(self, 'master dark is available',
                                self.master_dark_file is not None),
                  pre_condition(self, 'Image type is not DARK',
                                self.action.args.imtype != 'DARK'),
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
        self.log.info(f"  Found master dark file: {self.master_dark_file.name}")
        master_dark_ccddata = CCDData.read(self.master_dark_file, unit="adu")
        self.log.info(f"  Subtracting dark")
        self.action.args.ccddata = ccdproc.subtract_bias(self.action.args.ccddata,
                                                         master_dark_ccddata)

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: GainCorrect
##-----------------------------------------------------------------------------
class GainCorrect(BasePrimitive):
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

        gain = self.action.args.meta.get('GAIN', None)
        if gain is not None: self.log.debug(f'  Using gain = {gain}')
        if gain is None:
            gain = self.cfg['Telescope'].getfloat('gain', None)
            self.log.debug(f'  Got gain from config: {gain}')

        self.log.debug('  Gain correcting data')
        self.action.args.ccddata = ccdproc.gain_correct(self.action.args.ccddata,
                                           gain,
                                           gain_unit=u.electron/u.adu)

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: CreateDeviation
##-----------------------------------------------------------------------------
class CreateDeviation(BasePrimitive):
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

        read_noise = self.action.args.meta.get('read_noise', None)
        if read_noise is not None: self.log.debug(f'  Using read_noise = {read_noise}')
        if read_noise is None:
            read_noise = self.cfg['Telescope'].getfloat('read_noise', None)
            self.log.debug(f'  Got read_noise from config: {read_noise}')

        self.action.args.ccddata = ccdproc.create_deviation(self.action.args.ccddata,
                                           readnoise=read_noise*u.electron)

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: MakeMasterCalFrame
##-----------------------------------------------------------------------------
class MakeMasterCalFrame(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument
        self.context = context
        date_string = self.action.args.meta['UT date string']
#         if not hasattr(self.context, 'date_string'):
#             self.context[date_string] = {}
#         if self.action.args.imtype not in self.context[date_string].keys():
#             self.context[date_string][self.action.args.imtype] = []

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        imtype = self.action.args.imtype
        date_string = self.action.args.meta['UT date string']
        checks = [pre_condition(self, 'Skip image is not set',
                                not self.action.args.skip),
                  pre_condition(self, 'Image type is cal',
                                imtype in ['BIAS', 'DARK']),
#                   pre_condition(self, 'Connected to mongo',
#                                 self.mongo_iqmon is not None),
                  pre_condition(self, f'do_{imtype}_subtraction is True',
                                self.cfg['Calibrations'].getboolean(f'do_{imtype}_subtraction', True) is True),
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

        imtype = self.action.args.imtype
        date_string = self.action.args.meta['UT date string']
        self.context[date_string][self.action.args.imtype].append(self.action.args.ccddata)

        n_cals = len(self.context[date_string][imtype])
        self.log.info(f"Found {n_cals} {imtype} files for {date_string}")
        if n_cals >= self.cfg['Calibrations'].getint(f"min_{imtype}_frames"):
            self.log.info(f"Stacking {n_cals} {imtype} files")
            combined = ccdproc.combine(self.context[date_string][imtype],
                                       method='average',
                                       sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                       sigma_clip_func=np.ma.median, sigma_clip_dev_func=stats.mad_std,
                                      )
            self.log.info(f"  Combined.")
            combined_bias.meta['combined'] = True
            combined_bias.meta['ncomb'] = n_cals

            combined_filename = f'Master{imtype}_{date_string}.fits'
            combined_filepath = Path(self.cfg['Calibrations'].get('directory_for_masters'))
            combined_file = combined_filepath.joinpath(combined_filename)
            if combined_file.exists() is True:
                self.log.debug(f"  Deleting existing: {combined_file}")
                combined_file.unlink()
            self.log.info(f"  Saving: {combined_file}")
            combined.write(combined_file)

            return self.action.args

