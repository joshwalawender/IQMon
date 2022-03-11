from pathlib import Path
from datetime import datetime, timedelta
import pymongo
import subprocess

import numpy as np
from astropy import units as u
from astropy import coordinates as c
from astropy.io import fits
from astropy.nddata import CCDData
from astropy import stats
from astropy.time import Time

from keckdrpframework.primitives.base_primitive import BasePrimitive

from .utils import pre_condition, post_condition, get_memory_size, get_destination_dir


##-----------------------------------------------------------------------------
## Primitive: ReadFITS
##-----------------------------------------------------------------------------
class ReadFITS(BasePrimitive):
    """
    Read a FITS file in to a CCDData object as the data model for the pipeline.
    
    Initialize the following self.action.args properties:
    - fitsfilepath : a Path instance
    - skip : boolean
    - ccddata : a CCDData instance
    - meta : a dict of metadata about the image
    - mongodb : a pymongo collection for image data
    - start_time : a datetime instance recording when processing of this file began
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.log.info("")
        self.cfg = self.context.config.instrument
        # initialize values in the args for general use
        self.action.args.fitsfilepath = Path(self.action.args.name).expanduser().absolute()
        self.log.info(f"--> Input FITS file: {self.action.args.fitsfilepath}")
        self.action.args.skip = False

        if self.action.args.fitsfilepath.exists() != True:
            self.log.error('Skipping this file')
            self.action.args.skip = True

        self.action.args.ccddata = None
        self.action.args.destination_dir = None
        self.action.args.destination_file = None
        # initialize values in the args for use with science frames
        self.action.args.meta = {'telescope': self.cfg['Telescope'].get('name')}
        self.action.args.imtype = None
        self.action.args.header_pointing = None
        self.action.args.background = None
        self.action.args.wcs = None
        self.action.args.wcs_pointing = None
        self.action.args.objects = None
        self.action.args.calibration_catalog = None
        self.action.args.associated_calibrators = None

        # If we are reading a compressed file, use the uncompressed version of
        # the name for the database
        if self.action.args.fitsfilepath.suffix == '.fz':
            self.action.args.meta['fitsfile'] = '.'.join(self.action.args.fitsfilepath.name.split('.')[:-1])
        else:
            self.action.args.meta['fitsfile'] = self.action.args.fitsfilepath.name

        ## Connect to mongo
        try:
            self.log.debug('Connecting to mongo db')
            mongo_host = self.cfg['mongo'].get('host')
            mongo_port = self.cfg['mongo'].getint('port')
            mongo_db = self.cfg['mongo'].get('db')
            mongo_collection = 'iqmon'
            self.mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
            self.mongo_iqmon = self.mongoclient[mongo_db][mongo_collection]
        except Exception as e:
            self.log.error('Could not connect to mongo db')
            self.log.error(e)
            self.mongoclient = None
            self.mongo_iqmon = None

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'FITS file exists',
                                self.action.args.fitsfilepath.exists()),
                 ]
        if np.all(checks) is False:
            self.log.ERROR('Skipping this file')
            self.action.args.skip = True
        return np.all(checks)

    def _post_condition(self):
        """Check for conditions necessary to verify that the process run correctly"""
        checks = [post_condition(self, 'FITS file was read',
                                 self.action.args.ccddata is not None),
                  post_condition(self, 'FITS header was read',
                                 len(self.action.args.meta.keys()) > 1),
                  post_condition(self, 'Image type was set',
                                 self.action.args.imtype is not None),
                 ]
        return np.all(checks)

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation.
        """
        self.log.info(f"Running {self.__class__.__name__} action")
        self.action.args.start_time = datetime.now()

        # Read FITS file
        self.log.info(f'Reading: {self.action.args.fitsfilepath}')
        try:
            if self.action.args.fitsfilepath.suffix == '.fz':
                self.action.args.ccddata = CCDData.read(self.action.args.fitsfilepath,
                                                        unit="adu", hdu=1)
            else:
                self.action.args.ccddata = CCDData.read(self.action.args.fitsfilepath,
                                                        unit="adu")
        except Exception as err:
            self.log.error(f"Failed to read file!")
            self.log.error(err)
            self.action.args.skip = True
            return self.action.args

        # Read header metadata
        self.log.info(f'Reading FITS header: Memory Size of ccddata.meta: {get_memory_size(self.action.args.ccddata.meta):.1f} MB')
        hdr = self.action.args.ccddata.header
        for key, raw_read in self.cfg.items('Header'):
            raw_read = raw_read.split(',')
            if key != 'default_section' and key[:11] != 'type_string':
                hdr_key = raw_read[0]
                val = hdr.get(hdr_key, None)
                if len(raw_read) == 2 and val is not None:
                    type_string = raw_read[1]
                    if type_string == 'datetime':
                        val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
                    else:
                        val = eval(type_string)(val)
                elif len(raw_read) > 2:
                    raise TypeError(f'Could not parse "{raw_read}" in Header config')
                if val is not None:
                    self.action.args.meta[key] = val
                    self.log.debug(f"  {key} = {self.action.args.meta[key]} ({type(self.action.args.meta[key])})")
        # Set UT date string
        self.action.args.meta['UT date string'] = self.action.args.meta.get('date').strftime('%Y%m%dUT')

        # Set Image Type
        self.log.info('Determining image type')
        self.log.debug(f"  header imtype = {self.action.args.meta['imtype']}")
        self.action.args.imtype = f"unknown ({self.action.args.meta['imtype']})"
        for key, val in self.cfg.items('Header'):
            if key[:11] == 'type_string':
                if self.action.args.meta['imtype'] in val.split(','):
                    self.action.args.imtype = key[12:].upper()
        self.action.args.meta['imtype'] = self.action.args.imtype
        self.log.info(f"  Image type is {self.action.args.imtype}")

        # Build header pointing coordinate
        if self.action.args.imtype == 'OBJECT':
            self.log.info('Build astropy coordinate for header pointing')
            try:
                self.action.args.header_pointing = c.SkyCoord(self.action.args.meta.get('header_ra'),
                                                              self.action.args.meta.get('header_dec'),
                                                              frame='icrs',
                                                              unit=(u.hourangle, u.deg))
            except:
                self.action.args.header_pointing = None
            if self.action.args.header_pointing is not None:
                self.log.info(f'  {self.action.args.header_pointing.to_string("hmsdms", precision=1)}')
            else:
                self.log.warning(f"  Could not parse header RA and Dec")
                self.log.warning(f"  {self.action.args.meta.get('header_ra')}")
                self.log.warning(f"  {self.action.args.meta.get('header_dec')}")

        # Should we skip a previously processed file?
        if self.mongo_iqmon is not None:
            self.log.info('Check for previous analysis of this file')
            already_processed = [d for d in self.mongo_iqmon.find( {'fitsfile': self.action.args.meta['fitsfile']} )]
            
            # Set skip if image has analysis values
            if len(already_processed) != 0 and self.action.args.imtype == 'OBJECT':
                n_objects = already_processed[0].get('n_objects', 0)
                if n_objects > 0 and self.cfg['mongo'].getboolean('overwrite', False) is False:
                    self.log.info('  File is already in the database, skipping further processing')
                    self.action.args.skip = True

            # Read previously solved WCS if available
            if len(already_processed) != 0 and self.action.args.imtype == 'OBJECT':
                self.log.info(f'  Found {len(already_processed)} database entries for this file')
                db_wcs = already_processed[0].get('wcs_header', None)
                if db_wcs is None:
                    self.log.info('  Database entry does not contain WCS')
                else:
                    from astropy.wcs import WCS
                    self.log.info('  Found Previously Solved WCS')
                    self.action.args.wcs = WCS(db_wcs)
                    self.action.args.meta['wcs_header'] = self.action.args.wcs.to_header_string()
                    nx, ny = self.action.args.ccddata.data.shape
                    r, d = self.action.args.wcs.all_pix2world([nx/2.], [ny/2.], 1)
                    self.action.args.wcs_pointing = c.SkyCoord(r[0], d[0],
                                     frame='icrs', equinox='J2000',
                                     unit=(u.deg, u.deg),
                                     obstime=self.action.args.meta.get('date'))
                    self.action.args.meta['perr'] = self.action.args.wcs_pointing.separation(
                                                    self.action.args.header_pointing).to(u.arcmin).value
                    self.log.info(f'  Pointing error = {self.action.args.meta.get("perr"):.1f}')

        ## VYSOS Specific Code:
        if self.cfg['Telescope'].get('name') in ['V5', 'V20'] and self.action.args.skip is False:
            # Manually set filter for V5A
            if self.action.args.meta.get('filter', None) is None:
                self.action.args.meta['filter'] = 'PSr'
            # Pull focus info from status database
            try:
                self.log.info(f'  Reading focus data from DB')
                status_collection = self.mongoclient[self.cfg['mongo'].get('db')]['V5_focuser']
                date_obs = self.action.args.meta.get('date')
                querydict = {"date": {"$gt": date_obs,
                                      "$lt": date_obs+timedelta(minutes=2)}}
                focus_pos = [entry['position'] for entry in\
                             status_collection.find(querydict).sort(\
                             [('date', pymongo.ASCENDING)])]
                focus_temp = [entry['temperature'] for entry in\
                              status_collection.find(querydict).sort(\
                              [('date', pymongo.ASCENDING)])]
                pmean, pmed, pstd = stats.sigma_clipped_stats(focus_pos)
                tmean, tmed, tstd = stats.sigma_clipped_stats(focus_temp)
                if pstd < 2 and tstd < 0.7:
                    self.action.args.meta["focus_position"] = pmed
                    self.action.args.meta["focus_temperature"] = tmed
                    self.log.info(f"  Focus Position: {pmean:.0f} {pmed:.0f} {pstd:.2f}")
                    self.log.info(f"  Focus Temperature: {tmean:.1f} {tmed:.1f} {tstd:.3f}")
                else:
                    self.log.warning(f"  Focus std dev too high: {pstd:.1f} {tstd:.1f}")
            except Exception as e:
                self.log.warning(f'  Failed to read focus data from DB')
                self.log.warning(e)

        if self.mongoclient is not None:
            self.mongoclient.close()

        self.log.info(f"Memory Size after {self.__class__.__name__} action: {get_memory_size(self):.1f} MB")
#         self.log.info(f"Memory Size after {self.__class__.__name__} action on {self.action.args.meta['fitsfile']}: {get_memory_size(self):.1f} MB")
#         self.log.info(f"  Memory Size of self.action.args.ccddata.meta: {get_memory_size(self.action.args.ccddata.meta):.1f} MB")
#         self.log.info(f"  Memory Size of self.action.args.ccddata.header: {get_memory_size(self.action.args.ccddata.header):.1f} MB")

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: PopulateMetaData
##-----------------------------------------------------------------------------
class PopulateAdditionalMetaData(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'Skip image is not set',
                                self.action.args.skip != True),
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

        if self.action.args.imtype == 'OBJECT':
            if self.action.args.header_pointing is not None:
                self.log.info('Determine Moon info')
                site_lat = c.Latitude(self.cfg['Telescope'].getfloat('site_lat'), unit=u.degree)
                site_lon = c.Longitude(self.cfg['Telescope'].getfloat('site_lon'), unit=u.degree)
                site_elevation = self.cfg['Telescope'].getfloat('site_elevation') * u.meter
                loc = c.EarthLocation(site_lon, site_lat, site_elevation)
                pressure = self.cfg['Telescope'].getfloat('pressure', 700)*u.mbar
                obstime = Time(self.action.args.meta.get('date'), location=loc)
                altazframe = c.AltAz(location=loc, obstime=obstime, pressure=pressure)
                moon = c.get_moon(obstime)
                sun = c.get_sun(obstime)
                moon_alt = ((moon.transform_to(altazframe).alt).to(u.degree)).value
                moon_separation = (moon.separation(self.action.args.header_pointing).to(u.degree)).value\
                            if self.action.args.header_pointing is not None else None
                # Moon illumination formula from Meeus, ÒAstronomical 
                # Algorithms". Formulae 46.1 and 46.2 in the 1991 edition, 
                # using the approximation cos(psi) \approx -cos(i). Error 
                # should be no more than 0.0014 (p. 316). 
                moon_illum = 50*(1 - np.sin(sun.dec.radian)*np.sin(moon.dec.radian)\
                             - np.cos(sun.dec.radian)*np.cos(moon.dec.radian)\
                             * np.cos(sun.ra.radian-moon.ra.radian))
                self.action.args.meta['moon_alt'] = moon_alt
                self.action.args.meta['moon_separation'] = moon_separation
                self.action.args.meta['moon_illum'] = moon_illum
        elif self.action.args.imtype == 'BIAS':
            self.log.info('Determine image stats')
            mean, med, std = stats.sigma_clipped_stats(self.action.args.ccddata.data)
            self.log.info(f"  mean, med, std = {mean:.0f}, {med:.0f}, {std:.0f} (adu)")
            self.action.args.meta['mean adu'] = mean
            self.action.args.meta['median adu'] = med
            self.action.args.meta['std dev adu'] = std
        elif self.action.args.imtype == 'DARK':
            self.log.info('Determine image stats')
            mean, med, std = stats.sigma_clipped_stats(self.action.args.ccddata.data)
            self.log.info(f"  mean, med, std = {mean:.0f}, {med:.0f}, {std:.0f} (adu)")
            self.action.args.meta['mean adu'] = mean
            self.action.args.meta['median adu'] = med
            self.action.args.meta['std dev adu'] = std
        elif self.action.args.imtype in ['DOMEFLAT', 'TWIFLAT']:
            self.log.info('Determine image stats')
            mean, med, std = stats.sigma_clipped_stats(self.action.args.ccddata.data)
            self.log.info(f"  mean, med, std = {mean:.0f}, {med:.0f}, {std:.0f} (adu)")
            self.action.args.meta['mean adu'] = mean
            self.action.args.meta['median adu'] = med
            self.action.args.meta['std dev adu'] = std

#         self.log.info(f"Memory Size after {self.__class__.__name__} action: {get_memory_size(self):.1f} MB")
        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: CopyFile
##-----------------------------------------------------------------------------
class CopyFile(BasePrimitive):
    """
    Take the original raw image read by ReadFITS and write it out as a
    compressed FITS file to a new directory.
    
    If the system uses ACP, look for a similarly named log file and copy that
    as well.
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument
        # Set values for destination
        self.action.args.destination_dir = Path(get_destination_dir(self.cfg,
                                           date=self.action.args.meta['date']))
        self.action.args.destination_file = self.action.args.destination_dir / self.action.args.fitsfilepath.name

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'FITS file exists',
                                self.action.args.fitsfilepath.exists()),
                  pre_condition(self, 'Destination file does not exist',
                                self.action.args.destination_file.exists() is False),
                 ]
        return np.all(checks)

    def _post_condition(self):
        """Check for conditions necessary to verify that the process run correctly"""
        checks = [post_condition(self, 'File was copied',
                                 self.action.args.destination_file.exists() is True),
                  ]
        return np.all(checks)

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation.
        """
        self.log.info(f"Running {self.__class__.__name__} action")
        self.action.args.destination_dir.mkdir(parents=True, exist_ok=True)
        overwrite = self.cfg['FileHandling'].getboolean('overwrite', False)
        
        if self.action.args.fitsfilepath.suffix not in ['.fits', '.fts', '.fz']:
            msg = f'File extension not handled: {self.action.args.fitsfilepath.suffix}'
            self.log.error(msg)
            self.action.args.skip = True
            return self.action.args
        if self.action.args.destination_file.exists() is True and overwrite is True:
            self.action.args.destination_file.unlink()
        if self.action.args.destination_file.exists() is False:
            with fits.open(self.action.args.fitsfilepath, checksum=True) as hdul:
                hdul[0].add_checksum()
                hdul.writeto(self.action.args.destination_file, checksum=True)

            # Compress Destination File
            self.log.info('Compressing destination file')
            compressed_file = self.action.args.destination_file.parent / f"{self.action.args.destination_file.name}.fz"
            if compressed_file.exists() is True:
                compressed_file.unlink()
            self.log.debug(f'  Compressed file: {compressed_file}')
            try:
                self.log.debug(f'  Running fpack')
                subprocess.call(['fpack', self.action.args.destination_file])
            except Exception as e:
                self.log.error(f'  fpack failed')
                self.log.error(e)
            if compressed_file.exists() is True:
                self.action.args.destination_file.unlink()
                self.action.args.destination_file = compressed_file
                self.log.debug(f'  Destination file: {self.action.args.destination_file}')

        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: DeleteOriginal
##-----------------------------------------------------------------------------
class DeleteOriginal(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'delete_original flag is True',
                  self.cfg['FileHandling'].getboolean('delete_original', False) is True),
                  pre_condition(self, 'destination copy exists',
                  self.action.args.destination_file.exists() is True),
                  pre_condition(self, 'Original and destination are different',
                  self.action.args.destination_file != self.action.args.fitsfilepath),
                  ]
        return np.all(checks)

    def _post_condition(self):
        """Check for conditions necessary to verify that the process run correctly"""
        checks = [post_condition(self, 'Original file is gone',
                  self.action.args.fitsfilepath.exists() is False),
                  ]
        return np.all(checks)

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation.
        """
        self.log.info(f"Running {self.__class__.__name__} action")
        self.action.args.fitsfilepath.unlink()
        return self.action.args


##-----------------------------------------------------------------------------
## Primitive: ReleaseMemory
##-----------------------------------------------------------------------------
class ReleaseMemory(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        return True

    def _post_condition(self):
        """
        Check for conditions necessary to verify that the process ran
        correctly.
        """
        return True

    def _perform(self):
        """
        Returns an Argument() with the parameters that depend on this
        operation.
        """
        self.log.info(f"Running {self.__class__.__name__} action")
        self.action.args = None
        self.action.event.args = None
        self.action.config = None
        self.action.cfg = None
        self.action.mongoclient = None
        self.action.mongo_iqmon = None
#         self.log.info(f"Memory Size after {self.__class__.__name__} action: {get_memory_size(self):.1f} MB")
#         self.log.info(f"  Memory size of self.action.event.set_recurrent = {get_memory_size(self.action.event.set_recurrent):.1f} MB")
#         self.log.info(type(self.action.event.set_recurrent))
#         self.log.info(dir(self.action.event.set_recurrent))
        return self.action.args

