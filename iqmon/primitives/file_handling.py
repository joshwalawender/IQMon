from pathlib import Path
from datetime import datetime, timedelta

from astropy import units as u
from astropy import coordinates as c
from astropy.io import fits
from astropy.nddata import CCDData

from keckdrpframework.primitives.base_primitive import BasePrimitive

from .utils import pre_condition, post_condition, get_destination_dir


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
        self.action.args.skip = False
        self.action.args.ccddata = None
        self.action.args.destination_dir = None
        self.action.args.destination_file = None
        # initialize values in the args for use with science frames
        self.action.args.meta = {'telescope': self.cfg['Telescope'].get('name')}
        self.action.args.imtype = None
        self.action.args.header_pointing = None

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
            mongo_collection = self.cfg['mongo'].get('collection')
            mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
            self.action.args.mongodb = mongoclient[mongo_db][mongo_collection]
        except:
            self.log.error('Could not connect to mongo db')
            self.action.args.mongodb = None

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'FITS file exists',
                                self.action.args.fitsfilepath.exists()),
                 ]
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

        if self.action.args.mongodb is not None:
            already_processed = [d for d in self.action.args.mongodb.find( {'filename': self.action.args.meta['fitsfile']} )]
            if len(already_processed) != 0\
               and self.cfg['FileHandling'].getboolean('overwrite', False) is False:
                self.log.info(f"overwrite is {self.cfg['FileHandling'].getboolean('overwrite')}")
                self.log.info('  File is already in the database, skipping further processing')
                self.action.args.skip = True

        # Read FITS file
        self.log.info(f"  Reading: {self.action.args.meta['fitsfile']}")
        self.action.args.ccddata = CCDData.read(self.action.args.fitsfilepath,
                                                unit="adu")
        # Read header metadata
        hdr = self.action.args.ccddata.header
        for key, hdr_key in self.cfg.items('Header'):
            if hdr_key not in ['None', None]:
                raw_read = hdr.get(hdr_key, '').split(',')
                val = raw_read[0]
                if len(raw_read) == 2:
                    type_string = raw_read[1]
                    if type_string == 'datetime':
                        val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
                    else:
                        val = eval(type_string)(val)
                elif len(raw_read) > 2:
                    raise TypeError(f'Could not parse "{hdr_key}" in Header config')
                self.action.args.meta[key] = val
                self.log.debug(f"  {key} = {self.action.args.meta[key]} ({type(self.action.args.meta[key])})")

        # Set Image Type
        self.action.args.imtype = None
        for key, val in self.cfg.items('Header'):
            if key[:11] == 'type_string':
                if self.action.args.meta['imtype'] == self.cfg['Header'].get(key):
                    self.action.args.imtype = key[12:].upper()
        self.log.info(f"  Image type is {self.action.args.imtype}")

        ## VYSOS Specific Code:
        if self.cfg['telescope'].get('name') == 'V5':
            try:
                self.log.info(f'  Reading focus data from DB')
                status_collection = mongoclient[mongo_db][f'V5status']
                date_obs = datetime.strptime(self.action.args.ccddata.header.get('DATE-OBS'), '%Y-%m-%dT%H:%M:%S')
                querydict = {"date": {"$gt": date_obs, "$lt": date_obs+timedelta(minutes=2)},
                             "telescope": 'V5'}
                focus_pos = [entry['focuser_position'] for entry in\
                             status_collection.find(querydict).sort(\
                             [('date', pymongo.ASCENDING)])]
                focus_temp = [entry['focuser_temperature'] for entry in\
                              status_collection.find(querydict).sort(\
                              [('date', pymongo.ASCENDING)])]
                pmean, pmed, pstd = stats.sigma_clipped_stats(focus_pos)
                tmean, tmed, tstd = stats.sigma_clipped_stats(focus_temp)
                if pstd < 2 and tstd < 0.4:
                    self.action.args.meta["focus_position"] = pmed
                    self.action.args.meta["focus_temperature"] = tmed
                    self.log.info(f"  Focus Position: {pmean:.0f} {pmed:.0f} {pstd:.2f}")
                    self.log.info(f"  Focus Temperature: {tmean:.1f} {tmed:.1f} {tstd:.3f}")
                else:
                    self.log.warning(f"  Focus std dev too high: {pstd:.1f} {tstd:.1f}")
            except Exception as e:
                self.log.warning(f'  Failed to read focus data from DB')
                self.log.warning(e)

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

        if self.action.args.imtype == 'object':
            # Build header pointing
            try:
                self.action.args.header_pointing = c.SkyCoord(self.action.args.meta.get('RA'),
                                                              self.action.args.meta.get('DEC'),
                                                              frame='icrs',
                                                              unit=(u.hourangle, u.deg))
            except:
                self.action.args.header_pointing = None

            # Determine Moon Info
            if self.action.args.header_pointing is not None:
                site_lat = c.Latitude(self.cfg['Telescope'].getfloat('site_lat'), unit=u.degree)
                site_lon = c.Longitude(self.cfg['Telescope'].getfloat('site_lon'), unit=u.degree)
                site_elevation = self.cfg['Telescope'].getfloat('site_elevation') * u.meter
                loc = c.EarthLocation(site_lon, site_lat, site_elevation)
                pressure = self.cfg['Telescope'].getfloat('pressure', 700)*u.mbar
                obstime = Time(self.action.args.meta.get('date_obs'), location=loc)
                altazframe = c.AltAz(location=loc, obstime=obstime, pressure=pressure)
                moon = c.get_moon(obstime)
                sun = c.get_sun(obstime)
                moon_alt = ((moon.transform_to(altazframe).alt).to(u.degree)).value
                moon_separation = (moon.separation(self.action.args.header_pointing).to(u.degree)).value\
                            if self.action.args.header_pointing is not None else None
                # Moon illumination formula from Meeus, “Astronomical 
                # Algorithms". Formulae 46.1 and 46.2 in the 1991 edition, 
                # using the approximation cos(psi) \approx -cos(i). Error 
                # should be no more than 0.0014 (p. 316). 
                moon_illum = 50*(1 - np.sin(sun.dec.radian)*np.sin(moon.dec.radian)\
                             - np.cos(sun.dec.radian)*np.cos(moon.dec.radian)\
                             * np.cos(sun.ra.radian-moon.ra.radian))
                self.action.args.meta['moon_alt'] = moon_alt
                self.action.args.meta['moon_separation'] = moon_separation
                self.action.args.meta['moon_illum'] = moon_illum
        elif self.action.args.imtype == 'bias':
            pass
        elif self.action.args.imtype == 'dark':
            pass
        elif self.action.args.imtype in ['domeflat', 'twiflat']:
            pass

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

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'FITS file exists',
                                self.action.args.fitsfilepath.exists()),
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
        self.action.args.destination_dir = Path(get_destination_dir(self.cfg))
        destination_dir.mkdir(parents=True, exist_ok=True)
        
        overwrite = self.cfg['FileHandling'].getboolean('overwrite', False)
        
        if self.action.args.fitsfilepath.suffix in ['.fits', '.fts']:
            self.action.args.destination_file = destination_dir / f"{self.action.args.fitsfilepath.name}.fz"
        elif self.action.args.fitsfilepath.suffix == '.fz':
            self.action.args.destination_file = destination_dir / self.action.args.fitsfilepath.name
        else:
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

