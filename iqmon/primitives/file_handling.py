from pathlib import Path
from datetime import datetime, timedelta

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
        # initialize values in the args for use with science frames
        self.action.args.meta = {}

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

        hdr = self.action.args.ccddata.header
        for key, val in self.cfg.items('Header'):
            if (val is not 'None') and (hdr.get(val, None) is not None):
                self.action.args.meta[key] = hdr.get(val)
                self.log.debug(f"  {key} = {self.action.args.meta[val]}")
            else:
                self.log.warning(f"  Could not read {key} from header {val} keyword")

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
## Primitive: MoveFile
##-----------------------------------------------------------------------------
class MoveFile(BasePrimitive):
    """
    Take the original raw image read by ReadFITS and write it out as a
    compressed FITS file to a new directory.
    
    If the system uses ACP, look for a similarly named log file and move that
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
        checks = []
        return np.all(checks)

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation.
        """
        self.log.info(f"Running {self.__class__.__name__} action")

        return self.action.args

