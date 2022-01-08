import sys
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import pymongo

from keckdrpframework.primitives.base_primitive import BasePrimitive

from .utils import pre_condition, post_condition


##-----------------------------------------------------------------------------
## Primitive: RecordFile
##-----------------------------------------------------------------------------
class RecordFile(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

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
        checks = [pre_condition(self, 'Skip image is not set',
                                not self.action.args.skip),
                  pre_condition(self, 'norecord configuration is not set',
                                not self.cfg['Telescope'].getboolean('norecord', False)),
                  pre_condition(self, 'connected to mongoDB',
                                self.mongo_iqmon is not None),
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
        end_time = datetime.now()
        self.action.args.meta['processing_time'] = (end_time - self.action.args.start_time).total_seconds()
        self.log.info(f"  Processing time = {self.action.args.meta.get('processing_time'):.0f} s")

        self.log.info(f"Recording the following metadata:")
        for key in self.action.args.meta:
            val = self.action.args.meta.get(key)
            self.log.info(f"  {key:15s} : {val} ({type(val)})")

        # Remove old entries for this image file
        deletion = self.mongo_iqmon.delete_many( {'fitsfile': self.action.args.meta.get('fitsfile')} )
        self.log.info(f'  Deleted {deletion.deleted_count} previous entries for {self.action.args.meta.get("fitsfile")}')

        # Save new entry for this image file
        self.log.info('  Adding image info to mongo database')
        ## Save document
        try:
            inserted_id = self.mongo_iqmon.insert_one(self.action.args.meta).inserted_id
            self.log.debug(f"  Inserted document id: {inserted_id}")
        except:
            e = sys.exc_info()[0]
            self.log.error('Failed to add new document')
            self.log.error(e)
        self.mongoclient.close()

        return self.action.args


