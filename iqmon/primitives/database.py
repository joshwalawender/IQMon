from pathlib import Path
from datetime import datetime, timedelta
import numpy as np

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

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'FITS file was read',
                                self.action.args.ccddata is not None),
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

        self.log.info(f"Recording the following metadata:")
        for key in self.action.args.meta:
            val = self.action.args.meta.get(key)
            self.log.info(f"  {key:15s} : {val} ({type(val)})")

        return self.action.args


