from keckdrpframework.pipelines.base_pipeline import BasePipeline
from keckdrpframework.models.processing_context import ProcessingContext
from keckdrpframework.primitives.base_primitive import BasePrimitive
import numpy as np

# MODIFY THIS IMPORT to reflect the name of the module created in the primitives directory
from iqmon.primitives.file_handling import (ReadFITS,
                                            PopulateAdditionalMetaData,
                                            )
from iqmon.primitives.image_reduction import (SubtractBias,
                                              SubtractDark,
                                              GainCorrect,
                                              CreateDeviation,
                                              )
from iqmon.primitives.photometry import (MakeSourceMask,
                                         CreateBackground,
                                         ExtractStars,
                                         )
from iqmon.primitives.astrometry import SolveAstrometry
from iqmon.primitives.graphics import RenderJPEG
from iqmon.primitives.database import RecordFile


class AnalysisPipeline(BasePipeline):
    """
    This pipeline ingests files from their raw location on disk and rewrites
    them to the destination directory with a checksum.  It then (if spcified)
    deletes the original file.  Finally, some basic image information and
    statistics are recorded to the mongo database.
    
    This is meant to be a very quick sequence which moves the file and records
    the file's existence to the database.
    """
    event_table = {
        "next_file":         ("ReadFITS", "reading_file", "populate_metadata"),
        "populate_metadata": ("PopulateAdditionalMetaData", "populating_metadata", "subtract_bias"),
        "subtract_bias":     ("SubtractBias", "subtracting_bias", "subtract_dark"),
        "subtract_dark":     ("SubtractDark", "subtracting_dark", "gain_correct"),
        "gain_correct":      ("GainCorrect", "correcting_gain", "create_deviation"),
        "create_deviation":  ("CreateDeviation", "creating_deviation", "make_source_mask"),
        "make_source_mask":  ("MakeSourceMask", "making_source_mask", "create_background"),
        "create_background": ("CreateBackground", "creating_background", "extract"),
        "extract":           ("ExtractStars", "extracting_stars", "solve_astrometry"),
        "solve_astrometry":  ("SolveAstrometry", "solving", "render_jpeg"),
#         "get_calibrators": ("GetCalibrationStars", "getting_calibrators", "associate_calibrators"),
#         "associate_calibrators": ("AssociateCalibratorStars", "associating_calibrators", "render_jpeg"),
        "render_jpeg":       ("RenderJPEG", "rendering_jpeg", "record_analysis"),
        "record_analysis":   ("RecordFile", "recording", "release_memory"),
        "release_memory":    ("ReleaseMemory", "releaseing", None),
    }

    def __init__(self, context: ProcessingContext):
        BasePipeline.__init__(self, context)


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

        if hasattr(self.action.args, 'ccddata'):
            del(self.action.args.ccddata)
        if hasattr(self.action.args, 'background'):
            del(self.action.args.background)
        if hasattr(self.action.args, 'objects'):
            del(self.action.args.objects)

        return self.action.args

