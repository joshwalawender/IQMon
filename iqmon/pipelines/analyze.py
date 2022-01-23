from keckdrpframework.pipelines.base_pipeline import BasePipeline
from keckdrpframework.models.processing_context import ProcessingContext

# MODIFY THIS IMPORT to reflect the name of the module created in the primitives directory
from iqmon.primitives.file_handling import (ReadFITS,
                                            PopulateAdditionalMetaData,
                                            ReleaseMemory,
                                            )
from iqmon.primitives.image_reduction import (SubtractBias,
                                              SubtractDark,
                                              GainCorrect,
                                              CreateDeviation,
                                              MakeMasterCalFrame,
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
    """
    event_table = {
        "next_file":         ("ReadFITS",
                              "reading_file",
                              "populate_metadata"),
        "populate_metadata": ("PopulateAdditionalMetaData",
                              "populating_metadata",
                              "subtract_bias"),
        "subtract_bias":     ("SubtractBias",
                              "subtracting_bias",
                              "subtract_dark"),
        "subtract_dark":     ("SubtractDark",
                              "subtracting_dark",
                              "stack_cal_frames"),
        "stack_cal_frames":  ("MakeMasterCalFrame",
                              "making master cal frame",
                              "gain_correct"),
        "gain_correct":      ("GainCorrect",
                              "correcting_gain",
                              "create_deviation"),
        "create_deviation":  ("CreateDeviation",
                              "creating_deviation",
                              "make_source_mask"),
        "make_source_mask":  ("MakeSourceMask",
                              "making_source_mask",
                              "create_background"),
        "create_background": ("CreateBackground",
                              "creating_background",
                              "extract"),
        "extract":           ("ExtractStars",
                              "extracting_stars",
                              "solve_astrometry"),
        "solve_astrometry":  ("SolveAstrometry",
                              "solving",
                              "render_jpeg"),
#         "get_calibrators":   ("GetCalibrationStars",
#                               "getting_calibrators",
#                               "associate_calibrators"),
#         "associate_calibrators": ("AssociateCalibratorStars",
#                                   "associating_calibrators",
#                                   "render_jpeg"),
        "render_jpeg":       ("RenderJPEG",
                              "rendering_jpeg",
                              "record_analysis"),
        "record_analysis":   ("RecordFile",
                              "recording",
                              "release_memory"),
        "release_memory":    ("ReleaseMemory",
                              "releasing memory",
                              None),
    }

    def __init__(self, context: ProcessingContext):
        BasePipeline.__init__(self, context)
