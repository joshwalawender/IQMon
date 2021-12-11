from keckdrpframework.pipelines.base_pipeline import BasePipeline
from keckdrpframework.models.processing_context import ProcessingContext

# MODIFY THIS IMPORT to reflect the name of the module created in the primitives directory
from iqmon.primitives.file_handling import (ReadFITS,
                                            PopulateAdditionalMetaData,
                                            )
from iqmon.primitives.image_reduction import (SubtractBias,
                                              SubtractDark,
                                              GainCorrect,
                                              CreateDeviation,
                                              )


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
        "create_deviation":  ("CreateDeviation", "creating_deviation", None),

    }

    def __init__(self, context: ProcessingContext):
        BasePipeline.__init__(self, context)
