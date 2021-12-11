from keckdrpframework.pipelines.base_pipeline import BasePipeline
from keckdrpframework.models.processing_context import ProcessingContext

# MODIFY THIS IMPORT to reflect the name of the module created in the primitives directory
from iqmon.primitives.file_handling import (ReadFITS,
                                            PopulateAdditionalMetaData,
                                            CopyFile,
                                            DeleteOriginal)
from iqmon.primitives.database import RecordFile


class IngestPipeline(BasePipeline):
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
        "populate_metadata": ("PopulateAdditionalMetaData", "populating_metadata", "copy_file"),
        "copy_file":         ("CopyFile", "copying_file", "delete_original"),
        "delete_original":   ("DeleteOriginal", "deleting_original", "record_file"),
        "record_file":       ("RecordFile", "recording", None),
        
    }

    def __init__(self, context: ProcessingContext):
        BasePipeline.__init__(self, context)
