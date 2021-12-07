from keckdrpframework.pipelines.base_pipeline import BasePipeline
from keckdrpframework.models.processing_context import ProcessingContext

# MODIFY THIS IMPORT to reflect the name of the module created in the primitives directory
from iqmon.primitives.file_handling import (ReadFITS, MoveFile, DeleteOriginal,
                                            PopulateMetaData, RecordFile)


class Analyze(BasePipeline):
    """
    This pipeline ingests files from their raw location on disk and rewrites
    them to the destination directory with a checksum.  It then (if spcified)
    deletes the original file.  Finally, some basic image information and
    statistics are recorded to the mongo database.
    
    This is meant to be a very quick sequence which moves the file and records
    the file's existence to the database.
    """
    event_table = {
        "next_file":         ("ReadFITS", "reading_file", "science_metadata"),
        "populate_metadata": ("PopulateMetaData", "populating_metadata", None),
        
    }

    def __init__(self, context: ProcessingContext):
        BasePipeline.__init__(self, context)
    
    def template_action(self, action, context):
        self.logger.info("Running template_action on %s", action.args.name)
        return action.args
  