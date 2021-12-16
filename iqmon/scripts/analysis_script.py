from keckdrpframework.core.framework import Framework
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.utils.drpf_logger import getLogger
from keckdrpframework.tools.interface import FrameworkInterface, Arguments, Event, Framework, ConfigClass
from keckdrpframework.core import queues
import subprocess
import time
import argparse
import sys
import traceback
import pkg_resources
import logging.config
from pathlib import Path
from datetime import datetime
from glob import glob

from iqmon.pipelines.analyze import AnalysisPipeline


##-----------------------------------------------------------------------------
## Parse Arguments
##-----------------------------------------------------------------------------
def _parseArguments(in_args):
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                      description='')
    parser.add_argument("-v", "--verbose", dest="verbose",
           default=False, action="store_true",
           help="Be verbose.")
    parser.add_argument('-c', dest="config_file", type=str,
           help="Configuration file")
    parser.add_argument("-O", "--overwrite", dest="overwrite",
           default=False, action="store_true",
           help="Reprocess files if they already exist in database?  Only works for analyzeone.")
    parser.add_argument('input', type=str, nargs='?',
           help="input image file (full path)", default='')
    args = parser.parse_args(in_args[1:])

    return args


##-----------------------------------------------------------------------------
## Setup Framework
##-----------------------------------------------------------------------------
def setup_framework(args, pipeline=AnalysisPipeline):
    # START HANDLING OF CONFIGURATION FILES ##########
    pkg = 'iqmon'
    framework_config_file = "configs/framework_analysis.cfg"
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)

    framework_logcfg_file = 'configs/logger_analysis.cfg'
    framework_logcfg_fullpath = pkg_resources.resource_filename(pkg, framework_logcfg_file)

    # add PIPELINE specific config files
    if args.config_file is None:
        pipeline_config_file = 'configs/pipeline.cfg'
        pipeline_config_fullpath = pkg_resources.resource_filename(pkg, pipeline_config_file)
        pipeline_config = ConfigClass(pipeline_config_fullpath, default_section='DEFAULT')
    else:
        pipeline_config = ConfigClass(args.pipeline_config_file, default_section='DEFAULT')

    if args.overwrite is True:
        pipeline_config.set('mongo', 'overwrite', value='True')

    # END HANDLING OF CONFIGURATION FILES ##########

    try:
        framework = Framework(AnalysisPipeline, framework_config_fullpath)
        logging.config.fileConfig(framework_logcfg_fullpath)
        framework.config.instrument = pipeline_config
    except Exception as e:
        print("Failed to initialize framework, exiting ...", e)
        traceback.print_exc()
        sys.exit(1)

    # this part defines a specific logger for the pipeline, so that we can
    # separate the output of the pipeline from the output of the framework
    framework.context.pipeline_logger = getLogger(framework_logcfg_fullpath, name="pipeline")
    framework.logger = getLogger(framework_logcfg_fullpath, name="DRPF")
    framework.logger.info("Framework initialized")

    return framework


##-----------------------------------------------------------------------------
## Watch Directory
##-----------------------------------------------------------------------------
def watch_directory():
    args = _parseArguments(sys.argv)
    framework = setup_framework(args, pipeline=AnalysisPipeline)

    if args.input is not '':
        p = Path(args.input).expanduser()
        args.input = str(p)
        if p.exists() is False:
            framework.context.pipeline_logger.error(f'Could not find: {args.input}')
        else:
            base_path = [x for x in p.parents][-3]
            if p.is_file() is True:
                framework.context.pipeline_logger.info(f'Found file: {args.input}')
            elif p.is_dir() is True:
                framework.context.pipeline_logger.info(f'Found directory: {args.input}')
    else:
        path_str = framework.config.instrument.get('FileHandling', 'destination_dir')
        now = datetime.utcnow()
        path_str = path_str.replace('YYYY', f'{now.year:4d}')
        path_str = path_str.replace('MM', f'{now.month:02d}')
        path_str = path_str.replace('DD', f'{now.day:02d}')
        p = Path(path_str).expanduser().parent
        framework.logger.info(f'Setting data path: {p}')
        if p.exists() is False:
            p.mkdir(parents=True, exist_ok=True)
        args.input = str(p)

    framework.logger.info(f'Ingesting files from {args.input}')
    infiles = glob(f"{args.input}/{framework.config['DEFAULT']['file_type']}")
    framework.ingest_data(args.input, infiles, True)
    framework.start(False, False, False, True)


##-----------------------------------------------------------------------------
## Analyze One File
##-----------------------------------------------------------------------------
def analyze_one():
    args = _parseArguments(sys.argv)
    framework = setup_framework(args, pipeline=AnalysisPipeline)
    p = Path(args.input).expanduser().absolute()
    if p.exists() is False:
        print(f'Unable to find file: {p}')
        return

    print(f"Triggering ingest of {p}")
    framework.ingest_data(path=None, files=[str(p)], monitor=False)


##-----------------------------------------------------------------------------
## Change Watched Directory
##-----------------------------------------------------------------------------
def change_directory():
    args = _parseArguments(sys.argv)
    framework = setup_framework(args, pipeline=AnalysisPipeline)
    if args.input is not '':
        newdir = Path(args.input).expanduser().absolute()
    else:
        date_string = datetime.utcnow().strftime('%Y%m%dUT')
        newdir = Path(f'~/V5Data/Images/{date_string}').expanduser()

    args.input = str(newdir)
    if newdir.exists() is False:
        newdir.mkdir(parents=True)

    framework.ingest_data(path=args.input, files=None, monitor=True)


##-----------------------------------------------------------------------------
## List Queue
##-----------------------------------------------------------------------------
def list_queue():
    args = _parseArguments(sys.argv)
    pkg = 'iqmon'
    framework_config_file = "configs/framework_analysis.cfg"
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)
    drpif = FrameworkInterface(cfg)
    # Print pending Events
    if drpif.is_queue_ok():
        events = drpif.pending_events()
        print(f'Found {len(events)} in queue')
        if args.verbose is True:
            for event in events:
                print(event)
    else:
        print ("Pending events: Queue not available", drpif.queue)


##-----------------------------------------------------------------------------
## Clear Queue
##-----------------------------------------------------------------------------
def clear_queue():
    args = _parseArguments(sys.argv)

    pkg = 'iqmon'
    framework_config_file = "configs/framework_analysis.cfg"
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)

    drpif = FrameworkInterface(cfg)
    # Print pending Events
    if drpif.is_queue_ok():
        events = drpif.pending_events()
        print(f'Found {len(events)} in queue')
    else:
        print ("Pending events: Queue not available", drpif.queue)

    if drpif.is_queue_ok():
        drpif.stop_event_queue()
        print ("Queue manager stopped")
    else:
        print ("Queue manager already stopped")

if __name__ == "__main__":
    analyze_one()

