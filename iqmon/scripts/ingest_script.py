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

# the preferred way to import the pipeline is a direct import

from iqmon.pipelines.ingest import IngestPipeline


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
    parser.add_argument("--fz", "--fpack", dest="fpack",
           default=False, action="store_true",
           help="Accept fpack'd files (.fz extension) instead of .fts")
    parser.add_argument('input', type=str, nargs='?',
           help="input image file (full path)", default='')
    args = parser.parse_args(in_args[1:])

    return args


##-----------------------------------------------------------------------------
## Setup Framework
##-----------------------------------------------------------------------------
def setup_framework(args, pipeline=IngestPipeline):
    # START HANDLING OF CONFIGURATION FILES ##########
    pkg = 'iqmon'
    framework_config_file = "configs/framework.cfg"
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)

    framework_logcfg_file = 'configs/logger_ingest.cfg'
    framework_logcfg_fullpath = pkg_resources.resource_filename(pkg, framework_logcfg_file)

    # add PIPELINE specific config files
    if args.config_file is None:
        pipeline_config_file = 'configs/pipeline.cfg'
        pipeline_config_fullpath = pkg_resources.resource_filename(pkg, pipeline_config_file)
        pipeline_config = ConfigClass(pipeline_config_fullpath, default_section='DEFAULT')
    else:
        pipeline_config = ConfigClass(args.pipeline_config_file, default_section='DEFAULT')

    # END HANDLING OF CONFIGURATION FILES ##########

    try:
        framework = Framework(IngestPipeline, framework_config_fullpath)
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
## Analyze One File
##-----------------------------------------------------------------------------
def analyze_one():
    args = _parseArguments(sys.argv)
    p = Path(args.input).expanduser().absolute()
    if p.exists() is False:
        print(f'Unable to find file: {p}')
        return
    args.name = f"{p}"

    framework = setup_framework(args, pipeline=IngestPipeline)
    if args.fpack is True:
        framework.config['DEFAULT']['file_type'] = '*.fz'
    else:
        framework.config['DEFAULT']['file_type'] = '*.fts'

    pkg = 'iqmon'
    framework_config_file = "configs/framework.cfg"
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)
    queue = queues.get_event_queue(cfg.queue_manager_hostname,
                                   cfg.queue_manager_portnr,
                                   cfg.queue_manager_auth_code)
    if queue is None:
        print("Failed to connect to Queue Manager")
        return

    pending = queue.get_pending()
    event = Event("next_file", args)
    queue.put(event)
    framework.config['DEFAULT']['file_type'] = '*.fts'


##-----------------------------------------------------------------------------
## Analyze All Files in a Directory
##-----------------------------------------------------------------------------
def analyze_directory():
    args = _parseArguments(sys.argv)
    data_path = Path(args.input).expanduser().absolute()
    if data_path.exists() is False:
        print(f'Unable to find directory: {data_path}')
        return

    framework = setup_framework(args, pipeline=IngestPipeline)
    if args.fpack is True:
        framework.config['DEFAULT']['file_type'] = '*.fz'
    else:
        framework.config['DEFAULT']['file_type'] = '*.fts'

    pkg = 'iqmon'
    framework_config_file = "configs/framework.cfg"
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)
    queue = queues.get_event_queue(cfg.queue_manager_hostname,
                                   cfg.queue_manager_portnr,
                                   cfg.queue_manager_auth_code)
    if queue is None:
        print("Failed to connect to Queue Manager")
        return
    pending = queue.get_pending()

    infiles = [f for f in data_path.glob(framework.config['DEFAULT']['file_type'])]
    for infile in infiles:
        print(f"Ingesting {infile.name}")
        args.name = f"{infile}"
        event = Event("next_file", args)
        queue.put(event)
        framework.config['DEFAULT']['file_type'] = '*.fts'


##-----------------------------------------------------------------------------
## Watch Directory
##-----------------------------------------------------------------------------
def watch_directory():
    args = _parseArguments(sys.argv)
    framework = setup_framework(args, pipeline=IngestPipeline)
    if args.fpack is True:
        framework.config['DEFAULT']['file_type'] = '*.fz'
    else:
        framework.config['DEFAULT']['file_type'] = '*.fts'
    now = datetime.utcnow()
    data_path = framework.config.instrument.get('FileHandling', 'ingest_dir')
    data_path = data_path.replace('YYYY', f'{now.year:4d}')
    data_path = data_path.replace('MM', f'{now.month:02d}')
    data_path = data_path.replace('DD', f'{now.day:02d}')
    framework.logger.info(f'Setting data path: {data_path}')
    data_path = Path(data_path).expanduser()
    if data_path.exists() is False:
        data_path.mkdir(parents=True, exist_ok=True)

    framework.logger.info(f'Ingesting files from {data_path}')
    infiles = data_path.glob(framework.config['DEFAULT']['file_type'])
    framework.ingest_data(str(data_path), infiles, True)
    framework.start(False, False, False, True)


##-----------------------------------------------------------------------------
## List Queue
##-----------------------------------------------------------------------------
def list_queue():
    args = _parseArguments(sys.argv)
    pkg = 'iqmon'
    framework_config_file = "configs/framework.cfg"
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
    framework_config_file = "configs/framework.cfg"
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
