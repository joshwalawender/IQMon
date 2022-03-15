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
import gc

# the preferred way to import the pipeline is a direct import

from iqmon import get_webpage_config, get_all_configs
from iqmon.pipelines.ingest import IngestPipeline
from iqmon.pipelines.analyze import AnalysisPipeline


##-------------------------------------------------------------------------
## Function: get_memory_size
##-------------------------------------------------------------------------
def get_memory_size(input_obj):
    memory_size = 0
    ids = set()
    objects = [input_obj]
    while objects:
        new = []
        for obj in objects:
            if id(obj) not in ids:
                ids.add(id(obj))
                memory_size += sys.getsizeof(obj)
                new.append(obj)
        objects = gc.get_referents(*new)
    return memory_size/1024/1024


##-----------------------------------------------------------------------------
## Parse Arguments
##-----------------------------------------------------------------------------
def _parseArguments(in_args):
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                      description='')
    parser.add_argument("-v", "--verbose", dest="verbose",
           default=False, action="store_true",
           help="Be verbose.")
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
def setup_framework(args, pipeline=IngestPipeline,
                    framework_config_file="configs/framework.cfg",
                    framework_logcfg_file='configs/logger_ingest.cfg',
                    pipeline_config_file='configs/pipeline.cfg',
                    ):
    # START HANDLING OF CONFIGURATION FILES ##########
    pkg = 'iqmon'
    # Handle framework config
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    framework_logcfg_fullpath = pkg_resources.resource_filename(pkg, framework_logcfg_file)
    # add PIPELINE specific config files
    if pipeline_config_file[0] in ['~', '/']:
        # This is a full path specification
        pipeline_config_fullpath = str(Path(pipeline_config_file).expanduser())
    else:
        # This is a relative path specification
        pipeline_config_fullpath = pkg_resources.resource_filename(pkg, pipeline_config_file)
    pipeline_config = ConfigClass(pipeline_config_fullpath, default_section='DEFAULT')
    # END HANDLING OF CONFIGURATION FILES ##########
    try:
        framework = Framework(pipeline, framework_config_fullpath)
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
    framework.logger.debug(f"  Framework config file: {framework_config_file}")
    framework.logger.debug(f"  {framework_config_fullpath}")
    framework.logger.debug(f"  Pipeline config file: {pipeline_config_file}")
    framework.logger.debug(f"  {pipeline_config_fullpath}")

    return framework


##-----------------------------------------------------------------------------
## Analyze One File
##-----------------------------------------------------------------------------
def analyze_onefile(pipeline=IngestPipeline,
                    framework_config_file="configs/framework.cfg",
                    framework_logcfg_file='configs/logger_ingest.cfg',
                    pipeline_config_file = 'configs/pipeline.cfg'):
    args = _parseArguments(sys.argv)
    p = Path(args.input).expanduser().absolute()
    if p.exists() is False:
        print(f'Unable to find file: {p}')
        return
    args.name = f"{p}"

    framework = setup_framework(args, pipeline=pipeline,
                                framework_config_file=framework_config_file,
                                framework_logcfg_file=framework_logcfg_file,
                                pipeline_config_file=pipeline_config_file,
                                )
    if args.fpack is True:
        framework.config['DEFAULT']['file_type'] = '*.fz'
    else:
        framework.config['DEFAULT']['file_type'] = '*.fts'

    pkg = 'iqmon'
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)
    queue = queues.get_event_queue(cfg.queue_manager_hostname,
                                   cfg.queue_manager_portnr,
                                   cfg.queue_manager_auth_code)
    if queue is None:
        print("Failed to connect to Queue Manager")
    else:
        pending = queue.get_pending()
        event = Event("next_file", args)
        queue.put(event)
    if args.fpack is True:
        framework.config['DEFAULT']['file_type'] = '*.fts'


##-----------------------------------------------------------------------------
## Analyze All Files in a Directory
##-----------------------------------------------------------------------------
def analyze_directory(pipeline=IngestPipeline,
                      framework_config_file="configs/framework.cfg",
                      framework_logcfg_file='configs/logger_ingest.cfg',
                      pipeline_config_file = 'configs/pipeline.cfg',
                      ):
    args = _parseArguments(sys.argv)
    data_path = Path(args.input).expanduser().absolute()
    if data_path.exists() is False:
        print(f'Unable to find directory: {data_path}')
        return

    framework = setup_framework(args, pipeline=pipeline,
                                framework_config_file=framework_config_file,
                                framework_logcfg_file=framework_logcfg_file,
                                pipeline_config_file=pipeline_config_file,
                                )
    if args.fpack is True:
        framework.config['DEFAULT']['file_type'] = '*.fz'
    else:
        framework.config['DEFAULT']['file_type'] = '*.fts'

    pkg = 'iqmon'
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)
    queue = queues.get_event_queue(cfg.queue_manager_hostname,
                                   cfg.queue_manager_portnr,
                                   cfg.queue_manager_auth_code)
    if queue is None:
        print("Failed to connect to Queue Manager")
    else:
        pending = queue.get_pending()
        infiles = [f for f in data_path.glob(framework.config['DEFAULT']['file_type'])]
        for infile in infiles:
#             print(f"Ingesting {infile.name}")
            args.name = f"{infile}"
            event = Event("next_file", args)
            queue.put(event)
            framework.config['DEFAULT']['file_type'] = '*.fts'
    if args.fpack is True:
        framework.config['DEFAULT']['file_type'] = '*.fts'


##-----------------------------------------------------------------------------
## Watch Directory
##-----------------------------------------------------------------------------
def watch_directory(pipeline=IngestPipeline,
                    framework_config_file="configs/framework.cfg",
                    framework_logcfg_file='configs/logger_ingest.cfg',
                    pipeline_config_file = 'configs/pipeline.cfg',
                    ):
    args = _parseArguments(sys.argv)
    framework = setup_framework(args, pipeline=pipeline,
                                framework_config_file=framework_config_file,
                                framework_logcfg_file=framework_logcfg_file,
                                pipeline_config_file=pipeline_config_file,
                                )
    nevents1 = list_queue(framework_config_file=framework_config_file)
    framework.logger.info(f"  Found {nevents1} events in queue before ingesting input dir")
    if nevents1 > 0:
        framework.logger.info(f"  Clearing queue")
        clear_queue(pipeline=pipeline,
                    framework_config_file=framework_config_file,
                    framework_logcfg_file=framework_logcfg_file,
                    pipeline_config_file=pipeline_config_file,
                    )
        nevents2 = list_queue(framework_config_file=framework_config_file)
        framework.logger.info(f"  Found {nevents2} events in queue before ingesting input dir")

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
#     infiles = data_path.glob(framework.config['DEFAULT']['file_type'])
    framework.ingest_data(path=str(data_path), files=None, monitor=True)
    nevents3 = list_queue(framework_config_file=framework_config_file)
    framework.logger.info(f"  Found {nevents3} events in queue")
    framework.start(False, False, False, True)


##-----------------------------------------------------------------------------
## Change Watched Directory
##-----------------------------------------------------------------------------
def change_directory(newdir='~',
                     pipeline=IngestPipeline,
                     framework_config_file="configs/framework.cfg",
                     framework_logcfg_file='configs/logger_ingest.cfg',
                     pipeline_config_file = 'configs/pipeline.cfg',):
    args = _parseArguments(sys.argv)
    framework = setup_framework(args, pipeline=pipeline,
                                framework_config_file=framework_config_file,
                                framework_logcfg_file=framework_logcfg_file,
                                pipeline_config_file=pipeline_config_file,
                                )
    if args.input is not '':
        newdir = Path(args.input).expanduser().absolute()
    else:
        newdir = Path(newdir).expanduser()

    args.input = str(newdir)
    if newdir.exists() is False:
        newdir.mkdir(parents=True)

    framework.ingest_data(path=args.input, files=None, monitor=True)


##-----------------------------------------------------------------------------
## List Queue
##-----------------------------------------------------------------------------
def list_queue(framework_config_file="configs/framework.cfg"):
    args = _parseArguments(sys.argv)
    pkg = 'iqmon'
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
        return len(events)
    else:
        print ("Pending events: Queue not available", drpif.queue)
        return None


##-----------------------------------------------------------------------------
## Clear Queue
##-----------------------------------------------------------------------------
def clear_queue(pipeline=IngestPipeline,
                framework_config_file="configs/framework.cfg",
                framework_logcfg_file='configs/logger_ingest.cfg',
                pipeline_config_file = 'configs/pipeline.cfg',
                ):
    args = _parseArguments(sys.argv)
    pkg = 'iqmon'
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)
    drpif = FrameworkInterface(cfg)

    framework = setup_framework(args, pipeline=pipeline,
                                framework_config_file=framework_config_file,
                                framework_logcfg_file=framework_logcfg_file,
                                pipeline_config_file=pipeline_config_file,
                                )

    # Print pending Events
    if drpif.is_queue_ok():
        events = drpif.pending_events()
        print(f'Found {len(events)} in queue')
        for i in range(len(events)):
            e = framework.get_event()
        events = drpif.pending_events()
        print(f'Found {len(events)} in queue')

    else:
        print ("Pending events: Queue not available", drpif.queue)

#     if drpif.is_queue_ok():
#         drpif.stop_event_queue()
#         print ("Queue manager stopped")
#     else:
#         print ("Queue manager already stopped")


def examine_memory_use(pipeline=IngestPipeline,
                framework_config_file="configs/framework.cfg",
                framework_logcfg_file='configs/logger_ingest.cfg',
                pipeline_config_file = 'configs/pipeline.cfg',
                ):
    args = _parseArguments(sys.argv)
    pkg = 'iqmon'
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)
    cfg = ConfigClass(framework_config_fullpath)
    drpif = FrameworkInterface(cfg)

    framework = setup_framework(args, pipeline=pipeline,
                                framework_config_file=framework_config_file,
                                framework_logcfg_file=framework_logcfg_file,
                                pipeline_config_file=pipeline_config_file,
                                )

    print(f"framework = {get_memory_size(framework):.1f} MB")
    print(f"drpif = {get_memory_size(drpif):.1f} MB")



##-----------------------------------------------------------------------------
## Ingest Scripts
##-----------------------------------------------------------------------------
def ingest_start_queue():
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    watch_directory(pipeline=IngestPipeline,
                    framework_config_file="configs/framework_ingest.cfg",
                    framework_logcfg_file='configs/logger_ingest.cfg',
                    pipeline_config_file=pipeline_config_file,
                    )

def ingest_one():
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    analyze_onefile(pipeline=IngestPipeline,
                    framework_config_file="configs/framework_ingest.cfg",
                    framework_logcfg_file='configs/logger_ingest.cfg',
                    pipeline_config_file=pipeline_config_file,
                    )

def ingest_all():
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    analyze_directory(pipeline=IngestPipeline,
                      framework_config_file="configs/framework_ingest.cfg",
                      framework_logcfg_file='configs/logger_ingest.cfg',
                      pipeline_config_file=pipeline_config_file,
                      )

def ingest_list(framework_config_file="configs/framework_ingest.cfg"):
    list_queue(framework_config_file=framework_config_file)

def ingest_clear(framework_config_file="configs/framework_ingest.cfg"):
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    clear_queue(pipeline=IngestPipeline,
                framework_config_file="configs/framework_ingest.cfg",
                framework_logcfg_file='configs/logger_ingest.cfg',
                pipeline_config_file=pipeline_config_file,
                )



##-----------------------------------------------------------------------------
## Analysis Scripts
##-----------------------------------------------------------------------------
def analyze_start_queue():
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    watch_directory(pipeline=AnalysisPipeline,
                    framework_config_file="configs/framework_analysis.cfg",
                    framework_logcfg_file='configs/logger_analysis.cfg',
                    pipeline_config_file=pipeline_config_file,
                    )

def analyze_one():
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    analyze_onefile(pipeline=AnalysisPipeline,
                    framework_config_file="configs/framework_analysis.cfg",
                    framework_logcfg_file='configs/logger_analysis.cfg',
                    pipeline_config_file=pipeline_config_file,
                    )

def analyze_all():
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    analyze_directory(pipeline=AnalysisPipeline,
                      framework_config_file="configs/framework_analysis.cfg",
                      framework_logcfg_file='configs/logger_analysis.cfg',
                      pipeline_config_file=pipeline_config_file,
                      )

def analyze_list(framework_config_file="configs/framework_analysis.cfg"):
    list_queue(framework_config_file=framework_config_file)

def analyze_clear(framework_config_file="configs/framework_analysis.cfg"):
    clear_queue(framework_config_file=framework_config_file)

def analyze_cd():
    webcfg, cfgs = get_all_configs()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    if len(cfgs.keys()) != 0:
        if 'primary' in cfgs.keys():
            pipeline_cfg = cfgs['primary']
    args = _parseArguments(sys.argv)
    if args.input == '':
        now = datetime.utcnow()
        args.input = pipeline_cfg['FileHandling'].get('destination_dir')
        args.input = data_path.replace('YYYY', f'{now.year:4d}')
        args.input = data_path.replace('MM', f'{now.month:02d}')
        args.input = data_path.replace('DD', f'{now.day:02d}')
    change_directory(pipeline=AnalysisPipeline,
                     framework_config_file="configs/framework_analysis.cfg",
                     framework_logcfg_file='configs/logger_analysis.cfg',
                     pipeline_config_file=pipeline_config_file,
                     )


if __name__ == "__main__":
    webcfg = get_webpage_config()
    pipeline_config_file = webcfg['Telescopes'].get('pipeline_config_files').split(',')[0]
    examine_memory_use(pipeline=AnalysisPipeline,
                       framework_config_file="configs/framework_analysis.cfg",
                       framework_logcfg_file='configs/logger_analysis.cfg',
                       pipeline_config_file=pipeline_config_file,
                       )
