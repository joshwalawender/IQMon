from pathlib import Path
import configparser


##-------------------------------------------------------------------------
## Function: get_webpage_config
##-------------------------------------------------------------------------
def get_webpage_config():
    webcfg = configparser.ConfigParser()
    files_to_try = ['~/.iqmon_webpage.cfg', '~/iqmon_webpage.cfg',
                    Path(__file__).absolute().parent/'configs'/'webpage.cfg',
                    ]
    for file in files_to_try:
        webcfg_path = Path(file).expanduser()
        try:
            webcfg.read(webcfg_path, encoding='utf8')
        except:
            pass
        if webcfg.sections() != []:
            return webcfg
    return None


##-------------------------------------------------------------------------
## Function: get_pipeline_config
##-------------------------------------------------------------------------
def get_pipeline_config(config_file='~/.iqmon_pipeline.cfg'):
    cfg = configparser.ConfigParser()
    if config_file is not None:
        cfg_path = Path(config_file).expanduser()
        try:
            cfg.read(cfg_path, encoding='utf8')
        except:
            config_file = None
    if config_file is None:
        cfg_name='pipeline.cfg'
        cfg_path = Path(__file__).absolute().parent/'configs'/cfg_name
        cfg.read(cfg_path)
    return cfg


##-------------------------------------------------------------------------
## Function: get_all_configs
##-------------------------------------------------------------------------
def get_all_configs():
    webcfg = get_webpage_config(config_file=config_file)

    cfgs = {}
    pipeline_config_files = webcfg['Telescopes'].get('pipeline_config_files', None)
    if pipeline_config_files is not None:
        for i,pipeline_config_file in enumerate(pipeline_config_files.split(',')):
            pipeline_config_file = Path(pipeline_config_file).expanduser()
            if pipeline_config_file.exists() is True:
                cfg = get_pipeline_config(config_file=pipeline_config_file)
                name = cfg['Telescope'].get('name', None)
                if name is not None:
                    cfgs[name] = cfg
                    if i == 0:
                        cfgs['primary'] = name
    
    return webcfg, cfgs
