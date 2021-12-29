from pathlib import Path
import configparser


##-------------------------------------------------------------------------
## Function: get_config
##-------------------------------------------------------------------------
def get_config(config_file='~/iqmon.cfg'):
    cfg = configparser.ConfigParser()
    if config_file is not None:
        cfg_path = Path(config_file).expanduser()
        try:
            cfg.read(cfg_path)
        except:
            config_file = None
    if config_file is None:
        cfg_name='pipeline.cfg'
        cfg_path = Path(__file__).absolute().parent/'configs'/cfg_name
        cfg.read(cfg_path)
    return cfg
