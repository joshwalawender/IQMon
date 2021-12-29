from pathlib import Path
import sys
import os
import logging
import argparse
from datetime import datetime
from time import sleep
import pymongo
import requests
import configparser


##-------------------------------------------------------------------------
## Parse Command Line Arguments
##-------------------------------------------------------------------------
## create a parser object for understanding command-line arguments
p = argparse.ArgumentParser(description='''
''')
## add options
p.add_argument("-c", "--config", dest="config", type=str, default=None,
               help="The config file to use (full path).")
args = p.parse_args()


##-------------------------------------------------------------------------
## Query AAG Solo for Weather Data
##-------------------------------------------------------------------------
def get_aagsolo_once():
    devicename = 'AAGSolo'

    log = logging.getLogger('AAG')
    if len(log.handlers) < 1:
        log.setLevel(logging.DEBUG)
        ## Set up console output
        LogConsoleHandler = logging.StreamHandler()
        LogConsoleHandler.setLevel(logging.INFO)
        LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
        LogConsoleHandler.setFormatter(LogFormat)
        log.addHandler(LogConsoleHandler)

    # Read Config File
    cfg = configparser.ConfigParser()
    if args.config is not None:
        cfg_file = Path(args.config)
        try:
            cfg.read(cfg_path)
        except Exception as e:
            log.error(f"Could not read config file: {cfg_path}")
            log.error(e)
            args.config = None
    if args.config is None:
        cfg_name='pipeline.cfg'
        cfg_path = Path(__file__).absolute().parent.parent/'configs'/cfg_name
        cfg.read(cfg_path)

    log.info('Getting Weather status')
    IP = cfg[devicename].get('address', None)
    if IP is None:
        return
    querydate = datetime.utcnow()
    address = f'http://{IP}/cgi-bin/cgiLastData'
    try:
        r = requests.get(address)
    except:
        log.error(f'Failed to connect to {devicename}')
        return
    lines = r.text.splitlines()
    result = {}
    for line in lines:
        key, val = line.split('=')
        result[str(key)] = str(val)
        log.debug('  {} = {}'.format(key, val))
    log.info(f'  Got {len(result.keys())} results')
    mongodoc = {"date": datetime.strptime(result['dataGMTTime'], '%Y/%m/%d %H:%M:%S'),
                "querydate": querydate,
                "cwinfo": result['cwinfo'],
                "dew point": float(result['dewp']),
                "cloud value": float(result['clouds']),
                "outside temperature": float(result['temp']),
                "wind value": float(result['wind']),
                "gust value": float(result['gust']),
                "rain value": int(result['rain']),
                "light value": int(result['light']),
                "switch": int(result['switch']),
                "safe": {'1': True, '0': False}[result['safe']],
               }
    age = (mongodoc["querydate"] - mongodoc["date"]).total_seconds()
    if len(result.keys()) != len(mongodoc.keys()) or len(result.keys()) != 12:
        log.warning(f'Possible missing keys')
        log.info(result)
        log.info(f'  Prepared mongodoc with {len(mongodoc.keys())} keys')
        log.info(mongodoc)
    if age > 60:
        log.warning(f'Data age = {age:.1f} seconds')

    log.info(f'Connecting to mongoDB')
    mongo_host = cfg['mongo'].get('host')
    mongo_port = cfg['mongo'].getint('port')
    mongo_db = cfg['mongo'].get('db')
    mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
    collection = mongoclient[mongo_db][devicename]

    try:
        inserted_id = collection.insert_one(mongodoc).inserted_id
        log.info("  Inserted document with id: {}".format(inserted_id))
    except:
        e = sys.exc_info()[0]
        log.error('Failed to add new document')
        log.error(e)
    mongoclient.close()
    logging.shutdown()


def monitor_aag():
    devicename = 'AAGSolo'
    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)
    while True:
        get_aagsolo_once()
        sleep(cfg[devicename].getfloat('sleep', 60))


if __name__ == '__main__':
    monitor_aag()
