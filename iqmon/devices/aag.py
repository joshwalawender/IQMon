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
    log = logging.getLogger('AAG')
    if len(log.handlers) < 1:
        log.setLevel(logging.DEBUG)
        ## Set up console output
        LogConsoleHandler = logging.StreamHandler()
        LogConsoleHandler.setLevel(logging.INFO)
        LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                                      datefmt='%Y%m%d %H:%M:%S')
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
        cfg_name='pipeline_hokuula.cfg'
        cfg_path = Path(__file__).absolute().parent.parent/'configs'/cfg_name
        cfg.read(cfg_path)

    log.info('Getting Weather status')
    IP = cfg['AAGSolo'].get('address', None)
    if IP is None:
        return
    querydate = dt.utcnow()
    address = f'http://{IP}/cgi-bin/cgiLastData'
    try:
        r = requests.get(address)
    except:
        log.error('Failed to connect to AAG Solo')
        return
    lines = r.text.splitlines()
    result = {}
    for line in lines:
        key, val = line.split('=')
        result[str(key)] = str(val)
        log.debug('  {} = {}'.format(key, val))
    log.info('  Done.')
    mongodoc = {"AAG date": dt.strptime(result['dataGMTTime'], '%Y/%m/%d %H:%M:%S'),
                "AAG querydate": querydate,
                "AAG cloud value": float(result['clouds']),
                "AAG outside temperature": float(result['temp']),
                "AAG wind": float(result['wind']),
                "AAG gust": float(result['gust']),
                "AAG rain value": int(result['rain']),
                "AAG light value": int(result['light']),
                "AAG switch": int(result['switch']),
                "AAG safe": {'1': True, '0': False}[result['safe']],
               }
    age = (weatherdoc["querydate"] - weatherdoc["date"]).total_seconds()
    log.debug('  Data age = {:.1f} seconds'.format(age))

    log.debug(f'Connecting to mongoDB')

    mongo_host = cfg['mongo'].get('host')
    mongo_port = cfg['mongo'].getint('port')
    mongo_db = cfg['mongo'].get('db')
    collection_name = 'AAGSolo'
    mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
    collection = mongoclient[mongo_db][collection_name]

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
    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)
    while True:
        get_aagsolo_once()
        time.sleep(cfg['AAGSolo'].getfloat('sleep', 60))


if __name__ == '__main__':
    monitor_aag()
