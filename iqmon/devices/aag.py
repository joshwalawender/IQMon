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

from iqmon import get_webpage_config
from iqmon.devices import insert_mongodoc


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
    cfg = get_webpage_config()

    log = logging.getLogger(devicename)
    if len(log.handlers) < 1:
        log.setLevel(logging.DEBUG)
        ## Set up console output
        LogConsoleHandler = logging.StreamHandler()
        LogConsoleHandler.setLevel(logging.INFO)
        LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
        LogConsoleHandler.setFormatter(LogFormat)
        log.addHandler(LogConsoleHandler)

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
                "outside temperature": float(result['temp']),
                "cloud value": float(result['clouds']),
                "wind speed": float(result['wind']),
                "gust speed": float(result['gust']),
                "rain value": int(result['rain']),
                "light value": int(result['light']),
                "switch": int(result['switch']),
                "safe": {'1': True, '0': False}[result['safe']],
               }
    temperature_units = cfg[devicename].get('temperature_units', None)
    if temperature_units is not None:
        mongodoc['temperature units'] = temperature_units
    wind_speed_units = cfg[devicename].get('wind_speed_units', None)
    if wind_speed_units is not None:
        mongodoc['wind speed units'] = wind_speed_units

    age = (mongodoc["querydate"] - mongodoc["date"]).total_seconds()
    if len(result.keys()) != len(mongodoc.keys())-1 or len(result.keys()) != 12:
        log.warning(f'Possible missing keys')
        log.info(result)
        log.info(f'  Prepared mongodoc with {len(mongodoc.keys())} keys')
        log.info(mongodoc)
    if age > 60:
        log.warning(f'Data age = {age:.1f} seconds')

    insert_mongodoc(devicename, mongodoc, log=log)
    logging.shutdown()


def monitor_aag():
    devicename = 'AAGSolo'
    cfg = get_webpage_config()
    while True:
        get_aagsolo_once()
        sleep(cfg[devicename].getfloat('sleep', 60))


if __name__ == '__main__':
    monitor_aag()
