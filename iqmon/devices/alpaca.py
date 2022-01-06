import sys
import time
import requests
import json
import logging
from datetime import datetime, timedelta
import pymongo


from iqmon import get_webpage_config
from iqmon.devices import insert_mongodoc


##-------------------------------------------------------------------------
## get_alpaca
##-------------------------------------------------------------------------
def get_alpaca(address, devicetype, devicenumber=0):
    log = logging.getLogger(devicetype)
    if len(log.handlers) < 1:
        log.setLevel(logging.DEBUG)
        ## Set up console output
        LogConsoleHandler = logging.StreamHandler()
        LogConsoleHandler.setLevel(logging.DEBUG)
        LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
        LogConsoleHandler.setFormatter(LogFormat)
        log.addHandler(LogConsoleHandler)

    url = f"{address}/api/v1/{devicetype}/{devicenumber}"
    if url[0:7] != 'http://':
        url = 'http://' + url
    mongodoc = {'date': datetime.utcnow(),
                'source': f'ALPACA:{url}'}

    commands = {'dome': [('connected', bool),
                         ('shutterstatus', int),
                         ('atpark', bool),
                         ('athome', bool),
                         ('azimuth', float),
                         ('slaved', bool),
                         ('slewing', bool)],
                'telescope': [],
                'focuser': [],
                }

    log.info(f'Getting {devicetype} info')
    for command in commands[devicetype]:
        try:
            log.info(f'  Getting {url}/{command[0]}')
            r = requests.get(f"{url}/{command[0]}")
            j = json.loads(r.text)
            mongodoc[command[0]] = command[1](j['Value'])
        except:
            log.warning(f'  Failed to get {command[0]}')
        else:
            log.info(f'  Got {devicetype} {command[0]} = {mongodoc[command[0]]}')

    log.info(f'Got {len(mongodoc)} entries')
    insert_mongodoc(devicetype, mongodoc, log=log)
    logging.shutdown()


##-------------------------------------------------------------------------
## poll_ASCOM_devices
##-------------------------------------------------------------------------
def poll_ALPACA_devices():
    webcfg = get_webpage_config()
    sleeptime = webcfg['devices'].getint('polling_time', 30)
    while True:
        for device in ['telescope', 'focuser', 'dome']:
            deviceinfostring = webcfg['devices'].get(device, None)
            if deviceinfostring is not None:
                deviceinfo = deviceinfostring.split(',')
                if len(deviceinfo) < 2 and deviceinfo[0] == 'ALPACA':
                    raise Exception(f'Device specification incomplete: {deviceinfostring}')
                elif len(deviceinfo) == 2 and deviceinfo[0] == 'ALPACA':
                    address = deviceinfo[1]
                    get_alpaca(address, device, devicenumber=0)
                elif len(deviceinfo) > 2 and deviceinfo[0] == 'ALPACA':
                    for entry in deviceinfo[2:]:
                        deviceargs = {}
                        try:
                            key, val = entry.split(':')
                            deviceargs[key] = val
                        except:
                            raise Exception(f'Device arguments not parsed: {deviceinfostring}')
                    address = deviceinfo[1]
                    get_alpaca(address, device, **deviceargs)
        time.sleep(sleeptime)


if __name__ == '__main__':
    poll_ALPACA_devices()
