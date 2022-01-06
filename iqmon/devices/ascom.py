from pathlib import Path
import logging
import time
from datetime import datetime
import pymongo
import numpy as np

import win32com.client
import pywintypes

from iqmon import get_webpage_config
from iqmon.devices import insert_mongodoc


##-------------------------------------------------------------------------
## get_telescope
##-------------------------------------------------------------------------
def get_telescope(devicename):
    log = logging.getLogger('telescope')
    if len(log.handlers) < 1:
        log.setLevel(logging.DEBUG)
        ## Set up console output
        LogConsoleHandler = logging.StreamHandler()
        LogConsoleHandler.setLevel(logging.INFO)
        LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
        LogConsoleHandler.setFormatter(LogFormat)
        log.addHandler(LogConsoleHandler)

    mongodoc = {'date': datetime.utcnow(),
                'source': f'ASCOM:{devicename}'}

    log.info(f'Connecting to ASCOM object: {devicename}')
    try:
        ASCOMtel = win32com.client.Dispatch(devicename)
    except:
        log.error(f'Could not connect to ASCOM object: {devicename}')
        return

    try:
        mongodoc['connected'] = ASCOMtel.Connected,
        log.debug(f'  Telescope Connected = {mongodoc["connected"]}')
        if mongodoc['connected']:
            mongodoc['park'] = ASCOMtel.AtPark
            mongodoc['slewing'] = ASCOMtel.Slewing
            mongodoc['tracking'] = ASCOMtel.Tracking
            mongodoc['alt'] = float(ASCOMtel.Altitude)
            mongodoc['az'] = float(ASCOMtel.Azimuth)
    except pywintypes.com_error as err:
        log.warning('COM error:')
        log.warning(err)
        mongodoc['err'] = str(err)
    except:
        mongodoc['connected'] = False

    log.info(f'  Got {len(mongodoc)} entries')
    log.info(mongodoc)
    insert_mongodoc('telescope', mongodoc, log=log)
    logging.shutdown()


##-------------------------------------------------------------------------
## get_focuser
##-------------------------------------------------------------------------
def get_focuser(devicename, temperature_unit='C'):
    log = logging.getLogger('focuser')
    if len(log.handlers) < 1:
        log.setLevel(logging.DEBUG)
        ## Set up console output
        LogConsoleHandler = logging.StreamHandler()
        LogConsoleHandler.setLevel(logging.DEBUG)
        LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
        LogConsoleHandler.setFormatter(LogFormat)
        log.addHandler(LogConsoleHandler)

    mongodoc = {'date': datetime.utcnow(),
                'source': f'ASCOM:{devicename}'}

    log.info(f'Connecting to ASCOM object: {devicename}')
    try:
        focuser = win32com.client.Dispatch(devicename)
        if not focuser.Link:
            try:
                focuser.Link = True
            except:
                log.error('Could not start focuser ASCOM link.')
        log.debug('  Connected to focuser')
    except:
        log.error('Could not connect to focuser ASCOM object.')
        return

    ## Get Average of 3 Temperature Readings
    focuser_temps = []
    for i in range(0,3,1):
        try:
            newtemp = float(focuser.Temperature)
            log.debug('  Queried focuser temperature = {newtemp:.1f}')
            focuser_temps.append(newtemp)
        except Exception as e:
            log.warning('Failed to get focuser temperature')
            log.warning(e)
    if len(focuser_temps) > 0:
        ## Filter out bad values
        focuser_temp = np.median(focuser_temps)
        if (focuser_temp > -20) and (focuser_temp < 150):
            log.debug('  focuser temperature = {focuser_temp:.1f}')
            if temperature_unit == 'C':
                log.debug('  Converting C to F')
                mongodoc['temperature'] = focuser_temp*1.8+32
                mongodoc['temperature unit'] = 'F'
            else:
                mongodoc['temperature'] = focuser_temp
                mongodoc['temperature unit'] = temperature_unit

    ## Get Position
    try:
        mongodoc['position'] = int(focuser.Position)
        log.debug('  focuser position = {mongodoc["position"]:d}')
    except Exception as e:
        log.warning('Failed to get focuser position')
        log.warning(e)

    log.info(f'  Got {len(mongodoc)} entries')
    log.info(mongodoc)
    insert_mongodoc('focuser', mongodoc, log=log)
    logging.shutdown()


##-------------------------------------------------------------------------
## get_dome
##-------------------------------------------------------------------------
def get_dome(devicename):
    log.info(f'Connecting to ASCOM object: {devicename}')
    pass


##-------------------------------------------------------------------------
## poll_ASCOM_devices
##-------------------------------------------------------------------------
def poll_ASCOM_devices():
    webcfg = get_webpage_config()
    sleeptime = webcfg['devices'].getint('polling_time', 30)
    devfunctions = {'telescope': get_telescope,
                    'focuser': get_focuser,
                    'dome': get_dome,
                    }
    while True:
        for device in ['telescope', 'focuser', 'dome']:
            deviceinfostring = webcfg['devices'].get(device, None)
            if deviceinfostring is not None:
                deviceinfo = deviceinfostring.split(',')
                if len(deviceinfo) < 2 and deviceinfo[0] == 'ASCOM':
                    raise Exception(f'Device specification incomplete: {deviceinfostring}')
                elif len(deviceinfo) == 2 and deviceinfo[0] == 'ASCOM':
                    devfunctions[device](deviceinfo[1])
                elif len(deviceinfo) > 2 and deviceinfo[0] == 'ASCOM':
                    for entry in deviceinfo[2:]:
                        deviceargs = {}
                        try:
                            key, val = entry.split(':')
                            deviceargs[key] = val
                        except:
                            raise Exception(f'Device arguments not parsed: {deviceinfostring}')
                    devfunctions[device](deviceinfo[1], **deviceargs)
        time.sleep(sleeptime)


if __name__ == '__main__':
    poll_ASCOM_devices()
