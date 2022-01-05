import requests
import json


##-------------------------------------------------------------------------
## get_alpaca
##-------------------------------------------------------------------------
def get_alpaca(address, devicetype, devicenumber=0):
    log = logging.getLogger(devicetype)
    if len(log.handlers) < 1:
        log.setLevel(logging.DEBUG)
        ## Set up console output
        LogConsoleHandler = logging.StreamHandler()
        LogConsoleHandler.setLevel(logging.INFO)
        LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
        LogConsoleHandler.setFormatter(LogFormat)
        log.addHandler(LogConsoleHandler)

    url = f"{address}/api/v1/{devicetype}/{devicenumber}"
    mongodoc = {'date': datetime.utcnow(),
                'source': f'ALPACA:{url}'}

    commands = {'dome': [('connected', str),
                         ('shutterstatus', str),
                         ('atpark', str),
                         ('athome', str),
                         ('azimuth', float),
                         ('slaved', str),
                         ('slewing', str)],
                'telescope': [],
                'focuser': [],
                }

    commands = 
    log.info(f'Getting {devicetype} Status')
    for command in commands:
        try:
            r = requests.get(f"{url}/{command[0]}")
            j = json.loads(r.text)
            mongodoc[command] = command[1](j['Value'])
        except:
            logger.warning(f'  Failed to get {devicetype} status: {command}')
        else:
            logger.debug(f'  {command} = {mongodoc[command]}')

    log.info(f'  Got {len(mongodoc)} entries')
    insert_mongodoc(devicetype, mongodoc, log=log)
    logging.shutdown()


##-------------------------------------------------------------------------
## poll_ASCOM_devices
##-------------------------------------------------------------------------
def poll_ASCOM_devices():
    webcfg = get_webpage_config()
    
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


