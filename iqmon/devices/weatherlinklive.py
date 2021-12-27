## Import General Tools
from pathlib import Path
from datetime import datetime, timedelta
import configparser
import logging
import requests
import json
import pymongo

import socket
import struct
import time


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
## DavisWeatherLinkLive
##-------------------------------------------------------------------------
class DavisWeatherLinkLive():
    def __init__(self, IP='192.168.4.76',
                 mongoIP='192.168.4.49', mongoport=32768,
                 dbname='weather', collectionname='DavisWeatherLink'):

        self.log = logging.getLogger('DavisWeatherLink')
        if len(self.log.handlers) < 1:
            self.log.setLevel(logging.DEBUG)
            ## Set up console output
            LogConsoleHandler = logging.StreamHandler()
            LogConsoleHandler.setLevel(logging.DEBUG)
            LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                                          datefmt='%Y%m%d %H:%M:%S')
            LogConsoleHandler.setFormatter(LogFormat)
            self.log.addHandler(LogConsoleHandler)

        self.IP = IP
        self.name = 'DavisWeatherLinkLive'
        self.url = f'http://{self.IP}/v1/current_conditions'
        # Mongo Setup
        self.mongoIP = mongoIP
        self.mongoport = mongoport
        self.dbname = dbname
        self.collectionname = collectionname
        self.client = None
        self.db = None
        self.collection = None
        if mongoIP is not None:
            try:
                self.client = pymongo.MongoClient(mongoIP, mongoport)
                self.db = self.client[dbname]
                self.collection = self.db[collectionname]
                self.log.info(f'Connected to mongoDB at {self.mongoIP}:{self.mongoport}')
            except Exception as err:
                self.log.error(f'ERROR: failed to connect to mongoDB at {self.mongoIP}:{self.mongoport}')
                self.log.error(err)


    def get_data(self):
        try:
            r = requests.get(self.url)
        except Exception as err:
            self.log.error(f'Exception getting data from WeatherLinkLive:')
            self.log.error(err)
            return {}

        result = json.loads(r.text)
        data = result.get('data', [])
        error = result.get('error', None)
        conditions = data.get('conditions', None)
        if type(conditions) == list:
            conditions = conditions[0]
        if error is None:
            self.log.info(f"Got {len(conditions.keys())} data points from Weather Link Live")
        else:
            self.log.error(f'Got error from WeatherLinkLive at {data["timestamp"]}')
            self.log.error(error)

        mongodata = {'DWLL date': datetime.fromtimestamp(data.get('ts'))}
        keys = [('temp', 'DWLL outside temperature', float),
                ('hum', 'DWLL outside humidity', float),
                ('dew_point', 'DWLL outside dew point', float),
                ('wet_bulb', 'DWLL wet bulb', float),
                ('heat_index', 'DWLL heat index', float),
                ('wind_chill', 'DWLL wind chill', float),
                ('wind_speed_last', 'DWLL wind speed', float),
                ('wind_dir_last', 'DWLL wind direction', float),
                ('wind_speed_avg_last_10_min', 'DWLL wind speed (10 min avg)', float),
                ('wind_dir_scalar_avg_last_10_min', 'DWLL wind direction (10 min avg)', float),
                ('wind_speed_hi_last_10_min', 'DWLL wind speed (10 min high)', float),
                ('wind_dir_at_hi_speed_last_10_min', 'DWLL wind direction (at 10 min high)', float),
                ('rain_size', 'DWLL rain_size', float),
                ('rain_rate_last', 'DWLL rain rate', float),
                ('rain_rate_hi', 'DWLL rain rate hi', float),
                ('rainfall_last_15_min', 'DWLL rainfall last 15 min', float),
                ('rain_rate_hi_last_15_min', 'DWLL rain rate hi last 15 min', float),
                ('rainfall_last_60_min', 'DWLL rainfall last 1 hour', float),
                ('rainfall_last_24_hr', 'DWLL rainfall last 24 hours', float),
                ('rain_storm', 'DWLL rain storm total', float),
                ('rain_storm_start_at', 'DWLL rain storm start', datetime),
                ('temp_in', 'DWLL inside temperature', float),
                ('hum_in', 'DWLL inside humidity', float),
                ('bar_sea_level', 'DWLL bar sea level', float),
                ('bar_trend', 'DWLL bar trend', float),
                ('bar_absolute', 'DWLL bar absolute', float),
                ]
        for key, mongokey, datatype in keys:
            if conditions.get(key, None) is not None:
                mongodata[mongokey] = datatype(conditions.get(key, None))
        self.log.info(f"Parsed {len(mongodata.keys())} data points")

        return mongodata


    def poll(self, sleep=60):
        while True:
            self.log.info('Starting poll of weather data')
            try:
                data = self.get_data()
            except Exception as err:
                self.log.error(f'Error polling WeatherLinkLive:')
                self.log.error(err)
            if self.collection is not None:
                try:
                    inserted_id = self.collection.insert_one(data).inserted_id
                    self.log.info(f"  Inserted mongo document {inserted_id}")
                except Exception as e:
                    self.log.error(f"  Failed to insert mongo document")
                    self.log.error(e)
            self.log.info(f"Sleeping {sleep:.0f} s")
            time.sleep(sleep)


def get_davis_weather_link():

    # Read Config File
    cfg = configparser.ConfigParser()
    if args.config is not None:
        cfg_file = Path(args.config)
        try:
            cfg.read(cfg_path)
        except:
            print(f"Could not read config file: {cfg_path}")
            print(e)
            args.config = None
    if args.config is None:
        cfg_name='pipeline_hokuula.cfg'
        cfg_path = Path(__file__).absolute().parent.parent/'configs'/cfg_name
        cfg.read(cfg_path)

    sleeptime = cfg['DavisWeatherLink'].getfloat('sleep', 60)
    IP = cfg['DavisWeatherLink'].get('address', None)
    if IP is not None:
        d = DavisWeatherLinkLive(IP=IP,
                                 mongoIP=cfg['mongo'].get('host'),
                                 mongoport=cfg['mongo'].getint('port'),
                                 dbname=cfg['mongo'].get('db'),
                                 )
        d.poll(sleep=sleeptime)


if __name__ == '__main__':
    get_davis_weather_link()