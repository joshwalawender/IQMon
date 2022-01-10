## Import General Tools
from pathlib import Path
from datetime import datetime, timedelta
import configparser
import logging
import argparse
import requests
import json
import pymongo

import socket
import struct
import time

from iqmon import get_webpage_config


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
## DavisWeatherLink
##-------------------------------------------------------------------------
class DavisWeatherLink():
    def __init__(self, IP='192.168.4.76',
                 mongoIP='192.168.4.49', mongoport=49153,
                 dbname='weather', tzoffset=0,
                 temperature_units='F', wind_speed_units='mph'):
        self.IP = IP
        self.name = 'DavisWeatherLink'
        self.url = f'http://{self.IP}/v1/current_conditions'
        self.tzoffset = tzoffset
        self.temperature_units = temperature_units
        self.wind_speed_units = wind_speed_units

        self.log = logging.getLogger(self.name)
        if len(self.log.handlers) < 1:
            self.log.setLevel(logging.DEBUG)
            ## Set up console output
            LogConsoleHandler = logging.StreamHandler()
            LogConsoleHandler.setLevel(logging.DEBUG)
            LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                                          datefmt='%Y%m%d %H:%M:%S')
            LogConsoleHandler.setFormatter(LogFormat)
            self.log.addHandler(LogConsoleHandler)

        # Mongo Setup
        self.mongoIP = mongoIP
        self.mongoport = mongoport
        self.dbname = dbname
        self.client = None
        self.db = None
        self.collection = None
        if mongoIP is not None:
            try:
                self.client = pymongo.MongoClient(mongoIP, mongoport)
                self.db = self.client[dbname]
                self.collection = self.db[self.name]
                self.log.info(f'Connected to mongoDB at {self.mongoIP}:{self.mongoport}')
            except Exception as err:
                self.log.error(f'ERROR: failed to connect to mongoDB at {self.mongoIP}:{self.mongoport}')
                self.log.error(err)


    def get_data(self):
        try:
            r = requests.get(self.url)
        except Exception as err:
            self.log.error(f'Exception getting data from {self.name}')
            self.log.error(err)
            return {}

        result = json.loads(r.text)
        data = result.get('data', [])
        error = result.get('error', None)
        conditions_list = data.get('conditions', None)
        conditions = {}
        if type(conditions_list) == list:
            for entry in conditions_list:
                conditions.update(entry)
        else:
            conditions.update(conditions_list)

        if error is None:
            self.log.info(f"Got {len(conditions.keys())} data points from {self.name}")
        else:
            self.log.error(f'Got error from {self.name} at {data["timestamp"]}')
            self.log.error(error)

        mongodata = {'date': datetime.fromtimestamp(data.get('ts'))\
                             - timedelta(hours=self.tzoffset)}
        keys = [('temp', 'outside temperature', float),
                ('hum', 'outside humidity', float),
                ('dew_point', 'outside dew point', float),
                ('wet_bulb', 'wet bulb', float),
                ('heat_index', 'heat index', float),
                ('wind_chill', 'wind chill', float),
                ('wind_speed_last', 'wind speed', float),
                ('wind_dir_last', 'wind direction', float),
                ('wind_speed_avg_last_10_min', 'wind speed (10 min avg)', float),
                ('wind_dir_scalar_avg_last_10_min', 'wind direction (10 min avg)', float),
                ('wind_speed_hi_last_10_min', 'wind gust', float),
                ('wind_dir_at_hi_speed_last_10_min', 'wind gust direction', float),
                ('temp_in', 'inside temperature', float),
                ('hum_in', 'inside humidity', float),
                ('bar_sea_level', 'bar sea level', float),
                ('bar_trend', 'bar trend', float),
                ('bar_absolute', 'bar absolute', float),
                ]
        for key, mongokey, datatype in keys:
            if conditions.get(key, None) is not None:
                mongodata[mongokey] = datatype(conditions.get(key, None))
        # Handle rain units
        rain_size_dict = {1: 0.01, 2: 0.2, 3:  0.1, 4: 0.001}
        rain_size_units = {1: 'in', 2: 'mm', 3: 'mm', 4: 'in'}
        mongodata['rain size'] = rain_size_dict[int(conditions.get('rain_size'))]
        mongodata['rain units'] = rain_size_units[int(conditions.get('rain_size'))]
        mongodata['rain rate'] = conditions.get('rain_rate_last')*mongodata['rain size']
        mongodata['rain rate hi'] = conditions.get('rain_rate_hi')*mongodata['rain size']
        mongodata['rainfall last 15 min'] = conditions.get('rainfall_last_15_min')*mongodata['rain size']
        mongodata['rain rate hi last 15 min'] = conditions.get('rain_rate_hi_last_15_min')*mongodata['rain size']
        mongodata['rainfall last 1 hour'] = conditions.get('rainfall_last_60_min')*mongodata['rain size']
        mongodata['rainfall last 24 hours'] = conditions.get('rainfall_last_24_hr')*mongodata['rain size']
        mongodata['rain storm total'] = conditions.get('rain_storm')*mongodata['rain size']
        if conditions.get('rain_storm_start_at') is not None:
            mongodata['rain storm start'] = datetime.fromtimestamp(conditions.get('rain_storm_start_at'))\
                                            - timedelta(hours=self.tzoffset)

        self.log.info(f"Parsed {len(mongodata.keys())} data points")
        # Add units
        mongodata['temperature units'] = self.temperature_units
        mongodata['wind speed units'] = self.wind_speed_units

        return mongodata


    def poll(self, sleep=60):
        while True:
            self.log.info('Starting poll of weather data')
            try:
                data = self.get_data()
            except Exception as err:
                self.log.error(f'Error polling {self.name}')
                self.log.error(err)
                data = None
            if self.collection is not None and data is not None:
                try:
                    inserted_id = self.collection.insert_one(data).inserted_id
                    self.log.info(f"  Inserted mongo document {inserted_id}")
                except Exception as e:
                    self.log.error(f"  Failed to insert mongo document")
                    self.log.error(e)
            self.log.info(f"Sleeping {sleep:.0f} s")
            time.sleep(sleep)


def monitor_davis_weather_link():
    devicename = 'DavisWeatherLink'
    cfg = get_webpage_config()

    sleeptime = cfg[devicename].getfloat('sleep', 60)
    tzoffset = cfg[devicename].getfloat('TZoffset', 0)
    temperature_units = cfg[devicename].get('temperature_units', 'F')
    IP = cfg[devicename].get('address', None)
    if IP is not None:
        d = DavisWeatherLink(IP=IP,
                             mongoIP=cfg['mongo'].get('host'),
                             mongoport=cfg['mongo'].getint('port'),
                             dbname=cfg['mongo'].get('db'),
                             tzoffset=tzoffset,
                             temperature_units=temperature_units,
                             )
        d.poll(sleep=sleeptime)


if __name__ == '__main__':
    monitor_davis_weather_link()