from __future__ import division, print_function

## Import General Tools
import sys
import os
import argparse
import logging
from datetime import datetime as dt
from datetime import timedelta as tdelta
import re
import glob

import pymongo
from pymongo import MongoClient

from tornado.web import RequestHandler, Application, url, StaticFileHandler
import tornado.log as tlog

from astropy import units as u
from astropy.coordinates import SkyCoord
import ephem

import IQMon


##-------------------------------------------------------------------------
## Check for Images
##-------------------------------------------------------------------------
def get_nimages(telescope, date):
    path = os.path.join('/Volumes/Data_{}/Images/{}'.format(telescope, date))
    image_list = glob.glob(os.path.join(path, '{}*fts'.format(telescope)))
    return len(image_list)

##-------------------------------------------------------------------------
## Check Free Space on Drive
##-------------------------------------------------------------------------
def free_space(path):
    statvfs = os.statvfs(path)
    size_GB = statvfs.f_frsize * statvfs.f_blocks / 1024 / 1024 / 1024
    avail_GB = statvfs.f_frsize * statvfs.f_bfree / 1024 / 1024 / 1024
    pcnt_used = float(size_GB - avail_GB)/float(size_GB) * 100
    return (size_GB, avail_GB, pcnt_used)


##-----------------------------------------------------------------------------
## Handler for Status Page
##-----------------------------------------------------------------------------
class Status(RequestHandler):
    def get(self):
        tlog.app_log.info('Get request for Status recieved')
        nowut = dt.utcnow()
        now = nowut - tdelta(0,10*60*60)

        client = MongoClient('192.168.1.101', 27017)

        ##------------------------------------------------------------------------
        ## Use pyephem determine sunrise and sunset times
        ##------------------------------------------------------------------------
        Observatory = ephem.Observer()
        Observatory.lon = "-155:34:33.9"
        Observatory.lat = "+19:32:09.66"
        Observatory.elevation = 3400.0
        Observatory.temp = 10.0
        Observatory.pressure = 680.0
        Observatory.horizon = '0.0'

        Observatory.date = nowut
        TheSun = ephem.Sun()
        TheSun.compute(Observatory)
        sun = {}
        sun['alt'] = float(TheSun.alt) * 180. / ephem.pi
        sun['set'] = Observatory.next_setting(TheSun).datetime()
        sun['rise'] = Observatory.next_rising(TheSun).datetime()
        if sun['alt'] <= -18:
            sun['now'] = 'night'
        elif sun['alt'] > -18 and sun['alt'] <= -12:
            sun['now'] = 'astronomical twilight'
        elif sun['alt'] > -12 and sun['alt'] <= -6:
            sun['now'] = 'nautical twilight'
        elif sun['alt'] > -6 and sun['alt'] <= 0:
            sun['now'] = 'civil twilight'
        elif sun['alt'] > 0:
            sun['now'] = 'day'

        TheMoon = ephem.Moon()
        Observatory.date = nowut
        TheMoon.compute(Observatory)
        moon = {}
        moon['phase'] = TheMoon.phase
        moon['alt'] = TheMoon.alt * 180. / ephem.pi
        moon['set'] = Observatory.next_setting(TheMoon).datetime()
        moon['rise'] = Observatory.next_rising(TheMoon).datetime()
        if moon['alt'] > 0:
            moon['now'] = 'up'
        else:
            moon['now'] = 'down'

        tlog.app_log.info('  Ephem data calculated')

        ##---------------------------------------------------------------------
        ## Get Latest V20 Data
        ##---------------------------------------------------------------------
        v20status = client.vysos['V20.status']
#         v20entries = []
#         while (len(v20entries) < 1) and (nowut > dt(2015,1,1)):
#             v20entries = [entry for entry\
#                           in v20status.find(\
#                           {"UT date" : nowut.strftime('%Y%m%dUT')}\
#                           ).sort([('UT time', pymongo.ASCENDING)])]
#             if len(v20entries) > 0: v20data = v20entries[-1]
#             else: nowut = nowut - tdelta(1, 0)
#         nowut = dt.utcnow()

        v20data = v20status.find_one( {'current': True} )
        tlog.app_log.info('  v20data retrieved')

        try:
            try:
                v20clarity_time = v20data['boltwood timestamp']
            except:
                v20clarity_time = dt.strptime('{} {}'.format(\
                                  v20data['boltwood date'],\
                                  v20data['boltwood time'][:-3]),\
                                  '%Y-%m-%d %H:%M:%S')
            v20clarity_age = (now - v20clarity_time).total_seconds()/60.
            if v20clarity_age > 1: v20clarity_color = 'red'
            else: v20clarity_color = 'black'
        except:
            v20clarity_age = float('nan')
            v20clarity_color = 'red'
        tlog.app_log.info('  v20clarity_color determined')

        try:
            try:
                v20data_time = v20data['UT timestamp']
            except:
                v20data_time = dt.strptime('{} {}'.format(\
                               v20data['UT date'],\
                               v20data['UT time']),\
                               '%Y%m%dUT %H:%M:%S')
            v20data_age = (nowut - v20data_time).total_seconds()/60.
            if v20data_age > 1: v20data_color = 'red'
            else: v20data_color = 'black'
        except:
            v20data_age = float('nan')
            v20data_color = 'red'
        tlog.app_log.info('  v20data_color determined')

        ##---------------------------------------------------------------------
        ## Get Latest V5 Data
        ##---------------------------------------------------------------------
        v5status = client.vysos['V5.status']
#         v5entries = []
#         while (len(v5entries) < 1) and (nowut > dt(2015,1,1)):
#             v5entries = [entry for entry\
#                           in v5status.find( {"UT date" : nowut.strftime('%Y%m%dUT')} ).sort([('UT time', pymongo.ASCENDING)])]
#             if len(v5entries) > 0: v5data = v5entries[-1]
#             else: nowut = nowut - tdelta(1, 0)
#         nowut = dt.utcnow()

        v5data = v5status.find_one( {'current': True} )
        tlog.app_log.info('  v5data retrieved')

        try:
            try:
                v5clarity_time = v5data['boltwood timestamp']
            except:
                v5clarity_time = dt.strptime('{} {}'.format(\
                                  v5data['boltwood date'],\
                                  v5data['boltwood time'][:-3]),\
                                  '%Y-%m-%d %H:%M:%S')        
            v5clarity_age = (now - v5clarity_time).total_seconds()/60.
            if v5clarity_age > 1.: v5clarity_color = 'red'
            else: v5clarity_color = 'black'
        except:
            v5clarity_age = float('nan')
            v5clarity_color = 'red'
        tlog.app_log.info('  v5clarity_color determined')

        try:
            try:
                v5data_time = v5data['UT timestamp']
            except:
                v5data_time = dt.strptime('{} {}'.format(\
                               v5data['UT date'],\
                               v5data['UT time']),\
                               '%Y%m%dUT %H:%M:%S')
            v5data_age = (nowut - v5data_time).total_seconds()/60.
            if v5data_age > 1: v5data_color = 'red'
            else: v5data_color = 'black'
        except:
            v5data_age = float('nan')
            v5data_color = 'red'
        tlog.app_log.info('  v5clarity_color determined')

        ##---------------------------------------------------------------------
        ## Format and Color Code Boltwood Data
        ##---------------------------------------------------------------------
        wind_units = {'M': 'mph', 'K': 'kph', 'm': 'm/s'}
        rain_status = {0: 'Dry', 1: 'Recent Rain', 2: 'Raining'}
        wet_status = {0: 'Dry', 1: 'Recent Wet', 2: 'Wet'}
        cloud_condition = {0: 'Unknown', 1: 'Clear', 2: 'Cloudy', 3: 'Very Cloudy'}
        wind_condition = {0: 'Unknown', 1: 'Calm', 2: 'Windy', 3: 'Very Windy'}
        rain_condition = {0: 'Unknown', 1: 'Dry', 2: 'Wet', 3: 'Rain'}
        day_condition = {0: 'Unknown', 1: 'Dark', 2: 'Light', 3: 'Very Light'}
        roof_close = {0: 'Safe', 1: 'Unsafe'}

        if 'boltwood wind units' in v20data.keys():
            v20data['boltwood wind units'] = wind_units[v20data['boltwood wind units']]
        if 'boltwood rain status' in v20data.keys():
            if v20data['boltwood rain status'] == 0: v20data['boltwood rain status color'] = 'green'
            elif v20data['boltwood rain status'] == 1: v20data['boltwood rain status color'] = 'red'
            elif v20data['boltwood rain status'] == 2: v20data['boltwood rain status color'] = 'red'
            else: v20data['boltwood rain color'] = ''
            v20data['boltwood rain status string'] = rain_status[v20data['boltwood rain status']]
        else: v20data['boltwood rain status color'] = 'red'
        if 'boltwood wet status' in v20data.keys():
            if v20data['boltwood wet status'] == 0: v20data['boltwood wet status color'] = 'green'
            elif v20data['boltwood wet status'] == 1: v20data['boltwood wet status color'] = 'red'
            elif v20data['boltwood wet status'] == 2: v20data['boltwood wet status color'] = 'red'
            else: v20data['boltwood wet color'] = ''
            v20data['boltwood wet status string'] = wet_status[v20data['boltwood wet status']]
        else: v20data['boltwood wet status color'] = 'red'
        if 'boltwood cloud condition' in v20data.keys():
            if v20data['boltwood cloud condition'] == 0: v20data['boltwood cloud color'] = 'orange'
            elif v20data['boltwood cloud condition'] == 1: v20data['boltwood cloud color'] = 'green'
            elif v20data['boltwood cloud condition'] == 2: v20data['boltwood cloud color'] = 'orange'
            elif v20data['boltwood cloud condition'] == 3: v20data['boltwood cloud color'] = 'red'
            else: v20data['boltwood cloud color'] = ''
            v20data['boltwood cloud condition string'] = cloud_condition[v20data['boltwood cloud condition']]
        else: v20data['boltwood cloud color'] = 'red'
        if 'boltwood wind condition' in v20data.keys():
            if v20data['boltwood wind condition'] == 0: v20data['boltwood wind color'] = 'orange'
            elif v20data['boltwood wind condition'] == 1: v20data['boltwood wind color'] = 'green'
            elif v20data['boltwood wind condition'] == 2: v20data['boltwood wind color'] = 'orange'
            elif v20data['boltwood wind condition'] == 3: v20data['boltwood wind color'] = 'red'
            else: v20data['boltwood wind color'] = ''
            v20data['boltwood wind condition string'] = wind_condition[v20data['boltwood wind condition']]
        else: v20data['boltwood wind color'] = 'red'
        if 'boltwood rain condition' in v20data.keys():
            if v20data['boltwood rain condition'] == 0: v20data['boltwood rain color'] = 'orange'
            elif v20data['boltwood rain condition'] == 1: v20data['boltwood rain color'] = 'green'
            elif v20data['boltwood rain condition'] == 2: v20data['boltwood rain color'] = 'red'
            elif v20data['boltwood rain condition'] == 3: v20data['boltwood rain color'] = 'red'
            else: v20data['boltwood rain color'] = ''
            v20data['boltwood rain condition string'] = rain_condition[v20data['boltwood rain condition']]
        else: v20data['boltwood rain color'] = 'red'
        if 'boltwood day condition' in v20data.keys():
            if v20data['boltwood day condition'] == 0: v20data['boltwood day color'] = 'orange'
            elif v20data['boltwood day condition'] == 1: v20data['boltwood day color'] = 'green'
            elif v20data['boltwood day condition'] == 2: v20data['boltwood day color'] = 'orange'
            elif v20data['boltwood day condition'] == 3: v20data['boltwood day color'] = 'red'
            else: v20data['boltwood day color'] = ''
            v20data['boltwood day condition string'] = day_condition[v20data['boltwood day condition']]
        else: v20data['boltwood day color'] = 'red'
        if 'boltwood roof close' in v20data.keys():
            if v20data['boltwood roof close'] == 0: v20data['boltwood roof close color'] = 'green'
            elif v20data['boltwood roof close'] == 1: v20data['boltwood roof close color'] = 'red'
            else: v20data['boltwood roof close color'] = ''
            v20data['boltwood roof close string'] = roof_close[v20data['boltwood roof close']]
        else: v20data['boltwood roof close color'] = 'red'

        if 'boltwood wind units' in v5data.keys():
            v5data['boltwood wind units'] = wind_units[v5data['boltwood wind units']]
        if 'boltwood rain status' in v5data.keys():
            if v5data['boltwood rain status'] == 0: v5data['boltwood rain status color'] = 'green'
            elif v5data['boltwood rain status'] == 1: v5data['boltwood rain status color'] = 'red'
            elif v5data['boltwood rain status'] == 2: v5data['boltwood rain status color'] = 'red'
            else: v5data['boltwood rain color'] = ''
            v5data['boltwood rain status string'] = rain_status[v5data['boltwood rain status']]
        else: v5data['boltwood rain status color'] = 'red'
        if 'boltwood wet status' in v5data.keys():
            if v5data['boltwood wet status'] == 0: v5data['boltwood wet status color'] = 'green'
            elif v5data['boltwood wet status'] == 1: v5data['boltwood wet status color'] = 'red'
            elif v5data['boltwood wet status'] == 2: v5data['boltwood wet status color'] = 'red'
            else: v5data['boltwood wet color'] = ''
            v5data['boltwood wet status string'] = wet_status[v5data['boltwood wet status']]
        else: v5data['boltwood wet status color'] = 'red'
        if 'boltwood cloud condition' in v5data.keys():
            if v5data['boltwood cloud condition'] == 0: v5data['boltwood cloud color'] = 'orange'
            elif v5data['boltwood cloud condition'] == 1: v5data['boltwood cloud color'] = 'green'
            elif v5data['boltwood cloud condition'] == 2: v5data['boltwood cloud color'] = 'orange'
            elif v5data['boltwood cloud condition'] == 3: v5data['boltwood cloud color'] = 'red'
            else: v5data['boltwood cloud color'] = ''
            v5data['boltwood cloud condition string'] = cloud_condition[v5data['boltwood cloud condition']]
        else: v5data['boltwood cloud color'] = 'red'
        if 'boltwood wind condition' in v5data.keys():
            if v5data['boltwood wind condition'] == 0: v5data['boltwood wind color'] = 'orange'
            elif v5data['boltwood wind condition'] == 1: v5data['boltwood wind color'] = 'green'
            elif v5data['boltwood wind condition'] == 2: v5data['boltwood wind color'] = 'orange'
            elif v5data['boltwood wind condition'] == 3: v5data['boltwood wind color'] = 'red'
            else: v5data['boltwood wind color'] = ''
            v5data['boltwood wind condition string'] = wind_condition[v5data['boltwood wind condition']]
        else: v5data['boltwood wind color'] = 'red'
        if 'boltwood rain condition' in v5data.keys():
            if v5data['boltwood rain condition'] == 0: v5data['boltwood rain color'] = 'orange'
            elif v5data['boltwood rain condition'] == 1: v5data['boltwood rain color'] = 'green'
            elif v5data['boltwood rain condition'] == 2: v5data['boltwood rain color'] = 'red'
            elif v5data['boltwood rain condition'] == 3: v5data['boltwood rain color'] = 'red'
            else: v5data['boltwood rain color'] = ''
            v5data['boltwood rain condition string'] = rain_condition[v5data['boltwood rain condition']]
        else: v5data['boltwood rain color'] = 'red'
        if 'boltwood day condition' in v5data.keys():
            if v5data['boltwood day condition'] == 0: v5data['boltwood day color'] = 'orange'
            elif v5data['boltwood day condition'] == 1: v5data['boltwood day color'] = 'green'
            elif v5data['boltwood day condition'] == 2: v5data['boltwood day color'] = 'orange'
            elif v5data['boltwood day condition'] == 3: v5data['boltwood day color'] = 'red'
            else: v5data['boltwood day color'] = ''
            v5data['boltwood day condition string'] = day_condition[v5data['boltwood day condition']]
        else: v5data['boltwood day color'] = 'red'
        if 'boltwood roof close' in v5data.keys():
            if v5data['boltwood roof close'] == 0: v5data['boltwood roof close color'] = 'green'
            elif v5data['boltwood roof close'] == 1: v5data['boltwood roof close color'] = 'red'
            else: v5data['boltwood roof close color'] = ''
            v5data['boltwood roof close string'] = roof_close[v5data['boltwood roof close']]
        else: v5data['boltwood roof close color'] = 'red'
        tlog.app_log.info('  colors determined')

        if v20clarity_age > 10:
            if 'boltwood ambient temp' in v20data.keys():
                del v20data['boltwood ambient temp']
            if 'boltwood sky temp' in v20data.keys():
                del v20data['boltwood sky temp']
            if 'boltwood wind speed' in v20data.keys():
                del v20data['boltwood wind speed']
            if 'boltwood humidity' in v20data.keys():
                del v20data['boltwood humidity']
            if 'boltwood rain status string' in v20data.keys():
                del v20data['boltwood rain status string']
            if 'boltwood wet status string' in v20data.keys():
                del v20data['boltwood wet status string']
            if 'boltwood cloud condition string' in v20data.keys():
                del v20data['boltwood cloud condition string']
            if 'boltwood wind condition string' in v20data.keys():
                del v20data['boltwood wind condition string']
            if 'boltwood rain condition string' in v20data.keys():
                del v20data['boltwood rain condition string']
            if 'boltwood day condition string' in v20data.keys():
                del v20data['boltwood day condition string']
            if 'boltwood roof close string' in v20data.keys():
                del v20data['boltwood roof close string']
            tlog.app_log.info('  V20 Clarity data is old and has been removed.')

        if v5clarity_age > 10:
            if 'boltwood ambient temp' in v5data.keys():
                del v5data['boltwood ambient temp']
            if 'boltwood sky temp' in v5data.keys():
                del v5data['boltwood sky temp']
            if 'boltwood wind speed' in v5data.keys():
                del v5data['boltwood wind speed']
            if 'boltwood humidity' in v5data.keys():
                del v5data['boltwood humidity']
            if 'boltwood rain status string' in v5data.keys():
                del v5data['boltwood rain status string']
            if 'boltwood wet status string' in v5data.keys():
                del v5data['boltwood wet status string']
            if 'boltwood cloud condition string' in v5data.keys():
                del v5data['boltwood cloud condition string']
            if 'boltwood wind condition string' in v5data.keys():
                del v5data['boltwood wind condition string']
            if 'boltwood rain condition string' in v5data.keys():
                del v5data['boltwood rain condition string']
            if 'boltwood day condition string' in v5data.keys():
                del v5data['boltwood day condition string']
            if 'boltwood roof close string' in v5data.keys():
                del v5data['boltwood roof close string']
            tlog.app_log.info('  V5 Clarity data is old and has been removed.')


        ##---------------------------------------------------------------------
        ## Format and Color Code ACP Data
        ##---------------------------------------------------------------------
        ACP_connected = {True: 'Connected', False: 'Disconnected'}

        if 'ACP connected' in v20data.keys():
            v20data['ACP connected string'] = ACP_connected[v20data['ACP connected']]
            if (v20data['ACP connected']):
                v20data['ACP connected color'] = 'green'
                if ('ACP park status' in v20data.keys()) and\
                   ('ACP slewing status' in v20data.keys()) and\
                   ('ACP tracking status' in v20data.keys()):
                    P = v20data['ACP park status']
                    S = v20data['ACP slewing status']
                    T = v20data['ACP tracking status']
                    if P:
                        v20data['ACP status string'] = 'Parked'
                        v20data['ACP status color'] = ''
                    elif not P and not S and not T:
                        v20data['ACP status string'] = 'Stationary'
                        v20data['ACP status color'] = ''
                    elif not P and S and not T:
                        v20data['ACP status string'] = 'Slewing'
                        v20data['ACP status color'] = 'orange'
                    elif not P and S and T:
                        v20data['ACP status string'] = 'Slewing'
                        v20data['ACP status color'] = 'orange'
                    elif not P and not S and T:
                        v20data['ACP status string'] = 'Tracking'
                        v20data['ACP status color'] = 'green'
                    else:
                        v20data['ACP status string'] = '{}{}{}'.format(P,S,T)
                        v20data['ACP status color'] = 'red'
                if ('ACP target RA' in v20data.keys()) and ('ACP target Dec' in v20data.keys()):
                    v20c = SkyCoord(ra=v20data['ACP target RA']*u.degree,\
                                    dec=v20data['ACP target Dec']*u.degree,\
                                    frame='icrs')
                    v20coord = '{} {}'.format(\
                                              v20c.ra.to_string(sep=':', precision=1),\
                                              v20c.dec.to_string(sep=':', precision=1),\
                                             )
                else:
                    v20c = None
                    v20coord = ''
            else:
                v20data['ACP connected color'] = ''
                v20coord = ''
        else:
            v20data['ACP connected color'] = ''
            v20coord = ''
        tlog.app_log.info('  V20 ACP Connected Color determined')

        if 'ACP connected' in v5data.keys():
            v5data['ACP connected string'] = ACP_connected[v5data['ACP connected']]
            if (v5data['ACP connected']):
                v5data['ACP connected color'] = 'green'
                if ('ACP park status' in v5data.keys()) and\
                   ('ACP slewing status' in v5data.keys()) and\
                   ('ACP tracking status' in v5data.keys()):
                    P = v5data['ACP park status']
                    S = v5data['ACP slewing status']
                    T = v5data['ACP tracking status']
                    if P:
                        v5data['ACP status string'] = 'Parked'
                        v5data['ACP status color'] = ''
                    elif not P and not S and not T:
                        v5data['ACP status string'] = 'Stationary'
                        v5data['ACP status color'] = ''
                    elif not P and S and not T:
                        v5data['ACP status string'] = 'Slewing'
                        v5data['ACP status color'] = 'orange'
                    elif not P and S and T:
                        v5data['ACP status string'] = 'Slewing'
                        v5data['ACP status color'] = 'orange'
                    elif not P and not S and T:
                        v5data['ACP status string'] = 'Tracking'
                        v5data['ACP status color'] = 'green'
                    else:
                        v5data['ACP status string'] = '{}{}{}'.format(P,S,T)
                        v5data['ACP status color'] = 'red'
                if ('ACP target RA' in v5data.keys()) and ('ACP target Dec' in v5data.keys()):
                    v5c = SkyCoord(ra=v5data['ACP target RA']*u.degree,\
                                   dec=v5data['ACP target Dec']*u.degree,\
                                   frame='icrs')
                    v5coord = '{} {}'.format(\
                                             v5c.ra.to_string(sep=':', precision=1),\
                                             v5c.dec.to_string(sep=':', precision=1),\
                                            )
                else:
                    v5c = None
                    v5coord = ''
            else:
                v5data['ACP connected color'] = ''
                v5coord = ''
        else:
            v5data['ACP connected color'] = ''
            v5coord = ''
        tlog.app_log.info('  V5 ACP Connected Color determined')

        ##---------------------------------------------------------------------
        ## Get disk use info
        ##---------------------------------------------------------------------
        paths = {'Drobo': os.path.join('/', 'Volumes', 'Drobo'),\
                 'Data': os.path.expanduser('~'),\
                 'USB Drive B': os.path.join('/', 'Volumes', 'WD500B'),\
                 'USB Drive C': os.path.join('/', 'Volumes', 'WD500_C'),\
                 'Vega': os.path.join('/', 'Volumes', 'Data_V5'),\
                 'Black': os.path.join('/', 'Volumes', 'Data_V20'),\
                }

        disks = {}
        for disk in paths.keys():
            if os.path.exists(paths[disk]):
                size_GB, avail_GB, pcnt_used = free_space(paths[disk])
                if disk == 'Drobo':
                    size_GB -= 12226.56
                    avail_GB -= 12226.56
                    pcnt_used = float(size_GB - avail_GB)/float(size_GB) * 100
                disks[disk] = [size_GB, avail_GB, pcnt_used]

        tlog.app_log.info('  Disk use data determined')

        ##---------------------------------------------------------------------
        ## Render
        ##---------------------------------------------------------------------
        if nowut.hour >= 4:
            link_date_string = nowut.strftime('%Y%m%dUT')
        else:
            link_date_string = (nowut - tdelta(1,0)).strftime('%Y%m%dUT')

        tlog.app_log.info('  Rendering Status')
        self.render("status.html", title="VYSOS Status",\
                    now = now,\
                    nowut = nowut,\
                    link_date_string = link_date_string,\
                    v20clarity_age = v20clarity_age,\
                    v20clarity_color = v20clarity_color,\
                    v20data_time = v20data_time,\
                    v20data_age = v20data_age,\
                    v20data_color = v20data_color,\
                    v20data = v20data,\
                    v20coord = v20coord,\
                    v5clarity_age = v5clarity_age,\
                    v5clarity_color = v5clarity_color,\
                    v5data_time = v5data_time,\
                    v5data_age = v5data_age,\
                    v5data_color = v5data_color,\
                    v5data = v5data,\
                    v5coord = v5coord,\
                    moon = moon,\
                    sun = sun,\
                    disks = disks,\
                    v5_nimages = get_nimages('V5', link_date_string),\
                    v20_nimages = get_nimages('V20', link_date_string),\
                    )
        tlog.app_log.info('  Done')
