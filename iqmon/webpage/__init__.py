from pathlib import Path
import logging
import pymongo
from datetime import datetime, timedelta

from astroplan import Observer
from astropy import coordinates as c
from astropy.time import Time


##-------------------------------------------------------------------------
## Create logger object
##-------------------------------------------------------------------------
log = logging.getLogger('FlaskLogger')
log.setLevel(logging.DEBUG)
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
## Set up console output
## Set up file output
LogFileName = '/usr/local/var/log/flask.log'
LogFileHandler = logging.FileHandler(LogFileName)
LogFileHandler.setLevel(logging.DEBUG)
LogFileHandler.setFormatter(LogFormat)
log.addHandler(LogFileHandler)


##-------------------------------------------------------------------------
## Function: mongo_query
##-------------------------------------------------------------------------
def mongo_query(collection, query_dict, cfg, distinct=False, count=False,
                sort=[('date', pymongo.ASCENDING)]):
    log.debug(f'Connecting to mongo db, collection {collection}')
    mongo_host = cfg['mongo'].get('host')
    mongo_port = cfg['mongo'].getint('port')
    mongo_db = cfg['mongo'].get('db')
    mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
    mongo_iqmon = mongoclient[mongo_db][collection]
    if distinct is True:
        query_result = mongo_iqmon.distinct(query_dict)
    elif count is True:
        query_result = mongo_iqmon.find(query_dict).count()
    else:
        query_result = mongo_iqmon.find(query_dict, sort=sort)
    mongoclient.close()
    return query_result


##-------------------------------------------------------------------------
## Function: get_twilights
##-------------------------------------------------------------------------
def get_twilights(start, end, cfg, nsample=256):
    """ Determine sunrise and sunset times """
    location = c.EarthLocation(
        lat=cfg['Telescope'].getfloat('site_lat'),
        lon=cfg['Telescope'].getfloat('site_lon'),
        height=cfg['Telescope'].getfloat('site_elevation'),
    )
    obs = Observer(location=location, name=cfg['Telescope'].get('name'))
    sunset = obs.sun_set_time(Time(start), which='next').datetime
    sunrise = obs.sun_rise_time(Time(start), which='next').datetime

    # Calculate and order twilights and set plotting alpha for each
    twilights = [(start, 'start', 0.0),
                 (sunset, 'sunset', 0.0),
                 (obs.twilight_evening_civil(Time(start),
                                             which='next').datetime, 'ec', 0.1),
                 (obs.twilight_evening_nautical(Time(start),
                                                which='next').datetime, 'en', 0.2),
                 (obs.twilight_evening_astronomical(Time(start),
                                                    which='next').datetime, 'ea', 0.3),
                 (obs.twilight_morning_astronomical(Time(start),
                                                    which='next').datetime, 'ma', 0.5),
                 (obs.twilight_morning_nautical(Time(start),
                                                which='next').datetime, 'mn', 0.3),
                 (obs.twilight_morning_civil(Time(start),
                                             which='next').datetime, 'mc', 0.2),
                 (sunrise, 'sunrise', 0.1),
                 ]
    twilights.sort(key=lambda x: x[0])
    final = {'sunset': 0.1, 'ec': 0.2, 'en': 0.3, 'ea': 0.5,
             'ma': 0.3, 'mn': 0.2, 'mc': 0.1, 'sunrise': 0.0}
    twilights.append((end, 'end', final[twilights[-1][1]]))

    return twilights


##-------------------------------------------------------------------------
## Function: overplot_twilights
##-------------------------------------------------------------------------
def overplot_twilights(plot_list, plot_end, cfg, plot_ndays=1, log=None):
    for days in range(1,plot_ndays+1):
        if log is not None: log.info(f'Getting twilights for {days} days ago')
        twilights = get_twilights(plot_end-timedelta(days=days),
                                  plot_end-timedelta(days=days-1),
                                  cfg)
        for plot_info in plot_list:
            for j in range(len(twilights)-1):
                name, plot, top, bottom = plot_info
                plot.quad(top=[top], bottom=[bottom],
                           left=[twilights[j][0]], right=[twilights[j+1][0]],
                           color="blue", alpha=twilights[j+1][2])
        if log is not None: log.info(f'  Added twilights for {days} days ago')
    return

