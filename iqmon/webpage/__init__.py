from pathlib import Path
import logging
import pymongo
from datetime import datetime, timedelta
from time import sleep

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
def mongo_query(collection, query_dict, cfg,
                distinct=False, count=False, last=False,
                sort=[('date', pymongo.ASCENDING)]):
    log.debug(f'Connecting to mongo db, collection {collection}')
    mongo_host = cfg['mongo'].get('host')
    mongo_port = cfg['mongo'].getint('port')
    mongo_db = cfg['mongo'].get('db')
    mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
    mongo_iqmon = mongoclient[mongo_db][collection]
    if distinct is True:
        query_result = list(mongo_iqmon.distinct(query_dict))
    elif count is True:
        query_result = mongo_iqmon.find(query_dict).count()
    elif last is True:
        query_result = list(mongo_iqmon.find(query_dict,
                            sort=[('date', pymongo.DESCENDING)]).limit(1))
    else:
        query_result = list(mongo_iqmon.find(query_dict, sort=sort))
    mongoclient.close()
    return query_result


##-------------------------------------------------------------------------
## Function: get_twilights
##-------------------------------------------------------------------------
def get_sun_alt(now, loc):
    obj = c.get_sun(now)
    altaz = c.AltAz(location=loc, obstime=now,
                    pressure=0*u.mbar,
                    temperature=5*u.Celsius,
                    obswl=0.5*u.micron,
                    relative_humidity=10,
                   )
    alt = obj.transform_to(altaz).alt
    return alt


def get_twilights(start, end, webcfg, nsample=256):
    """ Determine sunrise and sunset times """
    location = c.EarthLocation(
        lat=webcfg['site'].getfloat('site_lat'),
        lon=webcfg['site'].getfloat('site_lon'),
        height=webcfg['site'].getfloat('site_elevation'),
    )
    UTdate_string = start.strftime('%Y%m%dUT')
    now = Time(datetime.strptime(UTdate_string, '%Y%m%dUT'))
    alt = get_sun_alt(now, loc)
    now += TimeDelta(np.floor(alt.value/15)*u.hour)
    alt = get_sun_alt(now, loc=loc)
    assert alt.value > 0
    # Sunset
    while alt.value > 0:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    sunset = now + TimeDelta(10*u.minute) # Add 10 minute fudge factor for refraction
    # Civil Dusk
    while alt.value > -6:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    civil_dusk = now
    # Nautical Dusk
    while alt.value > -12:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    nautical_dusk = now
    # Astronomincal Dusk
    while alt.value > -18:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    astronomical_dusk = now

    # Larger steps through night
    while alt.value < -18:
        now += TimeDelta(60*u.minute)
        alt = get_sun_alt(now, loc)
    now -= TimeDelta(60*u.minute)
    alt = get_sun_alt(now, loc=loc)
    
    # Astronomical Dawn
    while alt.value < -18:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    astronomical_dawn = now
    # Nautical Dawn
    while alt.value < -12:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    nautical_dawn = now
    # Civil Dawn
    while alt.value < -6:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    civil_dawn = now
    # Sunrise
    while alt.value < 0:
        now += TimeDelta(precision)
        alt = get_sun_alt(now, loc)
    sunrise = now - TimeDelta(10*u.minute) # Subtract 10 minute fudge factor for refraction

    return sunset, civil_dusk, nautical_dusk, astronomical_dusk, astronomical_dawn, nautical_dawn, civil_dawn, sunrise


def old_get_twilights(start, end, webcfg, nsample=256):
    """ Determine sunrise and sunset times """

    from astroplan import Observer

    location = c.EarthLocation(
        lat=webcfg['site'].getfloat('site_lat'),
        lon=webcfg['site'].getfloat('site_lon'),
        height=webcfg['site'].getfloat('site_elevation'),
    )
    obs = Observer(location=location, name=webcfg['site'].get('name', ''))
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
def overplot_twilights(plot_list, plot_end, webcfg, ndays=1, log=None):
    for days in range(1,ndays+1):
        if log is not None: log.info(f'Getting twilights for {days} days ago')
        try:
            twilights = get_twilights(plot_end-timedelta(days=days),
                                      plot_end-timedelta(days=days-1),
                                      webcfg)
            for plot_info in plot_list:
                for j in range(len(twilights)-1):
                    name, plot, top, bottom = plot_info
                    plot.quad(top=[top], bottom=[bottom],
                               left=[twilights[j][0]], right=[twilights[j+1][0]],
                               color="blue", alpha=twilights[j+1][2])
            if log is not None: log.info(f'  Added twilights for {days} days ago')
        except Exception as err:
            if log is not None: log.error(err)
            pass
    return

