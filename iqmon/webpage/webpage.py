from pathlib import Path
import configparser
import logging
import flask
import pymongo
from datetime import datetime, timedelta
import numpy as np

from iqmon import get_webpage_config, get_all_configs
from iqmon.webpage import mongo_query
from iqmon.webpage.weather_plot import generate_weather_plot
from iqmon.webpage.iqmon_plot import generate_iqmon_plot

app = flask.Flask(__name__)
log = logging.getLogger('FlaskLogger')


##-------------------------------------------------------------------------
## static_path: /static/plots/<string:telescope>/<string:date>/<string:filename>
##-------------------------------------------------------------------------
@app.route("/static/plots/<string:telescope>/<string:date>/<string:filename>")
def static_path(telescope, date, filename):
    log.info(f"Returning static path: /Users/vysosuser/plots/{telescope} {filename}")
    return flask.send_from_directory(f'/Users/vysosuser/plots/{telescope}/{date}', filename)


##-------------------------------------------------------------------------
## base: /
##-------------------------------------------------------------------------
@app.route("/")
def base():
    log.info(f'')
    log.info(f'Building {__name__} base')
    webcfg, cfgs = get_all_configs()
    
    if len(cfgs.keys()) == 0:
        log.info(f'Building {__name__} nightWeather')
        script, div = generate_weather_plot(webcfg, None, date=None)
        log.info(f"Rendering template")
        return flask.render_template('nightPlot.html',
                                     telescope=None,
                                     date=datetime.utcnow().strftime('%Y%m%dUT'),
                                     title='Weather',
                                     script=script,
                                     div=div,
                                     )
    else:
        if 'primary' in cfgs.keys():
            telescope = cfgs['primary']
            return status(telescope)
        else:
            return 'Hello World'


##-------------------------------------------------------------------------
## status: /<string:telescope>
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/")
def status(telescope):
    log.info(f'')
    log.info(f'Building {__name__} status')
    tick = datetime.utcnow()

    webcfg, cfgs = get_all_configs()
    if telescope not in cfgs.keys():
        return f'Could not find config for "{telescope}"'
    telcfg = cfgs[telescope]
    script, div = generate_weather_plot(webcfg, telcfg, plot_ndays=2, span_hours=12)

    ## Format currentweather
    log.info(f"Querying weather limits")
    query_result = mongo_query('weather_limits', {}, webcfg)
    weather_limits = query_result[0]

    ## Get Currentweather:
    ## temperature, clouds, wind, gust, rain, safe
    log.info(f'Get current weather')
    currentweather = {}
    devices = webcfg['Weather'].get('devices').split(',')
    query_dict = {'date': {'$gt': tick-timedelta(minutes=10), '$lt': tick}}
    for device in devices:
        deviceweather = mongo_query(device, query_dict, webcfg, last=True)
        if len(deviceweather) > 0:
            currentdeviceweather = deviceweather[-1]
        else:
            currentdeviceweather = {}
        # Look for outside temperature
        if 'outside temperature' in currentdeviceweather.keys()\
             and 'outside temperature' not in currentweather.keys():
            log.info(f"{device} outside temperature {currentdeviceweather['outside temperature']}")
            currentweather['outside temperature'] = currentdeviceweather['outside temperature']
        if 'temperature units' in currentdeviceweather.keys()\
             and 'temperature units' not in currentweather.keys():
            log.info(f"{device} temperature units {currentdeviceweather['temperature units']}")
            currentweather['temperature units'] = currentdeviceweather['temperature units']
        # Look for date
        if 'date' in currentdeviceweather.keys()\
             and 'date' not in currentweather.keys():
            log.info(f"{device} date {currentdeviceweather['date']}")
            currentweather['date'] = currentdeviceweather['date']
            currentweather['age'] = (tick - currentweather['date']).total_seconds()
        # Look for clouds
        if 'cloud value' in currentdeviceweather.keys()\
             and 'cloud value' not in currentweather.keys():
            log.info(f"{device} cloud value {currentdeviceweather['cloud value']}")
            currentweather['cloud value'] = currentdeviceweather['cloud value']
            if currentweather['cloud value'] > weather_limits['cloudy']:
                currentweather['cloud status'] = 'very cloudy'
            elif currentweather['cloud value'] > weather_limits['clear']:
                currentweather['cloud status'] = 'cloudy'
            else:
                currentweather['cloud status'] = 'clear'
        # Look for wind
        if 'wind speed' in currentdeviceweather.keys()\
             and 'wind speed' not in currentweather.keys():
            log.info(f"{device} wind speed {currentdeviceweather['wind speed']}")
            currentweather['wind speed'] = currentdeviceweather['wind speed']
            if currentweather['wind speed'] > weather_limits['windy']:
                currentweather['wind status'] = 'very windy'
            elif currentweather['wind speed'] > weather_limits['calm']:
                currentweather['wind status'] = 'windy'
            else:
                currentweather['wind status'] = 'calm'
        if 'wind speed units' in currentdeviceweather.keys()\
             and 'wind speed units' not in currentweather.keys():
            log.info(f"{device} wind speed units {currentdeviceweather['wind speed units']}")
            currentweather['wind speed units'] = currentdeviceweather['wind speed units']
        # Look for rain
        if 'rain value' in currentdeviceweather.keys()\
             and 'rain value' not in currentweather.keys():
            log.info(f"{device} rain value {currentdeviceweather['rain value']}")
            currentweather['rain value'] = currentdeviceweather['rain value']
            if currentweather['rain value'] > weather_limits['dry']:
                currentweather['rain status'] = 'dry'
            elif currentweather['rain value'] > weather_limits['wet']:
                currentweather['rain status'] = 'wet'
            else:
                currentweather['rain status'] = 'rain'

    ## Dome Status
    log.info(f"Querying dome status database")
    query_dict = {'date': {'$gt': tick-timedelta(minutes=5), '$lt': tick}}
    query_result = mongo_query(f'{telescope}_dome', query_dict, telcfg)
    domestatus = [d for d in query_result]
    log.info(f"  Got {len(domestatus)} data points")
    shutter_values = {0: 0, 1: 1, 2: 0, 3: 1, 4: 4}
    for i,d in enumerate(domestatus):
        if d['shutterstatus'] == 4 and i > 0:
            domestatus[i]['open_closed'] = domestatus[i-1]['open_closed']
        else:
            domestatus[i]['open_closed'] = shutter_values[d['shutterstatus']]
    ## Format currentstatus
    dome_string = {0: 'Open', 1: 'Closed', 2: 'Opening', 3: 'Closing', 4: 'Unknown'}
    dome_color = {0: 'green', 1: 'red', 2: 'orange', 3: 'orange', 4: 'black'}
    try:
        currentstatus = domestatus[-1]
    except:
        currentstatus = {}
    else:
        currentstatus['age'] = (datetime.utcnow() - currentstatus['date']).total_seconds()
        if currentstatus['shutterstatus'] == 4:
            query_dict = {'shutterstatus': {'$ne': 4}}
            query_result = mongo_query(f'{telescope}_dome', query_dict, telcfg,
                                       sort=[('date', pymongo.DESCENDING)])
            last_shutter = query_result.next()
            currentstatus['dome_string'] = dome_string[last_shutter['shutterstatus']]
            currentstatus['dome_color'] = dome_color[last_shutter['shutterstatus']]
        else:
            currentstatus['dome_string'] = dome_string[currentstatus['shutterstatus']]
            currentstatus['dome_color'] = dome_color[currentstatus['shutterstatus']]

    ## Telescope Status
    log.info(f"Querying telescope status database")
    query_dict = {'date': {'$gt': tick-timedelta(minutes=5), '$lt': tick}}
    query_result = mongo_query(f'{telescope}_telescope', query_dict, telcfg)
    telescopestatus = [d for d in query_result]
    log.info(f"  Got {len(telescopestatus)} data points")
    ## Format currentstatus
    if len(telescopestatus) > 0:
        log.info(currentstatus)
        currentstatus.update(telescopestatus[-1])
        log.info(currentstatus)
        if currentstatus['connected'] is True:
            currentstatus['alt'] = f"{currentstatus['alt']:.1f}"
            currentstatus['az'] = f"{currentstatus['az']:.1f}"
            if currentstatus['slewing'] is True:
                currentstatus['slew status'] = 'slewing'
            elif currentstatus['tracking'] is True:
                currentstatus['slew status'] = 'tracking'
            elif currentstatus['park'] is True:
                currentstatus['slew status'] = 'parked'
            else:
                currentstatus['slew status'] = 'stationary'
        else:
            currentstatus['slew status'] = ''
            currentstatus['alt'] = ''
            currentstatus['az'] = ''

    log.info(f"Rendering flask template")
    result = flask.render_template('status.html',
                                   date_string=datetime.utcnow().strftime('%Y%m%dUT'),
                                   telescope=telescope,
                                   currentweather=currentweather,
                                   currentstatus=currentstatus,
                                   now=datetime.now(),
                                   utcnow=datetime.utcnow(),
                                   script=script,
                                   div=div,
                                   image=webcfg['WebPage'].get('image', ''),
                                   image_link=webcfg['WebPage'].get('image_link', ''),
                                   image_title=webcfg['WebPage'].get('image_title', ''),
                                   )
    tock = datetime.utcnow()
    duration = (tock-tick).total_seconds()
    log.info(f'Page built in {duration:.2f} s')
    return result


##-------------------------------------------------------------------------
## nightWeather: /<string:telescope>/weather/<string:date>
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/weather/<string:date>")
def nightWeather(telescope, date):
    log.info(f'Building {__name__} nightWeather {telescope}')
    webcfg, cfgs = get_all_configs()
    if telescope not in cfgs.keys():
        return f'Could not find config for "{telescope}"'
    telcfg = cfgs[telescope]
    script, div = generate_weather_plot(webcfg, telcfg, date=date)

    log.info(f"Rendering template")
    return flask.render_template('nightPlot.html',
                                 telescope=telescope,
                                 date=date,
                                 title='Weather',
                                 script=script,
                                 div=div,
                                 )


##-------------------------------------------------------------------------
## nightReport: /<string:telescope>/report/<string:date>
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/report/<string:date>")
def nightReport(telescope, date):
    log.info(f'Building {__name__} iqmonReport {telescope}')
    webcfg, cfgs = get_all_configs()
    if telescope not in cfgs.keys():
        return f'Could not find config for "{telescope}"'
    script, div = generate_iqmon_plot(cfgs[telescope], date=date)

    log.info(f"Rendering template")
    return flask.render_template('nightPlot.html',
                                 telescope=telescope,
                                 date=date,
                                 title='IQMon Report',
                                 script=script,
                                 div=div,
                                 )


##-------------------------------------------------------------------------
## imageList: /<string:telescope>/images/<string:date>
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/images/<string:date>")
def imageList(telescope, date):
    log.info(f'Building {__name__} imageList {telescope} {date}')
    subject = date
    webcfg, cfgs = get_all_configs()
    if telescope not in cfgs.keys():
        return f'Could not find config for "{telescope}"'
    telcfg = cfgs[telescope]

    log.info(f"Querying image database")
    start = datetime.strptime(date, '%Y%m%dUT')
    end = start + timedelta(days=1)
    query_dict = {'telescope': telescope,
                  'date': {'$gt': start, '$lt': end}}
    image_list = mongo_query('iqmon', query_dict, telcfg)
    flat_count = len([d for d in image_list if d['imtype'] in ['FLAT', 'TWIFLAT']])
    cal_count = len([d for d in image_list if d['imtype'] in ['BIAS', 'DARK']])
    object_count = len([d for d in image_list if d['imtype'] in ['OBJECT']])
    total_processing_time = np.sum([d['processing_time'] for d in image_list])
    log.info(f"Got {len(image_list)} images")

    log.info(f"Getting IQMon limits")
    iqmon_limits = mongo_query('V5limits', {}, webcfg, sort=[])[0]
    for key in iqmon_limits:
        log.info(f"  {key} : {iqmon_limits[key]}")

    log.info(f"Assigning image flags")
    for i,image in enumerate(image_list):
        for key in iqmon_limits:
            if key in image.keys():
                if image[key] > iqmon_limits[key]:
                    image_list[i][f"{key} flag"] = True

    log.info(f"Rendering template")
    return flask.render_template('imageList.html',
                                 telescope=telescope,
                                 subject=subject,
                                 image_list=image_list,
                                 flat_count=flat_count,
                                 cal_count=cal_count,
                                 object_count=object_count,
                                 total_processing_time=total_processing_time,
                                 )


##-------------------------------------------------------------------------
## nightList: /<string:telescope>/nights/
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/nights/")
def nightList(telescope):
    log.info(f'Building {__name__} nightList {telescope}')
    webcfg, cfgs = get_all_configs()
    if telescope not in cfgs.keys():
        return f'Could not find config for "{telescope}"'
    telcfg = cfgs[telescope]

    log.info(f"Querying image database")
    query_dict = 'UT date string'
    query_result = mongo_query('iqmon', query_dict, telcfg, distinct=True)
    night_list = sorted(query_result, reverse=True)
    log.info(f"  Found {len(night_list)} nights")

    nights = []
    for datestr in night_list:
        log.info(f"Querying image database for {datestr}")
        night_info = {'date': datestr}
        start = datetime.strptime(datestr, '%Y%m%dUT')
        end = start + timedelta(days=1)
        query_dict = {'telescope': telescope,
                      'date': {'$gt': start, '$lt': end}}
        query_result = mongo_query('iqmon', query_dict, telcfg)
        image_list = [d for d in query_result]
        night_info['n images'] = len(image_list)
        nights.append(night_info)

    log.info(f"Rendering template")
    return flask.render_template('nightList.html',
                                 telescope=telescope,
                                 nights=nights,
                                 )


if __name__ == "__main__":
    app.run()
