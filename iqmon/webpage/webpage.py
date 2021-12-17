from pathlib import Path
import configparser
import logging
import flask
import pymongo
from datetime import datetime, timedelta

import numpy as np
from astroplan import Observer
from astropy import coordinates as c
from astropy.time import Time

from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.layouts import column
from bokeh.models import DatetimeTickFormatter, NumeralTickFormatter, DatetimeTicker

app = flask.Flask(__name__)


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
def mongo_query(collection, query_dict, distinct=False, count=False,
                sort=[('date', pymongo.ASCENDING)]):
    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)
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
def get_twilights(start, end, nsample=256):
    """ Determine sunrise and sunset times """
    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)

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
## Function: generate_weather_plot
##-------------------------------------------------------------------------
def generate_weather_plot(telescope, date=None, plot_ndays=1, span_hours=24):
    log.info(f"Querying weather limits")
    query_result = mongo_query('weather_limits', {})
    weather_limits = query_result[0]

    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)

    if date is None:
        end = datetime.utcnow()
        start = end - timedelta(days=plot_ndays)
    else:
        start = datetime.strptime(date, '%Y%m%dUT')
        end = start + timedelta(days=plot_ndays)

    ##-------------------------------------------------------------------------
    ## Weather Query
    log.info(f"Querying weather database")
    query_dict = {'date': {'$gt': start, '$lt': end}}
    query_result = mongo_query('weather', query_dict)
    weather = [d for d in query_result]
    log.info(f"  Got {len(weather)} data points")
    date = [w['date'] for w in weather]

    ##-------------------------------------------------------------------------
    ## Telescope Status Query
    log.info(f"Querying telescope status database")
    query_dict = {'date': {'$gt': start, '$lt': end}}
    query_result = mongo_query(f'{telescope}status', query_dict)
    telstatus = [d for d in query_result]
    log.info(f"  Got {len(telstatus)} data points")
    shutter_values = {0: 0, 1: 1, 2: 0, 3: 1, 4: 4}
    for i,d in enumerate(telstatus):
        if d['dome_shutterstatus'] == 4 and i > 0:
            telstatus[i]['dome_numerical_status'] = telstatus[i-1]['dome_numerical_status']
        else:
            telstatus[i]['dome_numerical_status'] = shutter_values[d['dome_shutterstatus']]

    ##-------------------------------------------------------------------------
    ## IQMon Query
    log.info(f"Querying IQMon results database")
    query_dict = {'telescope': telescope,
                  'date': {'$gt': start, '$lt': end}}
    query_result = mongo_query('iqmon', query_dict)
    iqmon = [d for d in query_result]
    log.info(f"  Got {len(iqmon)} data points")
    iqmon_obj_dates = [d['date'] for d in iqmon if d['imtype'] == 'OBJECT']
    iqmon_obj_alt = [d['alt']/90 for d in iqmon if d['imtype'] == 'OBJECT']
    iqmon_cal_dates = [d['date'] for d in iqmon if d['imtype'] in ['BIAS', 'DARK']]
    iqmon_cal_alt = [0.5 for d in iqmon if d['imtype'] in ['BIAS', 'DARK']]
    iqmon_flat_dates = [d['date'] for d in iqmon if d['imtype'] in ['FLAT', 'TWIFLAT', 'DOMEFLAT']]
    iqmon_flat_alt = [0.5 for d in iqmon if d['imtype'] in ['FLAT', 'TWIFLAT', 'DOMEFLAT']]

    markersize = 2


    ##-------------------------------------------------------------------------
    ## Temperature Plot
#     log.info('Build temperature plot')
#     temp = [w['temp']*1.8+32 for w in weather]
#     temperature_plot = figure(width=900, height=100, x_axis_type="datetime",
#                               y_range=(25,95),
#                               x_range=(end - timedelta(hours=12), end),
#                               )
#     temperature_plot.circle(date, temp,
#                             size=markersize, color="blue", alpha=0.8)
#     temperature_plot.yaxis.axis_label = 'Temp (F)'
#     temperature_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
#     temperature_plot.yaxis.ticker = [30, 50, 70, 90]
#     temperature_plot.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Cloudiness Plot
    log.info('Build cloudiness plot')
    clouds = [w['clouds'] for w in weather]
    cloudiness_plot = figure(width=900, height=100, x_axis_type="datetime",
                             y_range=(-45,5),
                             x_range=(end - timedelta(hours=span_hours), end),
                             )
    very_cloudy_x = [date[i] for i,val in enumerate(clouds)\
                     if val >= weather_limits['cloudy']]
    very_cloudy_y = [val for i,val in enumerate(clouds)\
                     if val >= weather_limits['cloudy']]
    cloudiness_plot.circle(very_cloudy_x, very_cloudy_y,
                           size=markersize, color="red", alpha=0.8)
    cloudy_x = [date[i] for i,val in enumerate(clouds)\
                if val < weather_limits['cloudy']\
                and val >= weather_limits['clear']]
    cloudy_y = [val for i,val in enumerate(clouds)\
                if val < weather_limits['cloudy']\
                and val >= weather_limits['clear']]
    cloudiness_plot.circle(cloudy_x, cloudy_y,
                           size=markersize, color="orange", alpha=0.8)
    clear_x = [date[i] for i,val in enumerate(clouds)\
               if val < weather_limits['clear']]
    clear_y = [val for i,val in enumerate(clouds)\
               if val < weather_limits['clear']]
    cloudiness_plot.circle(clear_x, clear_y,
                           size=markersize, color="green", alpha=0.8)
    cloudiness_plot.yaxis.axis_label = 'Cloudiness (C)'
    cloudiness_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    cloudiness_plot.xaxis.visible = False


    ##-------------------------------------------------------------------------
    ## Wind Plot
    log.info('Build wind plot')
    wind = [w['wind'] for w in weather]
    wind_plot = figure(width=900, height=100, x_axis_type="datetime",
                       y_range=(0,85), x_range=cloudiness_plot.x_range,
                       )
    very_windy_x = [date[i] for i,val in enumerate(wind)\
                     if val >= weather_limits['windy']]
    very_windy_y = [val for i,val in enumerate(wind)\
                     if val >= weather_limits['windy']]
    wind_plot.circle(very_windy_x, very_windy_y,
                     size=markersize, color="red", alpha=0.8)
    windy_x = [date[i] for i,val in enumerate(wind)\
                if val < weather_limits['windy']\
                and val >= weather_limits['calm']]
    windy_y = [val for i,val in enumerate(wind)\
                if val < weather_limits['windy']\
                and val >= weather_limits['calm']]
    wind_plot.circle(windy_x, windy_y,
                     size=markersize, color="orange", alpha=0.8)
    calm_x = [date[i] for i,val in enumerate(wind)\
               if val < weather_limits['calm']]
    calm_y = [val for i,val in enumerate(wind)\
               if val < weather_limits['calm']]
    wind_plot.circle(calm_x, calm_y,
                     size=markersize, color="green", alpha=0.8)
    wind_plot.yaxis.axis_label = 'Wind (kph)'
    wind_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    wind_plot.yaxis.ticker = [0, 20, 40, 60, 80]
    wind_plot.xaxis.visible = False


    ##-------------------------------------------------------------------------
    ## Rain Plot
    log.info('Build rain plot')
    rain = [w['rain'] for w in weather]
    rain_plot = figure(width=900, height=60, x_axis_type="datetime",
                       y_range=(1000,2800), x_range=cloudiness_plot.x_range,
                       )
    dry_x = [date[i] for i,val in enumerate(rain)\
             if val >= weather_limits['dry']]
    dry_y = [val for i,val in enumerate(rain)\
             if val >= weather_limits['dry']]
    rain_plot.circle(dry_x, dry_y,
                     size=markersize, color="green", alpha=0.8)
    wet_x = [date[i] for i,val in enumerate(rain)\
                if val < weather_limits['dry']\
                and val >= weather_limits['wet']]
    wet_y = [val for i,val in enumerate(rain)\
                if val < weather_limits['dry']\
                and val >= weather_limits['wet']]
    rain_plot.circle(wet_x, wet_y,
                     size=markersize, color="orange", alpha=0.8)
    rain_x = [date[i] for i,val in enumerate(rain)\
               if val < weather_limits['wet']]
    rain_y = [val for i,val in enumerate(rain)\
               if val < weather_limits['wet']]
    rain_plot.circle(rain_x, rain_y,
                     size=markersize, color="red", alpha=0.8)
    rain_plot.yaxis.axis_label = 'Rain'
    rain_plot.yaxis.formatter = NumeralTickFormatter(format="0.0a")
    rain_plot.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Safe Plot
    log.info('Build safe plot')
    safe = [w['safe'] for w in weather]
    safe_date = [w['date'] for w in weather]
    safe_plot = figure(width=900, height=50, x_axis_type="datetime",
                       y_range=(-0.2,1.2), x_range=cloudiness_plot.x_range,
                       )
    width = (max(date)-min(date))/len(date)/2
    safe_dates = [date for i,date in enumerate(date) if safe[i] == True]
    unsafe_dates = [date for i,date in enumerate(date) if safe[i] != True]
    safe_plot.circle(safe_dates, [1]*len(safe_dates),
                     size=markersize, color="green", alpha=0.8)
    safe_plot.circle(unsafe_dates, [0]*len(unsafe_dates),
                     size=markersize, color="red", alpha=0.8)
    safe_plot.yaxis.axis_label = 'Safe'
    safe_plot.xaxis.axis_label = 'Time (UT)'
    safe_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    safe_plot.yaxis.ticker = [0,1]
    safe_plot.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Telescope Status Plot
    log.info('Build Telescope Status plot')
    dome = [s['dome_numerical_status'] for s in telstatus]
    dome_date = [s['date'] for s in telstatus]
    dome_plot = figure(width=900, height=100, x_axis_type="datetime",
                       y_range=(-0.2,1.2), x_range=cloudiness_plot.x_range,
                       )
    dome_plot.line(dome_date, dome, line_width=2)
    dome_plot.circle(iqmon_obj_dates, iqmon_obj_alt,
                     size=markersize, color="blue", alpha=0.8)
    dome_plot.circle(iqmon_cal_dates, iqmon_cal_alt,
                     size=markersize, color="black", alpha=0.8)
    dome_plot.circle(iqmon_flat_dates, iqmon_flat_alt,
                     size=markersize, color="yellow", alpha=0.8)
    dome_plot.yaxis.axis_label = f'{cfg["Telescope"].get("name")}'
    dome_plot.xaxis.axis_label = 'Time (UT)'
    dome_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    dome_plot.yaxis.ticker = [0,1]
    dome_plot.xaxis.visible = True
    dome_plot.xaxis.formatter = DatetimeTickFormatter(hourmin=['%H:%M'])
    dome_plot.xaxis.ticker = DatetimeTicker(desired_num_ticks=24)
    dome_plot.xaxis.axis_label = 'UT Time'


    ##-------------------------------------------------------------------------
    ## Overplot Twilights
    for days in range(1,plot_ndays+1):
        twilights = get_twilights(end-timedelta(days=days), end-timedelta(days=days-1))
        for j in range(len(twilights)-1):
#            temperature_plot.quad(top=[95], bottom=[25],
#                                   left=[twilights[j][0]], right=[twilights[j+1][0]],
#                                   color="blue", alpha=twilights[j+1][2])
           cloudiness_plot.quad(top=[5], bottom=[-45],
                                left=[twilights[j][0]], right=[twilights[j+1][0]],
                                color="blue", alpha=twilights[j+1][2])
           wind_plot.quad(top=[100], bottom=[0],
                          left=[twilights[j][0]], right=[twilights[j+1][0]],
                          color="blue", alpha=twilights[j+1][2])
           rain_plot.quad(top=[2800], bottom=[1000],
                          left=[twilights[j][0]], right=[twilights[j+1][0]],
                          color="blue", alpha=twilights[j+1][2])
           safe_plot.quad(top=[1.2], bottom=[-0.2],
                          left=[twilights[j][0]], right=[twilights[j+1][0]],
                          color="blue", alpha=twilights[j+1][2])

    ##-------------------------------------------------------------------------
    ## Render
    log.info(f"Rendering bokeh plot")
    script, div = components(column(#temperature_plot,
                                    cloudiness_plot,
                                    rain_plot,
                                    wind_plot,
                                    safe_plot,
                                    dome_plot,
                                    ))

    return script, div, weather, telstatus


##-------------------------------------------------------------------------
## static_path: 
##-------------------------------------------------------------------------
@app.route("/static/plots/<string:telescope>/<string:date>/<string:filename>")
def static_path(telescope, date, filename):
    log.info(f"Returning static path: /Users/vysosuser/plots/{telescope} {filename}")
    return flask.send_from_directory(f'/Users/vysosuser/plots/{telescope}/{date}', filename)


##-------------------------------------------------------------------------
## status: /
##-------------------------------------------------------------------------
@app.route("/")
def hello():
    log.info(f'Building {__name__} hello')
    return status('V5')


##-------------------------------------------------------------------------
## status: /status
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/")
def status(telescope):
    log.info(f'')
    log.info(f'Building {__name__} status')
    tick = datetime.utcnow()

    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)

    script, div, weather, telstatus = generate_weather_plot(telescope, plot_ndays=2, span_hours=12)

    ## Format currentweather
    log.info(f"Querying weather limits")
    query_result = mongo_query('weather_limits', {})
    weather_limits = query_result[0]
    currentweather = weather[-1]
    currentweather['age'] = (datetime.utcnow() - currentweather['date']).total_seconds()
    currentweather['temp F'] = currentweather['temp']*1.8 + 32
    if currentweather['clouds'] > weather_limits['cloudy']:
        currentweather['cloud status'] = 'very cloudy'
    elif currentweather['clouds'] > weather_limits['clear']:
        currentweather['cloud status'] = 'cloudy'
    else:
        currentweather['cloud status'] = 'clear'
    if currentweather['wind'] > weather_limits['windy']:
        currentweather['wind status'] = 'very windy'
    elif currentweather['wind'] > weather_limits['calm']:
        currentweather['wind status'] = 'windy'
    else:
        currentweather['wind status'] = 'calm'
    if currentweather['rain'] > weather_limits['dry']:
        currentweather['rain status'] = 'dry'
    elif currentweather['rain'] > weather_limits['wet']:
        currentweather['rain status'] = 'wet'
    else:
        currentweather['rain status'] = 'rain'

    ## Format currentstatus
    dome_string = {0: 'Open', 1: 'Closed', 2: 'Opening', 3: 'Closing', 4: 'Unknown'}
    dome_color = {0: 'green', 1: 'red', 2: 'orange', 3: 'orange', 4: 'black'}
    currentstatus = telstatus[-1]
    currentstatus['age'] = (datetime.utcnow() - currentstatus['date']).total_seconds()
    if currentstatus['dome_shutterstatus'] == 4:
        query_dict = {'dome_shutterstatus': {'$ne': 4}}
        query_result = mongo_query(f'{telescope}status', query_dict,
                                   sort=[('date', pymongo.DESCENDING)])
        last_shutter = query_result.next()
        currentstatus['dome_string'] = dome_string[last_shutter['dome_shutterstatus']]
        currentstatus['dome_color'] = dome_color[last_shutter['dome_shutterstatus']]
    else:
        currentstatus['dome_string'] = dome_string[currentstatus['dome_shutterstatus']]
        currentstatus['dome_color'] = dome_color[currentstatus['dome_shutterstatus']]
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
        currentstatus['RA'] = ''
        currentstatus['DEC'] = ''

    log.info(f"Rendering flask template")
    result = flask.render_template('status.html',
                                   telescope=telescope,
                                   weather=weather,
                                   currentweather=currentweather,
                                   now=datetime.now(),
                                   utcnow=datetime.utcnow(),
                                   script=script,
                                   div=div,
                                   status=status,
                                   currentstatus=currentstatus,
                                   date_string=datetime.utcnow().strftime('%Y%m%dUT'),
                                   image=cfg['WebPage'].get('image', ''),
                                   image_link=cfg['WebPage'].get('image_link', ''),
                                   image_title=cfg['WebPage'].get('image_title', ''),
                                   )
    tock = datetime.utcnow()
    duration = (tock-tick).total_seconds()
    log.info(f'Page built in {duration:.2f} s')
    return result

##-------------------------------------------------------------------------
## imageList: /<string:telescope>/images/<string:date>
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/images/<string:date>")
def imageList(telescope, date):
    log.info(f'Building {__name__} imageList {telescope} {date}')
    subject = date

    log.info(f"Querying image database")
    start = datetime.strptime(date, '%Y%m%dUT')
    end = start + timedelta(days=1)
    query_dict = {'telescope': telescope,
                  'date': {'$gt': start, '$lt': end}}
    query_result = mongo_query('iqmon', query_dict)
    image_list = [d for d in query_result]
    flat_count = len([d for d in image_list if d['imtype'] in ['FLAT', 'TWIFLAT']])
    cal_count = len([d for d in image_list if d['imtype'] in ['BIAS', 'DARK']])
    object_count = len([d for d in image_list if d['imtype'] in ['OBJECT']])
    total_processing_time = np.sum([d['processing_time'] for d in image_list])
    log.info(f"Got {len(image_list)} images")

    log.info(f"Getting IQMon limits")
    query_result = mongo_query('V5limits', {})
    iqmon_limits = query_result[0]
    for key in iqmon_limits:
        log.info(f"  {key} : {iqmon_limits[key]}")

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

    log.info(f"Querying image database")
    query_dict = 'UT date string'
    query_result = mongo_query('iqmon', query_dict, distinct=True)
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
        query_result = mongo_query('iqmon', query_dict)
        image_list = [d for d in query_result]
        night_info['n images'] = len(image_list)
        nights.append(night_info)

    log.info(f"Rendering template")
    return flask.render_template('nightList.html',
                                 telescope=telescope,
                                 nights=nights,
                                 )


##-------------------------------------------------------------------------
## nightWeather: /<string:telescope>/nights/<string:date>
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/nights/<string:date>")
def nightWeather(telescope, date):
    log.info(f'Building {__name__} nightWeather {telescope}')

    script, div, weather, telstatus = generate_weather_plot(telescope, date=date)

    log.info(f"Rendering template")
    return flask.render_template('nightWeather.html',
                                 telescope=telescope,
                                 date=date,
                                 script=script,
                                 div=div,
                                 )


if __name__ == "__main__":
    app.run()
