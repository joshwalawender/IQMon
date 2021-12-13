from pathlib import Path
import configparser
import logging
import flask
import pymongo
from datetime import datetime, timedelta

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
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
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
def mongo_query(collection, query_dict,
                sort=[('date', pymongo.ASCENDING)]):
    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    log.debug(f'Reading config file {cfg_path}')
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)
    log.debug('Connecting to mongo db')
    mongo_host = cfg['mongo'].get('host')
    log.debug(f"  mongo_host : {mongo_host}")
    mongo_port = cfg['mongo'].getint('port')
    log.debug(f"  mongo_port : {mongo_port}")
    mongo_db = cfg['mongo'].get('db')
    log.debug(f"  mongo_db : {mongo_db}")
    mongoclient = pymongo.MongoClient(mongo_host, mongo_port)
    mongo_iqmon = mongoclient[mongo_db][collection]
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
## status: /
##-------------------------------------------------------------------------
@app.route("/")
def status():
    log.info(f'Building {__name__} status')

    log.info(f"Querying weather limits")
    query_result = mongo_query('weather_limits', {})
    weather_limits = query_result[0]
    for key in weather_limits:
        log.info(f"  {key} : {weather_limits[key]}")

    log.info(f"Querying weather database")
    end = datetime.now() 
    start = end - timedelta(days=0, minutes=10)
    query_dict = {'date': {'$gt': start, '$lt': end}}
    query_result = mongo_query('weather', query_dict)
    weather = [d for d in query_result]
    log.info(f"Got {len(weather)} data points")
    currentweather = weather[-1]
    currentweather['age'] = (datetime.now() - currentweather['date']).total_seconds()
    currentweather['temp F'] = currentweather['temp']*1.8 + 32

    log.info(f"Rendering template")
    return flask.render_template('status.html',
                                 weather=weather,
                                 currentweather=currentweather,
                                 now=datetime.now(),
                                 utcnow=datetime.utcnow(),
                                 )


##-------------------------------------------------------------------------
## weather_plot: /weather
##-------------------------------------------------------------------------
@app.route("/weather")
def weather():
    log.info(f'Building {__name__} weather')

    log.info(f"Querying weather limits")
    query_result = mongo_query('weather_limits', {})
    weather_limits = query_result[0]
    for key in weather_limits:
        log.info(f"  {key} : {weather_limits[key]}")

    log.info(f"Querying weather database")
    end = datetime.utcnow()
    oneday = timedelta(days=1)
    plot_ndays = 2
    start = end - timedelta(days=plot_ndays)
    query_dict = {'date': {'$gt': start, '$lt': end}}
    query_result = mongo_query('weather', query_dict)
    weather = [d for d in query_result]
    log.info(f"Got {len(weather)} data points")

    log.info(f'Examining current conditions / status')
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

    markersize = 2

    ##-------------------------------------------------------------------------
    ## Temperature Plot
    log.info('Build temperature plot')
    date = [w['date'] for w in weather]
    temp = [w['temp']*1.8+32 for w in weather]
    temperature_plot = figure(width=800, height=100, x_axis_type="datetime",
                              y_range=(25,95), x_range=(end - oneday, end),
                              )
    temperature_plot.circle(date, temp,
                            size=markersize, color="blue", alpha=0.8)
    temperature_plot.yaxis.axis_label = 'Temp (F)'
    temperature_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    temperature_plot.xaxis.formatter = DatetimeTickFormatter(hourmin=['%H:%M'])
    temperature_plot.xaxis.ticker = DatetimeTicker(desired_num_ticks=24)

    ##-------------------------------------------------------------------------
    ## Safe Plot
    log.info('Build safe plot')
    safe = [w['safe'] for w in weather]
    safe_plot = figure(width=800, height=70, x_axis_type="datetime",
                       y_range=(-0.2,1.2), x_range=temperature_plot.x_range,
                       )
    width = (max(date)-min(date))/len(date)/2
    safe_dates = [date for i,date in enumerate(date) if safe[i] == True]
    unsafe_dates = [date for i,date in enumerate(date) if safe[i] != True]
    safe_plot.circle(safe_dates, [1]*len(safe_dates),
                     size=markersize, color="green", alpha=0.8)
    safe_plot.circle(unsafe_dates, [0]*len(unsafe_dates),
                     size=markersize, color="red", alpha=0.8)
#     safe_plot.vbar(x=safe_dates, width=width, bottom=0, top=1,
#                    color="green", alpha=0.2)
#     safe_plot.vbar(x=unsafe_dates, width=width, bottom=0, top=1,
#                    color="red", alpha=0.2)
    safe_plot.yaxis.axis_label = 'Safe'
    safe_plot.xaxis.axis_label = 'Time (UT)'
    safe_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    safe_plot.yaxis.ticker = [0,1]
    safe_plot.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Cloudiness Plot
    log.info('Build cloudiness plot')
    clouds = [w['clouds'] for w in weather]
    cloudiness_plot = figure(width=800, height=120, x_axis_type="datetime",
                             y_range=(-45,5), x_range=temperature_plot.x_range,
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
#     cloudiness_plot.quad(top=[5],
#                          bottom=[weather_limits['cloudy']],
#                          left=[min(date)], right=[max(date)],
#                          color="red", alpha=0.15)
#     cloudiness_plot.quad(top=[weather_limits['cloudy']],
#                          bottom=[weather_limits['clear']],
#                          left=[min(date)], right=[max(date)],
#                          color="yellow", alpha=0.15)
#     cloudiness_plot.quad(top=[weather_limits['clear']],
#                          bottom=[-45],
#                          left=[min(date)], right=[max(date)],
#                          color="green", alpha=0.15)
#     cloudiness_plot.circle(date, clouds,
#                            size=markersize, color="#53777a", alpha=0.8)
    cloudiness_plot.yaxis.axis_label = 'Cloudiness (C)'
    cloudiness_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    cloudiness_plot.xaxis.visible = False


    ##-------------------------------------------------------------------------
    ## Wind Plot
    log.info('Build wind plot')
    wind = [w['wind'] for w in weather]
    wind_plot = figure(width=800, height=120, x_axis_type="datetime",
                       y_range=(0,80), x_range=temperature_plot.x_range,
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
#     wind_plot.quad(top=[80],
#                    bottom=[weather_limits['windy']],
#                    left=[min(date)], right=[max(date)],
#                    color="red", alpha=0.15)
#     wind_plot.quad(top=[weather_limits['windy']],
#                    bottom=[weather_limits['calm']],
#                    left=[min(date)], right=[max(date)],
#                    color="yellow", alpha=0.15)
#     wind_plot.quad(top=[weather_limits['calm']],
#                    bottom=[0],
#                    left=[min(date)], right=[max(date)],
#                    color="green", alpha=0.15)
#     wind_plot.circle(date, wind,
#                      size=markersize, color="#53777a", alpha=0.8)
    wind_plot.yaxis.axis_label = 'Wind Speed (kph)'
    wind_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    wind_plot.xaxis.visible = False


    ##-------------------------------------------------------------------------
    ## Rain Plot
    log.info('Build rain plot')
    rain = [w['rain'] for w in weather]
    rain_plot = figure(width=800, height=90, x_axis_type="datetime",
                       y_range=(1000,2800), x_range=temperature_plot.x_range,
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
#     rain_plot.quad(top=[2800],
#                    bottom=[weather_limits['dry']],
#                    left=[min(date)], right=[max(date)],
#                    color="green", alpha=0.15)
#     rain_plot.quad(top=[weather_limits['dry']],
#                    bottom=[weather_limits['wet']],
#                    left=[min(date)], right=[max(date)],
#                    color="yellow", alpha=0.15)
#     rain_plot.quad(top=[weather_limits['wet']],
#                    bottom=[1000],
#                    left=[min(date)], right=[max(date)],
#                    color="red", alpha=0.15)
#     rain_plot.circle(date, rain,
#                      size=markersize, color="#53777a", alpha=0.8)
    rain_plot.yaxis.axis_label = 'Rain'
    rain_plot.yaxis.formatter = NumeralTickFormatter(format="0.0a")
    rain_plot.xaxis.visible = False
#     rain_plot.yaxis.ticker = SingleIntervalTicker(desired_num_ticks=2)

    ##-------------------------------------------------------------------------
    ## Overplot Twilights
    for days in range(1,plot_ndays+1):
        twilights = get_twilights(end-timedelta(days=days), end-timedelta(days=days-1))
        for j in range(len(twilights)-1):
           temperature_plot.quad(top=[95], bottom=[25],
                                  left=[twilights[j][0]], right=[twilights[j+1][0]],
                                  color="blue", alpha=twilights[j+1][2])
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



    log.info(f"Rendering template")
    script, div = components(column(safe_plot,
                                    cloudiness_plot,
                                    wind_plot,
                                    rain_plot,
                                    temperature_plot,
                                    ))
    return flask.render_template('status.html',
                                 weather=weather,
                                 currentweather=currentweather,
                                 now=datetime.now(),
                                 utcnow=datetime.utcnow(),
                                 script=script,
                                 div=div,
                                 )


##-------------------------------------------------------------------------
## status: /<string:telescope>/images/<string:date>
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
    log.info(f"Got {len(image_list)} images")

    log.info(f"Rendering template")
    return flask.render_template('imageList.html',
                                 telescope=telescope,
                                 subject=subject,
                                 image_list=image_list,
                                 flat_count=flat_count,
                                 cal_count=cal_count,
                                 object_count=object_count,
                                 )


if __name__ == "__main__":
    app.run()
