from pathlib import Path
import configparser
import logging
import flask
import pymongo
from datetime import datetime, timedelta
import numpy as np

from iqmon.webpage import mongo_query, get_twilights, overplot_twilights
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
## status: /
##-------------------------------------------------------------------------
@app.route("/")
def hello():
    log.info(f'Building {__name__} hello')
    return status('V5')


##-------------------------------------------------------------------------
## status: /<string:telescope>
##-------------------------------------------------------------------------
@app.route("/<string:telescope>/")
def status(telescope):
    log.info(f'')
    log.info(f'Building {__name__} status')
    tick = datetime.utcnow()

    cfg_path = Path(__file__).parent.parent / 'configs' / 'pipeline.cfg'
    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)

    script, div, currentweather = generate_weather_plot(telescope, plot_ndays=2, span_hours=12)

    ## Format currentweather
    log.info(f"Querying weather limits")
    query_result = mongo_query('weather_limits', {})
    weather_limits = query_result[0]
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


    ## Telescope Status Query
    log.info(f"Querying telescope status database")
    query_dict = {'date': {'$gt': tick-timedelta(minutes=5), '$lt': tick}}
    query_result = mongo_query(f'{telescope}status', query_dict)
    telstatus = [d for d in query_result]
    log.info(f"  Got {len(telstatus)} data points")
    shutter_values = {0: 0, 1: 1, 2: 0, 3: 1, 4: 4}
    for i,d in enumerate(telstatus):
        if d['dome_shutterstatus'] == 4 and i > 0:
            telstatus[i]['dome_numerical_status'] = telstatus[i-1]['dome_numerical_status']
        else:
            telstatus[i]['dome_numerical_status'] = shutter_values[d['dome_shutterstatus']]
    ## Format currentstatus
    dome_string = {0: 'Open', 1: 'Closed', 2: 'Opening', 3: 'Closing', 4: 'Unknown'}
    dome_color = {0: 'green', 1: 'red', 2: 'orange', 3: 'orange', 4: 'black'}
    try:
        currentstatus = telstatus[-1]
    except:
        currentstatus = {}
    else:
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
                                   image=cfg['WebPage'].get('image', ''),
                                   image_link=cfg['WebPage'].get('image_link', ''),
                                   image_title=cfg['WebPage'].get('image_title', ''),
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

    script, div, currentweather = generate_weather_plot(telescope, date=date)

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

    script, div = generate_iqmon_plot(telescope, date=date)

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


if __name__ == "__main__":
    app.run()
