from pathlib import Path
import configparser
import logging
import pymongo
from datetime import datetime, timedelta
import numpy as np

from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.layouts import column
from bokeh.models import DatetimeTickFormatter, NumeralTickFormatter, DatetimeTicker

from iqmon.webpage import mongo_query, get_twilights, overplot_twilights

log = logging.getLogger('FlaskLogger')


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

    markersize = 2
    currentweather = {}

    ##-------------------------------------------------------------------------
    ## Temperature Plot
    log.info('Build temperature plot')
    plot_temperature = figure(width=900, height=100, x_axis_type="datetime",
                              y_range=(25,95),
                              x_range=(end - timedelta(hours=span_hours), end),
                              )

    plot_values = cfg['Weather'].get('plot_temperature').split(',')
    query_dict = {'date': {'$gt': start, '$lt': end}}
    for plot_value in plot_values:
        collection, name = plot_value.split(':')
        log.debug(f'  Querying mongo collection {collection}')
        query_result = mongo_query(collection, query_dict)
        plot_vals = np.array([(d['date'], d[name]) for d in query_result])
        log.debug(f'  Got {len(plot_vals)} entries')
        if 'temperature' not in currentweather.keys():
            currentweather['date'] = plot_vals[-1][0]
            currentweather['temp'] = plot_vals[-1][1]
        plot_temperature.circle(plot_vals[:,0],
                               plot_vals[:,1]*1.8+32,
                               size=markersize, color="blue", alpha=0.8)
    plot_temperature.yaxis.axis_label = 'Temp (F)'
    plot_temperature.yaxis.formatter = NumeralTickFormatter(format="0,0")
    plot_temperature.yaxis.ticker = [30, 50, 70, 90]
    plot_temperature.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Cloudiness Plot
    log.info('Build cloudiness plot')
    plot_cloudiness = figure(width=900, height=100, x_axis_type="datetime",
                             y_range=(-50,5),
                             x_range=(end - timedelta(hours=span_hours), end),
                             )
    plot_values = cfg['Weather'].get('plot_cloudiness').split(',')
    query_dict = {'date': {'$gt': start, '$lt': end}}
    for plot_value in plot_values:
        collection, name = plot_value.split(':')
        log.debug(f'  Querying mongo collection {collection}')
        query_result = mongo_query(collection, query_dict)
        plot_vals = np.array([(d['date'], d[name]) for d in query_result])
        log.debug(f'  Got {len(plot_vals)} entries')
        if 'clouds' not in currentweather.keys():
            currentweather['clouds'] = plot_vals[-1][1]
        where_vcloudy = np.where(plot_vals[:,1] >= weather_limits['cloudy'])
        plot_cloudiness.circle(plot_vals[where_vcloudy][:,0],
                               plot_vals[where_vcloudy][:,1],
                               size=markersize, color="red", alpha=0.8)

        where_cloudy = np.where((plot_vals[:,1] < weather_limits['cloudy'])\
                                & (plot_vals[:,1] <= weather_limits['clear']))
        plot_cloudiness.circle(plot_vals[where_cloudy][:,0],
                               plot_vals[where_cloudy][:,1],
                               size=markersize, color="orange", alpha=0.8)

        where_clear = np.where(plot_vals[:,1] < weather_limits['clear'])
        plot_cloudiness.circle(plot_vals[where_clear][:,0],
                               plot_vals[where_clear][:,1],
                               size=markersize, color="green", alpha=0.8)

    plot_cloudiness.yaxis.axis_label = 'Cloudiness (C)'
    plot_cloudiness.yaxis.formatter = NumeralTickFormatter(format="0,0")
    plot_cloudiness.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Wind Plot
    log.info('Build wind plot')
    plot_wind_speed = figure(width=900, height=100, x_axis_type="datetime",
                             y_range=(-3,85),
                             x_range=(end - timedelta(hours=span_hours), end),
                             )
    plot_values = cfg['Weather'].get('plot_wind_speed').split(',')
    query_dict = {'date': {'$gt': start, '$lt': end}}
    for plot_value in plot_values:
        collection, name = plot_value.split(':')
        log.debug(f'  Querying mongo collection {collection}')
        query_result = mongo_query(collection, query_dict)
        plot_vals = np.array([(d['date'], d[name]) for d in query_result])
        log.debug(f'  Got {len(plot_vals)} entries')
        if 'wind' not in currentweather.keys():
            currentweather['wind'] = plot_vals[-1][1]

        where_vwindy = np.where(plot_vals[:,1] >= weather_limits['windy'])
        plot_wind_speed.circle(plot_vals[where_vwindy][:,0],
                               plot_vals[where_vwindy][:,1],
                               size=markersize, color="red", alpha=0.8)

        where_windy = np.where((plot_vals[:,1] < weather_limits['windy'])\
                                & (plot_vals[:,1] <= weather_limits['calm']))
        plot_wind_speed.circle(plot_vals[where_windy][:,0],
                               plot_vals[where_windy][:,1],
                               size=markersize, color="orange", alpha=0.8)

        where_calm = np.where(plot_vals[:,1] < weather_limits['calm'])
        plot_wind_speed.circle(plot_vals[where_calm][:,0],
                               plot_vals[where_calm][:,1],
                               size=markersize, color="green", alpha=0.8)

    plot_wind_speed.yaxis.axis_label = 'Wind (kph)'
    plot_wind_speed.yaxis.formatter = NumeralTickFormatter(format="0,0")
    plot_wind_speed.yaxis.ticker = [0, 20, 40, 60, 80]
    plot_wind_speed.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Rain Plot
    log.info('Build rain plot')
    plot_rain = figure(width=900, height=60, x_axis_type="datetime",
                       y_range=(1000,2800),
                       x_range=(end - timedelta(hours=span_hours), end),
                       )
    plot_values = cfg['Weather'].get('plot_rain').split(',')
    query_dict = {'date': {'$gt': start, '$lt': end}}
    for plot_value in plot_values:
        collection, name = plot_value.split(':')
        log.debug(f'  Querying mongo collection {collection}')
        query_result = mongo_query(collection, query_dict)
        plot_vals = np.array([(d['date'], d[name]) for d in query_result])
        log.debug(f'  Got {len(plot_vals)} entries')
        if 'rain' not in currentweather.keys():
            currentweather['rain'] = plot_vals[-1][1]

        where_dry = np.where(plot_vals[:,1] >= weather_limits['dry'])
        plot_rain.circle(plot_vals[where_dry][:,0],
                         plot_vals[where_dry][:,1],
                         size=markersize, color="green", alpha=0.8)

        where_wet = np.where((plot_vals[:,1] < weather_limits['dry'])\
                             & (plot_vals[:,1] <= weather_limits['wet']))
        plot_rain.circle(plot_vals[where_wet][:,0],
                         plot_vals[where_wet][:,1],
                         size=markersize, color="orange", alpha=0.8)

        where_rain = np.where(plot_vals[:,1] < weather_limits['wet'])
        plot_rain.circle(plot_vals[where_rain][:,0],
                         plot_vals[where_rain][:,1],
                         size=markersize, color="red", alpha=0.8)

    plot_rain.yaxis.axis_label = 'Rain'
    plot_rain.yaxis.formatter = NumeralTickFormatter(format="0.0a")
    plot_rain.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Safe Plot
    log.info('Build safe plot')
    plot_safe = figure(width=900, height=50, x_axis_type="datetime",
                       y_range=(-0.2,1.2),
                       x_range=(end - timedelta(hours=span_hours), end),
                       )
    plot_values = cfg['Weather'].get('plot_safe').split(',')
    query_dict = {'date': {'$gt': start, '$lt': end}}

    for plot_value in plot_values:
        collection, name = plot_value.split(':')
        log.debug(f'  Querying mongo collection {collection}')
        query_result = mongo_query(collection, query_dict)
        plot_vals = np.array([(d['date'], d[name]) for d in query_result])
        log.debug(f'  Got {len(plot_vals)} entries')
        if 'safe' not in currentweather.keys():
            currentweather['safe'] = plot_vals[-1][1]
        where_safe = np.where(plot_vals[:,1] == True)
        plot_safe.circle(plot_vals[where_safe][:,0],
                         plot_vals[where_safe][:,1],
                         size=markersize, color="green", alpha=0.8)
        where_unsafe = np.where(plot_vals[:,1] != True)
        plot_safe.circle(plot_vals[where_unsafe][:,0],
                         plot_vals[where_unsafe][:,1],
                         size=markersize, color="red", alpha=0.8)
    plot_safe.yaxis.axis_label = 'Safe'
    plot_safe.xaxis.axis_label = 'Time (UT)'
    plot_safe.yaxis.formatter = NumeralTickFormatter(format="0,0")
    plot_safe.yaxis.ticker = [0,1]
    plot_safe.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Telescope Status Plot
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

    ## Build Telescope Status plot
    log.info('Build Telescope Status plot')
    dome = [s['dome_numerical_status'] for s in telstatus]
    dome_date = [s['date'] for s in telstatus]
    dome_plot = figure(width=900, height=100, x_axis_type="datetime",
                       y_range=(-0.2,1.2),
                       x_range=(end - timedelta(hours=span_hours), end),
                       )
    open_date = [dome_date[i] for i,d in enumerate(dome) if d < 0.5]
    open_dome = [d for i,d in enumerate(dome) if d < 0.5]
    closed_date = [dome_date[i] for i,d in enumerate(dome) if d >= 0.5]
    closed_dome = [d for i,d in enumerate(dome) if d >= 0.5]
    dome_plot.line(closed_date, closed_dome, line_width=4, color="black")
    dome_plot.line(open_date, open_dome, line_width=4, color="green")
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
    ## Render
    log.info(f"Overplotting twilights")
    plot_info_list = [('plot_temperature', plot_temperature, 25, 95),
                      ('plot_cloudiness', plot_cloudiness, 5, -50),
                      ('plot_wind_speed', plot_wind_speed, 100, -3),
                      ('plot_rain', plot_rain, 2800, 1000),
                      ('plot_safe', plot_safe, 1.2, -0.2),
                      ]
    plot_column_list = []
    plot_twilights_list = []
    for i,plot_info in enumerate(plot_info_list):
        if cfg['Weather'].get(plot_info[0], None) is not None:
            if i != 0:
                plot_info[1].x_range = plot_info_list[0][1].x_range
            plot_column_list.append(plot_info[1])
            plot_twilights_list.append(plot_info)
    overplot_twilights(plot_twilights_list, end, plot_ndays=plot_ndays, log=log)
    plot_column_list.append(dome_plot)

    log.info(f"Rendering bokeh plot")
    script, div = components(column(plot_column_list))

    return script, div, currentweather
