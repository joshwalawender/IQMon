from pathlib import Path
import configparser
import logging
import pymongo
from datetime import datetime, timedelta

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
    plot_info = [(cloudiness_plot, 5, -45),
                 (wind_plot, 100, 0),
                 (rain_plot, 2800, 1000),
                 (safe_plot, 1.2, -0.2),
                 ]
    overplot_twilights(plot_info, end, plot_ndays=plot_ndays)

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
