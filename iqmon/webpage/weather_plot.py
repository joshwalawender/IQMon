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

from iqmon import get_config
from iqmon.webpage import mongo_query, overplot_twilights

log = logging.getLogger('FlaskLogger')


##-------------------------------------------------------------------------
## Function: generate_weather_plot
##-------------------------------------------------------------------------
def generate_weather_plot(cfg, date=None, plot_ndays=1, span_hours=24):
    telescope = cfg['Telescope'].get('name')

    log.info(f"Querying weather limits")
    weather_limits = mongo_query('weather_limits', {}, cfg)[0]

    if date is None:
        end = datetime.utcnow()
        start = end - timedelta(days=plot_ndays)
    else:
        start = datetime.strptime(date, '%Y%m%dUT')
        end = start + timedelta(days=plot_ndays)

    markersize = 2

    ##-------------------------------------------------------------------------
    ## Temperature Plot
    if cfg['Weather'].get('plot_temperature', None) is not None:
        log.info('Build temperature plot')
        limit_string = cfg['Weather'].get(f'plot_temperature_limits', '25,95')
        ymin,ymax = limit_string.split(',')
        plot_temperature = figure(width=900, height=120, x_axis_type="datetime",
                                  y_range=(float(ymin),float(ymax)),
                                  x_range=(end - timedelta(hours=span_hours), end),
                                  )
        plot_values = cfg['Weather'].get('plot_temperature').split(',')
        query_dict = {'date': {'$gt': start, '$lt': end}}
        colors = ["blue", "black", "black"]
        alphas = [0.8, 0.4, 0.4]
        for i,plot_value in enumerate(plot_values):
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, cfg)
            plot_vals = np.array([(d['date'], d[name]) for d in query_result if name in d.keys()])
            if len(plot_vals) == 0:
                log.warning(f'Found 0 data points for {collection}:{name}')
            temperature_units = cfg[collection].get('temperature_units', 'F')
            if temperature_units == 'C':
                log.info('  Converting temperature plot from C to F')
                plot_vals[:,1] = plot_vals[:,1]*1.8+32

            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                plot_temperature.circle(plot_vals[:,0],
                                        plot_vals[:,1],
                                        legend_label=f"{name}",
                                        color=colors[i],
                                        fill_alpha=alphas[i],
                                        line_alpha=alphas[i],
                                        size=markersize)
        if len(plot_values) > 1:
            plot_temperature.legend.location = "top_left"
            plot_temperature.legend.margin = 0
            plot_temperature.legend.padding = 0
            plot_temperature.legend.spacing = 0
            plot_temperature.legend.label_text_font_size = '8pt'
        else:
            plot_temperature.legend.visible = False
        plot_temperature.yaxis.axis_label = 'Temp (F)'
        plot_temperature.yaxis.formatter = NumeralTickFormatter(format="0,0")
        plot_temperature.yaxis.ticker = [10,30, 50, 70, 90, 110]
        plot_temperature.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Humidity Plot
    if cfg['Weather'].get('plot_humidity', None) is not None:
        log.info('Build humidity plot')
        limit_string = cfg['Weather'].get(f'plot_humidity_limits', '40,100')
        ymin,ymax = limit_string.split(',')
        plot_humidity = figure(width=900, height=120, x_axis_type="datetime",
                                 y_range=(float(ymin),float(ymax)),
                                 x_range=(end - timedelta(hours=span_hours), end),
                                 )
        plot_values = cfg['Weather'].get('plot_humidity').split(',')
        query_dict = {'date': {'$gt': start, '$lt': end}}
        for i,plot_value in enumerate(plot_values):
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, cfg)
            plot_vals = np.array([(d['date'], d[name]) for d in query_result if name in d.keys()])
            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) == 0:
                log.warning(f'Found 0 data points for {collection}:{name}')
            else:
                if i == 0:
                    where_vhumid = np.where(plot_vals[:,1] >= weather_limits['very humid'])
                    plot_humidity.circle(plot_vals[where_vhumid][:,0],
                                           plot_vals[where_vhumid][:,1],
                                           size=markersize, color="red",
                                           line_alpha=0.8, fill_alpha=0.8)

                    where_humid = np.where((plot_vals[:,1] < weather_limits['very humid'])\
                                            & (plot_vals[:,1] >= weather_limits['humid']))
                    plot_humidity.circle(plot_vals[where_humid][:,0],
                                           plot_vals[where_humid][:,1],
                                           size=markersize, color="orange",
                                           line_alpha=0.8, fill_alpha=0.8)

                    where_nothumid = np.where(plot_vals[:,1] < weather_limits['not humid'])
                    plot_humidity.circle(plot_vals[where_nothumid][:,0],
                                           plot_vals[where_nothumid][:,1],
                                           legend_label=f"{name}",
                                           size=markersize, color="green",
                                           line_alpha=0.8, fill_alpha=0.8)
                else:
                    plot_humidity.circle(plot_vals[:,0],
                                         plot_vals[:,1],
                                         legend_label=f"{name}",
                                         color='black',
                                         line_alpha=0.4, fill_alpha=0.4,
                                         size=markersize)
        if len(plot_values) > 1:
            plot_humidity.legend.location = "top_left"
            plot_humidity.legend.margin = 0
            plot_humidity.legend.padding = 0
            plot_humidity.legend.spacing = 0
            plot_humidity.legend.label_text_font_size = '8pt'
        else:
            plot_humidity.legend.visible = False
        plot_humidity.yaxis.axis_label = 'Humidity (%)'
        plot_humidity.yaxis.formatter = NumeralTickFormatter(format="0,0")
        plot_humidity.yaxis.ticker = [0, 20, 40, 60, 80, 100]
        plot_humidity.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Cloudiness Plot
    if cfg['Weather'].get('plot_cloudiness', None) is not None:
        log.info('Build cloudiness plot')
        limit_string = cfg['Weather'].get(f'plot_cloudiness_limits', '-50,10')
        ymin,ymax = limit_string.split(',')
        plot_cloudiness = figure(width=900, height=120, x_axis_type="datetime",
                                 y_range=(float(ymin),float(ymax)),
                                 x_range=(end - timedelta(hours=span_hours), end),
                                 )
        plot_values = cfg['Weather'].get('plot_cloudiness').split(',')
        query_dict = {'date': {'$gt': start, '$lt': end}}
        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, cfg)
            plot_vals = np.array([(d['date'], d[name]) for d in query_result if name in d.keys()])
            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                where_vcloudy = np.where(plot_vals[:,1] >= weather_limits['cloudy'])
                plot_cloudiness.circle(plot_vals[where_vcloudy][:,0],
                                       plot_vals[where_vcloudy][:,1],
                                       size=markersize, color="red",
                                       line_alpha=0.8, fill_alpha=0.8)

                where_cloudy = np.where((plot_vals[:,1] < weather_limits['cloudy'])\
                                        & (plot_vals[:,1] >= weather_limits['clear']))
                plot_cloudiness.circle(plot_vals[where_cloudy][:,0],
                                       plot_vals[where_cloudy][:,1],
                                       size=markersize, color="orange",
                                       line_alpha=0.8, fill_alpha=0.8)

                where_clear = np.where(plot_vals[:,1] < weather_limits['clear'])
                plot_cloudiness.circle(plot_vals[where_clear][:,0],
                                       plot_vals[where_clear][:,1],
                                       size=markersize, color="green",
                                       line_alpha=0.8, fill_alpha=0.8)
        plot_cloudiness.yaxis.axis_label = 'Cloudiness (C)'
        plot_cloudiness.yaxis.formatter = NumeralTickFormatter(format="0,0")
        plot_cloudiness.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Wind Plot
    if cfg['Weather'].get('plot_wind_speed', None) is not None:
        log.info('Build wind plot')
        limit_string = cfg['Weather'].get(f'plot_wind_speed_limits', '-3,85')
        ymin,ymax = limit_string.split(',')
        plot_wind_speed = figure(width=900, height=120, x_axis_type="datetime",
                                 y_range=(float(ymin),float(ymax)),
                                 x_range=(end - timedelta(hours=span_hours), end),
                                 )
        plot_values = cfg['Weather'].get('plot_wind_speed').split(',')
        query_dict = {'date': {'$gt': start, '$lt': end}}
        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, cfg)
            plot_vals = np.array([(d['date'], d[name]) for d in query_result if name in d.keys()])
            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                where_vwindy = np.where(plot_vals[:,1] >= weather_limits['windy'])
                plot_wind_speed.circle(plot_vals[where_vwindy][:,0],
                                       plot_vals[where_vwindy][:,1],
                                       size=markersize, color="red",
                                       line_alpha=0.8, fill_alpha=0.8)

                where_windy = np.where((plot_vals[:,1] < weather_limits['windy'])\
                                        & (plot_vals[:,1] >= weather_limits['calm']))
                plot_wind_speed.circle(plot_vals[where_windy][:,0],
                                       plot_vals[where_windy][:,1],
                                       size=markersize, color="orange",
                                       line_alpha=0.8, fill_alpha=0.8)

                where_calm = np.where(plot_vals[:,1] < weather_limits['calm'])
                plot_wind_speed.circle(plot_vals[where_calm][:,0],
                                       plot_vals[where_calm][:,1],
                                       size=markersize, color="green",
                                       line_alpha=0.8, fill_alpha=0.8)
        plot_wind_speed.yaxis.axis_label = 'Wind (kph)'
        plot_wind_speed.yaxis.formatter = NumeralTickFormatter(format="0,0")
        plot_wind_speed.yaxis.ticker = [0, 20, 40, 60, 80]
        plot_wind_speed.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Rain Plot
    if cfg['Weather'].get('plot_rain', None) is not None:
        log.info('Build rain plot')
        limit_string = cfg['Weather'].get(f'plot_rain_limits', '500,2800')
        ymin,ymax = limit_string.split(',')
        plot_rain = figure(width=900, height=60, x_axis_type="datetime",
                                  y_range=(float(ymin),float(ymax)),
                           x_range=(end - timedelta(hours=span_hours), end),
                           )
        plot_values = cfg['Weather'].get('plot_rain').split(',')
        query_dict = {'date': {'$gt': start, '$lt': end}}
        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, cfg)
            plot_vals = np.array([(d['date'], d[name]) for d in query_result if name in d.keys()])
            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                where_dry = np.where(plot_vals[:,1] >= weather_limits['dry'])
                plot_rain.circle(plot_vals[where_dry][:,0],
                                 plot_vals[where_dry][:,1],
                                 size=markersize, color="green",
                                 line_alpha=0.8, fill_alpha=0.8)

                where_wet = np.where((plot_vals[:,1] < weather_limits['dry'])\
                                     & (plot_vals[:,1] >= weather_limits['wet']))
                plot_rain.circle(plot_vals[where_wet][:,0],
                                 plot_vals[where_wet][:,1],
                                 size=markersize, color="orange",
                                 line_alpha=0.8, fill_alpha=0.8)

                where_rain = np.where(plot_vals[:,1] < weather_limits['wet'])
                plot_rain.circle(plot_vals[where_rain][:,0],
                                 plot_vals[where_rain][:,1],
                                 size=markersize, color="red",
                                 line_alpha=0.8, fill_alpha=0.8)
        plot_rain.yaxis.axis_label = 'Rain'
        plot_rain.yaxis.formatter = NumeralTickFormatter(format="0.0a")
        plot_rain.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Safe Plot
    if cfg['Weather'].get('plot_safe', None) is not None:
        log.info('Build safe plot')
        plot_safe = figure(width=900, height=60, x_axis_type="datetime",
                           y_range=(-0.2,1.2),
                           x_range=(end - timedelta(hours=span_hours), end),
                           )
        plot_values = cfg['Weather'].get('plot_safe').split(',')
        query_dict = {'date': {'$gt': start, '$lt': end}}

        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, cfg)
            plot_vals = np.array([(d['date'], d[name]) for d in query_result if name in d.keys()])
            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                where_safe = np.where(plot_vals[:,1] == True)
                plot_safe.circle(plot_vals[where_safe][:,0],
                                 plot_vals[where_safe][:,1],
                                 size=markersize, color="green",
                                 line_alpha=0.8, fill_alpha=0.8)
                where_unsafe = np.where(plot_vals[:,1] != True)
                plot_safe.circle(plot_vals[where_unsafe][:,0],
                                 plot_vals[where_unsafe][:,1],
                                 size=markersize, color="red",
                                 line_alpha=0.8, fill_alpha=0.8)
        plot_safe.yaxis.axis_label = 'Safe'
        plot_safe.xaxis.axis_label = 'Time (UT)'
        plot_safe.yaxis.formatter = NumeralTickFormatter(format="0,0")
        plot_safe.yaxis.ticker = [0,1]
        plot_safe.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Telescope Status Plot
    ##-------------------------------------------------------------------------
    if cfg['Weather'].getboolean('plot_dome', False) is True:
        ## Telescope Status Query
        log.info(f"Querying telescope status database")
        query_dict = {'date': {'$gt': start, '$lt': end}}
        telstatus = mongo_query(f'{telescope}status', query_dict, cfg)
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
        iqmon = mongo_query('iqmon', query_dict, cfg)
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
        dome_plot = figure(width=900, height=60, x_axis_type="datetime",
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
                         size=markersize, color="blue",
                         line_alpha=0.8, fill_alpha=0.8)
        dome_plot.circle(iqmon_cal_dates, iqmon_cal_alt,
                         size=markersize, color="black",
                         line_alpha=0.8, fill_alpha=0.8)
        dome_plot.circle(iqmon_flat_dates, iqmon_flat_alt,
                         size=markersize, color="yellow",
                         line_alpha=0.8, fill_alpha=0.8)
        dome_plot.yaxis.axis_label = f'{cfg["Telescope"].get("name")}'
        dome_plot.xaxis.axis_label = 'Time (UT)'
        dome_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
        dome_plot.yaxis.ticker = [0,1]
        dome_plot.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Render
    log.info(f"Overplotting twilights")

    plot_names = ['plot_temperature', 'plot_humidity', 'plot_cloudiness',
                  'plot_wind_speed', 'plot_rain', 'plot_safe']
    plot_info_list = []
    for name in plot_names:
        if cfg['Weather'].get(name, None) is not None:
            log.info(f"  {name}")
            limit_string = cfg['Weather'].get(f'{name}_limits', None)
            if limit_string is not None:
                ymin,ymax = limit_string.split(',')
            else:
                ymin,ymax = 0,1
            plot_info_list.append([name, eval(name), float(ymax), float(ymin)])

    plot_column_list = []
    plot_twilights_list = []
    for i,plot_info in enumerate(plot_info_list):
        if cfg['Weather'].get(plot_info[0], None) is not None:
            if i != 0:
                plot_info[1].x_range = plot_info_list[0][1].x_range
            plot_column_list.append(plot_info[1])
            plot_twilights_list.append(plot_info)
    overplot_twilights(plot_twilights_list, end, cfg, plot_ndays=plot_ndays, log=log)

    if cfg['Weather'].getboolean('plot_dome', False) is True:
        plot_column_list.append(dome_plot)

    # Add time log
    plot_column_list[-1].plot_height += 50
    plot_column_list[-1].xaxis.visible = True
    plot_column_list[-1].xaxis.formatter = DatetimeTickFormatter(hourmin=['%H:%M'])
    plot_column_list[-1].xaxis.ticker = DatetimeTicker(desired_num_ticks=24)
    plot_column_list[-1].xaxis.axis_label = 'UT Time'

    log.info(f"Rendering bokeh plot for {plot_column_list}")
    script, div = components(column(plot_column_list))

    return script, div
