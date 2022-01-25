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

from iqmon.webpage import mongo_query, overplot_twilights

log = logging.getLogger('FlaskLogger')


##-------------------------------------------------------------------------
## Function: generate_weather_plot
##-------------------------------------------------------------------------
def generate_weather_plot(webcfg, telcfg, date=None, querybuffer_ndays=1, span_hours=24):
    log.info(f"Querying weather limits")
    weather_limits = mongo_query('weather_limits', {}, webcfg)[0]

    if date is None:
        if span_hours > 24*querybuffer_ndays:
            querybuffer_ndays = int(np.ceil(span_hours/24))
        end = datetime.utcnow()
        query_end = end
        start = end - timedelta(hours=span_hours)
        query_start = start - timedelta(days=querybuffer_ndays)
        query_days = querybuffer_ndays
    else:
        if span_hours > 24*querybuffer_ndays:
            querybuffer_ndays = int(np.ceil(span_hours/24))
        start = datetime.strptime(date, '%Y%m%dUT')
        end = start + timedelta(hours=span_hours)
        query_start = start - timedelta(days=querybuffer_ndays)
        query_end = end + timedelta(days=querybuffer_ndays)
        query_days = 2*querybuffer_ndays

    markersize = 2

    ##-------------------------------------------------------------------------
    ## Temperature Plot
    if webcfg['Weather'].get('plot_temperature', None) is not None:
        log.info('Build temperature plot')
        limit_string = webcfg['Weather'].get(f'plot_temperature_limits', '25,95')
        ymin,ymax = limit_string.split(',')
        height = webcfg['Weather'].getint('plot_temperature_height', 120)
        plot_temperature = figure(width=900, height=height, x_axis_type="datetime",
                                  y_range=(float(ymin),float(ymax)),
                                  x_range=(end - timedelta(hours=span_hours), end),
                                  )
        plot_values = webcfg['Weather'].get('plot_temperature').split(',')
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
        colors = ["blue", "black", "black"]
        alphas = [0.8, 0.4, 0.4]
        for i,plot_value in enumerate(plot_values):
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, webcfg)

            plot_vals = []
            for d in query_result:
                if name in d.keys():
                    if d.get(f"{name} units", None) == 'C':
                        plot_vals.append((d['date'], d[name]*1.8+32))
                    elif d.get(f"{name} unit", None) == 'C':
                        plot_vals.append((d['date'], d[name]*1.8+32))
                    elif d.get("temperature units", None) == 'C':
                        plot_vals.append((d['date'], d[name]*1.8+32))
                    else:
                        plot_vals.append((d['date'], d[name]))
            plot_vals = np.array(plot_vals)

            if len(plot_vals) == 0:
                log.warning(f'Found 0 data points for {collection}:{name}')

            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                plot_temperature.circle(plot_vals[:,0],
                                        plot_vals[:,1],
                                        legend_label=f"{plot_value}",
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
        plot_temperature.yaxis.ticker = [val for val in range(-50, 130, 20)]
        plot_temperature.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Humidity Plot
    if webcfg['Weather'].get('plot_humidity', None) is not None:
        log.info('Build humidity plot')
        limit_string = webcfg['Weather'].get(f'plot_humidity_limits', '40,100')
        ymin,ymax = limit_string.split(',')
        height = webcfg['Weather'].getint('plot_humidity_height', 120)
        plot_humidity = figure(width=900, height=height, x_axis_type="datetime",
                                 y_range=(float(ymin),float(ymax)),
                                 x_range=(end - timedelta(hours=span_hours), end),
                                 )
        plot_values = webcfg['Weather'].get('plot_humidity').split(',')
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
        for i,plot_value in enumerate(plot_values):
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, webcfg)
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

                    where_nothumid = np.where(plot_vals[:,1] < weather_limits['humid'])
                    plot_humidity.circle(plot_vals[where_nothumid][:,0],
                                           plot_vals[where_nothumid][:,1],
                                           legend_label=f"{plot_value}",
                                           size=markersize, color="green",
                                           line_alpha=0.8, fill_alpha=0.8)
                else:
                    plot_humidity.circle(plot_vals[:,0],
                                         plot_vals[:,1],
                                         legend_label=f"{plot_value}",
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
        plot_humidity.yaxis.ticker = [val for val in range(0, 120, 20)]
        plot_humidity.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Cloudiness Plot
    if webcfg['Weather'].get('plot_cloudiness', None) is not None:
        log.info('Build cloudiness plot')
        limit_string = webcfg['Weather'].get(f'plot_cloudiness_limits', '-50,10')
        ymin,ymax = limit_string.split(',')
        height = webcfg['Weather'].getint('plot_cloudiness_height', 120)
        plot_cloudiness = figure(width=900, height=height, x_axis_type="datetime",
                                 y_range=(float(ymin),float(ymax)),
                                 x_range=(end - timedelta(hours=span_hours), end),
                                 )
        plot_values = webcfg['Weather'].get('plot_cloudiness').split(',')
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, webcfg)
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
    if webcfg['Weather'].get('plot_wind_speed', None) is not None:
        log.info('Build wind plot')
        limit_string = webcfg['Weather'].get(f'plot_wind_speed_limits', '-3,85')
        ymin,ymax = limit_string.split(',')
        height = webcfg['Weather'].getint('plot_wind_speed_height', 120)
        plot_wind_speed = figure(width=900, height=height, x_axis_type="datetime",
                                 y_range=(float(ymin),float(ymax)),
                                 x_range=(end - timedelta(hours=span_hours), end),
                                 )
        plot_values = webcfg['Weather'].get('plot_wind_speed').split(',')
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, webcfg)

            plot_vals = []
            for d in query_result:
                if name in d.keys():
                    if d.get(f"{name} units", None) == 'kph':
                        plot_vals.append((d['date'], d[name]*1.61))
                    elif d.get(f"{name} unit", None) == 'kph':
                        plot_vals.append((d['date'], d[name]*1.61))
                    elif d.get("wind speed units", None) == 'kph':
                        plot_vals.append((d['date'], d[name]*1.61))
                    else:
                        plot_vals.append((d['date'], d[name]))
            plot_vals = np.array(plot_vals)

            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                where_vwindy = np.where(plot_vals[:,1] >= weather_limits['very windy'])
                plot_wind_speed.circle(plot_vals[where_vwindy][:,0],
                                       plot_vals[where_vwindy][:,1],
                                       size=markersize, color="red",
                                       line_alpha=0.8, fill_alpha=0.8)

                where_windy = np.where((plot_vals[:,1] < weather_limits['very windy'])\
                                        & (plot_vals[:,1] >= weather_limits['windy']))
                plot_wind_speed.circle(plot_vals[where_windy][:,0],
                                       plot_vals[where_windy][:,1],
                                       size=markersize, color="orange",
                                       line_alpha=0.8, fill_alpha=0.8)

                where_calm = np.where(plot_vals[:,1] < weather_limits['windy'])
                plot_wind_speed.circle(plot_vals[where_calm][:,0],
                                       plot_vals[where_calm][:,1],
                                       legend_label=f"{plot_value}",
                                       size=markersize, color="green",
                                       line_alpha=0.8, fill_alpha=0.8)

        if webcfg['Weather'].get('plot_wind_gust', None) is not None:
            plot_values = webcfg['Weather'].get('plot_wind_gust').split(',')
            query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
            for plot_value in plot_values:
                collection, name = plot_value.split(':')
                log.debug(f'  Querying mongo collection {collection}')
                query_result = mongo_query(collection, query_dict, webcfg)
                plot_vals = np.array([(d['date'], d[name]) for d in query_result if name in d.keys()])
                log.debug(f'  Got {len(plot_vals)} entries')
                if len(plot_vals) > 0:
                    plot_wind_speed.line(plot_vals[:,0],
                                         plot_vals[:,1],
                                         legend_label=f"{plot_value}",
                                         line_width=markersize, color="black",
                                         line_alpha=0.3)
        if len(plot_values) > 1 or webcfg['Weather'].get('plot_wind_gust', None) is not None:
            plot_wind_speed.legend.location = "top_left"
            plot_wind_speed.legend.margin = 0
            plot_wind_speed.legend.padding = 0
            plot_wind_speed.legend.spacing = 0
            plot_wind_speed.legend.label_text_font_size = '8pt'
        else:
            plot_wind_speed.legend.visible = False
        plot_wind_speed.yaxis.axis_label = 'Wind (mph)'
        plot_wind_speed.yaxis.formatter = NumeralTickFormatter(format="0,0")
        plot_wind_speed.yaxis.ticker = [val for val in range(0, 120, 20)]
        plot_wind_speed.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Rain Plot
    if webcfg['Weather'].get('plot_rain', None) is not None:
        log.info('Build rain plot')
        limit_string = webcfg['Weather'].get(f'plot_rain_limits', '500,2800')
        ymin,ymax = limit_string.split(',')
        height = webcfg['Weather'].getint('plot_rain_height', 60)
        plot_rain = figure(width=900, height=height, x_axis_type="datetime",
                                  y_range=(float(ymin),float(ymax)),
                           x_range=(end - timedelta(hours=span_hours), end),
                           )
        plot_values = webcfg['Weather'].get('plot_rain').split(',')
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, webcfg)
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
        plot_rain.yaxis.ticker = [val for val in range(500, 3000, 1000)]
        plot_rain.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Pressure Plot
    if webcfg['Weather'].get('plot_pressure', None) is not None:
        log.info('Build pressure plot')
        limit_string = webcfg['Weather'].get(f'plot_pressure_limits', '800,1050')
        ymin,ymax = limit_string.split(',')
        height = webcfg['Weather'].getint('plot_pressure_height', 120)
        plot_pressure = figure(width=900, height=height, x_axis_type="datetime",
                               y_range=(float(ymin),float(ymax)),
                               x_range=(end - timedelta(hours=span_hours), end),
                               )
        plot_values = webcfg['Weather'].get('plot_pressure').split(',')
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
        colors = ["blue", "black", "black"]
        alphas = [0.8, 0.4, 0.4]
        for i,plot_value in enumerate(plot_values):
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, webcfg)

            plot_vals = []
            for d in query_result:
                if name in d.keys():
                    if d.get(f"{name} units", None) == 'mbar':
                        plot_vals.append((d['date'], d[name]))
                    elif d.get(f"{name} units", None) == 'inHg':
                        plot_vals.append((d['date'], d[name]*33.86389))
                    else:
                        plot_vals.append((d['date'], d[name]*33.86389))
            plot_vals = np.array(plot_vals)

            if len(plot_vals) == 0:
                log.warning(f'Found 0 data points for {collection}:{name}')

            log.debug(f'  Got {len(plot_vals)} entries')
            if len(plot_vals) > 0:
                plot_pressure.circle(plot_vals[:,0],
                                     plot_vals[:,1],
                                     legend_label=f"{plot_value}",
                                     color=colors[i],
                                     fill_alpha=alphas[i],
                                     line_alpha=alphas[i],
                                     size=markersize)
        if len(plot_values) > 1:
            plot_pressure.legend.location = "top_left"
            plot_pressure.legend.margin = 0
            plot_pressure.legend.padding = 0
            plot_pressure.legend.spacing = 0
            plot_pressure.legend.label_text_font_size = '8pt'
        else:
            plot_pressure.legend.visible = False
        plot_pressure.yaxis.axis_label = 'Pressure (mbar)'
        plot_pressure.yaxis.formatter = NumeralTickFormatter(format="0,0")
        plot_pressure.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Safe Plot
    if webcfg['Weather'].get('plot_safe', None) is not None:
        log.info('Build safe plot')
        height = webcfg['Weather'].getint('plot_safe_height', 60)
        plot_safe = figure(width=900, height=height, x_axis_type="datetime",
                           y_range=(-0.2,1.2),
                           x_range=(end - timedelta(hours=span_hours), end),
                           )
        plot_values = webcfg['Weather'].get('plot_safe').split(',')
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}

        for plot_value in plot_values:
            collection, name = plot_value.split(':')
            log.debug(f'  Querying mongo collection {collection}')
            query_result = mongo_query(collection, query_dict, webcfg)
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
    ## Dome Status Plot
    ##-------------------------------------------------------------------------
    if webcfg['Weather'].getboolean('plot_dome', False) is True and telcfg is not None:
        telescope = telcfg['Telescope'].get('name')
        ## Telescope Status Query
        log.info(f"Querying dome status database")
        query_dict = {'date': {'$gt': query_start, '$lt': query_end}}
        domestatus = mongo_query(f'{telescope}_dome', query_dict, telcfg)
        log.info(f"  Got {len(domestatus)} data points")
        # 0=open, 1=closed, 2=opening, 3=closing
        shutter_values = {0: 0, 1: 1, 2: 0.25, 3: 0.75, 4: 4}
        for i,d in enumerate(domestatus):
            if d['shutterstatus'] == 4 and i > 0:
                domestatus[i]['open_closed'] = domestatus[i-1]['open_closed']
            else:
                domestatus[i]['open_closed'] = shutter_values[d['shutterstatus']]

        ## IQMon Query
        log.info(f"Querying IQMon results database")
        query_dict = {'telescope': telescope,
                      'date': {'$gt': query_start, '$lt': query_end}}
        iqmon = mongo_query('iqmon', query_dict, telcfg)
        log.info(f"  Got {len(iqmon)} data points")
        iqmon_obj_dates = [d['date'] for d in iqmon if d['imtype'] == 'OBJECT']
        iqmon_obj_alt = [d['alt']/90 for d in iqmon if d['imtype'] == 'OBJECT']
        iqmon_cal_dates = [d['date'] for d in iqmon if d['imtype'] in ['BIAS', 'DARK']]
        iqmon_cal_alt = [0.5 for d in iqmon if d['imtype'] in ['BIAS', 'DARK']]
        iqmon_flat_dates = [d['date'] for d in iqmon if d['imtype'] in ['FLAT', 'TWIFLAT', 'DOMEFLAT']]
        iqmon_flat_alt = [0.5 for d in iqmon if d['imtype'] in ['FLAT', 'TWIFLAT', 'DOMEFLAT']]

        ## Build Telescope Status plot
        log.info('Build Telescope Status plot')
        dome = [s['open_closed'] for s in domestatus]
        dome_date = [s['date'] for s in domestatus]
        height = webcfg['Weather'].getint('plot_dome_height', 60)
        dome_plot = figure(width=900, height=height, x_axis_type="datetime",
                           y_range=(-0.2,1.2),
                           x_range=(end - timedelta(hours=span_hours), end),
                           )
#         open_date = [dome_date[i] for i,d in enumerate(dome) if d < 0.5]
#         open_dome = [d for i,d in enumerate(dome) if d < 0.5]
#         closed_date = [dome_date[i] for i,d in enumerate(dome) if d >= 0.5]
#         closed_dome = [d for i,d in enumerate(dome) if d >= 0.5]
#         dome_plot.line(closed_date, closed_dome, line_width=4, color="black")
#         dome_plot.line(open_date, open_dome, line_width=4, color="green")

        dome_plot.line(dome_date, dome, line_width=4, color="black", line_alpha=0.8)

        # IQMon Files
        dome_plot.circle(iqmon_obj_dates, iqmon_obj_alt,
                         size=markersize, color="blue",
                         line_alpha=0.8, fill_alpha=0.8)
        dome_plot.circle(iqmon_cal_dates, iqmon_cal_alt,
                         size=markersize, color="black",
                         line_alpha=0.8, fill_alpha=0.8)
        dome_plot.circle(iqmon_flat_dates, iqmon_flat_alt,
                         size=markersize, color="yellow",
                         line_alpha=0.8, fill_alpha=0.8)
        dome_plot.yaxis.axis_label = f'{telcfg["Telescope"].get("name")}'
        dome_plot.xaxis.axis_label = 'Time (UT)'
        dome_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
        dome_plot.yaxis.ticker = [0,1]
        dome_plot.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Render
    log.info(f"Overplotting twilights")
    plot_names = ['plot_temperature', 'plot_humidity', 'plot_cloudiness',
                  'plot_wind_speed', 'plot_rain', 'plot_pressure', 'plot_safe']
    plot_info_list = []
    for name in plot_names:
        if webcfg['Weather'].get(name, None) is not None:
            log.info(f"  {name}")
            limit_string = webcfg['Weather'].get(f'{name}_limits', None)
            if limit_string is not None:
                ymin,ymax = limit_string.split(',')
            else:
                ymin,ymax = 0,1
            yrange = float(ymax)-float(ymin)
            plot_info_list.append([name, eval(name), float(ymax)+yrange, float(ymin)-yrange])

    plot_column_list = []
    plot_twilights_list = []
    for i,plot_info in enumerate(plot_info_list):
        if webcfg['Weather'].get(plot_info[0], None) is not None:
            if i != 0:
                plot_info[1].x_range = plot_info_list[0][1].x_range
            plot_column_list.append(plot_info[1])
            plot_twilights_list.append(plot_info)
    overplot_twilights(plot_twilights_list, query_end, webcfg, ndays=query_days, log=log)

    if webcfg['Weather'].getboolean('plot_dome', False) is True and telcfg is not None:
        dome_plot.x_range = plot_info_list[0][1].x_range
        plot_column_list.append(dome_plot)

    # Add time log
    plot_column_list[-1].plot_height += 40
    plot_column_list[-1].xaxis.visible = True
    plot_column_list[-1].xaxis.formatter = DatetimeTickFormatter(hourmin=['%H:%M'])
    plot_column_list[-1].xaxis.ticker = DatetimeTicker(desired_num_ticks=24)
    plot_column_list[-1].xaxis.axis_label = 'UT Time'

    log.info(f"Rendering bokeh plot for {plot_column_list}")
    script, div = components(column(plot_column_list))

    return script, div
