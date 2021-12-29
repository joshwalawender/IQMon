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

from iqmon.webpage import mongo_query

log = logging.getLogger('FlaskLogger')


##-------------------------------------------------------------------------
## Function: generate_iqmon_plot
##-------------------------------------------------------------------------
def generate_iqmon_plot(telescope, date=None, plot_ndays=1, span_hours=24):
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
    ## IQMon Query
    log.info(f"Querying IQMon results database")
    query_dict = {'telescope': telescope,
                  'date': {'$gt': start, '$lt': end}}
    query_result = mongo_query('iqmon', query_dict)
    iqmon = [d for d in query_result]
    log.info(f"  Got {len(iqmon)} data points")

    dates_fwhm = [d['date'] for d in iqmon if 'fwhm' in d.keys()]
    fwhm = [d['fwhm'] for d in iqmon if 'fwhm' in d.keys()]
    dates_elip = [d['date'] for d in iqmon if 'ellipticity' in d.keys()]
    elip = [d['ellipticity'] for d in iqmon if 'ellipticity' in d.keys()]

    markersize = 2

    ##-------------------------------------------------------------------------
    ## FWHM Plot
    log.info('Build FWHM plot')
    fwhm_plot = figure(width=900, height=200, x_axis_type="datetime",
                       y_range=(0,10),
                       x_range=(end - timedelta(hours=span_hours), end),
                       )
    fwhm_plot.circle(dates_fwhm, fwhm,
                     size=markersize, color="blue", alpha=0.8)
    fwhm_plot.yaxis.axis_label = 'FWHM (pix)'
    fwhm_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
    fwhm_plot.xaxis.visible = False

    ##-------------------------------------------------------------------------
    ## Ellipticity Plot
    log.info('Build ellipticity plot')
    ellipticity_plot = figure(width=900, height=200, x_axis_type="datetime",
                       y_range=(1,2),
                       x_range=fwhm_plot.x_range,
                       )
    ellipticity_plot.circle(dates_elip, elip,
                     size=markersize, color="blue", alpha=0.8)
    ellipticity_plot.yaxis.axis_label = 'ellipticity'
    ellipticity_plot.yaxis.formatter = NumeralTickFormatter(format="0.0")
    ellipticity_plot.xaxis.visible = True

    ellipticity_plot.xaxis.formatter = DatetimeTickFormatter(hourmin=['%H:%M'])
    ellipticity_plot.xaxis.ticker = DatetimeTicker(desired_num_ticks=24)
    ellipticity_plot.xaxis.axis_label = 'UT Time'

    ##-------------------------------------------------------------------------
    ## Render
    log.info(f"Rendering bokeh plot")
    script, div = components(column(fwhm_plot,
                                    ellipticity_plot,
                                    ))

    return script, div
