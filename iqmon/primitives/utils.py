from pathlib import Path
from datetime import datetime, timedelta
import logging

import numpy as np
import ephem
from astropy.table import Table, Column
from astropy import stats
from astroquery.vizier import Vizier


##-----------------------------------------------------------------------------
## Function: get_destination_dir
##-----------------------------------------------------------------------------
def get_destination_dir(cfg):
    raw_string = cfg['FileHandling'].get('destination_dir')
    nowut = datetime.utcnow()
    y = f'{nowut.year:4d}'
    m = f'{nowut.month:02d}'
    d = f'{nowut.day:02d}'
    result = raw_string.replace('YYYY', y).replace('MM', m).replace('DD', d)
    return result


##-----------------------------------------------------------------------------
## Function: get_destination_dir
##-----------------------------------------------------------------------------
def get_jpeg_dir(cfg, obstime):
    raw_string = cfg['jpeg'].get('directory')
    y = f'{obstime.year:4d}'
    m = f'{obstime.month:02d}'
    d = f'{obstime.day:02d}'
    result = raw_string.replace('YYYY', y).replace('MM', m).replace('DD', d)
    return result


##-----------------------------------------------------------------------------
## Function: get_sunrise_sunset
##-----------------------------------------------------------------------------
def get_sunrise_sunset(start):
    obs = ephem.Observer()
    obs.lon = "-155:34:33.9"
    obs.lat = "+19:32:09.66"
    obs.elevation = 3400.0
    obs.temp = 10.0
    obs.pressure = 680.0
    obs.date = start.strftime('%Y/%m/%d 10:00:00')

    obs.horizon = '0.0'
    result = {'sunset': obs.previous_setting(ephem.Sun()).datetime(),
              'sunrise': obs.next_rising(ephem.Sun()).datetime(),
             }
    obs.horizon = '-6.0'
    result['evening_civil_twilight'] = obs.previous_setting(ephem.Sun(),
                                           use_center=True).datetime()
    result['morning_civil_twilight'] = obs.next_rising(ephem.Sun(),
                                           use_center=True).datetime()
    obs.horizon = '-12.0'
    result['evening_nautical_twilight'] = obs.previous_setting(ephem.Sun(),
                                              use_center=True).datetime()
    result['morning_nautical_twilight'] = obs.next_rising(ephem.Sun(),
                                              use_center=True).datetime()
    obs.horizon = '-18.0'
    result['evening_astronomical_twilight'] = obs.previous_setting(ephem.Sun(),
                                                  use_center=True).datetime()
    result['morning_astronomical_twilight'] = obs.next_rising(ephem.Sun(),
                                                  use_center=True).datetime()
    return result


##-----------------------------------------------------------------------------
## Function: Evaluate pre and post conditions
##-----------------------------------------------------------------------------
def pre_condition(primitive, name, condition,
                  fail_level=logging.DEBUG,
                  success_level=logging.DEBUG):
    if condition is True:
        primitive.log.log(success_level,
            f'Precondition for {primitive.__class__.__name__} "{name}" satisfied')
    else:
        primitive.log.log(fail_level,
            f'Precondition for {primitive.__class__.__name__} "{name}" failed')
    return condition


def post_condition(primitive, name, condition,
                   fail_level=logging.WARNING,
                   success_level=logging.DEBUG):
    if condition is True:
        primitive.log.log(success_level,
            f'Postcondition for {primitive.__class__.__name__} "{name}" satisfied')
    else:
        primitive.log.log(fail_level,
            f'Postcondition for {primitive.__class__.__name__} "{name}" failed')
    return condition


##-----------------------------------------------------------------------------
## Function: download_vizier
##-----------------------------------------------------------------------------
def download_vizier(pointing, radius, catalog='UCAC4', band='i', maglimit=None):
    catalogs = {'UCAC4': 'I/322A', 'Gaia': 'I/345/gaia2'}
    if catalog not in catalogs.keys():
        print(f'{catalog} not in {catalogs.keys()}')
        raise NotImplementedError
    if band not in ['r', 'i']:
        print(f'Band {band} not supported')
        raise NotImplementedError

    columns = {'UCAC4': ['_RAJ2000', '_DEJ2000', 'rmag', 'imag'],
               'Gaia': ['RA_ICRS', 'DE_ICRS', 'Gmag', 'RPmag']}
    ra_colname = {'UCAC4': '_RAJ2000',
                  'Gaia': 'RA_ICRS'}
    dec_colname = {'UCAC4': '_DEJ2000',
                   'Gaia': 'DE_ICRS'}
    mag_colname = {'UCAC4': f'{band}mag',
                   'Gaia': 'RPmag'}
    filter_string = '>0' if maglimit is None else f"<{maglimit}"
    column_filter = {mag_colname[catalog]: filter_string}

    v = Vizier(columns=columns[catalog],
               column_filters=column_filter)
    v.ROW_LIMIT = 2e4

    try:
        stars = Table(v.query_region(pointing, catalog=catalogs[catalog],
                                     radius=c.Angle(radius, "deg"))[0])
        stars.add_column( Column(data=stars[ra_colname[catalog]], name='RA') )
        stars.add_column( Column(data=stars[dec_colname[catalog]], name='DEC') )
        stars.add_column( Column(data=stars[mag_colname[catalog]], name='mag') )
    except:
        stars = None
    return stars


##-----------------------------------------------------------------------------
## Function: get_panstarrs
##-----------------------------------------------------------------------------
def get_panstarrs(cfg, field_name, pointing, filter, radius=0.40,
                  maglimit=None, log=None):
    catalogname = cfg['Photometry'].get('calibration_catalog')
    band = {'PSi': 'i', 'PSr': 'r'}[filter]
    if maglimit is None: maglimit = 25

    ## First check if we have a pre-downloaded catalog for this field
    local_catalog_path = Path(cfg['Photometry'].get('local_catalog_path', '.'))
    local_catalog_file = local_catalog_path.joinpath(f'{field_name}_{band}{maglimit*10:03.0f}.cat')
    if local_catalog_file.exists() is True:
        ## Read local file
        if log: log.debug(f'  Reading {local_catalog_file}')
        pscat = Table.read(local_catalog_file, format='ascii.csv')
    else:
        ## Download
        if log: log.debug(f'  Downloading from Mast')
        from astroquery.mast import Catalogs
#         cols = ['objName', 'objID', 'objInfoFlag', 'qualityFlag', 'raMean',
#                 'decMean', 'raMeanErr', 'decMeanErr', 'epochMean', 'nDetections',
#                 'ng', 'nr', 'ni', 'gMeanApMag', 'gMeanApMagErr', 'gMeanApMagStd',
#                 'gMeanApMagNpt', 'gFlags', 'rMeanApMag', 'rMeanApMagErr',
#                 'rMeanApMagStd', 'rMeanApMagNpt', 'rFlags', 'iMeanApMag',
#                 'iMeanApMagErr', 'iMeanApMagStd', 'iMeanApMagNpt', 'iFlags']
        cols = ['objName', 'objID', 'raMean', 'decMean', 'raMeanErr', 'decMeanErr',
                'gMeanApMag', 'gMeanApMagErr', 'gMeanApMagStd',
                'gMeanApMagNpt', 'gFlags', 'rMeanApMag', 'rMeanApMagErr',
                'rMeanApMagStd', 'rMeanApMagNpt', 'rFlags', 'iMeanApMag',
                'iMeanApMagErr', 'iMeanApMagStd', 'iMeanApMagNpt', 'iFlags']

        if band in ['i', 'r', 'g']:
            pscat = Catalogs.query_region(pointing, radius=radius,
                             catalog="Panstarrs", table="mean", data_release="dr2",
                             sort_by=[("desc", f"{band}MeanApMag")], columns=cols,
                             iMeanApMag=[("gte", 0), ("lte", maglimit)],
                             )
        else:
            pscat = Catalogs.query_region(pointing, radius=radius,
                             catalog="Panstarrs", table="mean", data_release="dr2",
                             columns=cols,
                             )

#         if band == 'i':
#             pscat = Catalogs.query_region(pointing, radius=radius,
#                              catalog="Panstarrs", table="mean", data_release="dr2",
#                              sort_by=[("desc", f"{band}MeanApMag")], columns=cols,
#                              iMeanApMag=[("gte", 0), ("lte", maglimit)],
#                              )
#         elif band == 'r':
#             pscat = Catalogs.query_region(pointing, radius=radius,
#                              catalog="Panstarrs", table="mean", data_release="dr2",
#                              sort_by=[("desc", f"{band}MeanApMag")], columns=cols,
#                              rMeanApMag=[("gte", 0), ("lte", maglimit)],
#                              )
#         elif band == 'g':
#             pscat = Catalogs.query_region(pointing, radius=radius,
#                              catalog="Panstarrs", table="mean", data_release="dr2",
#                              sort_by=[("desc", f"{band}MeanApMag")], columns=cols,
#                              gMeanApMag=[("gte", 0), ("lte", maglimit)],
#                              )
#         else:
#             pscat = Catalogs.query_region(pointing, radius=radius,
#                              catalog="Panstarrs", table="mean", data_release="dr2",
#                              columns=cols,
#                              )
        if log: log.debug(f'  Got {len(pscat)} entries total')
        if log: log.debug(f'  Got {len(pscat)} entries with {band}-band magnitudes')
        if log: log.debug(f'  Writing {local_catalog_file}')
        pscat.write(local_catalog_file, format='ascii.csv')

    # Filter based on magnitude
    if maglimit is not None:
        pscat = pscat[pscat[f'{band}MeanApMag'] <= maglimit]

    return pscat


##-----------------------------------------------------------------------------
## Function: sigma_clipping_line_fit
##-----------------------------------------------------------------------------
def sigma_clipping_line_fit(xdata, ydata, nsigma=3, maxiter=7, maxcleanfrac=0.3,
                            intercept_fixed=False, intercept0=0, slope0=1,
                            log=None):
        if log: log.debug('  Running sigma_clipping_line_fit')
        npoints = len(xdata)
        if log: log.debug(f'  npoints = {npoints}')
        fit = fitting.LinearLSQFitter()
        line_init = models.Linear1D(slope=slope0, intercept=intercept0)
        line_init.intercept.fixed = intercept_fixed
        fitted_line = fit(line_init, xdata, ydata)
        deltas = ydata - fitted_line(xdata)
        mean, median, std = stats.sigma_clipped_stats(deltas)
        cleaned = np.array(abs(deltas) < nsigma*std)
#         if log: log.debug(cleaned)
        if log: log.debug(f'  fitted slope = {fitted_line.slope.value:3g}')
        if log: log.debug(f'  std = {std:4g}')
        if log: log.debug(f'  n_cleaned = {np.sum(cleaned)}')
        for iteration in range(1, maxiter+1):
            last_std = std
            new_fit = fit(line_init, xdata[cleaned], ydata[cleaned])
            deltas = ydata - new_fit(xdata)
            mean, median, std = stats.sigma_clipped_stats(deltas)
            cleaned = cleaned | np.array(abs(deltas) < nsigma*std)
            if np.sum(~cleaned)/npoints > maxcleanfrac:
                if log: log.debug(f'  Exceeded maxcleanfrac of {maxcleanfrac}')
                return fitted_line
            if std > last_std:
                if log: log.debug(f'  StdDev increased')
                return fitted_line
            else:
                fitted_line = new_fit
#             if log: log.debug(cleaned)
            if log: log.debug(f'  {iteration} fitted slope = {fitted_line.slope.value:3g}')
            if log: log.debug(f'  {iteration} std = {std:4g}')
            if log: log.debug(f'  {iteration} n_cleaned = {np.sum(cleaned)}')

        return fitted_line


##-----------------------------------------------------------------------------
## Function: estimate_f0
##-----------------------------------------------------------------------------
def estimate_f0(A, band='i'):
    '''
    1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1
    https://archive.is/20121204144725/http://www.astro.utoronto.ca/~patton/astro/mags.html#selection-587.2-587.19
    band cent    dl/l    Flux0   Reference
    U    0.36    0.15    1810    Bessel (1979)
    B    0.44    0.22    4260    Bessel (1979)
    V    0.55    0.16    3640    Bessel (1979)
    R    0.64    0.23    3080    Bessel (1979)
    I    0.79    0.19    2550    Bessel (1979)
    J    1.26    0.16    1600    Campins, Reike, & Lebovsky (1985)
    H    1.60    0.23    1080    Campins, Reike, & Lebovsky (1985)
    K    2.22    0.23    670     Campins, Reike, & Lebovsky (1985)
    g    0.52    0.14    3730    Schneider, Gunn, & Hoessel (1983)
    r    0.67    0.14    4490    Schneider, Gunn, & Hoessel (1983)
    i    0.79    0.16    4760    Schneider, Gunn, & Hoessel (1983)
    z    0.91    0.13    4810    Schneider, Gunn, & Hoessel (1983)
    '''
    tabledata = {'band': ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'g', 'r', 'i', 'z'],
                 'cent': [0.36, 0.44, 0.55, 0.64, 0.79, 1.26, 1.60, 2.22, 0.52, 0.67, 0.79, 0.91],
                 'dl/l': [0.15, 0.22, 0.16, 0.23, 0.19, 0.16, 0.23, 0.23, 0.14, 0.14, 0.16, 0.13],
                 'Flux0': [1810, 4260, 3640, 3080, 2550, 1600, 1080, 670 , 3730, 4490, 4760, 4810],
                }
    t = Table(tabledata)
    band = t[t['band'] == band]
    dl = 0.16 # dl/l (for i band)
    f0 = band['Flux0'] * 1.51e7 * A * band['dl/l'] # photons / sec
    return f0[0]


##-----------------------------------------------------------------------------
## Function: mode
##-----------------------------------------------------------------------------
def mode(data):
    '''
    Return mode of image.  Assumes int values (ADU), so uses binsize of one.
    '''
    bmin = np.floor(min(data.ravel())) - 1./2.
    bmax = np.ceil(max(data.ravel())) + 1./2.
    bins = np.arange(bmin,bmax,1)
    hist, bins = np.histogram(data.ravel(), bins=bins)
    centers = (bins[:-1] + bins[1:]) / 2
    w = np.argmax(hist)
    mode = int(centers[w])
    return mode


##-----------------------------------------------------------------------------
## Function: find_master
##-----------------------------------------------------------------------------
def build_master_file_name(meta, master_type, date_string):
    if master_type in ['BIAS']:
        master_file_name = f"MasterBIAS_{date_string}.fits"
    elif master_type in ['DARK']:
        exptime = int(meta.get('exptime'))
        master_file_name = f"MasterDARK_{exptime:03d}s_{date_string}.fits"
    elif master_type in ['FLAT']:
        master_file_name = f"MasterFLAT_{meta.get('filter')}_{date_string}.fits"
    else:
        master_file_name = None
    return master_file_name


def find_master(master_directory, master_type, meta):
    # Find master bias file
    if master_directory is not None:
        master_directory = Path(master_directory)
    else:
        return None
    if master_directory.exists() is False:
        return None

    # Build expected file name
    date_string = meta.get('date').strftime('%Y%m%dUT')
    master_file_name = build_master_file_name(meta, master_type, date_string)
    master_file = master_directory.joinpath(master_file_name)

    # Go hunting for the files
    if master_file.exists() is True:
        return master_file
    else:
        # Look for bias within 10 days
        count = 0
        while master_file.exists() is False and count <= 10:
            count += 1
            # Days before
            date_string = (meta.get('date')-timedelta(count)).strftime('%Y%m%dUT')
            master_file_name = build_master_file_name(meta, master_type, date_string)
            master_file = master_directory.joinpath(master_file_name)
            if master_file.exists() is True:
                return master_file
            # Days after
            date_string = (meta.get('date')+timedelta(count)).strftime('%Y%m%dUT')
            master_file_name = build_master_file_name(meta, master_type, date_string)
            master_file = master_directory.joinpath(master_file_name)
            if master_file.exists() is True:
                return master_file
        if master_file.exists() is False:
            return None
        return master_file
