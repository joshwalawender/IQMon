#!/usr/bin/env python
# encoding: utf-8
"""
IQMon.py

Created by Josh Walawender on 2013-07-22.
Copyright (c) 2013 . All rights reserved.
"""
from __future__ import division, print_function

## Import General Tools
import sys
import os
import re
import stat
import shutil
import datetime
import subprocess
import logging
import yaml
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pyplot


## Import Astronomy Specific Tools
import ephem
import astropy.units as u
import astropy.io.fits as fits
import astropy.coordinates as coords
import astropy.table as table
import astropy.wcs as wcs
import astropy.io.ascii as ascii


__version__ = '1.3.3'


##-----------------------------------------------------------------------------
## Mode Function
##-----------------------------------------------------------------------------
def mode(data, binsize):
    bmin = math.floor(min(data)/binsize)*binsize - binsize/2.
    bmax = math.ceil(max(data)/binsize)*binsize + binsize/2.
    bins = np.arange(bmin,bmax,binsize)
    hist, bins = np.histogram(data, bins=bins)
    centers = (bins[:-1] + bins[1:]) / 2
    foo = zip(hist, centers)
    return max(foo)[1]



##-----------------------------------------------------------------------------
## Define Telescope object to hold telescope information
##-----------------------------------------------------------------------------
class Telescope(object):
    '''
    Contains information about the telescope that can be passed to methods and
    functions.  The concept for operation is that the user will write a simple
    script which creates a telescope object and assigned values to all its
    properties (or sets them to None).
    '''
    def __init__(self, config_file):
        ## System Properties
        self.temp_file_path = None
        self.plot_file_path = None
        self.logs_file_path = None
        ## Telescope Properties
        self.name = None
        self.focal_length = None
        self.pixel_size = None
        self.aperture = None
        self.gain = None
        self.saturation = None
        self.site = ephem.Observer()
        ## Preferences
        self.threshold_FWHM = None
        self.threshold_pointing_err = None
        self.threshold_ellipticity = None
        self.threshold_zeropoint = None
        self.SCAMP_aheader = None
        self.units_for_FWHM = None
        self.ROI = None
        self.PSF_measurement_radius = None
        self.pointing_marker_size = None
        self.SExtractor_params = None
        self.SCAMP_params = None
        self.catalog_info = None
        ## Derived Properties
        self.nXPix = None
        self.nYPix = None
        self.pixel_scale = None
        self.f_ratio = None

        ## Read YAML Config File
        if not os.path.exists(config_file):
            print('WARNING: Configuration file {} not found'.format(config_file))
            raise IOError
        with open(config_file, 'r') as yaml_string:
            config = yaml.load(yaml_string)
        if not isinstance(config, dict):
            print('Configuration file {} not parsed as dict'.format(config_file))
            sys.exit(1)
        self.config = config
        ## Populate Configured Properties
        if 'name' in config.keys(): self.name = str(config['name'])
        if 'temp_file_path' in config.keys():
            self.temp_file_path = os.path.expanduser(config['temp_file_path'])
        if 'plot_file_path' in config.keys():
            self.plot_file_path = os.path.expanduser(config['plot_file_path'])
        if 'logs_file_path' in config.keys():
            self.logs_file_path = os.path.expanduser(config['logs_file_path'])
        if 'focal_length' in config.keys(): self.focal_length = config['focal_length'] * u.mm
        if 'pixel_size' in config.keys(): self.pixel_size = config['pixel_size'] * u.um
        if 'aperture' in config.keys(): self.aperture = config['aperture'] * u.mm
        if 'gain' in config.keys(): self.gain = config['gain'] / u.adu
        if 'saturation' in config.keys(): self.saturation = config['saturation'] * u.adu
        if 'threshold_FWHM' in config.keys(): self.threshold_FWHM = config['threshold_FWHM'] * u.pix
        if 'threshold_pointing_err' in config.keys(): self.threshold_pointing_err = config['threshold_pointing_err'] * u.arcmin
        if 'threshold_ellipticity' in config.keys(): self.threshold_ellipticity = config['threshold_ellipticity']
        if 'threshold_zeropoint' in config.keys(): self.threshold_zeropoint = config['threshold_zeropoint']
        if 'SCAMP_aheader' in config.keys():
            self.SCAMP_aheader = os.path.expanduser(config['SCAMP_aheader'])
        if 'units_for_FWHM' in config.keys(): self.units_for_FWHM = getattr(u, config['units_for_FWHM'])
        if 'ROI' in config.keys(): self.ROI = str(config['ROI'])
        if 'PSF_measurement_radius' in config.keys(): self.PSF_measurement_radius = config['PSF_measurement_radius'] * u.pix
        if 'pointing_marker_size' in config.keys(): self.pointing_marker_size = config['pointing_marker_size'] * u.arcmin
        if 'SExtractor_params' in config.keys(): self.SExtractor_params = config['SExtractor_params']
        if 'SCAMP_params' in config.keys(): self.SCAMP_params = config['SCAMP_params']
        if 'catalog' in config.keys(): self.catalog_info = config['catalog']

        ## Determine Pixel Scale and F-Ratio
        assert self.pixel_size
        assert self.focal_length
        self.pixel_scale = self.pixel_size.to(u.mm)/self.focal_length.to(u.mm)*u.radian.to(u.arcsec)*u.arcsec/u.pix
        self.f_ratio = self.focal_length.to(u.mm)/self.aperture.to(u.mm)

        ## Define astropy.units Equivalency for Arcseconds and Pixels
        self.pixel_scale_equivalency = [(u.pix, u.arcsec,
             lambda pix: (pix*u.radian.to(u.arcsec) * self.pixel_size / self.focal_length).decompose().value,
             lambda arcsec: (arcsec/u.radian.to(u.arcsec) * self.focal_length / self.pixel_size).decompose().value
             )]

        ## create paths
        assert self.temp_file_path
        assert self.plot_file_path
        assert self.logs_file_path
        paths_to_check = [self.temp_file_path, self.plot_file_path, self.logs_file_path]
        paths_to_create = []
        for path in paths_to_check:
            while not os.path.exists(path):
                paths_to_create.append(path)
                path = os.path.split(path)[0]
        while len(paths_to_create) > 0:
            os.mkdir(paths_to_create.pop())

        ## Assert required properties
        assert self.name



    def __del__(self):
        pass
#         print('Deleted telescope object')

    def __enter__(self):
        return self

    def __exit__(self ,type, value, traceback):
        self.__del__()


##-----------------------------------------------------------------------------
## Define Image object which holds information and methods for analysis
##-----------------------------------------------------------------------------
class Image(object):
    '''
    The Image object represents a single input image to the IQMon process.

    When defined, the image objects requires a filename to a valid fits file.
    '''
    def __init__(self, input, tel):
        self.start_process_time = datetime.datetime.now()
        if os.path.exists(input):
            fits_file_directory, fits_filename = os.path.split(input)
            self.raw_file = input
            self.raw_file_name = fits_filename
            self.raw_file_directory = fits_file_directory
            self.raw_file_basename, self.file_ext = os.path.splitext(fits_filename)
        else:
            self.raw_file = None
            self.raw_file_name = None
            self.raw_file_directory = None
            raise IOError("File {0} does not exist".format(input))
        ## Confirm that input tel is an IQMon.Telescope object
        assert isinstance(tel, Telescope)
        self.tel = tel
        ## Initialize values to None
        self.logger = None
        self.working_file = None
        self.header = None
        self.exptime = None
        self.catalog_filter = None
        self.object_name = None
        self.image_WCS = None
        self.astrometry_solved = None
        self.coordinate_of_center_pixel = None
        self.coordinate_from_header = None
        self.n_stars_SExtracted = None
        self.SExtractor_background = None
        self.SExtractor_background_RMS = None
        self.temp_files = []
        self.SExtractor_catalogfile = None
        self.SExtractor_results = None
        self.position_angle = None
        self.zero_point = None
        self.zero_point_mode = None
        self.zero_point_median = None
        self.zero_point_average = None
        self.zero_point_correlation = None
        self.zero_point_plotfile = None
        self.total_process_time = None
        self.FWHM = None
        self.FWHM_median = None
        self.FWHM_mode = None
        self.FWHM_average = None
        self.ellipticity = None
        self.ellipticity_median = None
        self.ellipticity_mode = None
        self.ellipticity_average = None
        self.PSF_plot_file = None
        self.pointing_error = None
        self.image_flipped = None
        self.jpeg_file_names = []
        self.cropped = False
        self.crop_x1 = None
        self.crop_x2 = None
        self.crop_y1 = None
        self.crop_y2 = None
        self.original_nXPix = None
        self.original_nYPix = None
        self.SCAMP_catalog = None
        self.SCAMP_successful = False
        self.catalog_name = None
        self.catalog_data = None
        self.flags = {
                      'FWHM': False,\
                      'ellipticity': False,\
                      'pointing error': False,\
                      'zero point': False,\
                      'other': False,\
                     }

    def __del__(self):
#         print('Deleting object referring to {}'.format(self.raw_file_name))
        if self.logger:
            self.logger = None
        if self.temp_files:
            for item in self.temp_files:
                if os.path.exists(item):
#                     print("  Deleting {0}".format(item))
                    os.remove(item)

    def __enter__(self):
        return self

    def __exit__(self ,type, value, traceback):
        self.__del__()


    ##-------------------------------------------------------------------------
    ## Make Logger Object
    ##-------------------------------------------------------------------------
    def get_logger(self, logger):
        '''
        If calling from another program which has its own logger object, pass
        that logger to IQMon with this method.
        '''
        self.logger = logger

        ## Print Configuration to Log
        self.logger.debug('Using configuration:')
        for entry in self.tel.config.keys():
            self.logger.debug('  {} = {}'.format(entry, self.tel.config[entry]))


    def make_logger(self, logfile=None, clobber=False, verbose=False):
        '''
        Create the logger object to use when processing.  Takes as input the
        full path to the file to write the log to and verboase, a boolean value
        which will increase the verbosity of the concole log (the file log will
        always be at debug level).
        '''
        if not logfile:
            logfile = os.path.join(self.tel.logs_file_path, '{}_IQMon.log'.format(self.raw_file_basename))
        self.logfile = logfile
        self.logfilename = os.path.split(self.logfile)[1]
        if clobber:
            if os.path.exists(logfile): os.remove(logfile)
        self.logger = logging.getLogger(self.raw_file_basename)
        if len(self.logger.handlers) < 1:
            self.logger.setLevel(logging.DEBUG)
            LogFileHandler = logging.FileHandler(logfile)
            LogFileHandler.setLevel(logging.DEBUG)
            LogConsoleHandler = logging.StreamHandler()
            if verbose:
                LogConsoleHandler.setLevel(logging.DEBUG)
            else:
                LogConsoleHandler.setLevel(logging.INFO)
            LogFormat = logging.Formatter('%(asctime)23s %(levelname)8s: %(message)s')
            LogFileHandler.setFormatter(LogFormat)
            LogConsoleHandler.setFormatter(LogFormat)
            self.logger.addHandler(LogConsoleHandler)
            self.logger.addHandler(LogFileHandler)

        ## Put initial lines in log
        self.logger.info('')
        self.logger.info("###### Processing Image {} ######".format(self.raw_file_name))
        self.logger.info('IQMon version = {}'.format(__version__))
        self.logger.info('')

        ## Print Configuration to Log
        if 'name' in self.tel.config.keys():
            self.logger.debug('Using configuration for telescope: {}'.format(self.tel.config['name']))
        else:
            self.logger.debug('Using configuration:')
        for entry in self.tel.config.keys():
            self.logger.debug('  {} = {}'.format(entry, self.tel.config[entry]))


    ##-------------------------------------------------------------------------
    ## Read Header
    ##-------------------------------------------------------------------------
    def read_header(self):
        '''
        Read information from the image fits header.
        '''
        start_time = datetime.datetime.now()
        self.logger.info("Reading image header.")
        if not self.working_file:
            with fits.open(self.raw_file, ignore_missing_end=True) as hdulist:
                self.header = hdulist[0].header
                self.nYPix, self.nXPix = hdulist[0].data.shape
                self.logger.debug('  Image size is: {},{}'.format(\
                                                           self.nXPix, self.nYPix))
        else:
            with fits.open(self.working_file, ignore_missing_end=True) as hdulist:
                self.header = hdulist[0].header
                self.nYPix, self.nXPix = hdulist[0].data.shape
                self.logger.debug('  Image size is: {},{}'.format(\
                                                           self.nXPix, self.nYPix))



        ## Get exposure time from header (assumes seconds)
        try:
            self.exptime = float(self.header['EXPTIME']) * u.s
        except:
            self.exptime = None
            self.logger.debug("  No exposure time value found in header")
        else:
            self.logger.debug("  Exposure time = {0:.1f} s".format(\
                                                   self.exptime.to(u.s).value))
        ## Get filter from header
        try:
            self.filter = str(self.header['FILTER'])
        except:
            self.filter = None
            self.logger.debug("  No filter value found in header")
        else:
            self.logger.debug("  filter = {}".format(self.filter))
        ## Get object name from header
        try:
            self.object_name = self.header["OBJECT"]
        except:
            self.object_name = None
            self.logger.debug("  No object value found in header")
        else:
            self.logger.debug("  Header object name = {0}".format(
                                                             self.object_name))
        ## Get Observation Date and Time from header
        ## (assumes YYYY-MM-DDTHH:MM:SS format)
        try:
            self.observation_date = self.header["DATE-OBS"]
        except:
            self.observation_date = None
            self.logger.debug("  No date value found in header")
        else:
            self.logger.debug("  Header date = {0}".format(self.observation_date))
        ## Get Site Latitude from header (assumes decimal degrees)
        try:
            self.latitude = self.header["LAT-OBS"] * u.deg
        except:
            self.latitude = None
            self.logger.debug("  No latitude value found in header")
        else:
            self.logger.debug("  Header latitude = {0:.4f} deg".format(\
                                                 self.latitude.to(u.deg).value))
        ## Get Site Longitude from header (assumes decimal degrees)
        try:
            self.longitude = self.header["LONG-OBS"] * u.deg
        except:
            self.longitude = None
            self.logger.debug("  No longitiude value found in header")
        else:
            self.logger.debug("  Header longitiude = {0:.4f} deg".format(\
                                               self.longitude.to(u.deg).value))
        ## Get Site Altitude from header (assumes meters)
        try:
            self.altitude = self.header["ALT-OBS"] * u.meter
        except:
            self.altitude = None
            self.logger.debug("  No altitude value found in header")
        else:
            self.logger.debug("  Header altitude = {0:.0f} meters".format(\
                                              self.altitude.to(u.meter).value))

        ## Read Header Coordinates in to astropy coordinates object
        self.coordinate_from_header = None
        if ('RA' in self.header.keys()) and ('DEC' in self.header.keys()):
            ## Header RA is : separated
            coord_string = '{} {}'.format(self.header['RA'], self.header['DEC'])
            self.logger.debug('  Parsing: "{}" as hours and degrees'.format(coord_string))
            try:
                self.coordinate_from_header = coords.SkyCoord(coord_string,\
                                                unit=(u.hour, u.degree),\
                                                frame='icrs')
            except:
                self.logger.debug('  Parsing: "{}" as hours and degrees failed'.format(coord_string))
            
            if not self.coordinate_from_header:
                self.logger.info('  Could not parse coordinate strings from header')
                self.logger.info('  RA = {}'.format(self.header['RA']))
                self.logger.info('  DEC = {}'.format(self.header['DEC']))

        ## Read WCS
        try:
            self.image_WCS = wcs.WCS(self.header)
        except:
            self.image_WCS = None
            self.logger.info("  No WCS found in image header")
        else:
            self.logger.debug("  Found WCS in image header.")
            for item in self.image_WCS.to_header().cards:
                self.logger.debug('    {}'.format(item))

        ## Determine PA of Image
        if self.image_WCS:
            self.orientation_from_wcs()
            if self.position_angle:
                self.logger.debug("  Position angle of WCS is {0:.1f} deg".format(\
                                              self.position_angle.to(u.deg).value))
                if self.image_flipped:
                    self.logger.debug("  Image is mirrored.")

        ## Determine Alt, Az, Moon Sep, Moon Illum using ephem module
        if self.observation_date and self.latitude and self.longitude and self.coordinate_from_header:
            ## Populate site object properties
            SiteDate = "/".join(self.observation_date[0:10].split("-"))
            SiteTime = self.observation_date[11:]        
            self.tel.site.date = ephem.Date(SiteDate+" "+SiteTime)
            self.tel.site.lat = str(self.latitude.to(u.deg).value)
            self.tel.site.lon = str(self.longitude.to(u.deg).value)
            if self.altitude:
                self.tel.site.elevation = self.altitude.to(u.meter).value
            ## Do calculations using ephem
            RADEC_string = ','.join(self.coordinate_from_header.to_string('hmsdms', sep=':').split())
            TargetObject = ephem.readdb("Target,f|M|F7,{},2.02,2000".format(RADEC_string))
            TargetObject.compute(self.tel.site)
            self.target_alt = TargetObject.alt * 180./ephem.pi * u.deg
            self.target_az = TargetObject.az * 180./ephem.pi * u.deg
            self.logger.debug("  Target Alt, Az = {0:.1f}, {1:.1f}".format(\
                                             self.target_alt.to(u.deg).value,\
                                             self.target_az.to(u.deg).value))
            self.target_zenith_angle = 90.*u.deg - self.target_alt
            self.airmass = 1.0/math.cos(self.target_zenith_angle.to(u.radian).value)*(1.0 - 0.0012*(1.0/(math.cos(self.target_zenith_angle.to(u.radian).value)**2 - 1.0)))
            self.logger.debug("  Target airmass (calculated) = {0:.2f}".format(\
                                                                 self.airmass))
            ## Calculate Moon Position and Illumination
            TheMoon = ephem.Moon()
            TheMoon.compute(self.tel.site)
            self.moon_phase = TheMoon.phase
            self.moon_sep = ephem.separation(TargetObject, TheMoon)
            self.moon_sep = self.moon_sep * 180./ephem.pi * u.deg
            self.moon_alt = TheMoon.alt * 180./ephem.pi * u.deg
            if self.moon_alt > 0:
                self.logger.debug("  A {0:.0f} percent illuminated Moon is {1:.0f} deg from target.".format(\
                                       self.moon_phase,\
                                       self.moon_sep.to(u.deg).value))
            else:
                self.logger.debug("  A {0:.0f} percent illuminated Moon is down.".format(\
                                       self.moon_phase))
        else:
            self.target_alt = None
            self.target_az = None
            self.moon_phase = None
            self.moon_sep = None
            self.moon_alt = None
            self.target_zenith_angle = None
            self.airmass = None
            self.logger.warning("Object and Moon positions not calculated.")

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done reading image header in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Edit Header
    ##-------------------------------------------------------------------------
    def edit_header(self, keyword, value):
        '''
        Edit a single keyword in the image fits header.
        '''
        self.logger.info('Editing image header: {} = {}'.format(keyword, value))
#         hdulist = fits.open(self.working_file,\
#                             ignore_missing_end=True, mode='update')
        with fits.open(self.working_file, ignore_missing_end=True, mode='update') as hdulist:
            hdulist[0].header[keyword] = value
            hdulist.flush()


    ##-------------------------------------------------------------------------
    ## Read Image
    ##-------------------------------------------------------------------------
    def read_image(self):
        '''
        Read the raw image and write out a working image in the IQMon temporary
        directory.
        
        - For the moment, this only copies a fits file from the original
          location to the IQMon tmp directory.
        - Later implement file format conversion from CRW, CR2, DNG, etc to
          fits using dcraw.
          
        Notes on reading raw files:
        * dcraw converts to ppm file, then pamtofits converts to fits (pamtofits 
          sends output to STDOUT, so must redirect to a fits file).
        * pamtofits is part of netpbm, but need to install netpbm-bin on fink
          to get the command line programs
        
        dcraw options:
        * -4 makes linear file
        * -D makes totally raw file (Without -D option, color interpolation is
             done.  Without -D option, get raw pixel values).
        '''
        start_time = datetime.datetime.now()
        chmod_code = stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH | stat.S_IWOTH
        if self.working_file:
            if os.path.exists(self.working_file): os.remove(self.working_file)
        ## fts extension:  make working copy and rename to .fits
        if self.file_ext == '.fts':
            self.logger.info('Making working copy of raw image: {}'.format(\
                                                            self.raw_file_name))
            self.working_file = os.path.join(self.tel.temp_file_path,\
                                             self.raw_file_basename+'.fits')
            shutil.copy2(self.raw_file, self.working_file)
            os.chmod(self.working_file, chmod_code)
            self.temp_files.append(self.working_file)
            self.file_ext = '.fits'
        ## fits extension:  make working copy
        elif self.file_ext == '.fits':
            self.logger.info('Making working copy of raw image: {}'.format(\
                                                            self.raw_file_name))
            self.working_file = os.path.join(self.tel.temp_file_path,\
                                             self.raw_file_name)
            shutil.copy2(self.raw_file, self.working_file)
            os.chmod(self.working_file, chmod_code)
            self.temp_files.append(self.working_file)
            self.file_ext = '.fits'
        ## DSLR file:  convert to fits
        elif self.file_ext in ['.dng', '.DNG', '.cr2', '.CR2']:
            self.logger.info('Converting {} to fits format'.format(\
                                                            self.raw_file_name))
            ## Make copy of raw file
            self.working_file = os.path.join(self.tel.temp_file_path, self.raw_file_name)
            self.logger.debug('Copying {} to {}'.format(self.raw_file, self.working_file))
            shutil.copy2(self.raw_file, self.working_file)
            self.logger.debug('Setting working file permissions for {}'.format(self.working_file))
            os.chmod(self.working_file, chmod_code)
            self.temp_files.append(self.working_file)
            ## Use dcraw to convert to ppm file
            command = ['dcraw', '-t', '2', '-4', self.working_file]
            self.logger.debug('Executing dcraw: {}'.format(repr(command)))
            subprocess.call(command)
            ppm_file = os.path.join(self.tel.temp_file_path, self.raw_file_basename+'.ppm')
            if os.path.exists(ppm_file):
                self.working_file = ppm_file
                self.temp_files.append(self.working_file)
            else:
                self.logger.critical('dcraw failed.  Could not find ppm file.')
            ## Use pamtofits to convert to fits file
            fits_file = os.path.join(self.tel.temp_file_path, self.raw_file_basename+'.fits')
            if os.path.exists(fits_file): os.remove(fits_file)
            conversion_tools = ['pamtofits', 'pnmtofits']
            for conversion_tool in conversion_tools:
                if not os.path.exists(fits_file):
                    command = '{} {} > {}'.format(conversion_tool, self.working_file, fits_file)
                    self.logger.debug('Trying {}: {}'.format(conversion_tool, command))
                    try:
                        subprocess.call(command, shell=True)
                    except:
                        pass
            if os.path.exists(fits_file):
                self.working_file = fits_file
                self.file_ext = self.file_ext = os.path.splitext(self.working_file)[1]
                self.temp_files.append(self.working_file)
            else:
                self.logger.critical('PPM to fits conversion failed.  Could not find fits file.')
            ## Write new fits file with only green image
            self.logger.debug('Only keeping green channel for analysis')
            with fits.open(self.working_file, 'update') as hdulist:
                if len(hdulist) == 1:
                    data = hdulist[0].data
                    green_data = data[1]
                    hdulist[0].data = green_data
                    hdulist.flush()
        else:
            self.logger.warning('Unrecognixed file extension: {}'.format(\
                                                                self.file_ext))
            self.working_file = os.path.join(self.tel.temp_file_path,\
                                             self.raw_file_name)
            sys.exit(1)

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done making working copy of image in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Dark Subtract Image
    ##-------------------------------------------------------------------------
    def dark_subtract(self, Darks):
        '''
        Create master dark and subtract from image.
        
        Input the filename of the appropriate master dark or a list of filenames
        which will be median combined to make the master dark.
        
        '''
        start_time = datetime.datetime.now()
        self.logger.info("Dark subtracting image.")
        self.logger.debug("  Opening image data.")
        hdulist_image = fits.open(self.working_file, mode='update')
        ## Load master dark if provided, but if multiple files input, combine
        ## them in to master dark, then load combined master dark.
        if len(Darks) == 1:
            self.logger.debug("  Found master dark.  Opening master dark data.")
            hdulist_dark = fits.open(Darks[0])
            MasterDarkData = hdulist_dark[0].data
            hdulist_dark.close()
        elif len(Darks) > 1:
            self.logger.info("  Median combining {0} darks.".format(len(Darks)))
            ## Combine multiple darks frames
            DarkData = []
            for Dark in Darks:
                hdulist = fits.open(Dark)
                DarkData.append(hdulist[0].data)
            DarkData = np.array(DarkData)
            MasterDarkData = np.median(DarkData, axis=0)
            ## Save Master Dark to Fits File
            DataPath = os.path.split(self.raw_file)[0]
            DataNightString = os.path.split(DataPath)[1]
#             MasterDarkFilename = "MasterDark_"+self.tel.name+"_"+DataNightString+"_"+str(int(math.floor(self.exptime.to(u.s).value)))+".fits"
            MasterDarkFilename = 'MasterDark_{}_{}_{}.fits'.format(\
                                      self.tel.name,\
                                      DataNightString,\
                                      str(int(math.floor(self.exptime.to(u.s).value)))
                                      )
            MasterDarkFile  = os.path.join(self.tel.temp_file_path,\
                                           MasterDarkFilename)
            hdu_MasterDark = fits.PrimaryHDU(MasterDarkData)
            hdulist_MasterDark = fits.HDUList([hdu_MasterDark])
            hdulist_MasterDark.header = hdulist[0].header
            hdulist_MasterDark.header['history'] = \
                       "Combined {0} images to make this master dark.".format(\
                                                                    len(Darks))
            self.logger.debug("  Writing master dark file: {0}".format(\
                                                               MasterDarkFile))
            hdulist_MasterDark.writeto(MasterDarkFile)
        else:
            self.logger.error("No input dark files detected.")
        ## Now Subtract MasterDark from Image
        self.logger.debug("  Subtracting dark from image.")
        ImageData = hdulist_image[0].data
        DifferenceImage = ImageData - MasterDarkData
        hdulist_image[0].data = DifferenceImage
        hdulist_image.flush()
        self.logger.debug("  Median level of image = {0}".format(
                                                         np.median(ImageData)))
        self.logger.debug("  Median level of dark = {0}".format(\
                                                    np.median(MasterDarkData)))
        self.logger.debug("  Median level of dark subtracted = {0}".format(\
                                                   np.median(DifferenceImage)))
        hdulist_image.close()
        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done with dark subtraction in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Crop Image
    ##-------------------------------------------------------------------------
    def crop(self):
        '''
        Crop working image to region of interest.
        '''
        assert self.tel.ROI
        
        self.logger.info('Cropping image to {}'.format(self.tel.ROI))
        ## Parse ROI String
        try:
            MatchROI = re.match("\[?(\d{1,5}):(\d{1,5}),(\d{1,5}):(\d{1,5})\]?",\
                                self.tel.ROI)
        except:
            self.logger.warning("Could not parse ROI string: {}".format(\
                                                                 self.tel.ROI))
        else:
            self.crop_x1 = int(MatchROI.group(1))
            self.crop_x2 = int(MatchROI.group(2))
            self.crop_y1 = int(MatchROI.group(3))
            self.crop_y2 = int(MatchROI.group(4))
            self.logger.debug("  Cropping Image To [{0}:{1},{2}:{3}]".format(\
                                                    self.crop_x1, self.crop_x2,\
                                                    self.crop_y1, self.crop_y2))
#             hdulist = fits.open(self.working_file, mode="update")
            with fits.open(self.working_file, mode="update") as hdulist:
                hdulist[0].data = hdulist[0].data[self.crop_y1:self.crop_y2,self.crop_x1:self.crop_x2]
                hdulist.flush()
            self.cropped = True
            self.original_nXPix = self.nXPix
            self.original_nYPix = self.nYPix


    ##-------------------------------------------------------------------------
    ## Solve Astrometry Using astrometry.net
    ##-------------------------------------------------------------------------
    def solve_astrometry(self, timeout=None, downsample=4):
        '''
        Solve astrometry in the working image using the astrometry.net solver.
        '''
        start_time = datetime.datetime.now()
        self.logger.info("Attempting to solve WCS using Astrometry.net solver.")
        AstrometryCommand = ["solve-field", "-l", "5", "-O", "-p", "-T",
                             "-L", str(self.tel.pixel_scale.value*0.75),
                             "-H", str(self.tel.pixel_scale.value*1.25),
                             "-u", "arcsecperpix", "-z", str(downsample), self.working_file]
        AstrometrySTDOUT = open(os.path.join(self.tel.temp_file_path, 'astrometry_output.txt'), 'w')
        self.temp_files.append(os.path.join(self.tel.temp_file_path, 'astrometry_output.txt'))

        self.logger.debug('  Calling astrometry.net with: {}'.format(' '.join(AstrometryCommand)))
        if timeout:
            StartTime = datetime.datetime.now()
            astrometry_process = subprocess.Popen(AstrometryCommand, stdout=AstrometrySTDOUT, stderr=AstrometrySTDOUT)
            EndTime = datetime.datetime.now()
            duration = EndTime - StartTime

            wait = True
            while wait:
                if astrometry_process.poll() == None:
                    ## Process has not finished
                    EndTime = datetime.datetime.now()
                    duration = EndTime - StartTime
                    if duration.seconds > timeout:
                        astrometry_process.terminate()
                        wait = False
                        self.logger.warning("Astrometry.net timed out.")
                else:
                    ## Process finished
                    wait = False
            rtncode = astrometry_process.returncode
        else:
            StartTime = datetime.datetime.now()
            rtncode = subprocess.call(AstrometryCommand, stdout=AstrometrySTDOUT, stderr=AstrometrySTDOUT)
            EndTime = datetime.datetime.now()
            duration = EndTime - StartTime

        with open(os.path.join(self.tel.temp_file_path, 'astrometry_output.txt'), 'r') as AstrometrySTDOUT:
            output = AstrometrySTDOUT.readlines()

        if rtncode != 0:
            self.logger.warning("Astrometry.net failed.")
            for line in output:
                self.logger.warning('  astrometry.net output: {}'.format(line.strip('\n')))
        else:
            for line in output:
                self.logger.debug('  Astrometry.net Output: {}'.format(line.strip('\n')))
            total_process_time = (EndTime - StartTime).total_seconds()
            self.logger.debug("  Astrometry.net Processing Time: {:.1f} s".format(\
                                                           total_process_time))
            IsFieldCenter = re.search("Field center:\s\(RA\sH:M:S,\sDec D:M:S\)\s=\s\((\d{1,2}:\d{2}:\d{2}\.\d+,\s[+-]?\d{1,2}:\d{2}:\d{2}\.\d+)\)", ''.join(output))
            if IsFieldCenter:
                self.logger.info("  Astrometry.net field center is: {}".format(\
                                                    IsFieldCenter.group(1)))
            else:
                self.logger.warning("Could not parse field center from astrometry.net output.")
                for line in output:
                    self.logger.warning('  astrometry.net output: {}'.format(line.strip('\n')))

            NewFile = self.working_file.replace(self.file_ext, ".new")
            NewFitsFile = self.working_file.replace(self.file_ext, ".new.fits")
            if os.path.exists(NewFile):
                self.logger.debug("  Found {}".format(NewFile))
                self.logger.debug("  Astrometry.net succeeded")
                if os.path.exists(NewFitsFile): os.remove(NewFitsFile)
                os.rename(NewFile, NewFitsFile)
                self.astrometry_solved = True
                self.working_file = NewFitsFile
            else:
                self.logger.warning("No new file created by astrometry.net")
                self.astrometry_solved = False
            ## Add files created by astrometry.net to temp_files list
            self.temp_files.append(os.path.join(self.tel.temp_file_path,\
                                      self.raw_file_basename+".axy"))
            self.temp_files.append(os.path.join(self.tel.temp_file_path,\
                                      self.raw_file_basename+".wcs"))
            self.temp_files.append(os.path.join(self.tel.temp_file_path,\
                                      self.raw_file_basename+".solved"))
            self.temp_files.append(os.path.join(self.tel.temp_file_path,\
                                      self.raw_file_basename+".rdls"))
            self.temp_files.append(os.path.join(self.tel.temp_file_path,\
                                      self.raw_file_basename+".match"))
            self.temp_files.append(os.path.join(self.tel.temp_file_path,
                                      self.raw_file_basename+".corr"))
            self.temp_files.append(os.path.join(self.tel.temp_file_path,\
                                      self.raw_file_basename+".new.fits"))
            self.temp_files.append(os.path.join(self.tel.temp_file_path,\
                                      self.raw_file_basename+"-indx.xyls"))

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done with astrometry.net in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-----------------------------------------------------------------------------
    ## Determine Orientation from WCS
    ##-----------------------------------------------------------------------------
    def orientation_from_wcs(self):
        '''
        Given an astropy.wcs.WCS object, return a tuple containing the pixel scale,
        position angle (in degrees), and the flipped state (a boolean) of the image
        calculated from the PCn_m matrix (no distortions considered).
        '''
        ## Check that the WCS exists in the image_WCS object
        if not self.image_WCS:
            self.read_header()
        if not self.image_WCS:
            self.logger.warning('No WCS found in header.  No orientation calculated.')
        else:
            ## Check if the image_WCS is actually an astropy.wcs.WCS object
            if isinstance(self.image_WCS, wcs.WCS):
                ## By using the wcs to_header to make a new WCS object, we 
                ## ensure that the CD matrix, if it exists, is converted to PC
                header = wcs.WCS(self.image_WCS.to_header()).to_header()
                if ('CTYPE1' in header.keys()) and ('CTYPE2' in header.keys()) and\
                   ('WCSAXES' in header.keys()):
                    if (header['CTYPE1'][0:4] == 'RA--') or (header['CTYPE1'][0:4] == 'DEC-') and\
                       (header['CTYPE2'][0:4] == 'RA--') or (header['CTYPE2'][0:4] == 'DEC-') and\
                       (int(header['WCSAXES']) == 2) and\
                       ('PC1_1' in header.keys()) and\
                       ('PC1_2' in header.keys()) and\
                       ('PC2_1' in header.keys()) and\
                       ('PC2_2' in header.keys()) and\
                       ('CDELT1' in header.keys()) and\
                       ('CDELT2' in header.keys()):
                        ## If the wcs in header format meets all of the above
                        ## assumptions, do nothing and proceed to header analysis.
                        self.logger.debug('  Header has expected keywords')
#                         self.logger.debug('  {}'.format(header))
                    else:
                        self.logger.warning('WCS does not match expected contents.')
                        for key in header.keys():
                            self.logger.debug('  {:8s} = {}'.format(key, header[key]))
                        header = None
                        self.image_WCS = None
                else:
                    self.logger.warning('WCS does not have expected keywords.')
                    for key in header.keys():
                        self.logger.debug('  {:8s} = {}'.format(key, header[key]))
                    header = None
                    self.image_WCS = None

        if header:
            ## By using the wcs to_header to make a new WCS object, we convert CD to PC
            PC = wcs.WCS(self.image_WCS.to_header()).wcs.pc
            cdelt1 = float(header['CDELT1'])
            cdelt2 = float(header['CDELT2'])

            ## Determine Pixel Scale
            result1 = PC.dot(np.array([[0], [1]]))
            pixel_scale1 = cdelt1*(math.sqrt(result1[0][0]**2 + result1[1][0]**2))*3600.
            result2 = PC.dot(np.array([[1], [0]]))
            pixel_scale2 = cdelt2*(math.sqrt(result2[0][0]**2 + result2[1][0]**2))*3600.
            ## Just average the pixel scale values for each direction
            pixel_scale = np.mean([pixel_scale1, pixel_scale2]) * u.arcsec/u.pix
            self.wcs_pixel_scale = pixel_scale

            ## Determine Position Angle
            PCnorm = pixel_scale.to(u.deg/u.pix).value
            angles = np.array([90*u.deg.to(u.radian) - np.arccos(PC[0][0]/PCnorm),\
                               90*u.deg.to(u.radian) - np.arcsin(-1*PC[0][1]/PCnorm),\
                               90*u.deg.to(u.radian) - np.arcsin(PC[1][0]/PCnorm),\
                               90*u.deg.to(u.radian) - np.arccos(PC[1][1]/PCnorm),\
                              ]) * u.radian
            self.position_angle = angles[~np.isnan(angles)].mean().to(u.deg)

            ## Determine Flip State
            flipped = np.linalg.det(PC) > 0
            self.image_flipped = flipped
        else:
            self.wcs_pixel_scale = None
            self.position_angle = None
            self.image_flipped = None


    ##-------------------------------------------------------------------------
    ## Determine Pointing Error
    ##-------------------------------------------------------------------------
    def determine_pointing_error(self):
        '''
        Determine pointing error (difference between object's coordinates and
        solved WCS).
        '''
        self.logger.info("Detemining pointing error based on WCS solution")
        try:
            center_from_WCS = self.image_WCS.wcs_pix2world([[self.nXPix/2, self.nYPix/2]], 1)
            self.logger.debug("  Using coordinates of center point: {0} {1}".format(\
                                 center_from_WCS[0][0], center_from_WCS[0][1]))
            self.coordinate_of_center_pixel = coords.SkyCoord(\
                                              ra=center_from_WCS[0][0],\
                                              dec=center_from_WCS[0][1],\
                                              unit=(u.degree, u.degree),\
                                              frame='icrs')
            self.pointing_error = self.coordinate_of_center_pixel.separation(\
                                                   self.coordinate_from_header)
            self.logger.debug("  Header Coordinate: {}".format(self.coordinate_from_header.to_string(style='hmsdms', precision=1)))
            self.logger.debug("  Center Coordinate: {}".format(self.coordinate_of_center_pixel.to_string(style='hmsdms', precision=1)))
            self.logger.info("  Pointing Error is {:.2f} arcmin".format(\
                                                self.pointing_error.arcminute))
        except:
            self.logger.warning("Pointing error not calculated.")
        ## Flag pointing error
        try:
            if self.pointing_error.arcminute > self.tel.threshold_pointing_err.to(u.arcmin).value:
                self.flags['pointing error'] = True
            else:
                self.flags['pointing error'] = False
        except:
            pass


    ##-------------------------------------------------------------------------
    ## Run SExtractor
    ##-------------------------------------------------------------------------
    def run_SExtractor(self, assoc=False):
        '''
        Run SExtractor on image.
        '''
        start_time = datetime.datetime.now()
        assert type(self.tel.gain) == u.quantity.Quantity
        assert type(self.tel.pixel_scale) == u.quantity.Quantity

        if assoc and self.catalog_data:
            self.catalog_filter = None
            if self.header['FILTER']:
                if self.header['FILTER'] in self.catalog_data.keys():
                    self.catalog_filter = self.header['FILTER']
                elif self.tel.config['catalog'][self.filter] in self.catalog_data.keys():
                    self.catalog_filter = self.tel.config['catalog'][self.filter]
                else:
                    self.logger.warning('  Filter in header ({}), not found in catalog table.'.format(\
                                                        self.header['FILTER']))
                    self.catalog_filter = None
            else:
                self.catalog_filter = None
            ## Find Filter to Use
            if not self.catalog_filter:
                filters = ['r', 'R2mag', 'i', 'Imag', 'g', 'V', 'B', 'B2mag']
                for filt in filters:
                    if not self.catalog_filter and (filt in self.catalog_data.keys()):
                        self.logger.info('  Using {} filter for catalog magnitudes.'.format(filt))
                        self.catalog_filter = filt
                if not self.catalog_filter:
                    ## Choose whatever is in catalog
                    self.catalog_data.keys().remove('ID')
                    self.catalog_data.keys().remove('RA')
                    self.catalog_data.keys().remove('Dec')
                    self.catalog_filter = self.catalog_data.keys()[0]
            self.logger.info('  Using {} filter for catalog magnitudes.'.format(self.catalog_filter))

        ## Set up file names
        self.SExtractor_catalogfile = os.path.join(self.tel.temp_file_path,\
                                               self.raw_file_basename+".cat")
        self.temp_files.append(self.SExtractor_catalogfile)

        ## Remove catalog file from previous run of SExtractor (if it exists)
        if os.path.exists(self.SExtractor_catalogfile):
            os.remove(self.SExtractor_catalogfile)

        sextractor_output_param_file = os.path.join(self.tel.temp_file_path,\
                                                   '{}.param'.format(self.raw_file_basename))
        if os.path.exists(sextractor_output_param_file):
            os.remove(sextractor_output_param_file)
        defaultparamsFO = open(sextractor_output_param_file, 'w')
        params = [
                  'XWIN_IMAGE', 'YWIN_IMAGE', 
                  'AWIN_IMAGE', 'BWIN_IMAGE', 'FWHM_IMAGE', 'THETAWIN_IMAGE',
                  'ERRAWIN_IMAGE', 'ERRBWIN_IMAGE', 'ERRTHETAWIN_IMAGE',
                  'ELONGATION', 'ELLIPTICITY',
                  'FLUX_AUTO', 'FLUXERR_AUTO', 'MAG_AUTO', 'MAGERR_AUTO',
                  'FLAGS', 'FLAGS_WEIGHT', 'FLUX_RADIUS'
                 ]
        if assoc: params.append('VECTOR_ASSOC(3)')
        for param in params:
            defaultparamsFO.write(param+'\n')
        defaultparamsFO.close()
        self.temp_files.append(sextractor_output_param_file)

        ## Compare input parameters dict to default
        SExtractor_default = {
                             'CATALOG_NAME': self.SExtractor_catalogfile,
                             'CATALOG_TYPE': 'FITS_LDAC',
                             'PARAMETERS_NAME': sextractor_output_param_file,
                             'GAIN': self.tel.gain.value,
                             'GAIN_KEY': 'GAIN',
                             'PIXEL_SCALE': '{:.3f}'.format(self.tel.pixel_scale.value),
                             'CHECKIMAGE_TYPE': 'NONE',
                            }

        ## Use optional sextractor params
        if not self.tel.SExtractor_params:
            SExtractor_params = SExtractor_default
        else:
            SExtractor_params = self.tel.SExtractor_params
            for key in SExtractor_default.keys():
                if not key in self.tel.SExtractor_params.keys():
                    SExtractor_params[key] = SExtractor_default[key]

        if assoc:
            ## Create Assoc file with pixel coordinates of catalog stars
            assoc_file = os.path.join(self.tel.temp_file_path, self.raw_file_basename+'_assoc.txt')
            self.temp_files.append(assoc_file)
            if os.path.exists(assoc_file): os.remove(assoc_file)

            with open(assoc_file, 'w') as assocFO:
                for star in self.catalog_data:
                    pix = self.image_WCS.wcs_world2pix([[star['RA'], star['Dec']]], 1)
                    try:
                        assocFO.write('{:8.1f} {:8.1f} {:8.1f}\n'.format(\
                                                        pix[0][0], pix[0][1],\
                                                        star[self.catalog_filter],\
                                                        ))
                    except ValueError:
                        assocFO.write('{:8.1f} {:8.1f} {:8.1f}\n'.format(\
                                                        pix[0][0], pix[0][1],\
                                                        float('nan'),\
                                                        ))
                    except:
                        print('Skipped: {} {} {}'.format(pix[0][0], pix[0][1], star[self.catalog_filter]))

            ## Add ASSOC parameters
            original_params = self.tel.SExtractor_params
            self.tel.SExtractor_params['ASSOC_NAME'] = assoc_file
            self.tel.SExtractor_params['ASSOC_DATA'] = '0'
            self.tel.SExtractor_params['ASSOC_PARAMS'] = '1,2,3'
            self.tel.SExtractor_params['ASSOC_RADIUS'] = '2.0'
            self.tel.SExtractor_params['ASSOC_TYPE'] = 'NEAREST'
            self.tel.SExtractor_params['ASSOCSELEC_TYPE'] = 'MATCHED'

        ## Run SExtractor
        SExtractorCommand = ["sex", self.working_file]
        for key in SExtractor_params.keys():
            SExtractorCommand.append('-{}'.format(key))
            SExtractorCommand.append('{}'.format(SExtractor_params[key]))
        self.logger.info("Invoking SExtractor")
        self.logger.debug("  SExtractor command: {}".format(' '.join(SExtractorCommand)))
        try:
            SExSTDOUT = subprocess.check_output(SExtractorCommand,\
                             stderr=subprocess.STDOUT, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            self.logger.error("SExtractor failed.  Command: {}".format(e.cmd))
            self.logger.error("SExtractor failed.  Returncode: {}".format(e.returncode))
            self.logger.error("SExtractor failed.  Output: {}".format(e.output))
            self.SExtractor_results = None
            self.SExtractor_background = None
            self.SExtractor_background_RMS = None
        except:
            self.logger.error("SExtractor process failed: {0} {1} {2}".format(\
                      sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
            self.SExtractor_results = None
            self.SExtractor_background = None
            self.SExtractor_background_RMS = None
        else:
            for line in SExSTDOUT.splitlines():
                line.replace("[1A", "")
                line.replace("[1M>", "")
                MatchVersion = re.search('SExtractor (\d+\.\d+\.\d+) started on', line)
                if MatchVersion:
                    self.logger.info('  SExtractor version = {}'.format(MatchVersion.group(1)))
                if not re.match(".*Setting up background map.*", line) and\
                   not re.match(".*Line:\s[0-9]*.*", line):
                    self.logger.debug("  SExtractor Output: {}".format(line))
            ## Extract Number of Stars from SExtractor Output
            pos = SExSTDOUT.find("sextracted ")
            IsSExCount = re.match("\s*([0-9]+)\s+", SExSTDOUT[pos+11:pos+21])
            if IsSExCount:
                nSExtracted = int(IsSExCount.group(1))
                self.logger.info("  SExtractor found {0} sources.".format(nSExtracted))
            ## Extract Background Level from SExtractor Output
            pos = SExSTDOUT.find("Background: ")
            IsSExBkgnd = re.match("\s*([0-9\.]+)\s*", SExSTDOUT[pos+11:pos+21])
            if IsSExBkgnd:
                self.SExtractor_background = float(IsSExBkgnd.group(1))
                self.logger.info("  SExtractor background is {0:.1f}".format(\
                                                   self.SExtractor_background))
            else:
                self.SExtractor_background = None
            ## Extract Background RMS from SExtractor Output
            IsSExtractor_background_RMS = re.match("\s*RMS:\s([0-9\.]+)\s*",\
                                                   SExSTDOUT[pos+21:pos+37])
            if IsSExtractor_background_RMS:
                self.SExtractor_background_RMS = float(IsSExtractor_background_RMS.group(1))
                self.logger.info("  SExtractor background RMS is {0:.1f}".format(\
                                               self.SExtractor_background_RMS))
            else:
                self.SExtractor_background_RMS = None

            ## If No Output Catalog Created ...
            if not os.path.exists(self.SExtractor_catalogfile):
                self.logger.warning("SExtractor failed to create catalog.")
                self.SExtractor_catalogfile = None

            ## Read FITS_LDAC SExtractor Catalog
            self.logger.debug("  Reading SExtractor output catalog.")
            hdu = fits.open(self.SExtractor_catalogfile)
            results = table.Table(hdu[2].data)

            rows_to_remove = []
            for i in range(0,len(results)):
                if results['FLAGS'][i] != 0:
                    rows_to_remove.append(i)
            if len(rows_to_remove) > 0:
                results.remove_rows(rows_to_remove)

            self.SExtractor_results = results
            SExImageRadius = []
            SExAngleInImage = []
            assoc_x = []
            assoc_y = []
            assoc_catmag = []
            for star in self.SExtractor_results:
                SExImageRadius.append(math.sqrt((self.nXPix/2-star['XWIN_IMAGE'])**2 +\
                                                (self.nYPix/2-star['YWIN_IMAGE'])**2))
                SExAngleInImage.append(math.atan((star['XWIN_IMAGE']-self.nXPix/2) /\
                                                 (self.nYPix/2-star['YWIN_IMAGE']))*180.0/math.pi)
                if assoc:
                    assoc_x.append(star['VECTOR_ASSOC'][0])
                    assoc_y.append(star['VECTOR_ASSOC'][1])
                    assoc_catmag.append(star['VECTOR_ASSOC'][2])
            self.SExtractor_results.add_column(table.Column(\
                                    data=SExImageRadius, name='ImageRadius'))
            self.SExtractor_results.add_column(table.Column(\
                                    data=SExAngleInImage, name='AngleInImage'))
            self.n_stars_SExtracted = len(self.SExtractor_results)
            self.logger.info("  Read in {0} stars from SExtractor catalog (after filtering).".format(\
                                                      self.n_stars_SExtracted))
            if assoc:
                self.SExtractor_results.add_column(table.Column(\
                                                   data=assoc_x, name='assoc_x'))
                self.SExtractor_results.add_column(table.Column(\
                                                   data=assoc_y, name='assoc_y'))
                self.SExtractor_results.add_column(table.Column(\
                                                   data=assoc_catmag, name='assoc_catmag'))
                self.tel.SExtractor_params = original_params

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done running SExtractor in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Determine Image FWHM from SExtractor Catalog
    ##-------------------------------------------------------------------------
    def determine_FWHM(self):
        '''
        Determine typical FWHM of image from SExtractor results.
        '''
        if self.n_stars_SExtracted > 1:
            self.logger.info('Analyzing SExtractor results to determine typical image quality.')
            if not self.tel.PSF_measurement_radius:
                IQRadiusFactor = 1.0
                DiagonalRadius = math.sqrt((self.nXPix/2)**2+(self.nYPix/2)**2)
                self.tel.PSF_measurement_radius = DiagonalRadius*IQRadiusFactor
                self.logger.info('  Using all stars in image.')

            CentralFWHMs = [star['FWHM_IMAGE']\
                            for star in self.SExtractor_results\
                            if (star['ImageRadius'] <= self.tel.PSF_measurement_radius.to(u.pix).value)]
            CentralEllipticities = [star['ELLIPTICITY']\
                                    for star in self.SExtractor_results\
                                    if (star['ImageRadius'] <= self.tel.PSF_measurement_radius.to(u.pix).value)]
            CentralSNRs = [star['FLUX_AUTO'] / star['FLUXERR_AUTO']\
                           for star in self.SExtractor_results\
                           if (star['ImageRadius'] <= self.tel.PSF_measurement_radius.to(u.pix).value)]

            if len(CentralFWHMs) > 3:
                self.FWHM_mode = mode(CentralFWHMs, 0.2) * u.pix
                self.FWHM_median = np.median(CentralFWHMs) * u.pix
                self.FWHM_average = np.average(CentralFWHMs, weights=CentralSNRs) * u.pix
                self.FWHM = self.FWHM_average
                self.ellipticity_mode = mode(CentralEllipticities, 0.05) 
                self.ellipticity_median = np.median(CentralEllipticities)
                self.ellipticity_average = np.average(CentralEllipticities, weights=CentralSNRs)
                self.ellipticity = self.ellipticity_average
                self.logger.debug("  Using {0} stars in central {1} to determine PSF quality.".format(\
                                                                len(CentralFWHMs),\
                                                                self.tel.PSF_measurement_radius))
                self.logger.info("  Mode FWHM in inner region is {0:.2f} pixels".format(\
                                                        self.FWHM_mode.to(u.pix).value))
                self.logger.info("  Median FWHM in inner region is {0:.2f} pixels".format(\
                                                        self.FWHM_median.to(u.pix).value))
                self.logger.info("  Average FWHM in inner region is {0:.2f} pixels".format(\
                                                        self.FWHM_average.to(u.pix).value))

                self.logger.info("  Mode Ellipticity in inner region is {0:.2f}".format(\
                                                                 self.ellipticity_mode))
                self.logger.info("  Median Ellipticity in inner region is {0:.2f}".format(\
                                                                 self.ellipticity_median))
                self.logger.info("  Average Ellipticity in inner region is {0:.2f}".format(\
                                                                 self.ellipticity_average))
            else:
                self.logger.warning("  Only detected {} stars in central region.  No FWHM or ellipticity calculated.".format(len(CentralFWHMs)))
                self.FWHM_mode = None
                self.FWHM_median = None
                self.FWHM_average = None
                self.FWHM = None
                self.ellipticity_mode = None
                self.ellipticity_median = None
                self.ellipticity_average = None
                self.ellipticity = None
        else:
            self.FWHM_mode = None
            self.FWHM_median = None
            self.FWHM_average = None
            self.FWHM = None
            self.ellipticity_mode = None
            self.ellipticity_median = None
            self.ellipticity_average = None
            self.ellipticity = None
        ## Flag FWHM
        try:
            if self.FWHM > self.tel.threshold_FWHM.to(u.pix, equivalencies=self.tel.pixel_scale_equivalency):
                self.flags['FWHM'] = True
            else:
                self.flags['FWHM'] = False
        except:
            pass
        ## Check ellipticity
        try:
            if self.ellipticity > self.tel.threshold_ellipticity:
                self.flags['ellipticity'] = True
            else:
                self.flags['ellipticity'] = False
        except:
            pass


    ##-------------------------------------------------------------------------
    ## Make PSF Statistics Plots
    ##-------------------------------------------------------------------------
    def make_PSF_plot(self, filename=None):
        '''
        Make various plots for analysis of image quality.
        '''
        start_time = datetime.datetime.now()

        if not self.FWHM:
            self.logger.warning('No FWHM statistics found.  Skipping PSF plot creation.')
            self.PSF_plot_filename = None
            self.PSF_plot_file = None
            return
        else:
            if filename:
                self.PSF_plot_filename = filename
            else:
                self.PSF_plot_filename = self.raw_file_basename+'_PSFinfo.png'
            self.logger.info('Generating plots of PSF statistics: {}'.format(self.PSF_plot_filename))
            self.PSF_plot_file = os.path.join(self.tel.plot_file_path, self.PSF_plot_filename)

            ellip_threshold = 0.15
            star_angles = [star['THETAWIN_IMAGE'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
            image_angles = [star['AngleInImage'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
            star_x = [star['XWIN_IMAGE'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
            star_y = [star['YWIN_IMAGE'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
            uncorrected_diffs = [star['THETAWIN_IMAGE']-star['AngleInImage'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]

            CentralFWHMs = [star['FWHM_IMAGE'] for star in self.SExtractor_results if (star['ImageRadius'] <= self.tel.PSF_measurement_radius.to(u.pix).value)]
            CentralEllipticities = [star['ELLIPTICITY'] for star in self.SExtractor_results if (star['ImageRadius'] <= self.tel.PSF_measurement_radius.to(u.pix).value)]

            nstars = len(star_angles)
            self.logger.debug('  Found {} stars with ellipticity greater than {:.2f}.'.format(\
                                                          nstars, ellip_threshold))

            angle_diffs = []
            for angle in uncorrected_diffs:
                if angle < -90:
                    angle_diffs.append(angle + 90.)
                elif angle > 90:
                    angle_diffs.append(angle - 90.)
                else:
                    angle_diffs.append(angle)
            angle_binsize = 10
            diff_hist, diff_bins = np.histogram(angle_diffs, bins=angle_binsize*(np.arange(37)-18))
            angle_hist, angle_bins = np.histogram(star_angles, bins=angle_binsize*(np.arange(37)-18))
            angle_centers = (diff_bins[:-1] + diff_bins[1:]) / 2

            ellip_binsize = 0.05
            ellip_bmin = math.floor(min(CentralEllipticities)/ellip_binsize)*ellip_binsize - ellip_binsize/2.
            ellip_bmax = math.ceil(max(CentralEllipticities)/ellip_binsize)*ellip_binsize + ellip_binsize/2.
            ellip_bins = np.arange(ellip_bmin,ellip_bmax,ellip_binsize)
            ellip_hist, ellip_bins = np.histogram(CentralEllipticities, bins=ellip_bins)
            ellip_centers = (ellip_bins[:-1] + ellip_bins[1:]) / 2

            fwhm_binsize = 0.2
            fwhm_bmin = math.floor(min(CentralFWHMs)/fwhm_binsize)*fwhm_binsize - fwhm_binsize/2.
            fwhm_bmax = math.ceil(max(CentralFWHMs)/fwhm_binsize)*fwhm_binsize + fwhm_binsize/2.
            fwhm_bins = np.arange(fwhm_bmin,fwhm_bmax,fwhm_binsize)
            fwhm_hist, fwhm_bins = np.histogram(CentralFWHMs, bins=fwhm_bins)
            fwhm_centers = (fwhm_bins[:-1] + fwhm_bins[1:]) / 2

            star_angle_mean = np.mean(star_angles)
            star_angle_median = np.median(star_angles)
            angle_diff_mean = np.mean(angle_diffs)
            angle_diff_median = np.median(angle_diffs)
            self.logger.debug('  Mean Stellar PA = {:.0f}'.format(star_angle_mean))
            self.logger.debug('  Median Stellar PA = {:.0f}'.format(star_angle_median))
            self.logger.debug('  Mean Difference Angle = {:.0f}'.format(angle_diff_mean))
            self.logger.debug('  Median Difference Angle = {:.0f}'.format(angle_diff_median))

            if self.PSF_plot_file:
                self.logger.debug('  Generating figure {}'.format(self.PSF_plot_file))

                pyplot.ioff()
                fig = pyplot.figure(figsize=(10,11), dpi=100)

                TopLeft = pyplot.axes([0.000, 0.750, 0.465, 0.235])
                pyplot.title('Histogram of FWHM Values for {}'.format(self.raw_file_name), size=10)
                pyplot.bar(fwhm_centers, fwhm_hist, align='center', width=0.7*fwhm_binsize)
                pyplot.plot([self.FWHM_mode.to(u.pix).value, self.FWHM_mode.to(u.pix).value], [0, 1.1*max(fwhm_hist)],\
                            'ro-', linewidth=2, label='Mode FWHM', alpha=0.4)
                pyplot.plot([self.FWHM_median.to(u.pix).value, self.FWHM_median.to(u.pix).value], [0, 1.1*max(fwhm_hist)],\
                            'go-', linewidth=2, label='Median FWHM', alpha=0.4)
                pyplot.plot([self.FWHM_average.to(u.pix).value, self.FWHM_average.to(u.pix).value], [0, 1.1*max(fwhm_hist)],\
                            'bo-', linewidth=2, label='Mode FWHM', alpha=1.0)
                pyplot.xlabel('FWHM (pixels)', size=10)
                pyplot.ylabel('N Stars', size=10)
                pyplot.xlim(0,np.percentile(CentralFWHMs, 95)+1)
                pyplot.xticks(size=10)
                pyplot.yticks(size=10)

                TopRight = pyplot.axes([0.535, 0.750, 0.465, 0.235])
                pyplot.title('Histogram of Elliptiticty Values for {}'.format(self.raw_file_name), size=10)
                pyplot.bar(ellip_centers, ellip_hist, align='center', width=0.7*ellip_binsize)
                pyplot.plot([self.ellipticity_mode, self.ellipticity_mode], [0, 1.1*max(ellip_hist)],\
                            'ro-', linewidth=2, label='Mode Ellipticity', alpha=0.4)
                pyplot.plot([self.ellipticity_median, self.ellipticity_median], [0, 1.1*max(ellip_hist)],\
                            'go-', linewidth=2, label='Median Ellipticity', alpha=0.4)
                pyplot.plot([self.ellipticity_average, self.ellipticity_average], [0, 1.1*max(ellip_hist)],\
                            'bo-', linewidth=2, label='Mode Ellipticity', alpha=1.0)
                pyplot.xlabel('Ellipticity', size=10)
                pyplot.ylabel('N Stars', size=10)
                pyplot.xlim(0,1)
                pyplot.xticks(0.1*np.arange(11), size=10)
                pyplot.yticks(size=10)

                MiddleLeft = pyplot.axes([0.000, 0.375, 0.465, 0.320])
                MiddleLeft.set_aspect('equal')
                pyplot.title('Average FWHM scaled from {:.1f} pix to {:.1f} pix'.format(0.8*self.FWHM.to(u.pix).value, 2.0*self.FWHM.to(u.pix).value), size=10)
                if self.n_stars_SExtracted > 20000:
                    gridsize = 20
                else:
                    gridsize = 10
                pyplot.hexbin(self.SExtractor_results['XWIN_IMAGE'].data,\
                              self.SExtractor_results['YWIN_IMAGE'].data,\
                              self.SExtractor_results['FWHM_IMAGE'].data,\
                              gridsize=gridsize,\
                              mincnt=5,\
                              vmin=0.8*self.FWHM.to(u.pix).value,\
                              vmax=2.0*self.FWHM.to(u.pix).value,\
                              alpha=0.5,\
                              cmap='Reds')
    #             center_region = pyplot.Circle((self.nXPix/2, self.nYPix/2),\
    #                                    radius=self.tel.PSF_measurement_radius/self.nXPix,\
    #                                    color='k')
    #             MiddleLeft.add_artist(center_region)
                pyplot.xlabel('X Pixels', size=10)
                pyplot.ylabel('Y Pixels', size=10)
                pyplot.xlim(0,self.nXPix)
                pyplot.ylim(0,self.nYPix)
                pyplot.xticks(size=10)
                pyplot.yticks(size=10)

                MiddleRight = pyplot.axes([0.535, 0.375, 0.465, 0.320])
                MiddleRight.set_aspect('equal')
                pyplot.title('Average Ellipticity scaled from 0.25 to 0.75', size=10)
                if self.n_stars_SExtracted > 20000:
                    gridsize = 20
                else:
                    gridsize = 10
                pyplot.hexbin(self.SExtractor_results['XWIN_IMAGE'].data,\
                              self.SExtractor_results['YWIN_IMAGE'].data,\
                              self.SExtractor_results['ELLIPTICITY'].data,\
                              gridsize=gridsize,\
                              mincnt=5,\
                              vmin=0.25, vmax=0.75,\
                              alpha=0.5,\
                              cmap='Reds')
    #             MiddleRight.add_artist(center_region)
                pyplot.xlabel('X Pixels', size=10)
                pyplot.ylabel('Y Pixels', size=10)
                pyplot.xlim(0,self.nXPix)
                pyplot.ylim(0,self.nYPix)
                pyplot.xticks(size=10)
                pyplot.yticks(size=10)

                BottomLeft = pyplot.axes([0.000, 0.0, 0.465, 0.320])
                pyplot.title('Correlation of Ellipticity with Image Radius', size=10)
                if self.n_stars_SExtracted > 4000:
                    bins = 40
                elif self.n_stars_SExtracted > 2000:
                    bins = 30
                else:
                    bins = 20
                pyplot.hist2d(self.SExtractor_results['ImageRadius'],\
                              self.SExtractor_results['ELLIPTICITY'],\
                              bins=bins, cmap='binary')
                pyplot.xlabel('r (pixels)', size=10)
                pyplot.ylabel('Ellipticity', size=10)
                pyplot.xlim(0, math.sqrt(self.nXPix**2 + self.nYPix**2)/2.)
                pyplot.ylim(0, 1.0)
                pyplot.xticks(size=10)
                pyplot.yticks(size=10)

                BottomRight = pyplot.axes([0.535, 0.0, 0.465, 0.320])
                BottomRight.set_aspect('equal')
                pyplot.title('Correlation Between PSF Angle and Position in Image', size=10)
                pyplot.hist2d(star_angles, image_angles, bins=bins, cmap='binary')
                pyplot.xlabel('Stellar PSF PA', size=10)
                pyplot.ylabel('Image PA', size=10)
                pyplot.xlim(-100,100)
                pyplot.xticks(30*(np.arange(7)-3), size=10)
                pyplot.ylim(-100,100)
                pyplot.yticks(30*(np.arange(7)-3), size=10)

                pyplot.savefig(self.PSF_plot_file, dpi=100, bbox_inches='tight', pad_inches=0.10)
                pyplot.close(fig)

            end_time = datetime.datetime.now()
            elapzed_time = end_time - start_time
            self.logger.info('  Done making PSF plot in {:.1f} s'.format(elapzed_time.total_seconds()))

    ##-------------------------------------------------------------------------
    ## Is the Image Blank
    ##-------------------------------------------------------------------------
    def is_blank(self, threshold=None, area=None):
        '''
        '''
        self.logger.info('Checking if image is blank')
        nstars_threshold = 5

        ## Save original SExtractor parameters
        if 'DETECT_THRESH' in self.tel.SExtractor_params.keys():
            dt = self.tel.SExtractor_params['DETECT_THRESH']
        else:
            dt = None
        if 'ANALYSIS_THRESH' in self.tel.SExtractor_params.keys():
            at = self.tel.SExtractor_params['ANALYSIS_THRESH']
        else:
            at = None
        if 'DETECT_MINAREA' in self.tel.SExtractor_params.keys():
            da = self.tel.SExtractor_params['DETECT_MINAREA']
        else:
            da = None
        ## Set new (temporary) parmaters
        if threshold:
            self.tel.SExtractor_params['DETECT_THRESH'] = threshold
            self.tel.SExtractor_params['ANALYSIS_THRESH'] = threshold
        if area:
            self.tel.SExtractor_params['DETECT_MINAREA'] = area
        
        ## Run SExtractor
        self.run_SExtractor()
        stars = [entry for entry in self.SExtractor_results if entry['FLAGS'] == 0]
#        filtered_stars = [star for star in stars if star['BWIN_IMAGE'] > 1.0]

        ## Edit SExtractor parameters back to original state
        if dt:
            self.tel.SExtractor_params['DETECT_THRESH'] = dt
        else:
            if 'DETECT_THRESH' in self.tel.SExtractor_params: del self.tel.SExtractor_params['DETECT_THRESH']
        if at:
            self.tel.SExtractor_params['ANALYSIS_THRESH'] = at
        else:
            if 'ANALYSIS_THRESH' in self.tel.SExtractor_params: del self.tel.SExtractor_params['ANALYSIS_THRESH']
        if da:
            self.tel.SExtractor_params['DETECT_MINAREA'] = da
        else:
            if 'DETECT_MINAREA' in self.tel.SExtractor_params: del self.tel.SExtractor_params['DETECT_MINAREA']

        ## Reset all SExtractor results to None, so we don't confuse these results with meaningful ones
        self.SExtractor_catalogfile = None
        self.SExtractor_results = None
        self.n_stars_SExtracted = None
        self.SExtractor_background = None
        self.SExtractor_background_RMS = None

        ## If few stars found, the image is blank
        if len(stars) < nstars_threshold:
            self.logger.warning('  Only {} bright stars detected.  Image appears blank'.format(len(filtered_stars)))
            self.flags['other'] = True
            return True
        else:
            self.logger.info('  Found {} bright stars.'.format(len(stars)))
            self.flags['other'] = False
            return False


    ##-------------------------------------------------------------------------
    ## Run SCAMP
    ##-------------------------------------------------------------------------
    def run_SCAMP(self):
        '''
        Run SCAMP on SExtractor output catalog.
        '''
        start_time = datetime.datetime.now()
        ## Parameters for SCAMP
        if self.tel.SCAMP_aheader:
            SCAMP_aheader = self.tel.SCAMP_aheader
        else:
            SCAMP_aheader = 'scamp.ahead'

        SCAMP_default = {
                        'AHEADER_GLOBAL': SCAMP_aheader,
                        'SAVE_REFCATALOG': 'N',
                        'REFOUT_CATPATH': self.tel.temp_file_path,
                        'MERGEDOUTCAT_NAME': os.path.join(self.tel.temp_file_path, 'scamp.cat'),
                        'MERGEDOUTCAT_TYPE': 'FITS_LDAC',
                        'CHECKPLOT_RES': '1200,1200',
                        'CHECKPLOT_TYPE': 'NONE',
                        'SOLVE_PHOTOM': 'Y',
                        'ASTRINSTRU_KEY': 'QRUNID',
                        'WRITE_XML': 'N',
                        'XML_NAME': os.path.join(self.tel.temp_file_path, 'scamp.xml'),
                        }

        if not self.tel.SCAMP_params:
            SCAMP_params = SCAMP_default
        else:
            SCAMP_params = self.tel.SCAMP_params
            for key in SCAMP_default.keys():
                if not key in self.tel.SCAMP_params.keys():
                    SCAMP_params[key] = SCAMP_default[key]

        SCAMPCommand = ["scamp", self.SExtractor_catalogfile]
        for key in SCAMP_params.keys():
            SCAMPCommand.append('-{}'.format(key))
            SCAMPCommand.append('{}'.format(SCAMP_params[key]))
        self.logger.info("Running SCAMP")
        if SCAMP_aheader:
            self.logger.info("  Using SCAMP aheader file: {}".format(SCAMP_aheader))
        self.logger.debug("  SCAMP command: {}".format(' '.join(SCAMPCommand)))
        try:
            SCAMP_STDOUT = subprocess.check_output(SCAMPCommand,\
                             stderr=subprocess.STDOUT, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            self.logger.error("SCAMP failed.  Command: {}".format(e.cmd))
            self.logger.error("SCAMP failed.  Returncode: {}".format(e.returncode))
            self.logger.error("SCAMP failed.  Output: {}".format(e.output))
        except:
            self.logger.error("SCAMP process failed: {0}".format(sys.exc_info()[0]))
            self.logger.error("SCAMP process failed: {0}".format(sys.exc_info()[1]))
            self.logger.error("SCAMP process failed: {0}".format(sys.exc_info()[2]))
        else:
            StartAstrometricStats = False
            EndAstrometricStats = False
            for line in SCAMP_STDOUT.splitlines():
                MatchVersion = re.search('SCAMP (\d+\.\d+\.\d+) started on', line)
                if MatchVersion:
                    self.logger.info('  SCAMP version = {}'.format(MatchVersion.group(1)))
                if re.search('Astrometric stats \(external\)', line):
                    StartAstrometricStats = True
                if re.search('Generating astrometric plots', line):
                    EndAstrometricStats = True
                if StartAstrometricStats and not EndAstrometricStats:
                    self.logger.info("  SCAMP Output: "+line)
                else:
                    self.logger.debug("  SCAMP Output: "+line)

        ## Populate FITS header with SCAMP derived header values in .head file
        head_filename = '{}.head'.format(self.raw_file_basename)
        head_file = os.path.join(self.tel.temp_file_path, head_filename)
        if os.path.exists(head_file):
            self.temp_files.append(head_file)
            self.logger.info('  Writing SCAMP results to fits header on {}'.format(\
                                                            self.working_file))
            missfits_cmd = 'missfits -SAVE_TYPE REPLACE -WRITE_XML N {}'.format(\
                                                             self.working_file)
            self.logger.debug('  Running: {}'.format(missfits_cmd))
            output = subprocess.check_output(missfits_cmd, shell=True,\
                             stderr=subprocess.STDOUT, universal_newlines=True)
            output = str(output)
            for line in output.splitlines():
                self.logger.debug(line)
            self.SCAMP_successful = True
        else:
            self.logger.critical('No .head file found from SCAMP.  SCAMP failed.')
            self.SCAMP_successful = False

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done running SCAMP in {:.1f} s'.format(elapzed_time.total_seconds()))
        return self.SCAMP_successful


    ##-------------------------------------------------------------------------
    ## Run SWarp
    ##-------------------------------------------------------------------------
    '''
    Run SWarp on the image (after SCAMP distortion solution) to de-distort it.
    '''
    def run_SWarp(self):
        start_time = datetime.datetime.now()
        ## Parameters for SWarp
        swarp_file = os.path.join(self.tel.temp_file_path, 'swarpped.fits')
        if os.path.exists(swarp_file): os.remove(swarp_file)
        SWarp_params = {'IMAGEOUT_NAME': swarp_file,
                        'COPY_KEYWORDS': 'FILTER,OBJECT,AIRMASS,DATE-OBS,LAT-OBS,LONG-OBS,ALT-OBS,RA,DEC',
                        'WRITE_XML': 'N',
                        'XML_NAME': os.path.join(self.tel.temp_file_path, 'swarp.xml'),
                        'FSCALASTRO_TYPE': 'NONE',
                        'SUBTRACT_BACK': 'N',
                       }
        SWarpCommand = ["swarp", self.working_file]
        for key in SWarp_params.keys():
            SWarpCommand.append('-{}'.format(key))
            SWarpCommand.append('{}'.format(SWarp_params[key]))
        self.logger.info("Running SWarp.")
        self.logger.debug("  SWarp command: {}".format(' '.join(SWarpCommand)))
        try:
            SWarp_STDOUT = subprocess.check_output(SWarpCommand,\
                                                   stderr=subprocess.STDOUT,\
                                                   universal_newlines=True)
        except subprocess.CalledProcessError as e:
            self.logger.error("SWarp failed.  Command: {}".format(e.cmd))
            self.logger.error("SWarp failed.  Returncode: {}".format(e.returncode))
            self.logger.error("SWarp failed.  Output: {}".format(e.output))
        except:
            self.logger.error("SWarp process failed: {0}".format(sys.exc_info()[0]))
            self.logger.error("SWarp process failed: {0}".format(sys.exc_info()[1]))
            self.logger.error("SWarp process failed: {0}".format(sys.exc_info()[2]))
        else:
            for line in SWarp_STDOUT.splitlines():
                MatchVersion = re.search('SWarp (\d+\.\d+\.\d+) started on', line)
                if MatchVersion:
                    self.logger.info('  SWarp version = {}'.format(MatchVersion.group(1)))
                if not re.search('Resampling line', line) and not re.search('Setting up background map at', line):
                    self.logger.debug("  SWarp Output: "+line)
        ## Replace working_file with SWarp output file
        if os.path.exists(swarp_file):
            self.logger.debug('  SWarp process succeeded.')
            self.logger.debug('  Moving SWarpped file to working file.')
            if os.path.exists(self.working_file): os.remove(self.working_file)
            os.rename(swarp_file, self.working_file)
            assert os.path.exists(self.working_file)

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done running SWarp in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Get Vizier Catalog
    ##-------------------------------------------------------------------------
    def get_catalog(self, max_stars=50000):
        '''
        Get a catalog using astroquery
        '''
        start_time = datetime.datetime.now()
        if self.image_WCS:
            import astroquery
            import astroquery.vizier
            catalog = self.tel.catalog_info['name']
            self.logger.info("Querying Vizier for {} catalog.".format(catalog))
            if catalog == 'USNO-B1.0':
                columns = ['_RAJ2000', '_DEJ2000', 'USNO-B1.0', 'B1mag', 'B2mag', 'R1mag', 'R2mag', 'Imag']
                catfilt = str(self.tel.catalog_info[self.filter])
                upperlimit = '<{:.1f}'.format(self.tel.catalog_info['magmax'])
                column_filters = {catfilt:upperlimit}
            elif catalog == 'UCAC4':
                columns = ['_RAJ2000', '_DEJ2000', 'UCAC4', 'Bmag', 'Vmag', 'gmag', 'rmag', 'imag']
                catfilt = str(self.tel.catalog_info[self.filter])
                upperlimit = '<{:.1f}'.format(self.tel.catalog_info['magmax'])
                column_filters = {catfilt:upperlimit}

            self.logger.debug('  Getting columns: {}'.format(columns))
            self.logger.debug('  Using column_filters: {}'.format(column_filters))

            viz = astroquery.vizier.Vizier(catalog=catalog,\
                                           columns=columns,\
                                           column_filters=column_filters)
            viz.ROW_LIMIT = max_stars

            center_from_WCS = self.image_WCS.wcs_pix2world([[self.nXPix/2, self.nYPix/2]], 1)
            self.coordinate_of_center_pixel = coords.SkyCoord(\
                                              ra=center_from_WCS[0][0],\
                                              dec=center_from_WCS[0][1],\
                                              unit=(u.degree, u.degree),\
                                              frame='icrs')
            footprint = self.image_WCS.calc_footprint()
            RAs = [val[0] for val in footprint]
            DECs = [val[1] for val in footprint]
            dRA = (max(RAs) - min(RAs))
            if dRA > 180:
                dRA = (min(RAs)+360. - max(RAs))
            dRA = dRA*math.cos(center_from_WCS[0][1]*u.deg.to(u.radian))
            dDEC = (max(DECs) - min(DECs))
            self.logger.debug("  Center Coordinate: {}".format(self.coordinate_of_center_pixel.to_string(style='hmsdms', precision=1)))

            vizier_data = viz.query_region(coordinates=self.coordinate_of_center_pixel,\
                                           width=dRA*u.deg, height=dDEC*u.deg,\
                                           catalog=catalog)
            n_stars = len(vizier_data[0])
            self.logger.info("  Retrieved {} lines from {} catalog.".format(n_stars, catalog))
            self.catalog_name = catalog
            self.catalog_data = vizier_data[0]

            ## Standardize Column Names
            if catalog == 'USNO-B1.0':
                self.catalog_data.rename_column('USNO-B1.0', 'ID')
                self.catalog_data.rename_column('_RAJ2000', 'RA')
                self.catalog_data.rename_column('_DEJ2000', 'Dec')
            if catalog == 'UCAC4':
                self.catalog_data.rename_column('UCAC4', 'ID')
                self.catalog_data.rename_column('_RAJ2000', 'RA')
                self.catalog_data.rename_column('_DEJ2000', 'Dec')
        else:
            self.logger.info("No image WCS, so catalog query skipped")
            self.catalog_name = None
            self.catalog_data = None

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done retrieving Vizier catalog in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Get UCAC4 Catalog for Image from Local File
    ##-------------------------------------------------------------------------
    def get_local_UCAC4(self,\
                      local_UCAC_command="/Volumes/Data/UCAC4/access/u4test",\
                      local_UCAC_data="/Volumes/Data/UCAC4/u4b"):
        '''
        Get a list of stars which are in the image from a local UCAC catalog.
        '''
        start_time = datetime.datetime.now()
        assert type(self.coordinate_of_center_pixel) == coords.SkyCoord

        if not os.path.exists(local_UCAC_command):
            self.logger.warning('Cannot find local UCAC command: {}'.format(\
                                                           local_UCAC_command))
        elif not os.path.exists(local_UCAC_data):
            self.logger.warning('Cannot find local UCAC data: {}'.format(local_UCAC_data))
        else:
            center_from_WCS = self.image_WCS.wcs_pix2world([[self.nXPix/2, self.nYPix/2]], 1)
            footprint = self.image_WCS.calc_footprint()
            RAs = [val[0] for val in footprint]
            DECs = [val[1] for val in footprint]
            dRA = (max(RAs) - min(RAs))
            if dRA > 180:
                dRA = (min(RAs)+360. - max(RAs))
            dDEC = (max(DECs) - min(DECs))
#             corners = self.image_WCS.wcs_pix2world([[0, 0], [self.nXPix, 0],\
#                                  [0, self.nYPix], [self.nXPix, self.nYPix]], 1)
#             field_size_RA = (max(corners[:,0]) - min(corners[:,0]))*math.cos(self.coordinate_of_center_pixel.dec.degree*u.deg.to(u.radian))
#             field_size_DEC = max(corners[:,1]) - min(corners[:,1])
            self.logger.info("Getting stars from local UCAC4 catalog.")
            UCACcommand = [local_UCAC_command,\
                           "{:.4f}".format(self.coordinate_of_center_pixel.ra.degree),\
                           "{:.4f}".format(self.coordinate_of_center_pixel.dec.degree),\
                           "{:.2f}".format(dRA),\
                           "{:.2f}".format(dDEC),\
                           local_UCAC_data]
            self.logger.debug("  Using command: {}".format(' '.join(UCACcommand)))
            if os.path.exists("ucac4.txt"): os.remove("ucac4.txt")
            result = subprocess.call(UCACcommand)
            if os.path.exists('ucac4.txt'):
                catalog_file_path = os.path.join(self.tel.temp_file_path,\
                                                      'ucac4.txt')
                shutil.move('ucac4.txt', catalog_file_path)
                self.temp_files.append(catalog_file_path)
            else:
                self.logger.warning('  No ucac4.txt output file found.  Trying again.')
                result = subprocess.call(UCACcommand)
                if os.path.exists('ucac4.txt'):
                    catalog_file_path = os.path.join(self.tel.temp_file_path,\
                                                          'ucac4.txt')
                    shutil.move('ucac4.txt', catalog_file_path)
                    self.temp_files.append(catalog_file_path)
                else:
                    self.logger.warning('  No ucac4.txt output file found.  Trying again.')
                    return False

            ## Read in UCAC catalog
            colnames = ('id', 'RA', 'Dec', 'mag1', 'mag2', 'smag', 'ot', 'dsf',\
                        'RAepoch', 'Decepoch', 'dRA', 'dde', 'nt', 'nu', 'nc',\
                        'pmRA', 'pmDec', 'sRA', 'sDec', '2mass', 'j', 'h', 'k',\
                        'e2mphos', 'icq_flag', 'B', 'V', 'g', 'r', 'i')
            colstarts = (0, 10, 24, 36, 43, 50, 54, 57, 60, 68, 76, 80, 84, 87,\
                         90, 93, 100, 107, 111, 115, 126, 133, 140, 147, 159,\
                         168, 175, 182, 189, 196)
            colends =   (9, 22, 35, 42, 49, 53, 56, 59, 67, 75, 79, 83, 86, 89,\
                         92, 99, 106, 110, 114, 125, 132, 139, 146, 158, 167,\
                         174, 181, 188, 195, 202)
            self.catalog_name = 'UCAC4'
            self.catalog_data = ascii.read(catalog_file_path,\
                                           Reader=ascii.FixedWidthNoHeader,\
                                           data_start=1, guess=False,\
                                           names=colnames,\
                                           col_starts=colstarts,\
                                           col_ends=colends,\
                                          )
            ## Standardize Column Names
            self.catalog_data.remove_columns(['smag', 'ot', 'dsf',\
                                              'mag1', 'mag2', 'dRA', 'dde',\
                                              'nt', 'nu', 'nc',\
                                              'pmRA', 'pmDec', 'sRA', 'sDec',\
                                              '2mass', 'RAepoch', 'Decepoch',\
                                              'e2mphos', 'icq_flag'])
            self.catalog_data.rename_column('id', 'ID')
            self.catalog_data.rename_column('B', 'Bmag')
            self.catalog_data.rename_column('V', 'Vmag')
            self.catalog_data.rename_column('g', 'gmag')
            self.catalog_data.rename_column('r', 'rmag')
            self.catalog_data.rename_column('i', 'imag')

            ## Filter catalog for magnitude limits
            faint_stars_to_remove = []
            for i in range(0,len(self.catalog_data)):
                entry = self.catalog_data[i]
                if entry[self.tel.catalog_info[self.filter]] > self.tel.catalog_info['magmax']:
                    faint_stars_to_remove.append(i)
            if len(faint_stars_to_remove) > 0:
                self.logger.info('  Removing {} faint stars ({} > {}) from catalog.'.format(\
                                  len(faint_stars_to_remove),\
                                  self.tel.catalog_info[self.filter],\
                                  self.tel.catalog_info['magmax'],\
                                  ))
                self.catalog_data.remove_rows(faint_stars_to_remove)

            nUCACStars = len(self.catalog_data)
            self.logger.info("  Retrieved {} stars from UCAC catalog.".format(nUCACStars))

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done retrieving local UCAC4 catalog in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Measure Zero Point
    ##-------------------------------------------------------------------------
    def measure_zero_point(self, plot=False):
        '''
        Estimate the zero point of the image by comparing the instrumental
        magnitudes as determined by SExtractor to the catalog magnitues.
        '''
        start_time = datetime.datetime.now()
        self.logger.info('Analyzing SExtractor results to determine photometric zero point')

        if self.SExtractor_results and ('assoc_catmag' in self.SExtractor_results.keys()) and ('MAG_AUTO' in self.SExtractor_results.keys()):
            min_stars = 10

            zero_points = [entry['assoc_catmag'] - entry['MAG_AUTO']\
                           for entry in self.SExtractor_results\
                           if (entry['FLAGS'] == 0)\
                           and not np.isnan(entry['assoc_catmag'])\
                           and not (float(entry['assoc_catmag']) == 0.0)]
            SNR = [entry['FLUX_AUTO'] / entry['FLUXERR_AUTO']\
                   for entry in self.SExtractor_results\
                   if (entry['FLAGS'] == 0)\
                   and not np.isnan(entry['assoc_catmag'])\
                   and not (float(entry['assoc_catmag']) == 0.0)]

            self.logger.info('  Using {} stars with {} catalog magnitude'.format(len(zero_points), self.catalog_filter))

    #         self.zero_point_correlation = np.corrcoef(zip(list(entry['assoc_catmag']), list(entry['MAG_AUTO'])))
    #         print(self.zero_point_correlation)

            if len(zero_points) < min_stars:
                self.logger.warning('  Zero point not calculated.  Only {} catalog stars found.'.format(\
                                 len(zero_points)))
            else:
                self.zero_point_mode = mode(zero_points, 0.1)
                self.zero_point_median = np.median(zero_points)
                self.zero_point_average = np.average(zero_points, weights=SNR)
                self.logger.info('  Mode Zero Point = {:.2f}'.format(self.zero_point_mode))
                self.logger.info('  Median Zero Point = {:.2f}'.format(self.zero_point_median))
                self.logger.info('  Weighted Average Zero Point = {:.2f}'.format(self.zero_point_average))
                self.zero_point = self.zero_point_average

                ## Check zero point
                if self.tel.threshold_zeropoint and self.zero_point:
                    if self.zero_point < self.tel.threshold_zeropoint:
                        self.flags['zero point'] = True
                    else:
                        self.flags['zero point'] = False
                else:
                    self.flags['zero point'] = False

                end_time = datetime.datetime.now()
                elapzed_time = end_time - start_time
                self.logger.info('  Done measuring zero point in {:.1f} s'.format(elapzed_time.total_seconds()))

                ## Make Plot if Requested
                if plot:
                    self.make_zero_point_plot()


    ##-------------------------------------------------------------------------
    ## Make Zero Point Plot
    ##-------------------------------------------------------------------------
    def make_zero_point_plot(self):
        start_time = datetime.datetime.now()
        self.logger.info('Making ZeroPoint Plot')
        self.zero_point_plotfilename = self.raw_file_basename+'_ZeroPoint.png'
        self.zero_point_plotfile = os.path.join(self.tel.plot_file_path,\
                                                self.zero_point_plotfilename)

        catalog_mags = [entry['assoc_catmag']\
                        for entry in self.SExtractor_results\
                        if (entry['FLAGS'] == 0)\
                        and not np.isnan(entry['assoc_catmag'])
                        and not (float(entry['assoc_catmag']) == 0.0)]
        instrumental_mags = [entry['MAG_AUTO']\
                             for entry in self.SExtractor_results\
                             if (entry['FLAGS'] == 0)\
                             and not np.isnan(entry['assoc_catmag'])
                             and not (float(entry['assoc_catmag']) == 0.0)]
        zero_points = [entry['assoc_catmag'] - entry['MAG_AUTO']\
                       for entry in self.SExtractor_results\
                       if (entry['FLAGS'] == 0)\
                       and not np.isnan(entry['assoc_catmag'])
                       and not (float(entry['assoc_catmag']) == 0.0)]
        xpix = [entry['XWIN_IMAGE']\
                for entry in self.SExtractor_results\
                if (entry['FLAGS'] == 0)\
                and not np.isnan(entry['assoc_catmag'])
                and not (float(entry['assoc_catmag']) == 0.0)]
        ypix = [entry['YWIN_IMAGE']\
                for entry in self.SExtractor_results
                if (entry['FLAGS'] == 0)\
                and not np.isnan(entry['assoc_catmag'])
                and not (float(entry['assoc_catmag']) == 0.0)]
        residuals = [entry['assoc_catmag'] - entry['MAG_AUTO'] - self.zero_point\
                     for entry in self.SExtractor_results
                     if (entry['FLAGS'] == 0)\
                     and not np.isnan(entry['assoc_catmag'])
                     and not (float(entry['assoc_catmag']) == 0.0)]

        zp_binsize = 0.1
        bmin = math.floor(min(zero_points)/zp_binsize)*zp_binsize - zp_binsize/2.
        bmax = math.ceil(max(zero_points)/zp_binsize)*zp_binsize + zp_binsize/2.
        zp_bins = np.arange(bmin,bmax,zp_binsize)
        zp_hist, zp_bins = np.histogram(zero_points, bins=zp_bins)
        zp_centers = (zp_bins[:-1] + zp_bins[1:]) / 2

        pyplot.ioff()
        fig = pyplot.figure(figsize=(10,11), dpi=100)

        reject_percent = 3.0
        padding = 0.5
        hist2d_binsize = 0.1

        ## Correlation of Instrumental Magnitude with Catalog Magnitude
        TopLeft = pyplot.axes([0.000, 0.650, 0.465, 0.335])
        TopLeft.set_aspect('equal')
        xmin = math.floor( ( (np.percentile(catalog_mags, reject_percent)-padding)*2))/2.
        xmax = math.ceil( ( (np.percentile(catalog_mags, 100.-reject_percent)+padding)*2))/2.
        ymin = math.floor( ( (np.percentile(instrumental_mags, reject_percent)-padding)*2))/2.
        ymax = math.ceil( ( (np.percentile(instrumental_mags, 100.-reject_percent)+padding)*2))/2.
        xbins = list(np.arange(xmin,xmax,hist2d_binsize))
        ybins = list(np.arange(ymin,ymax,hist2d_binsize))
        pyplot.title('Correlation of Instrumental and Calalog Magnitudes', size=10)
        pyplot.hist2d(catalog_mags, instrumental_mags, bins=[xbins, ybins], cmap='binary')
        pyplot.xlabel('{} {} Magnitude'.format(self.catalog_name, self.catalog_filter), size=10)
        pyplot.ylabel('Instrumental Magnitude', size=10)
        pyplot.xticks(np.arange(-5,25,1))
        pyplot.yticks(np.arange(-25,25,1))
        pyplot.grid()
        pyplot.ylim(ymin,ymax)
        pyplot.xlim(xmin,xmax)
        ## Overplot Line of Zero Point
        catmag = [-5,30]
        fitmag = [(val-self.zero_point) for val in catmag]
        pyplot.plot(catmag, fitmag, 'k-')


        ## Plot Histogram of Zero Point Values
        TopRight = pyplot.axes([0.535, 0.650, 0.465, 0.335])
        pyplot.title('Histogram of Zero Point Values for {}'.format(self.raw_file_name), size=10)
        pyplot.plot([self.zero_point_mode, self.zero_point_mode], [0, 1.1*max(zp_hist)],\
                    'ro-', linewidth=2, label='Mode Zero Point', alpha=0.4)
        pyplot.plot([self.zero_point_median, self.zero_point_median], [0, 1.1*max(zp_hist)],\
                    'go-', linewidth=2, label='Median Zero Point', alpha=0.4)
        pyplot.plot([self.zero_point_average, self.zero_point_average], [0, 1.1*max(zp_hist)],\
                    'bo-', linewidth=2, label='Mode Zero Point', alpha=1.0)
        pyplot.bar(zp_centers, zp_hist, align='center', width=0.7*zp_binsize)
        pyplot.xlabel('Zero Point', size=10)
        pyplot.ylabel('N Stars', size=10)
        pyplot.xlim(np.percentile(zero_points, reject_percent)-padding,\
                    np.percentile(zero_points, 100.-reject_percent)+padding)
        pyplot.yticks(size=10)

        ## Plot Residuals
        MiddleLeft = pyplot.axes([0.000, 0.275, 0.465, 0.320])
        xmin = math.floor( ( (np.percentile(catalog_mags, reject_percent)-padding)*4))/4.
        xmax = math.ceil( ( (np.percentile(catalog_mags, 100.-reject_percent)+padding)*4))/4.
        ymin = math.floor( ( (np.percentile(residuals, reject_percent)-padding)*4))/4.
        ymax = math.ceil( ( (np.percentile(residuals, 100.-reject_percent)+padding)*4))/4.
        xbins = list(np.arange(xmin,xmax,hist2d_binsize))
        ybins = list(np.arange(ymin,ymax,hist2d_binsize))
        pyplot.hist2d(catalog_mags, residuals, bins=[xbins, ybins], cmap='binary')
        pyplot.xlabel('{} {} Magnitude'.format(self.catalog_name, self.catalog_filter), size=10)
        pyplot.ylabel('Magnitude Residuals', size=10)
        pyplot.grid()
        pyplot.ylim(ymin,ymax)
        pyplot.xlim(xmin,xmax)
        ## Overplot Line of Zero Point
        catmag = [-5,30]
        fitmag = [0, 0]
        pyplot.plot(catmag, fitmag, 'k-')

        ## Plot Spatial Distribution of Residuals
        range = [-0.5, 0.5]
        MiddleRight = pyplot.axes([0.535, 0.275, 0.465, 0.320])
        MiddleRight.set_aspect('equal')
        pyplot.title('Residuals scaled from {:+.1f} to {:+.1f}'.format(range[0], range[1]), size=10)
        if len(residuals) > 20000:
            gridsize = 20
        else:
            gridsize = 10
        pyplot.hexbin(xpix, ypix, residuals,\
                      gridsize=gridsize,\
                      mincnt=5,\
                      vmin=range[0], vmax=range[1],\
                      alpha=0.5,\
                      cmap='Reds')
        pyplot.xlabel('X Pixels', size=10)
        pyplot.ylabel('Y Pixels', size=10)
        pyplot.xlim(0,self.nXPix)
        pyplot.ylim(0,self.nYPix)
        pyplot.xticks(size=10)
        pyplot.yticks(size=10)

        pyplot.savefig(self.zero_point_plotfile, dpi=100,\
                       bbox_inches='tight', pad_inches=0.10)
        pyplot.close(fig)

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done making zero point plot in {:.1f} s'.format(elapzed_time.total_seconds()))



    ##-------------------------------------------------------------------------
    ## Make JPEG of Image (using matplotlib)
    ##-------------------------------------------------------------------------
    def make_JPEG(self, jpeg_file_name, binning=1, p1=0.15, p2=0.5,\
                      mark_pointing=False,\
                      mark_detected_stars=False,\
                      mark_catalog_stars=False,\
                      mark_saturated=False,\
                      make_hist=False,\
                      transform=None,
                      crop=None,
                      quality=70,
                     ):
        '''
        Make jpegs of image.
        '''
        start_time = datetime.datetime.now()
        self.logger.info('Making jpeg: {}'.format(jpeg_file_name))
        jpeg_file = os.path.join(self.tel.plot_file_path, jpeg_file_name)

        from PIL import Image, ImageDraw
        import skimage.exposure as skiex

        self.logger.debug('  Opening working file')
        with fits.open(self.working_file, ignore_missing_end=True) as hdulist:
            data = hdulist[0].data
        data_zero = np.ma.masked_equal(data, 0)
        data_nonzero = data_zero[~data_zero.mask]

        ## Make exposure histogram (of unscaled data)
        if make_hist:
            self.logger.info('  Make histogram of unscaled data.')
            histogram_plot_file = os.path.join(self.tel.plot_file_path, '{}_hist.png'.format(self.raw_file_basename))
            hist_low = np.percentile(data_nonzero.ravel(), p1)
            hist_high = np.percentile(data_nonzero.ravel(), 100.-p2)
            hist_nbins = 128
            hist_binsize = (hist_high-hist_low)/128
            hist_bins = np.arange(hist_low,hist_high,hist_binsize)
            self.logger.debug('  Histogram range: {} {}.'.format(hist_low, hist_high))
            pyplot.ioff()
            fig = pyplot.figure()
            pyplot.hist(data.ravel(), bins=hist_bins, label='binsize = {:4f}'.format(hist_binsize))
            pyplot.xlim(hist_low,hist_high)
            pyplot.legend(loc='best')
            pyplot.xlabel('Pixel value')
            pyplot.ylabel('Number of Pixels')
            self.logger.info('  Saving histogram to {}.'.format(histogram_plot_file))
            pyplot.savefig(histogram_plot_file)
            pyplot.close(fig)

        ## Rescale data using arcsinh transform for jpeg
        self.logger.debug('  Rescaling image data using arcsinh')
        rescaled_data = np.arcsinh(data_zero)
        rescaled_data = rescaled_data / rescaled_data.max()
        rescaled_data_nonzero = rescaled_data[~rescaled_data.mask]
        low = np.percentile(rescaled_data_nonzero, p1)
        high = np.percentile(rescaled_data_nonzero, 100.-p2)
        self.logger.debug('  Clipping data using {} and {} percentiles.'.format(p1, 100.-p2))
        self.logger.debug('  Clipping data using {} and {} rescaled values.'.format(low, high))
        opt_img = skiex.exposure.rescale_intensity(rescaled_data, in_range=(low,high))
        jpegdata = (opt_img * 255.).astype('uint8')

        ## Create PIL Image object
        im = Image.fromarray(jpegdata).convert('RGB')
        draw = ImageDraw.Draw(im)

        ## If mark_pointing is set
        if mark_pointing and self.coordinate_from_header and self.image_WCS:
            xy = self.image_WCS.wcs_world2pix([[self.coordinate_from_header.ra.degree,\
                                                self.coordinate_from_header.dec.degree]], 1)[0]
            x = int(xy[0])
            y = int(xy[1])
            self.logger.debug('  Marking crosshairs at (x, y) = ({}, {})'.format(im.size[0]/2, im.size[1]/2))
            line_color = 'yellow'
            draw.line((im.size[0]/2+0, 0, im.size[0]/2+0, im.size[1]), fill=line_color)
            draw.line((0, im.size[1]/2+0, im.size[0], im.size[1]/2+0), fill=line_color)

            ## Draw Crosshair Over Pointing Location from Header
            self.logger.debug('  Marking pointing at (x, y) = ({}, {})'.format(x, y))
            crosshair_color = 'cyan'
            ms = int((self.tel.pointing_marker_size.to(u.arcsec) / self.tel.pixel_scale).to(u.pix).value)/2
            self.logger.debug('  Pointing marker diameter is {} = {} pix'.format(self.tel.pointing_marker_size.to(u.arcmin), ms*2))
            thickness = 5
            for i in range(-1*int((thickness-1)/2),int((thickness+1)/2),1):
                draw.line((x-1.5*ms, y+i, x-0.5*ms, y+i), fill=crosshair_color)
                draw.line((x+1.5*ms, y+i, x+0.5*ms, y+i), fill=crosshair_color)
                draw.line((x+i, y-1.5*ms, x+i, y-0.5*ms), fill=crosshair_color)
                draw.line((x+i, y+1.5*ms, x+i, y+0.5*ms), fill=crosshair_color)
            radii = np.linspace(ms, ms+thickness, thickness)
            for r in radii:
                draw.ellipse((x-r, y-r, x+r, y+r), outline=crosshair_color)

        ## Mark Catalog Stars
        if mark_catalog_stars and self.catalog_data:
            if self.FWHM:
                ms = max([7, 2.1*math.ceil(self.FWHM.to(u.pix).value)])
            else:
                ms = 7
            circle_color = 'red'
            self.logger.debug('  Marking catalog stars with {} radius {} circles'.format(ms, circle_color))

            for star in self.catalog_data:
                xy = self.image_WCS.wcs_world2pix([[float(star['RA']), float(star['Dec'])]], 1)[0]
                x = int(xy[0])
                y = int(xy[1])
                thickness = 2
                radii = np.linspace(ms, ms+thickness, thickness)
                for r in radii:
                    draw.ellipse((x-r, y-r, x+r, y+r), outline=circle_color)

        ## Mark Detected Stars
        if mark_detected_stars and self.SExtractor_results:
            if self.FWHM:
                ms = max([6, 2*math.ceil(self.FWHM.to(u.pix).value)])
            else:
                ms = 6
            circle_color = 'green'
            self.logger.debug('  Marking detected stars with {} radius {} circles'.format(ms, circle_color))
            for star in self.SExtractor_results:
                x = star['XWIN_IMAGE']
                y = star['YWIN_IMAGE']
                thickness = 2
                radii = np.linspace(ms, ms+thickness, thickness)
                for r in radii:
                    draw.ellipse((x-r, y-r, x+r, y+r), outline=circle_color)

        ## Flag Saturated Pixels
        if mark_saturated and self.tel.saturation:
            saturated_color = 'red'
            with fits.open(self.working_file, ignore_missing_end=True) as hdulist:
                data_raw = hdulist[0].data
            data_saturated = np.ma.masked_greater(data_raw, self.tel.saturation)
            indices = np.where(data_saturated.mask == 1)
            if len(indices) > 4:
                xy = zip(indices[1], indices[0])
                draw.point(xy, fill=saturated_color)

        ## Flip jpeg
        if transform:
            self.logger.debug('  Transforming (flipping/rotating) jpeg: {}'.format(transform))
            if transform == 'flip_vertical':
                im = im.transpose(Image.FLIP_TOP_BOTTOM)
            elif transform == 'flip_horizontal':
                im = im.transpose(Image.FLIP_LEFT_RIGHT)
            elif transform == 'rotate90':
                im = im.transpose(Image.ROTATE_90)
            elif transform == 'rotate180':
                im = im.transpose(Image.ROTATE_180)
            elif transform == 'rotate270':
                im = im.transpose(Image.ROTATE_270)
            else:
                self.logger.warning('  Transform "{}" not understood.'.format(transform))
                self.logger.warning('  No transform performed.'.format(transform))

        ## If crop is set
        if crop:
            if len(crop) == 4:
                self.logger.debug('  Cropping image to {}'.format(crop))
                im = im.crop(crop)

        ## If binning is set create thumbnail
        if binning > 1:
            size = (int(data.shape[0]/binning), int(data.shape[1]/binning))
            self.logger.debug('  Resizing image by binning factor of {}'.format(binning))
            im.thumbnail(size, Image.ANTIALIAS)

        ## Save to JPEG
        self.logger.info('  Saving jpeg (p1={:.1f}, p2={:.1f}), binning={}, quality={:.0f}) to: {}'.format(p1, p2, binning, quality, jpeg_file_name))
        im.save(jpeg_file, 'JPEG', quality=quality)
        self.jpeg_file_names.append(jpeg_file_name)

        end_time = datetime.datetime.now()
        elapzed_time = end_time - start_time
        self.logger.info('  Done making JPEG in {:.1f} s'.format(elapzed_time.total_seconds()))


    ##-------------------------------------------------------------------------
    ## Clean Up by Deleting Temporary Files
    ##-------------------------------------------------------------------------
    def clean_up(self):
        '''
        Clean up by deleting temporary files.
        '''
        self.logger.info("Cleaning Up Temporary Files.")
        for item in self.temp_files:
            if os.path.exists(item):
                self.logger.debug("  Deleting {0}".format(item))
                os.remove(item)


    ##-------------------------------------------------------------------------
    ## Append Line With Image Info to HTML File List
    ##-------------------------------------------------------------------------
    def add_web_log_entry(self, htmlImageList, fields=None):
        '''
        This function adds one line to the HTML table of images.  The line
        contains the image info extracted by IQMon.
        '''
        if not fields: fields=["Date and Time", "Filename", "Alt", "Az",\
                               "Airmass", "MoonSep", "MoonIllum", "FWHM",\
                               "ellipticity", "Background", "PErr", "PosAng",\
                               "ZeroPoint", "nStars", "ProcessTime"]
        ## If HTML file does not yet exist, create it and insert header
        ## from template file.
        self.logger.info('Adding results to HTML table: {}'.format(htmlImageList))
        if not os.path.exists(htmlImageList):
            self.logger.debug("  HTML file does not exist.  Creating it.")
            HTML = open(htmlImageList, 'w')
            header = ['<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">',
                      '<html lang="en">',
                      '<head>',
                      '    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">',
                      '    <title>IQMon Results</title>',
                      '    <style>',
                      '        table',
                      '        {',
                      '            border-collapse:collapse;',
                      '        }',
                      '        table,th,td',
                      '        {',
                      '            border:1px solid black;',
                      '            text-align:center;',
                      '        }',
                      '    </style>',
                      '</head>',
                      '<body>',
                      '    <h2>IQMon Results for {}</h2>'.format(self.tel.name),
                      '    <table>',
                      '        <tr>']
            if "Date and Time" in fields:
                header.append('        <th style="width:150px">Exposure Start<br>(Date and Time UT)</th>')
            if "Filename" in fields:
                header.append('        <th style="width:420px">Image File Name</th>')
            if "Target" in fields:
                header.append('        <th style="width:120px">Target Name</th>')
            if "ExpTime" in fields:
                header.append('        <th style="width:50px">Exp Time (s)</th>')
            if "Alt" in fields:
                header.append('        <th style="width:50px">Alt (deg)</th>')
            if "Az" in fields:
                header.append('        <th style="width:50px">Az (deg)</th>')
            if "Airmass" in fields:
                header.append('        <th style="width:50px">Airmass</th>')
            if "MoonSep" in fields:
                header.append('        <th style="width:50px">Moon Sep (deg)</th>')
            if "MoonIllum" in fields:
                header.append('        <th style="width:50px">Moon Illum. (%)</th>')
            if "FWHM" in fields:
                header.append('        <th style="width:60px">FWHM ({})</th>'.format(str(self.tel.units_for_FWHM)))
            if "ellipticity" in fields:
                header.append('        <th style="width:50px">Ellip.</th>')
            if "Background" in fields:
                header.append('        <th style="width:70px">Background<br>[RMS]</th>')
            if "PErr" in fields:
                header.append('        <th style="width:70px">Pointing Error (arcmin)</th>')
            if "PosAng" in fields:
                header.append('        <th style="width:50px">WCS Pos. Angle</th>')
            if "ZeroPoint" in fields:
                header.append('        <th style="width:50px">Zero Point (mag)</th>')
            if "nStars" in fields:
                header.append('        <th style="width:50px">N Stars</th>')
            if "ProcessTime" in fields:
                header.append('        <th style="width:50px">Process Time (sec)</th>')
            header.append('        </tr>')
            header.append('    </body>')
            header.append('</html>')
            for headerline in header:
                HTML.write(headerline)
            HTML.close()
        ## If HTML file does exist, we need to strip off the lines which
        ## end the file, so we can append more data to the table.
        else:
            self.logger.debug("HTML file exists.  Copying contents.")
            HTML = open(htmlImageList, 'r')
            existingContent = HTML.read().split("\n")
            HTML.close()
            HTML = open(htmlImageList, 'w')
            for line in existingContent:
                IsEndTable = re.match("\s*</table><br><br>\s*", line)
                IsEndBody = re.match("\s*</body>\s*", line)
                IsEndHTML = re.match("\s*</html>\s*", line)
                if not IsEndTable and not IsEndBody and not IsEndHTML:
                    HTML.write(line+"\n")
        ## Write Lines for this Image to HTML File
        HTML = open(htmlImageList, 'a')
        HTML.write("    <tr>\n")
        ## Write Observation Date and Time
        if "Date and Time" in fields:
            HTML.write("      <td style='color:black;text-align:left'>{0}</td>\n".format(self.observation_date))
        ## Write Filename (and links to jpegs)
        if "Filename" in fields:
            if len(self.jpeg_file_names) == 0:
                JPEG1_html = ""
                JPEG2_html = ""
                JPEG3_html = ""
            elif len(self.jpeg_file_names) == 1:
                JPEG1_html = "<a href='{}'>".format(os.path.join("..", "..", "Plots", self.jpeg_file_names[0]))
                JPEG2_html = ""
                JPEG3_html = ""
            elif len(self.jpeg_file_names) == 2:
                JPEG1_html = "<a href='{}'>".format(os.path.join("..", "..", "Plots", self.jpeg_file_names[0]))
                JPEG2_html = " (<a href='{}'>JPEG2</a>)".format(os.path.join("..", "..", "Plots", self.jpeg_file_names[1]))
                JPEG3_html = ""
            elif len(self.jpeg_file_names) >= 3:
                JPEG1_html = "<a href='{}'>".format(os.path.join("..", "..", "Plots", self.jpeg_file_names[0]))
                JPEG2_html = " (<a href='{}'>JPEG2</a>)".format(os.path.join("..", "..", "Plots", self.jpeg_file_names[1]))
                JPEG3_html = " (<a href='{}'>JPEG3</a>)".format(os.path.join("..", "..", "Plots", self.jpeg_file_names[2]))
            if self.PSF_plot_file:
                PSFplot_html = " (<a href='{}'>PSF</a>)".format(os.path.join("..", "..", "Plots", self.PSF_plot_filename))
            else:
                PSFplot_html = ""
            if self.zero_point_plotfile:
                ZPplot_html = " (<a href='{}'>ZP</a>)".format(os.path.join("..", "..", "Plots", self.zero_point_plotfilename))
            else:
                ZPplot_html = ""
            log_html = " (<a href='{}'>log</a>)".format(os.path.join("..", "..", "Logs", self.tel.config['name'], self.logfilename))
            htmlline = "      <td style='color:black;text-align:left'>" + JPEG1_html + "{}</a>".format(self.raw_file_basename) + JPEG2_html + JPEG3_html + PSFplot_html + ZPplot_html + log_html + "</td>\n"
            HTML.write(htmlline)
        ## Write Target Name
        if "Target" in fields:
            if self.object_name:
                HTML.write("      <td style='color:black'>{0:}</td>\n".format(self.object_name))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        ## Write Exposure Time
        if "ExpTime" in fields:
            if self.exptime:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.exptime.to(u.s).value))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        ## Write Alt, Az, airmass, moon separation, and moon phase
        if "Alt" in fields:
            if self.target_alt:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.target_alt.to(u.deg).value))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        if "Az" in fields:
            if self.target_az:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.target_az.to(u.deg).value))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        if "Airmass" in fields:
            if self.target_az:
                HTML.write("      <td style='color:{0}'>{1:.2f}</td>\n".format("black", self.airmass))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        if "MoonSep" in fields:
            if self.moon_sep:
                if self.moon_alt > 0:
                    HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.moon_sep.to(u.deg).value))
                else:
                    HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", "down"))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        if "MoonIllum" in fields:
            if self.moon_phase:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.moon_phase))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        ## Write FWHM and ellipticity
        if "FWHM" in fields:
            if self.FWHM:
                ## Decide whether to flag FWHM value with red color
                if self.flags['FWHM']:
                    colorFWHM = "#FF5C33"
                else:
                    colorFWHM = "#70DB70"
                ## Convert FWHM value to appropriate units for HTML output
                if self.tel.units_for_FWHM == u.arcsec:
                    FWHM_for_HTML = (self.FWHM * u.radian.to(u.arcsec)*self.tel.pixel_size.to(u.mm)/self.tel.focal_length.to(u.mm)).value
                else:
                    FWHM_for_HTML = self.FWHM.value
                HTML.write("      <td style='background-color:{0}'>{1:.2f}</td>\n".format(colorFWHM, FWHM_for_HTML))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("#FF5C33", ""))
        if "ellipticity" in fields:
            if self.ellipticity:
                ## Decide whether to flag ellipticity value with red color
                if self.flags['ellipticity']:
                    colorEllipticity = "#FF5C33"
                else:
                    colorEllipticity = "#70DB70"
                HTML.write("      <td style='background-color:{0}'>{1:.2f}</td>\n".format(colorEllipticity, self.ellipticity))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("#FF5C33", ""))
        ## Write SExtractor background and background RMS
        if "Background" in fields:
            if self.SExtractor_background and self.SExtractor_background_RMS:
                HTML.write("      <td style='color:{0}'>{1:.0f} [{2:.0f}]</td>\n".format("black", self.SExtractor_background, self.SExtractor_background_RMS))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Write pointing error
        if "PErr" in fields:
            if self.pointing_error:
                ## Decide whether to flag pointing error value with red color
                if self.flags['pointing error']:
                    colorpointing_error = "#FF5C33"
                else:
                    colorpointing_error = "#70DB70"
                ## Write HTML
                HTML.write("      <td style='background-color:{0}'>{1:.1f}</td>\n".format(colorpointing_error, self.pointing_error.arcminute))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("#FF5C33", ""))
        ## Write WCS position angle
        if "PosAng" in fields:
            if self.position_angle:
                HTML.write("      <td style='color:{}'>{:.1f}</td>\n".format("black", self.position_angle.to(u.deg).value))
            else:
                HTML.write("      <td style='color:{}'>{}</td>\n".format("black", ""))
        ## Write zero point
        if "ZeroPoint" in fields:
            if self.zero_point:
                ## Decide whether to flag pointing error value with red color
                if self.flags['zero point']:
                    colorzero_point = "#FF5C33"
                else:
                    colorzero_point = "#70DB70"
                HTML.write("      <td style='background-color:{}'>{:.2f}</td>\n".format(colorzero_point, self.zero_point))
            else:
                HTML.write("      <td style='background-color:{}'>{}</td>\n".format("white", ""))
        ## Write number of stars detected by SExtractor
        if "nStars" in fields:
            if self.n_stars_SExtracted:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", self.n_stars_SExtracted))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Write process time
        if "ProcessTime" in fields:
            if self.total_process_time:
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.total_process_time))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Complete Table
        HTML.write("    </tr>\n")
        HTML.write("  </table><br><br>\n")
        HTML.write("</body>\n")
        HTML.write("</html>\n")
        HTML.close()


    ##-------------------------------------------------------------------------
    ## Append Line With Image Info to YAML Text File
    ##-------------------------------------------------------------------------
    def add_yaml_entry(self, summary_file):
        self.logger.info("Writing YAML Summary File: {}".format(summary_file))
        result_list = []
        if os.path.exists(summary_file):
            self.logger.debug('  Reading existing summary file.')
            with open(summary_file, 'r') as yaml_string:
                result_list = yaml.load(yaml_string)
        ## Form dictionary with new result info
        try:
            FWHM_median_pix = self.FWHM_median.to(u.pix).value
        except:
            FWHM_median_pix = None
        try:
            FWHM_mode_pix = self.FWHM_mode.to(u.pix).value
        except:
            FWHM_mode_pix = None
        try:
            FWHM_pix = self.FWHM.to(u.pix).value
        except:
            FWHM_pix = None
        try:
            pointing_error_arcmin = self.pointing_error.arcminute
        except:
            pointing_error_arcmin = None
        try:
            alt = self.target_alt.to(u.deg).value
        except:
            alt = None
        try:
            az = self.target_az.to(u.deg).value
        except:
            az = None
        try:
            moon_sep = self.moon_sep.to(u.deg).value
        except:
            moon_sep = None
        try:
            posang = self.position_angle.to(u.deg).value
        except:
            posang = None
        new_result = {
                      'filename': self.raw_file_name,\
                      'exposure_start': self.observation_date,\
                      'FWHM_median_pix': str(FWHM_median_pix),\
                      'FWHM_mode_pix': str(FWHM_mode_pix),\
                      'FWHM_pix': str(FWHM_pix),\
                      'ellipticity_median': str(self.ellipticity_median),\
                      'ellipticity_mode': str(self.ellipticity_mode),\
                      'ellipticity': str(self.ellipticity),\
                      'n_stars': str(self.n_stars_SExtracted),\
                      'background': str(self.SExtractor_background),\
                      'background_rms': str(self.SExtractor_background_RMS),\
                      'pointing_error_arcmin': str(pointing_error_arcmin),\
                      'zero_point': str(self.zero_point),\
                      'alt': str(alt),\
                      'az': str(az),\
                      'airmass': str(self.airmass),\
                      'moon_separation': str(moon_sep),\
                      'moon_illumination': str(self.moon_phase),\
                      'WCS_position_angle': str(posang),\
                      'process_time': str(self.total_process_time),\
                      'flags': str(self.flags),\
                      'IQMon Version': str(__version__),\
                     }
        result_list.append(new_result)
        yaml_string = yaml.dump(result_list)
        with open(summary_file, 'w') as output:
            output.write(yaml_string)


    ##-------------------------------------------------------------------------
    ## Append Line With Image Info to Summary Text File
    ##-------------------------------------------------------------------------
    def add_summary_entry(self, summaryFile):
        self.logger.info("Writing Summary File Entry.")
        self.logger.debug("  Summary File: {0}".format(summaryFile))
        ## Read in previous data
        if not os.path.exists(summaryFile):
            self.logger.debug("  Making new astropy table object")
            SummaryTable = table.Table(names=("ExpStart", "File", "FWHM (pix)", "Ellipticity",\
                                       "Alt (deg)", "Az (deg)", "Airmass", "pointing_error (arcmin)", \
                                       "ZeroPoint", "nStars", "Background", "Background RMS"),\
                                 dtype=('S22', 'S100', 'f4', 'f4', 'f4', 'f4',\
                                        'f4', 'f4', 'f4', 'i4', 'f4', 'f4'),\
                                 masked=True)
        else:
            self.logger.debug("  Reading astropy table object from file: {0}".format(\
                                                                  summaryFile))
            try:
                SummaryTable = ascii.read(summaryFile, guess=False,
                                          header_start=0, data_start=1,
                                          Reader=ascii.basic.Basic,
                                          delimiter="\s",
                                          fill_values=('--', '0'),
                                          converters={
                                          'ExpStart': [ascii.convert_numpy('S22')],
                                          'File': [ascii.convert_numpy('S100')],
                                          'FWHM (pix)': [ascii.convert_numpy('f4')],
                                          'Ellipticity': [ascii.convert_numpy('f4')],
                                          'Alt (deg)': [ascii.convert_numpy('f4')],
                                          'Az (deg)': [ascii.convert_numpy('f4')],
                                          'Airmass': [ascii.convert_numpy('f4')],
                                          'pointing_error (arcmin)': [ascii.convert_numpy('f4')],
                                          'ZeroPoint': [ascii.convert_numpy('f4')],
                                          'nStars': [ascii.convert_numpy('i4')],
                                          'Background': [ascii.convert_numpy('f4')],
                                          'Background RMS': [ascii.convert_numpy('f4')]
                                          })
            except:
                self.logger.critical("Failed to read summary file: {0} {1} {2}".format(\
                                     sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
        ## Astropy table writer can not write None to table initialized
        ## with type.  If any outputs are None, change to some value.
        tableMask = np.zeros(12)
        ## observation_date
        if self.observation_date: observation_date = self.observation_date
        else: 
            observation_date = ""
            tableMask[0] = True
        ## FileName
        if self.raw_file_name: raw_file_name = self.raw_file_name
        else: 
            raw_file_name = ""
            tableMask[1] = True
        ## FWHM
        if self.FWHM: FWHM = self.FWHM.to(u.pix).value
        else:
            FWHM = 0.
            tableMask[2] = True
        ## Ellipticity
        if self.ellipticity: ellipticity = self.ellipticity
        else:
            ellipticity = 0.
            tableMask[3] = True
        ## Target Alt
        if self.target_alt: target_alt = self.target_alt.to(u.deg).value
        else:
            target_alt = 0.
            tableMask[4] = True
        ## Target Az
        if self.target_az: target_az = self.target_az.to(u.deg).value
        else:
            target_az = 0.
            tableMask[5] = True
        ## Airmass
        if self.airmass: airmass = self.airmass
        else:
            airmass = 0.
            tableMask[6] = True
        ## Pointing Error
        if self.pointing_error: pointing_error = self.pointing_error.arcminute
        else:
            pointing_error = 0.
            tableMask[7] = True
        ## Zero Point
        if self.zero_point: zeroPoint = self.zero_point
        else:
            zeroPoint = 0.
            tableMask[8] = True
        ## n_stars_SExtracted
        if self.n_stars_SExtracted: n_stars_SExtracted = self.n_stars_SExtracted
        else: 
            n_stars_SExtracted = 0.
            tableMask[9] = True
        ## SExtractor Background
        if self.SExtractor_background: SExtractor_background = self.SExtractor_background
        else:
            SExtractor_background = 0.
            tableMask[10] = True
        ## SExtractor Background RMS
        if self.SExtractor_background_RMS: SExtractor_background_RMS = self.SExtractor_background_RMS
        else:
            SExtractor_background_RMS = 0.
            tableMask[11] = True
        ## Add row to table
        self.logger.debug("  Writing new row to log table.  Filename: {0}".format(raw_file_name))
        SummaryTable.add_row((observation_date, raw_file_name,
                              FWHM, ellipticity,
                              target_alt, target_az,
                              airmass, pointing_error,
                              zeroPoint, n_stars_SExtracted,
                              SExtractor_background, SExtractor_background_RMS),
                              mask=tableMask)
        ## Write Table to File
        self.logger.debug("  Writing new summary file.")
        ascii.write(SummaryTable, summaryFile,
                    Writer=ascii.basic.Basic)


    ##-------------------------------------------------------------------------
    ## Calcualte Process Time
    ##-------------------------------------------------------------------------
    def calculate_process_time(self):
        '''
        Determine how long it took for IQMon to process this image.  Determined
        by subtracting the starting time (determined on the initialization of 
        the image object) to the ending time (determined by this method).
        '''
        self.end_process_time = datetime.datetime.now()
        self.total_process_time = (self.end_process_time - self.start_process_time).total_seconds()
        self.logger.info("IQMon processing time = {0:.1f} seconds".format(self.total_process_time))

