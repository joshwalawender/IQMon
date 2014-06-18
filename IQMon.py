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
import shutil
import datetime
import subprocess
import logging
import math
import numpy as np
import matplotlib.pyplot as pyplot


## Import Astronomy Specific Tools
import ephem
import astropy.units as u
from astropy.io import fits
import astropy.coordinates as coords
from astropy import table
from astropy import wcs
from astropy.io import ascii



##-----------------------------------------------------------------------------
## Function to Determine Orientation from WCS
##-----------------------------------------------------------------------------
def image_orientation(WCS):
    assert type(WCS) == astropy.wcs.WCS
    header = WCS.to_header()
    PC = np.array([ [float(header['PC11']), float(header['PC12'])],\
                    [float(header['PC21']), float(header['PC22'])] ])

    ## Determine Pixel Scale
    result1 = PC.dot(np.array([[0], [1]]))
    pixel_scale1 = (math.sqrt(result1[0][0]**2 + result1[1][0]**2))*3600.
    result2 = PC.dot(np.array([[1], [0]]))
    pixel_scale2 = (math.sqrt(result2[0][0]**2 + result2[1][0]**2))*3600.
    pixel_scale = np.mean([pixel_scale1, pixel_scale2]) * u.arcsec/u.pix

    ## Determine Position Angle
    ang1 = math.acos(PC[0][0])
    ang2 = math.acos(PC[0][1])
    ang3 = math.acos(PC[1][0])
    ang4 = math.acos(PC[1][1])
    PA = (270 - np.mean([ang1, ang2, ang3, ang4])*180/math.pi) * u.deg

    ## Determine Flip State
    flipped = np.linalg.det(PC) > 0

    return pixel_scale, PA, flipped


##-----------------------------------------------------------------------------
## Define Config object to hold IQMon configuration information
##-----------------------------------------------------------------------------
class Config(object):
    '''
    Contains configuration information for IQMon package that can be passed to
    methods and functions.

    Properties:
    - pathIQMonExec: path to the IQMon folder where sextractor default files
                     are stored.  Typically this is the directory where IQMon
                     is installed.
    - pathLog:       path where IQMon's log files should be written
                     (i.e. ~/IQMon/Logs/)
    - pathPlots:     path where plots and jpegs should be written
                     (i.e. ~/IQMon/Plots/)
    - pathTemp:      path where temporary files should be written
                     (i.e. ~/IQMon/tmp/)
    - pathConfig:    path for config files.  Typically the base path for the
                     above paths
    '''
    _singletons = dict()

    def __new__(cls):
        if not cls in cls._singletons:
            cls._singletons[cls] = object.__new__(cls)
        return cls._singletons[cls]

    def __init__(self):
        '''
        Read and parse configuration file.
        - Assumes that file is .IQMonConfig in the user's home directory.
        - No defaults set, if config file not read, then values default to
          None.
        '''

        ## Look for configuration file
        HomePath = os.path.expandvars("$HOME")
        ConfigFilePath = os.path.join(HomePath, ".IQMonConfig")
        if os.path.exists(ConfigFilePath):
            ConfigFile = open(ConfigFilePath, 'r')
            ConfigFileLines = ConfigFile.readlines()
            ConfigFile.close()
        else:
            ConfigFileLines = None

        ## read configuration file
        for line in ConfigFileLines:
            IsIQMonExecPath = re.match("IQMONPATH\s=\s([\w/\-\.]+)", line)
            if IsIQMonExecPath:
                self.pathIQMonExec = os.path.abspath(IsIQMonExecPath.group(1))
            IsLogPath = re.match("IQMONLOGS\s=\s([\w/\-\.]+)", line)
            if IsLogPath:
                self.pathLog = os.path.abspath(IsLogPath.group(1))
            IsPlotsPath = re.match("IQMONPLOTS\s=\s([\w/\-\.]+)", line)
            if IsPlotsPath:
                self.pathPlots = os.path.abspath(IsPlotsPath.group(1))
            IstmpPath = re.match("IQMONTMP\s=\s([\w/\-\.]+)", line)
            if IstmpPath:
                self.pathTemp = os.path.abspath(IstmpPath.group(1))
            IsConfigPath = re.match("CONFIGPATH\s=\s([\w/\-\.]+)", line)
            if IsConfigPath:
                self.pathConfig = os.path.abspath(IsConfigPath.group(1))

        ## Create Log Path if it doesn't exist
        SplitPath = [self.pathLog]
        CreatePaths = []
        while not os.path.exists(SplitPath[0]):
            CreatePaths.append(SplitPath[0])
            SplitPath = os.path.split(SplitPath[0])
        while len(CreatePaths) > 0:
            os.mkdir(CreatePaths.pop())
        ## Create Plots Path if it doesn't exist
        SplitPath = [self.pathPlots]
        CreatePaths = []
        while not os.path.exists(SplitPath[0]):
            CreatePaths.append(SplitPath[0])
            SplitPath = os.path.split(SplitPath[0])
        while len(CreatePaths) > 0:
            os.mkdir(CreatePaths.pop())
        ## Create temp Path if it doesn't exist
        SplitPath = [self.pathTemp]
        CreatePaths = []
        while not os.path.exists(SplitPath[0]):
            CreatePaths.append(SplitPath[0])
            SplitPath = os.path.split(SplitPath[0])
        while len(CreatePaths) > 0:
            os.mkdir(CreatePaths.pop())



##-----------------------------------------------------------------------------
## Define Telescope object to hold telescope information
##-----------------------------------------------------------------------------
class Telescope(object):
    '''
    Contains information about the telescope that can be passed to methods and
    functions.  The concept for operation is that the user will write a simple
    script which creates a telescope object and assigned values to all it's
    properties (or sets them to None).  The processing is then just a sequence
    of calls to the IQMon.Image object and it's methods.

    Properties:
      name:          A short name for the telescope
      long_name:      A long name for the telescope
      focal_length:   The focal length of the telescope
      pixel_size:     The size of the pixels in the CCD camera
      pixel_scale:    The estimated pixel scale of the telescope in arcsec per
                     pixel.  This should be an astropy Unit object with a value
                     and units.
      fRatio:        The F/ratio of the telescope.
      aperture:      The aperture of the telescope.
      gain:          The gain (in electrons per ADU) of the CCD.
      nXPix:         The number of pixels in the X direction of the CCD.
      nYPix:         The number of pixels in the Y direction of the CCD.
      site:          A pyephem site object defining the site of the telescope
      units_for_FWHM:  A quantity with the units for the FWHM output.  The units
                     should be reducible to either pixels or arcseconds.
      ROI:           The "region of interest" to crop the image to.  Format
                     should be '[x1:x2,y1:y2]' or 'x1:x2,y1:y2'.
      threshold_FWHM: The FWHM above which the image should be flagged in the
                     HTML output log.
      threshold_pointing_err:   
      threshold_ellipticity:   
      SExtractorPhotAperture: 
      SExtractorSeeing:       The initial seeing estimate for SExtractor.
    '''
    _singletons = dict()

    def __new__(cls):
        if not cls in cls._singletons:
            cls._singletons[cls] = object.__new__(cls)
        return cls._singletons[cls]

    def __init__(self):
        self.name = None
        self.long_name = None
        self.SCAMP_aheader = None
        self.focal_length = None
        self.pixel_size = None
        self.aperture = None
        self.gain = None
        self.nXPix = None
        self.nYPix = None
        self.units_for_FWHM = None
        self.ROI = None
        self.threshold_FWHM = None
        self.threshold_pointing_err = None
        self.threshold_ellipticity = None
        self.pixel_scale = None
        self.fRatio = None
        self.SExtractor_params = None
        self.distortionOrder = 1
        self.site = None
        self.pointing_marker_size = 1*u.arcmin
        self.PSF_measurement_radius = None
        
    def check_units(self):
        '''
        Checks whether the telescope properties have the right type.  If a unit
        is expected, checks whether the input has units and whether it is
        reducible to the expected unit.  If input has no units and they are
        expected, then adds the default unit.
        '''
        ## name is a string
        assert type(self.name) == str
        ## long_name is a string
        assert type(self.long_name) == str
        ## Default focal_length units to mm
        if type(self.focal_length) == u.quantity.Quantity:
            assert self.focal_length.to(u.mm)
        else:
            self.focal_length *= u.mm
        ## Default pixel_size units to microns
        if type(self.pixel_size) == u.quantity.Quantity:
            assert self.pixel_size.to(u.micron)
        else:
            self.pixel_size *= u.micron
        ## Default aperture to units of mm
        if type(self.aperture) == u.quantity.Quantity:
            assert self.aperture.to(u.mm)
        else:
            self.aperture *= u.mm
        ## Default gain to units of 1/ADU
        if type(self.gain) == u.quantity.Quantity:
            assert self.gain.to(1/u.adu)
        else:
            self.gain *= 1./u.adu
        ## Default units_for_FWHM to units of arcsec
        if type(self.units_for_FWHM) == u.quantity.Quantity:
            assert self.units_for_FWHM.unit in [u.arcsec, u.pix]
        else:
            self.units_for_FWHM *= u.pix
        ## ROI is string
        assert type(self.ROI) == str
        ## Default threshold_FWHM to same units as units_for_FWHM
        if type(self.threshold_FWHM) == u.quantity.Quantity:
            assert self.threshold_FWHM.unit in [u.arcsec, u.pix]
        else:
            self.threshold_FWHM *= u.pix
        ## Default threshold_pointing_err to units of arcmin
        if type(self.threshold_pointing_err) == u.quantity.Quantity:
            assert self.threshold_pointing_err.to(u.arcmin)
        else:
            self.threshold_pointing_err *= u.arcmin
        ## Default threshold_ellipticity to dimensionless
        if type(self.threshold_ellipticity) == u.quantity.Quantity:
            assert self.threshold_ellipticity.to(u.dimensionless_unscaled)
        else:
            assert float(self.threshold_ellipticity) >= 0
            assert float(self.threshold_ellipticity) <= 1.
        ## Default pixel_scale to units of arcsec per pixel
        if type(self.pixel_scale) == u.quantity.Quantity:
            assert self.pixel_scale.to(u.arcsec / u.pix)
        else:
            self.pixel_scale *= u.arcsec / u.pix
        ## Default fRatio to dimensionless
        if type(self.fRatio) == u.quantity.Quantity:
            assert self.fRatio.to(u.dimensionless_unscaled)
        else:
            assert float(self.fRatio)


    ##-------------------------------------------------------------------------
    ## Define astropy.units Equivalency for Arcseconds and Pixels
    ##-------------------------------------------------------------------------
    def define_pixel_scale(self):
        self.pixel_scale_equivalency = [(u.pix, u.arcsec,
             lambda pix: (pix*u.radian.to(u.arcsec) * self.pixel_size / self.focal_length).decompose().value,
             lambda arcsec: (arcsec/u.radian.to(u.arcsec) * self.focal_length / self.pixel_size).decompose().value
             )]


##-----------------------------------------------------------------------------
## Define Image object which holds information and methods for analysis
##-----------------------------------------------------------------------------
class Image(object):
    '''
    The IQMon.Image object represents a single input image to the IQMon
    process.

    When defined, the image objects requires both a filename to a valid fits
    file and an IQMon.Config object.

    Properties:

    Methods:
    '''
    def __init__(self, input, tel=None, config=None):
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
        if tel:
            assert type(tel) == Telescope
            self.tel = tel
        ## Confirm that input config is an IQMon.Config object
        if config:
            assert type(config) == Config
            self.config = config
        ## Initialize values to None
        self.logger = None
        self.working_file = None
        self.header = None
        self.exptime = None
        self.catalog_filter = None
        self.object_name = None
        self.astrometry_solved = None
        self.coordinate_of_center_pixel = None
        self.coordinate_from_header = None
        self.n_stars_SExtracted = None
        self.SExtractor_background = None
        self.SExtractor_background_RMS = None
        self.temp_files = []
        self.SExtractor_catalog = None
        self.SExtractor_results = None
        self.position_angle = None
        self.zeroPoint = None
        self.zeroPoint_plotfile = None
        self.total_process_time = None
        self.FWHM = None
        self.ellipticity = None
        self.PSF_plotfile = None
        self.pointing_error = None
        self.image_flipped = None
        self.jpeg_file_names = []
        self.check_image_file = None
        self.cropped = False
        self.crop_x1 = None
        self.crop_x2 = None
        self.crop_y1 = None
        self.crop_y2 = None
        self.original_nXPix = None
        self.original_nYPix = None
        self.SCAMP_catalog = None
        self.catalog_file_path = None

    ##-------------------------------------------------------------------------
    ## Make Logger Object
    ##-------------------------------------------------------------------------
    def GetLogger(self, logger):
        '''
        If calling from another prohram which has its own logger object, pass
        that logger to IQMon with this method.
        '''
        self.logger = logger


    def MakeLogger(self, IQMonLogFileName, verbose):
        '''
        Create the logger object to use when processing.  Takes as input the
        full path to the file to write the log to and verboase, a boolean value
        which will increase the verbosity of the concole log (the file log will
        always be at debug level).
        '''
        self.logger = logging.getLogger('IQMonLogger')
        self.logger.setLevel(logging.DEBUG)
        LogFileHandler = logging.FileHandler(IQMonLogFileName)
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


    ##-------------------------------------------------------------------------
    ## Read Header
    ##-------------------------------------------------------------------------
    def ReadHeader(self):
        '''
        Read information from the image fits header.
        '''
        hdulist = fits.open(self.working_file, ignore_missing_end=True)
        self.header = hdulist[0].header
        self.image = hdulist[0].data
        hdulist.close()
        self.logger.info("Reading image header.")
        
        ## Get exposure time from header (assumes seconds)
        try:
            self.exptime = float(self.header['EXPTIME']) * u.s
        except:
            self.exptime = None
            self.logger.debug("  No exposure time value found in header")
        else:
            self.logger.debug("  Exposure time = {0:.1f} s".format(\
                                                   self.exptime.to(u.s).value))
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
            self.dateObs = self.header["DATE-OBS"]
        except:
            self.dateObs = None
            self.logger.debug("  No date value found in header")
        else:
            self.logger.debug("  Header date = {0}".format(self.dateObs))
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


        ## Determine Image Size in Pixels
        self.nYPix, self.nXPix = self.image.shape
        self.logger.debug('  Image size is: {},{}'.format(\
                                                       self.nXPix, self.nYPix))

        ## Read Header Coordinates in to astropy coordinates object
        if (len(self.header['RA'].split(":")) != 3) and (len(self.header['DEC'].split(":")) != 3):
            ## Header RA is : separated
            self.coordinate_from_header = coords.ICRS('{} {}'.format(\
                                            self.header['RA'],\
                                            self.header['DEC']),\
                                            unit=(u.hour, u.degree))
        elif (len(self.header['RA'].split(" ")) == 3) and (len(self.header['DEC'].split(" ")) == 3):
            ## Header RA is space separated
            self.header['RA'] = ":".join(self.header['RA'].split(" "))
            self.header['DEC'] = ":".join(self.header['DEC'].split(" "))
            self.coordinate_from_header = coords.ICRS('{} {}'.format(\
                                            self.header['RA'],\
                                            self.header['DEC']),\
                                            unit=(u.hour, u.degree))
        elif float(self.header['RA']) and float(self.header['DEC']):
            ## Header RA is decimal.  Assume degrees.
            self.coordinate_from_header = coords.ICRS('{} {}'.format(\
                                            self.header['RA'],\
                                            self.header['DEC']),\
                                            unit=(u.degree, u.degree))
        else:
            self.logger.info('  Could not parse coordinate strings from header')


        ## Read WCS
        try:
            self.image_WCS = wcs.WCS(self.header)
        except:
            self.image_WCS = None
            self.logger.info("  No WCS found in image header")
        else:
            self.logger.debug("  Found WCS in image header.")

        ## Determine PA of Image
        if self.image_WCS:
            self.wcs_pixel_scale, self.position_angle, self.image_flipped = \
                                              image_orientation(self.image_WCS)
            self.logger.debug("  Position angle of WCS is {0:.1f} deg".format(\
                                          self.position_angle.to(u.deg).value))
            if self.image_flipped:
                self.logger.debug("  Image is mirrored.")


        ## Determine Alt, Az, Moon Sep, Moon Illum using ephem module
        if self.dateObs and self.latitude and self.longitude:
            ## Populate site object properties
            SiteDate = "/".join(self.dateObs[0:10].split("-"))
            SiteTime = self.dateObs[11:]        
            self.tel.site.date = ephem.Date(SiteDate+" "+SiteTime)
            self.tel.site.lat = str(self.latitude.to(u.deg).value)
            self.tel.site.lon = str(self.longitude.to(u.deg).value)
            if self.altitude:
                self.tel.site.elevation = self.altitude.to(u.meter).value
            ## Do calculations using ephem
            TargetObject = ephem.readdb("Target,f|M|F7,"+ImageRA+","+ImageDEC+",2.02,2000")
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


    ##-------------------------------------------------------------------------
    ## Add To Header
    ##-------------------------------------------------------------------------
    def EditHeader(self, keyword, value):
        '''
        Edit information in the image fits header.
        '''
        self.logger.info('Editing image header: {} = {}'.format(keyword, value))
        hdulist = fits.open(self.working_file,\
                            ignore_missing_end=True, mode='update')
        hdulist[0].header[keyword] = value
        hdulist.flush()


    ##-------------------------------------------------------------------------
    ## Read Image
    ##-------------------------------------------------------------------------
    def ReadImage(self):
        '''
        Read the raw image and write out a working image in the IQMon temporary
        directory.
        
        - For the moment, this only copies a fits file from the original
          location to the IQMon tmp directory.
        - Later implement file format conversion from CRW, CR2, DNG, etc to
          fits using dcraw.
        '''
        if os.path.exists(self.working_file): os.remove(self.working_file)
        if self.file_ext == '.fts':
            self.logger.info('Making working copy of raw image: {}'.format(\
                                                            self.working_file))
            self.working_file = os.path.join(self.config.pathTemp,\
                                             self.raw_file_basename+'.fits')
            shutil.copy2(self.raw_file, self.working_file)
            os.chmod(self.working_file, 0666)
            self.temp_files.append(self.working_file)
            self.file_ext = '.fits'
        elif self.file_ext == '.fits':
            self.logger.info('Making working copy of raw image: {}'.format(\
                                                            self.working_file))
            self.working_file = os.path.join(self.config.pathTemp,\
                                             self.raw_file_name)
            shutil.copy2(self.raw_file, self.working_file)
            os.chmod(self.working_file, 0666)
            self.temp_files.append(self.working_file)
            self.file_ext = '.fits'
        elif self.file_ext in ['.dng', '.DNG', '.cr2', '.CR2']:
            self.logger.info('Converting {} to fits format'.format(\
                                                            self.working_file))
            sys.exit(1)
        else:
            self.logger.warning('Unrecognixed file extension: {}'.format(\
                                                                self.file_ext))
            self.working_file = os.path.join(self.config.pathTemp,\
                                             self.raw_file_name)
            sys.exit(1)




    ##-------------------------------------------------------------------------
    ## Dark Subtract Image
    ##-------------------------------------------------------------------------
    def DarkSubtract(self, Darks):
        '''
        Create master dark and subtract from image.
        
        Input the filename of the appropriate master dark.  May want to write
        own function to make the master dark given input file data.
        '''
        self.logger.info("Dark subtracting image.")
        self.logger.debug("  Opening image data.")
        hdulist_image = fits.open(self.working_file, mode='update')
        ## Load master dark if provided, but if multiple files input, combine
        ## them in to master dark, then load combined master dark.
        if len(Darks) == 1:
            self.logger.debug("  Found master dark.  Opening master dark data.")
            hdulist_dark = fits.open(Darks[0])
            MasterDarkData = hdulist_dark[0].data
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
            MasterDarkFile  = os.path.join(self.config.pathTemp,\
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


    ##-------------------------------------------------------------------------
    ## Crop Image
    ##-------------------------------------------------------------------------
    def Crop(self):
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
            hdulist = fits.open(self.working_file, mode="update")
            hdulist[0].data = hdulist[0].data[self.crop_y1:self.crop_y2,self.crop_x1:self.crop_x2]
            hdulist.flush()
            hdulist.close()
            self.cropped = True
            self.original_nXPix = self.nXPix
            self.original_nYPix = self.nYPix


    ##-------------------------------------------------------------------------
    ## Solve Astrometry Using astrometry.net
    ##-------------------------------------------------------------------------
    def SolveAstrometry(self):
        '''
        Solve astrometry in the working image using the astrometry.net solver.
        '''
        self.logger.info("Attempting to solve WCS using Astrometry.net solver.")
        AstrometryCommand = ["solve-field", "-l", "5", "-O", "-p", "-T",
                             "-L", str(self.tel.pixel_scale.value*0.90),
                             "-H", str(self.tel.pixel_scale.value*1.10),
                             "-u", "arcsecperpix", "-z", "4", self.working_file]
        AstrometrySTDOUT = ""

        try:
            StartTime = datetime.datetime.now()
            AstrometrySTDOUT = subprocess.check_output(AstrometryCommand, 
                               stderr=subprocess.STDOUT, universal_newlines=True)
            EndTime = datetime.datetime.now()
        except subprocess.CalledProcessError as e:
            self.logger.warning("Astrometry.net failed.")
            for line in e.output.split("\n"):
                self.logger.error(line)
            self.astrometry_solved = False
        except:
            self.logger.error("solve-field process failed: {0} {1} {2}".format(\
                       sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
        else:
            total_process_time = (EndTime - StartTime).total_seconds()
            self.logger.debug("  Astrometry.net Processing Time: {.1f} s".format(\
                                                           total_process_time))
            pos = AstrometrySTDOUT.find("Field center: (RA H:M:S, Dec D:M:S) = ")
            if pos != -1:
                IsFieldCenter = re.match("\s*(\d{1,2}:\d{2}:\d{2}\.\d+,\s-?\d{1,2}:\d{2}:\d{2}\.\d+).*", 
                                         AstrometrySTDOUT[pos+40:pos+75])
                if IsFieldCenter:
                    self.logger.info("  Astrometry.net field center is: {}".format(\
                                                        IsFieldCenter.group(1)))
                else:
                    self.logger.warning("Could not parse field center from astrometry.net output.")
                    for line in AstrometrySTDOUT.split("\n"):
                        self.logger.warning("  %s" % line)
            else:
                for line in AstrometrySTDOUT.split("\n"):
                    self.logger.warning("  %s" % line)
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
            self.temp_files.append(os.path.join(self.config.pathTemp,\
                                      self.raw_file_basename+".axy"))
            self.temp_files.append(os.path.join(self.config.pathTemp,\
                                      self.raw_file_basename+".wcs"))
            self.temp_files.append(os.path.join(self.config.pathTemp,\
                                      self.raw_file_basename+".solved"))
            self.temp_files.append(os.path.join(self.config.pathTemp,\
                                      self.raw_file_basename+".rdls"))
            self.temp_files.append(os.path.join(self.config.pathTemp,\
                                      self.raw_file_basename+".match"))
            self.temp_files.append(os.path.join(self.config.pathTemp,
                                      self.raw_file_basename+".corr"))
            self.temp_files.append(os.path.join(self.config.pathTemp,\
                                      self.raw_file_basename+".new.fits"))
            self.temp_files.append(os.path.join(self.config.pathTemp,\
                                      self.raw_file_basename+"-indx.xyls"))


    ##-------------------------------------------------------------------------
    ## Determine Pointing Error
    ##-------------------------------------------------------------------------
    def DeterminePointingError(self):
        '''
        Determine pointing error (difference between objects coordinates and
        solved WCS).
        '''
        self.logger.info("Detemining pointing error based on WCS solution")
        if self.image_WCS:
            center_from_WCS = self.image_WCS.wcs_pix2world([[self.nXPix/2, self.nYPix/2]], 1)
            self.logger.debug("  Using coordinates of center point: {0} {1}".format(\
                                 center_from_WCS[0][0], center_from_WCS[0][1]))
            self.coordinate_of_center_pixel = coords.ICRS(\
                                              ra=center_from_WCS[0][0],\
                                              dec=center_from_WCS[0][1],\
                                              unit=(u.degree, u.degree))
            self.pointing_error = self.coordinate_of_center_pixel.separation(\
                                                   self.coordinate_from_header)
            self.logger.debug("  Target Coordinates are:  {}".format(
                              self.coordinate_from_header.to_string(sep=":", precision=1, alwayssign=True))),
            self.logger.debug("  WCS of Central Pixel is: {}".format(
                              self.coordinate_of_center_pixel.to_string(sep=":", precision=1, alwayssign=True)))
            self.logger.info("  Pointing Error is {:.2f} arcmin".format(\
                                                self.pointing_error.arcminute))
        else:
            self.logger.warning("Pointing error not calculated.")


    ##-------------------------------------------------------------------------
    ## Run SExtractor
    ##-------------------------------------------------------------------------
    def RunSExtractor(self, assoc=False):
        '''
        Run SExtractor on image.
        '''
        assert type(self.tel.gain) == u.quantity.Quantity
        assert type(self.tel.pixel_scale) == u.quantity.Quantity
        
        if assoc:
            if self.header['FILTER']:
                if self.header['FILTER'] in self.catalog.keys():
                    self.catalog_filter = self.header['FILTER']
                else:
                    self.logger.warning('  Filter from header ({}), not found in UCAC catalog table.'.format(\
                                                        self.header['FILTER']))
                    self.logger.info('  Using r filter for catalog magnitudes.')
                    self.catalog_filter = 'r'
            else:
                self.logger.warning('  Filter from header ({}), not found in UCAC catalog table.'.format(\
                                                        self.header['FILTER']))
                self.logger.info('  Using r filter for catalog magnitudes.')
                self.catalog_filter = 'r'

        ## Set up file names
        self.SExtractor_catalog = os.path.join(self.config.pathTemp,\
                                               self.raw_file_basename+".cat")
        self.temp_files.append(self.SExtractor_catalog)

        sextractor_output_param_file = os.path.join(self.config.pathTemp,\
                                                   'default.param')
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

        self.check_image_file = os.path.join(self.config.pathTemp,\
                                        self.raw_file_basename+"_bksub.fits")
        self.temp_files.append(self.check_image_file)
        ## Compare input parameters dict to default
        SExtractor_default = {
                             'CATALOG_NAME': self.SExtractor_catalog,
                             'CATALOG_TYPE': 'FITS_LDAC',
                             'PARAMETERS_NAME': sextractor_output_param_file,
                             'GAIN': self.tel.gain.value,
                             'GAIN_KEY': 'GAIN',
                             'PIXEL_SCALE': '{:.3f}'.format(self.tel.pixel_scale.value),
                             'CHECKIMAGE_TYPE': '-BACKGROUND',
                             'CHECKIMAGE_NAME': self.check_image_file,
                            }

        ## Use command line sextractor params
        if not self.tel.SExtractor_params:
            SExtractor_params = SExtractor_default
        else:
            SExtractor_params = self.tel.SExtractor_params
            for key in SExtractor_default.keys():
                if not key in self.tel.SExtractor_params.keys():
                    SExtractor_params[key] = SExtractor_default[key]

        if assoc:
            assert os.path.exists(self.catalog_file_path)
            assert os.path.exists(os.path.join(self.config.pathTemp, 'scamp.xml'))
            assert self.catalog_filter in self.catalog.keys()

            ## Create Assoc file with pixel coordinates of catalog stars
            assoc_file = os.path.join(self.config.pathTemp, 'assoc.txt')
            self.temp_files.append(assoc_file)
            if os.path.exists(assoc_file): os.remove(assoc_file)
            assocFO = open(assoc_file, 'w')
            for star in self.catalog:
                pix = self.image_WCS.wcs_world2pix([[star['RA'], star['Dec']]], 1)
                try:
                    assocFO.write('{:8.1f} {:8.1f} {:8.1f}\n'.format(\
                                                    pix[0][0], pix[0][1],\
                                                    star[self.catalog_filter]))
                except:
                    pass
            assocFO.close()

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
        self.logger.debug("  SExtractor command: {}".format(repr(SExtractorCommand)))
        try:
            SExSTDOUT = subprocess.check_output(SExtractorCommand,\
                             stderr=subprocess.STDOUT, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            self.logger.error("SExtractor failed.  Command: {}".format(e.cmd))
            self.logger.error("SExtractor failed.  Returncode: {}".format(e.returncode))
            self.logger.error("SExtractor failed.  Output: {}".format(e.output))
        except:
            self.logger.error("SExtractor process failed: {0} {1} {2}".format(\
                      sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
        else:
            for line in SExSTDOUT.splitlines():
                line.replace("[1A", "")
                line.replace("[1M>", "")
                if not re.match(".*Setting up background map.*", line) and\
                   not re.match(".*Line:\s[0-9]*.*", line):
                    self.logger.info("  SExtractor Output: {}".format(line))
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
            if not os.path.exists(self.SExtractor_catalog):
                self.logger.warning("SExtractor failed to create catalog.")
                self.SExtractor_catalog = None

            ## Read FITS_LDAC SExtractor Catalog
            self.logger.debug("  Reading SExtractor output catalog.")
            hdu = fits.open(self.SExtractor_catalog)
            self.SExtractor_results = table.Table(hdu[2].data)
            SExImageRadius = []
            SExAngleInImage = []
            zp_diff = []
            for star in self.SExtractor_results:
                SExImageRadius.append(math.sqrt((self.nXPix/2-star['XWIN_IMAGE'])**2 +\
                                                (self.nYPix/2-star['YWIN_IMAGE'])**2))
                SExAngleInImage.append(math.atan((star['XWIN_IMAGE']-self.nXPix/2) /\
                                                 (self.nYPix/2-star['YWIN_IMAGE']))*180.0/math.pi)
                if assoc:
                    zp_diff.append(star['VECTOR_ASSOC'][2] - star['MAG_AUTO'])
            self.SExtractor_results.add_column(table.Column(\
                                    data=SExImageRadius, name='ImageRadius'))
            self.SExtractor_results.add_column(table.Column(\
                                    data=SExAngleInImage, name='AngleInImage'))
            self.n_stars_SExtracted = len(self.SExtractor_results)
            self.logger.info("  Read in {0} stars from SExtractor catalog.".format(\
                                                      self.n_stars_SExtracted))
        if assoc:
            self.SExtractor_results.add_column(table.Column(\
                                               data=zp_diff, name='MagDiff'))
            self.tel.SExtractor_params = original_params


    ##-------------------------------------------------------------------------
    ## Determine Image FWHM from SExtractor Catalog
    ##-------------------------------------------------------------------------
    def DetermineFWHM(self):
        '''
        Determine typical FWHM of image from SExtractor results.
        '''
        if self.n_stars_SExtracted > 1:
            self.logger.info('Analyzing SExtractor results to determine typical image quality.')
            if self.tel.PSF_measurement_radius:
                self.logger.info('  Using stars in the inner {} pixels.'.format(\
                                              self.tel.PSF_measurement_radius))
                IQRadius = self.tel.PSF_measurement_radius
            else:
                IQRadiusFactor = 1.0
                DiagonalRadius = math.sqrt((self.nXPix/2)**2+(self.nYPix/2)**2)
                IQRadius = DiagonalRadius*IQRadiusFactor
            CentralFWHMs = [star['FWHM_IMAGE'] for star in self.SExtractor_results if star['ImageRadius'] <= IQRadius]
            CentralEllipticities = [star['ELLIPTICITY'] for star in self.SExtractor_results if star['ImageRadius'] <= IQRadius]
            CentralAs = [star['AWIN_IMAGE'] for star in self.SExtractor_results if star['ImageRadius'] <= IQRadius]
            CentralBs = [star['BWIN_IMAGE'] for star in self.SExtractor_results if star['ImageRadius'] <= IQRadius]
            if len(CentralFWHMs) > 3:
                self.FWHM = np.median(CentralFWHMs) * u.pix
                self.ellipticity = np.median(CentralEllipticities)
                self.major_axis = np.median(CentralAs) * u.pix
                self.minor_axis = np.median(CentralBs) * u.pix
            else:
                self.logger.warning("  Not enough stars detected in central region of image to form median FWHM.")
            self.logger.debug("  Using {0} stars in central region to determine FWHM and ellipticity.".format(\
                                                            len(CentralFWHMs)))
            self.logger.info("  Median FWHM in inner region is {0:.2f} pixels".format(\
                                                    self.FWHM.to(u.pix).value))
            self.logger.info("  Median Minor Axis in inner region is {0:.2f}".format(\
                                        2.355*self.minor_axis.to(u.pix).value))
            self.logger.info("  Median Major Axis in inner region is {0:.2f}".format(\
                                        2.355*self.major_axis.to(u.pix).value))
            self.logger.info("  Median Ellipticity in inner region is {0:.2f}".format(\
                                                             self.ellipticity))
        else:
            self.FWHM = None
            self.ellipticity = None


    ##-------------------------------------------------------------------------
    ## Make Ellipticity Plot
    ##-------------------------------------------------------------------------
    def MakePSFplot(self, filename=None):
        '''
        Plot ellipticity vectors of stars in image.
        Plot histogram of theta - image_angle values.
        '''
        self.logger.info('Calculating histogram of PSF angles.')
        if filename:
            self.PSF_plot_filename = filename
        else:
            self.PSF_plot_filename = self.raw_file_basename+'_PSFinfo.png'
        self.PSF_plotfile = os.path.join(self.config.pathPlots, self.PSF_plot_filename)

        ellip_threshold = 0.15
        star_angles = [star['THETAWIN_IMAGE'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
        image_angles = [star['AngleInImage'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
        star_x = [star['XWIN_IMAGE'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
        star_y = [star['YWIN_IMAGE'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
        uncorrected_diffs = [star['THETAWIN_IMAGE']-star['AngleInImage'] for star in self.SExtractor_results if star['ELLIPTICITY'] >= ellip_threshold]
        
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
        ellip_hist, ellip_bins = np.histogram(self.SExtractor_results['ELLIPTICITY'],\
                                              bins=ellip_binsize*np.arange(21))
        ellip_centers = (ellip_bins[:-1] + ellip_bins[1:]) / 2

        fwhm_binsize = 0.1
        fwhm_hist, fwhm_bins = np.histogram(self.SExtractor_results['FWHM_IMAGE'],\
                                               bins=fwhm_binsize*np.arange(91))
        fwhm_centers = (fwhm_bins[:-1] + fwhm_bins[1:]) / 2

        star_angle_mean = np.mean(star_angles)
        star_angle_median = np.median(star_angles)
        angle_diff_mean = np.mean(angle_diffs)
        angle_diff_median = np.median(angle_diffs)
        self.logger.debug('  Mean Stellar PA = {:.0f}'.format(star_angle_mean))
        self.logger.debug('  Median Stellar PA = {:.0f}'.format(star_angle_median))
        self.logger.debug('  Mean Difference Angle = {:.0f}'.format(angle_diff_mean))
        self.logger.debug('  Median Difference Angle = {:.0f}'.format(angle_diff_median))

        if self.PSF_plotfile:
            self.logger.debug('  Generating figure {}'.format(self.PSF_plotfile))

            pyplot.ioff()
            pyplot.figure(figsize=(12,11), dpi=100)
            Left1 = pyplot.axes([0.000, 0.750, 0.465, 0.240])
            pyplot.title('PSF Statistics for {}'.format(self.raw_file_name), size=10)
            pyplot.bar(ellip_centers, ellip_hist, align='center', width=0.7*ellip_binsize)
            pyplot.xlabel('Ellipticity', size=10)
            pyplot.ylabel('N Stars', size=10)
            pyplot.xlim(0,1)
            pyplot.xticks(0.1*np.arange(11), size=10)
            pyplot.yticks(size=10)

            Left2 = pyplot.axes([0.000, 0.460, 0.465, 0.240])
            pyplot.bar(fwhm_centers, fwhm_hist, align='center', width=0.7*fwhm_binsize)
            pyplot.xlabel('FWHM (pixels)', size=10)
            pyplot.ylabel('N Stars', size=10)
            if self.FWHM:
                pyplot.xlim(0,self.FWHM.to(u.pix).value + 3)
            else:
                pyplot.xlim(0,6)
            pyplot.xticks(size=10)
            pyplot.yticks(size=10)

            Left3 = pyplot.axes([0.000, 0.0, 0.465, 0.400])
            pyplot.plot(self.SExtractor_results['ImageRadius'],\
                        self.SExtractor_results['ELLIPTICITY'],\
                        'k,')
            pyplot.title('Correlation of Ellipticity with Image Radius')
            pyplot.xlabel('r (pixels)', size=10)
            pyplot.ylabel('Ellipticity', size=10)
            pyplot.xlim(0, math.sqrt(self.nXPix**2 + self.nYPix**2)/2.)
            pyplot.ylim(0, 1.0)
            pyplot.xticks(size=10)
            pyplot.yticks(size=10)

#             Left3 = pyplot.axes([0.000, 0.0, 0.465, 0.400])
#             Left3.set_aspect('equal')
#             pyplot.plot(star_x, star_y, 'k,')
#             pyplot.title('Positions of {:d}/{:d} stars with ellipticity > {:.2f}'.format(nstars, self.n_stars_SExtracted, ellip_threshold), size=10)
#             pyplot.xlabel('X (pixels)', size=10)
#             pyplot.ylabel('Y (pixels)', size=10)
#             pyplot.xlim(0, self.nXPix)
#             pyplot.ylim(0, self.nYPix)
#             pyplot.xticks(size=10)
#             pyplot.yticks(size=10)

            Right1 = pyplot.axes([0.535, 0.750, 0.465, 0.240])
            pyplot.title('PSF Angles for {:d}/{:d} stars with ellipticity > {:.2f}'.format(\
                            nstars, self.n_stars_SExtracted, ellip_threshold),\
                            size=10)
            pyplot.bar(angle_centers, angle_hist, align='center', width=0.7*angle_binsize)
            pyplot.xlabel('Stellar PSF PA', size=10)
            pyplot.ylabel('N Stars', size=10)
            pyplot.xlim(-90,90)
            pyplot.xticks(30*(np.arange(7)-3), size=10)
            pyplot.yticks(size=10)

            Right2 = pyplot.axes([0.535, 0.460, 0.465, 0.240])
            pyplot.bar(angle_centers, diff_hist, align='center', width=0.7*angle_binsize)
            pyplot.xlabel('Stellar PSF PA - Image PA', size=10)
            pyplot.ylabel('N Stars', size=10)
            pyplot.xlim(-90,90)
            pyplot.xticks(30*(np.arange(7)-3), size=10)
            pyplot.yticks(size=10)
            
            Right3 = pyplot.axes([0.535, 0.0, 0.465, 0.400])
            Right3.set_aspect('equal')
            pyplot.plot(star_angles, image_angles, 'k.')
            pyplot.title('Correlation Between PSF Angle and Position in Image', size=10)
            pyplot.xlabel('Stellar PSF PA', size=10)
            pyplot.ylabel('Image PA', size=10)
            pyplot.xlim(-100,100)
            pyplot.xticks(30*(np.arange(7)-3), size=10)
            pyplot.ylim(-100,100)
            pyplot.yticks(30*(np.arange(7)-3), size=10)

            pyplot.savefig(self.PSF_plotfile, dpi=100, bbox_inches='tight', pad_inches=0.10)


    ##-------------------------------------------------------------------------
    ## Run SCAMP
    ##-------------------------------------------------------------------------
    def RunSCAMP(self, catalog='USNO-B1', mergedcat_name='scamp.cat', mergedcat_type='ASCII_HEAD'):
        '''
        Run SCAMP on SExtractor output catalog.
        '''
        ## Parameters for SCAMP
        if self.tel.SCAMP_aheader:
            SCAMP_aheader = self.tel.SCAMP_aheader
        else:
            SCAMP_aheader = 'scamp.ahead'
        SCAMP_params = {
                        'DISTORT_DEGREES': self.tel.distortionOrder,
                        'SCAMP_aheader_GLOBAL': SCAMP_aheader,
                        'ASTREF_CATALOG': catalog,
                        'SAVE_REFCATALOG': 'N',
                        'REFOUT_CATPATH': self.config.pathTemp,
                        'MERGEDOUTCAT_NAME': mergedcat_name,
                        'MERGEDOUTCAT_TYPE': mergedcat_type,
                        'CHECKPLOT_RES': '1200,1200',
                        'CHECKPLOT_TYPE': 'FGROUPS,DISTORTION,ASTR_REFERROR2D,ASTR_REFERROR1D,PHOT_ZPCORR,ASTR_REFSYSMAP',
                        'CHECKPLOT_NAME': 'fgroups,distortion,astr_referror2d,astr_referror1d,phot_zpcorr,astr_refsysmap',
                        'CROSSID_RADIUS': 6.0,
                        'SOLVE_PHOTOM': 'Y',
                        'ASTRINSTRU_KEY': 'QRUNID',
                        'WRITE_XML': 'Y',
                        'XML_NAME': os.path.join(self.config.pathTemp, 'scamp.xml'),
                        }
        SCAMPCommand = ["scamp", self.SExtractor_catalog]
        for key in SCAMP_params.keys():
            SCAMPCommand.append('-{}'.format(key))
            SCAMPCommand.append('{}'.format(SCAMP_params[key]))
        self.logger.info("Running SCAMP on {} catalog with distortion order {}.".format(\
                                            catalog, self.tel.distortionOrder))
        if SCAMP_aheader:
            self.logger.info("  Using SCAMP_aheader file: {}".format(SCAMP_aheader))
        self.logger.debug("  SCAMP command: {}".format(SCAMPCommand))
        try:
            SCAMP_STDOUT = subprocess.check_output(SCAMPCommand,\
                             stderr=subprocess.STDOUT, universal_newlines=True)
            self.temp_files.append(os.path.join(self.config.pathTemp, 'scamp.xml'))
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
                if re.search('Astrometric stats \(external\)', line):
                    StartAstrometricStats = True
                if re.search('Generating astrometric plots', line):
                    EndAstrometricStats = True
                if StartAstrometricStats and not EndAstrometricStats:
                    self.logger.info("  SCAMP Output: "+line)
                else:
                    self.logger.debug("  SCAMP Output: "+line)
        ## Store Output Catalog Name
        if os.path.exists(mergedcat_name):
            self.SCAMP_catalog = mergedcat_name

        ## Populate FITS header with SCAMP derived header values in .head file
        head_file = os.path.splitext(self.working_file)[0]+'.head'
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
        else:
            self.logger.critical('No .head file found from SCAMP.')
            sys.exit(1)


    ##-------------------------------------------------------------------------
    ## Run SWarp
    ##-------------------------------------------------------------------------
    '''
    Run SWarp on the image (after SCAMP distortion solution) to de-distort it.
    '''
    def RunSWarp(self):
        ## Parameters for SWarp
        swarp_file = os.path.join(self.config.pathTemp, 'swarpped.fits')
        if os.path.exists(swarp_file): os.remove(swarp_file)
        SWarp_params = {'IMAGEOUT_NAME': swarp_file,
                        'COPY_KEYWORDS': 'FILTER,OBJECT,AIRMASS,DATE-OBS,LAT-OBS,LONG-OBS,ALT-OBS,RA,DEC',
                        'WRITE_XML': 'Y',
                        'XML_NAME': os.path.join(self.config.pathTemp, 'swarp.xml'),
                       }
        SWarpCommand = ["swarp", self.working_file]
        for key in SWarp_params.keys():
            SWarpCommand.append('-{}'.format(key))
            SWarpCommand.append('{}'.format(SWarp_params[key]))
        self.logger.info("Running SWarp.")
        self.logger.debug("  SWarp command: {}".format(SWarpCommand))
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
                self.logger.debug("  SWarp Output: "+line)
        ## Replace working_file with SWarp output file
        if os.path.exists(swarp_file):
            self.temp_files.append(os.path.join(self.config.pathTemp, 'swarp.xml'))
            self.logger.debug('  SWarp process succeeded.')
            self.logger.debug('  Moving SWarpped file to working file.')
            if os.path.exists(self.working_file): os.remove(self.working_file)
            os.rename(swarp_file, self.working_file)
            assert os.path.exists(self.working_file)


    ##-------------------------------------------------------------------------
    ## Get UCAC4 Catalog for Image from Local File
    ##-------------------------------------------------------------------------
    def GetLocalUCAC4(self,\
                      local_UCAC_command="/Volumes/Data/UCAC4/access/u4test",\
                      local_UCAC_data="/Volumes/Data/UCAC4/u4b"):
        '''
        Determine zero point by comparing measured magnitudes with catalog
        magnitudes.
        '''
        assert type(self.coordinate_of_center_pixel) == coords.builtin_systems.ICRS

        if not os.path.exists(local_UCAC_command):
            self.logger.warning('Cannot find local UCAC command: {}'.format(\
                                                           local_UCAC_command))
        elif not os.path.exists(local_UCAC_data):
            self.logger.warning('Cannot find local UCAC data: {}'.format(local_UCAC_data))
        else:
            corners = self.image_WCS.wcs_pix2world([[0, 0], [self.nXPix, 0], [0, self.nYPix], [self.nXPix, self.nYPix]], 1)
            field_size_RA = max(corners[:,0]) - min(corners[:,0])
            field_size_DEC = max(corners[:,1]) - min(corners[:,1])
            self.logger.info("Getting stars from local UCAC4 catalog.")
            UCACcommand = [local_UCAC_command,
                           "{:.4f}".format(self.coordinate_of_center_pixel.ra.degree),
                           "{:.4f}".format(self.coordinate_of_center_pixel.dec.degree),
                           "{:.2f}".format(field_size_RA),
                           "{:.2f}".format(field_size_DEC),
                           local_UCAC_data]
            self.logger.debug("  Using command: {}".format(UCACcommand))
            if os.path.exists("ucac4.txt"): os.remove("ucac4.txt")
            result = subprocess.call(UCACcommand)
            if os.path.exists('ucac4.txt'):
                self.catalog_file_path = os.path.join(self.config.pathTemp, 'ucac4.txt')
                shutil.move('ucac4.txt', self.catalog_file_path)
                self.temp_files.append(self.catalog_file_path)

            ## Read in UCAC catalog
            self.catalog = ascii.read(self.catalog_file_path, Reader=ascii.FixedWidthNoHeader,
                                      data_start=1, guess=False,
                                      names=('id', 'RA', 'Dec', 'mag1', 'mag2', 'smag', 'ot', 'dsf', 'RAepoch', 'Decepoch', 'dRA', 'dde', 'nt', 'nu', 'nc', 'pmRA', 'pmDec', 'sRA', 'sDec', '2mass', 'j', 'h', 'k', 'e2mphos', 'icq_flag', 'B', 'V', 'g', 'r', 'i'),
                                      col_starts=(0, 10, 24, 36, 43, 50, 54, 57, 60, 68, 76, 80, 84, 87, 90, 93, 100, 107, 111, 115, 126, 133, 140, 147, 159, 168, 175, 182, 189, 196),
                                      col_ends=(  9, 22, 35, 42, 49, 53, 56, 59, 67, 75, 79, 83, 86, 89, 92, 99, 106, 110, 114, 125, 132, 139, 146, 158, 167, 174, 181, 188, 195, 202),
                                     )
            nUCACStars = len(self.catalog)
            self.logger.info("  Retrieved {} lines from UCAC catalog.".format(nUCACStars))


    ##-------------------------------------------------------------------------
    ## Measure Zero Point
    ##-------------------------------------------------------------------------
    def MeasureZeroPoint(self, plot=False):
        '''
        '''
        assert 'VECTOR_ASSOC' in self.SExtractor_results.keys()
        assert 'MagDiff' in self.SExtractor_results.keys()
        ZeroPoint_mean = np.mean(self.SExtractor_results['MagDiff'])
        ZeroPoint_median = np.median(self.SExtractor_results['MagDiff'])
        self.logger.debug('Mean Zero Point = {:.2f}'.format(ZeroPoint_mean))
        self.logger.info('Median Zero Point = {:.2f}'.format(ZeroPoint_median))
        self.zeroPoint = ZeroPoint_median
        ## Make Plot if Requested
        if plot:
            self.logger.info('Making ZeroPoint Plot')
            self.zeroPoint_plotfile = os.path.join(self.config.pathPlots, self.raw_file_basename+'_ZeroPoint.png')
            pyplot.figure(figsize=(9,11), dpi=100)

            Fig1 = pyplot.axes([0.0, 0.5, 1.0, 0.4])
            pyplot.title('Instrumental Magnitudes vs. Calalog Magnitudes (Zero Point = {:.2f})'.format(self.zeroPoint))
            pyplot.plot(self.SExtractor_results['VECTOR_ASSOC'].data[:,2], self.SExtractor_results['MAG_AUTO'],\
                        'bo', markersize=4, markeredgewidth=0)
            pyplot.xlabel('UCAC4 {} Magnitude'.format(self.catalog_filter))
            pyplot.ylabel('Instrumental Magnitude')
            pyplot.grid()
            reject_fraction = 0.01
            ## Set Limits to XXth percentile of magnitudes in Y axis
            sorted_inst_mag = sorted(self.SExtractor_results['MAG_AUTO'])
            minmax_idx_inst_mag = [int(reject_fraction*len(sorted_inst_mag)), int((1.0-reject_fraction)*len(sorted_inst_mag))]
            pyplot.ylim(math.floor(sorted_inst_mag[minmax_idx_inst_mag[0]]), math.ceil(sorted_inst_mag[minmax_idx_inst_mag[1]]))
            ## Set Limits to XXth percentile of magnitudes in X axis
            sorted_cat_mag = sorted(self.SExtractor_results['VECTOR_ASSOC'].data[:,2])
            minmax_idx_cat_mag = [int(reject_fraction*len(sorted_cat_mag)), int((1.0-reject_fraction)*len(sorted_cat_mag))]
            pyplot.xlim(math.floor(sorted_cat_mag[minmax_idx_cat_mag[0]]), math.ceil(sorted_cat_mag[minmax_idx_cat_mag[1]]))
            ## Plot Fitted Line
            fit_mags_cat = [math.floor(sorted_cat_mag[minmax_idx_cat_mag[0]]), math.ceil(sorted_cat_mag[minmax_idx_cat_mag[1]])]
            fit_mags_inst = fit_mags_cat - self.zeroPoint
            pyplot.plot(fit_mags_cat, fit_mags_inst, 'k-', alpha=0.5, label='Zero Point = {:.2f}'.format(self.zeroPoint))

            Fig2 = pyplot.axes([0.0, 0.0, 1.0, 0.4])
            pyplot.title('Magnitude Residuals (Zero Point = {:.2f})'.format(self.zeroPoint))
            residuals = (self.SExtractor_results['MAG_AUTO'].data + self.zeroPoint) - self.SExtractor_results['VECTOR_ASSOC'].data[:,2]
            pyplot.plot(self.SExtractor_results['VECTOR_ASSOC'].data[:,2], residuals, \
                        'bo', markersize=4, markeredgewidth=0)
            pyplot.xlabel('UCAC4 {} Magnitude'.format(self.catalog_filter))
            pyplot.ylabel('Magnitude Residual')
            pyplot.grid()
            ## Set Limits to XXth percentile of magnitudes in X axis
            reject_fraction = 0.01
            sorted_cat_mag = sorted(self.SExtractor_results['VECTOR_ASSOC'].data[:,2])
            minmax_idx_cat_mag = [int(reject_fraction*len(sorted_cat_mag)), int((1.0-reject_fraction)*len(sorted_cat_mag))]
            pyplot.xlim(math.floor(sorted_cat_mag[minmax_idx_cat_mag[0]]), math.ceil(sorted_cat_mag[minmax_idx_cat_mag[1]]))
            ## Set Limits to XXth percentile of magnitudes in Y axis
            sorted_residuals = sorted(residuals)
            minmax_idx_residuals = [int(reject_fraction*len(sorted_residuals)), int((1.0-reject_fraction)*len(sorted_residuals))]
            pyplot.ylim(sorted_residuals[minmax_idx_residuals[0]]-0.1, sorted_residuals[minmax_idx_residuals[1]]+0.1)
            ## Plot Zero Line
            pyplot.plot(fit_mags_cat, [0, 0], 'k-', alpha=0.5, label='Zero Point = {:.2f}'.format(self.zeroPoint))

            pyplot.savefig(self.zeroPoint_plotfile, dpi=100, bbox_inches='tight', pad_inches=0.10)


    ##-------------------------------------------------------------------------
    ## Make JPEG of Image
    ##-------------------------------------------------------------------------
    def MakeJPEG(self, jpegFileName, binning=1, markCatalogStars=False, markDetectedStars=False, markPointing=False, backgroundSubtracted=False, p1=0.2, p2=1.0):
        '''
        Make jpegs of image.
        '''
        nStarsLimit = 5000
        jpegFile = os.path.join(self.config.pathPlots, jpegFileName)
        self.logger.info("Making jpeg (binning = {0}): {1}.".format(binning, jpegFileName))
        if os.path.exists(jpegFile): os.remove(jpegFile)
        binningString = str(1./binning*100)+"%"
        JPEGcommand = ["convert", "-contrast-stretch", "{}%,{}%".format(p1, p2), "-compress", "JPEG", "-quality", "70", "-resize", binningString]
        self.logger.debug('  Base convert command: {}'.format(' '.join(JPEGcommand)))
        ## Mark Intended Pointing Coordinates as read from header
        if markPointing and self.image_WCS:
            self.logger.debug("  Marking target pointing in jpeg.")
            markSize = (self.tel.pointing_marker_size.to(u.arcsec)/self.tel.pixel_scale).value/binning
            ## Mark Central Pixel with a White Cross
            JPEGcommand.append("-stroke")
            JPEGcommand.append("white")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("3")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            pixelCenter = [self.nXPix/2/binning, self.nYPix/2/binning]
            self.logger.debug("  Marking central pixel of JPEG: {:.1f},{:.1f}".format(pixelCenter[0], pixelCenter[1]))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0], pixelCenter[1]+markSize,
                               pixelCenter[0], pixelCenter[1]+markSize*0.3))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0]+markSize, pixelCenter[1],
                               pixelCenter[0]+markSize*0.3, pixelCenter[1]))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0], pixelCenter[1]-markSize,
                               pixelCenter[0], pixelCenter[1]-markSize*0.3))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0]-markSize, pixelCenter[1],
                               pixelCenter[0]-markSize*0.3, pixelCenter[1]))
            ## Mark Coordinates of Target with a blue Circle
            JPEGcommand.append("-stroke")
            JPEGcommand.append("blue")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("3")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            ## This next block of code seems to make the call to wcs_world2pix
            ## happy, but I'm not sure I understand why.
            foo = np.array([[self.coordinate_from_header.ra.hour*15., self.coordinate_from_header.dec.radian*180./math.pi], 
                            [self.coordinate_from_header.ra.hour*15., self.coordinate_from_header.dec.radian*180./math.pi]])
            targetPixel = (self.image_WCS.wcs_world2pix(foo, 1)[0])
            self.logger.debug("  Pixel of target on raw image: {:.1f},{:.1f}".format(targetPixel[0], targetPixel[1]))
            ## Adjust target pixel value for different origin in ImageMagick
            TargetXPos = targetPixel[0]
            if not self.cropped:
                TargetYPos = self.nYPix - targetPixel[1]
            else:
                TargetYPos = self.original_nYPix - targetPixel[1]
            ## Adjust target pixel value for cropping
            if self.cropped:
                TargetXPos = TargetXPos - self.crop_x1
                TargetYPos = TargetYPos - self.crop_y1
            ## Adjust target pixel value for binning
            TargetXPos = TargetXPos/binning
            TargetYPos = TargetYPos/binning
            self.logger.debug("  Marking pixel of target on JPEG: {:.1f},{:.1f}".format(TargetXPos, TargetYPos))
            JPEGcommand.append('-draw')
            JPEGcommand.append("circle %d,%d %d,%d" % (TargetXPos, TargetYPos,
                               TargetXPos+markSize/2, TargetYPos))
            ## Write label describing marking of pointing
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            JPEGcommand.append("text 200,40 'Blue circle centered on target is {:.1f} arcmin diameter.'".format(self.tel.pointing_marker_size.to(u.arcmin).value))


        ## Mark Stars Detected by SExtractor
        if markDetectedStars and self.SExtractor_results:
            nStarsMarked = 0
            self.logger.debug("  Marking stars found by SExtractor in jpeg.")
            JPEGcommand.append("-stroke")
            JPEGcommand.append("red")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("1")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            if self.FWHM:
                MarkRadius=max([6, 2*math.ceil(self.FWHM.value)])/binning
            else:
                MarkRadius = 6
            sortedSExtractor_results = np.sort(self.SExtractor_results, order=['MAG_AUTO'])
            for star in sortedSExtractor_results:
                nStarsMarked += 1
                if nStarsMarked <= nStarsLimit:
                    MarkXPos = star['XWIN_IMAGE']/binning
                    MarkYPos = (self.nYPix - star['YWIN_IMAGE'])/binning
                    JPEGcommand.append('-draw')
                    JPEGcommand.append("circle %d,%d %d,%d" % (MarkXPos, MarkYPos, MarkXPos+MarkRadius, MarkYPos))
                else:
                    self.logger.info("  Only marked brightest {} stars found in image.".format(nStarsLimit))
                    break
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            if nStarsMarked > nStarsLimit:
                JPEGcommand.append("text 200,80 'Red circles indicate SExtractor detected stars.  Marked {} brightest stars out of {} detected.'".format(nStarsLimit, self.n_stars_SExtracted))
            else:
                JPEGcommand.append("text 200,80 'Red circles indicate SExtractor detected stars.  Marked {} detected stars.'".format(self.n_stars_SExtracted))
        ## Mark Catalog Stars
        if markCatalogStars and self.image_WCS:
            ## Need to check if header includes distortion terms
            nStarsMarked = 0
            self.logger.debug("  Marking stars from catalog in jpeg.")
            JPEGcommand.append("-stroke")
            JPEGcommand.append("green")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("1")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            if self.FWHM:
                MarkRadius=max([6, 2*math.ceil(self.FWHM.value)])/binning
            else:
                MarkRadius = 6
            sorted_catalog = np.sort(self.catalog, order=['mag1'])
            for star in sorted_catalog:
                nStarsMarked += 1
                if nStarsMarked <= nStarsLimit:
                    pix = self.image_WCS.wcs_world2pix([[star['RA'], star['Dec']]], 1)
                    MarkXPos = pix[0][0]/binning
                    MarkYPos = (self.nYPix - pix[0][1])/binning
                    JPEGcommand.append('-draw')
                    JPEGcommand.append("circle %d,%d %d,%d" % (MarkXPos, MarkYPos, MarkXPos+MarkRadius, MarkYPos))
                else:
                    self.logger.info("  Only marked brightest {} stars found in image.".format(nStarsLimit))
                    break
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            if nStarsMarked > nStarsLimit:
                JPEGcommand.append("text 200,120 'Green circles indicate catalog stars.  Marked {} brightest stars out of {} in catalog.'".format(nStarsLimit, len(self.catalog)))
            else:
                JPEGcommand.append("text 200,120 'Green circles indicate catalog stars.  Marked {} catalog stars.'".format(self.n_stars_SExtracted))
        ## Use background subtracted image generated by SExtractor
        if not backgroundSubtracted:
            JPEGcommand.append(self.working_file)
        else:
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            JPEGcommand.append("text 200,120 'Background Subtracted Image'")
            JPEGcommand.append(self.check_image_file)
        JPEGcommand.append(jpegFile)
        self.logger.debug("  Issuing convert command to create jpeg from {}.".format(self.working_file))
        try:
            ConvertSTDOUT = subprocess.check_output(JPEGcommand, stderr=subprocess.STDOUT, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            self.logger.error("Failed to create jpeg.")
            for line in e.output.split("\n"):
                self.logger.error(line)
        except OSError as e:
            self.logger.error("Failed to create jpeg.")
            for line in e.strerror.split("\n"):
                self.logger.error(line)
        except:
            self.logger.error("Convert process failed: {0} {1} {2}".format(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
        else:
            for line in ConvertSTDOUT.split("\n"):
                if len(line) > 0:
                    self.logger.debug(line)
            self.jpeg_file_names.append(jpegFileName)


    ##-------------------------------------------------------------------------
    ## Make JPEG of Image
    ##-------------------------------------------------------------------------
    def old_MakeJPEG(self, jpegFileName, markDetectedStars=False, markPointing=False, rotate=False, binning=1, backgroundSubtracted=False):
        '''
        Make jpegs of image.
        '''
        nStarsMarked = 0
        nStarsLimit = 5000
        jpegFile = os.path.join(self.config.pathPlots, jpegFileName)
        self.logger.info("Making jpeg (binning = {0}): {1}.".format(binning, jpegFileName))
        if os.path.exists(jpegFile): os.remove(jpegFile)
        binningString = str(1./binning*100)+"%"
        JPEGcommand = ["convert", "-contrast-stretch", "0.9%,1%", "-compress", "JPEG", "-quality", "70", "-resize", binningString]
        ## Mark Target Coordinates as read from header
        if markPointing and self.image_WCS and self.coordinate_from_header:
            self.logger.debug("Marking target pointing in jpeg.")
            markSize = (self.tel.pointing_marker_size.to(u.arcsec)/self.tel.pixel_scale).value/binning
            ## Mark Central Pixel with a White Cross
            JPEGcommand.append("-stroke")
            JPEGcommand.append("white")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("3")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            pixelCenter = [self.nXPix/2/binning, self.nYPix/2/binning]
            self.logger.debug("Marking central pixel of JPEG: {:.1f},{:.1f}".format(pixelCenter[0], pixelCenter[1]))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0], pixelCenter[1]+markSize,
                               pixelCenter[0], pixelCenter[1]+markSize*0.3))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0]+markSize, pixelCenter[1],
                               pixelCenter[0]+markSize*0.3, pixelCenter[1]))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0], pixelCenter[1]-markSize,
                               pixelCenter[0], pixelCenter[1]-markSize*0.3))
            JPEGcommand.append('-draw')
            JPEGcommand.append("line %d,%d %d,%d" % (pixelCenter[0]-markSize, pixelCenter[1],
                               pixelCenter[0]-markSize*0.3, pixelCenter[1]))
            ## Mark Coordinates of Target with a blue Circle
            JPEGcommand.append("-stroke")
            JPEGcommand.append("blue")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("3")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            ## This next block of code seems to make the call to wcs_world2pix
            ## happy, but I'm not sure I understand why.
            foo = np.array([[self.coordinate_from_header.ra.hour*15., self.coordinate_from_header.dec.radian*180./math.pi], 
                            [self.coordinate_from_header.ra.hour*15., self.coordinate_from_header.dec.radian*180./math.pi]])
            targetPixel = (self.image_WCS.wcs_world2pix(foo, 1)[0])
            self.logger.debug("Pixel of target on raw image: {:.1f},{:.1f}".format(targetPixel[0], targetPixel[1]))
            ## Adjust target pixel value for different origin in ImageMagick
            TargetXPos = targetPixel[0]
            if not self.cropped:
                TargetYPos = self.nYPix - targetPixel[1]
            else:
                TargetYPos = self.original_nYPix - targetPixel[1]
            ## Adjust target pixel value for cropping
            if self.cropped:
                TargetXPos = TargetXPos - self.crop_x1
                TargetYPos = TargetYPos - self.crop_y1
            ## Adjust target pixel value for binning
            TargetXPos = TargetXPos/binning
            TargetYPos = TargetYPos/binning
            self.logger.debug("Marking pixel of target on JPEG: {:.1f},{:.1f}".format(TargetXPos, TargetYPos))
            JPEGcommand.append('-draw')
            JPEGcommand.append("circle %d,%d %d,%d" % (TargetXPos, TargetYPos,
                               TargetXPos+markSize/2, TargetYPos))
        ## Mark Stars Detected by SExtractor
        if markDetectedStars and self.SExtractor_results:
            self.logger.debug("Marking stars found by SExtractor in jpeg.")
            JPEGcommand.append("-stroke")
            JPEGcommand.append("red")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("1")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            if self.FWHM:
                MarkRadius=max([4, 2*math.ceil(self.FWHM.value)])
            else:
                MarkRadius = 4
            sortedSExtractor_results = np.sort(self.SExtractor_results, order=['MAG_AUTO'])
            for star in sortedSExtractor_results:
                nStarsMarked += 1
                if nStarsMarked <= nStarsLimit:
                    MarkXPos = star['XWIN_IMAGE']
                    MarkYPos = self.nYPix - star['YWIN_IMAGE']
                    JPEGcommand.append('-draw')
                    JPEGcommand.append("circle %d,%d %d,%d" % (MarkXPos, MarkYPos, MarkXPos+MarkRadius, MarkYPos))
                else:
                    self.logger.info("  Only marked brightest {} stars found in image.".format(nStarsLimit))
                    break
        ## Rotate jpeg according to WCS (for images which have not had SWarp applied)
        if rotate and self.position_angle:
            self.logger.debug("Rotating jpeg by {0:.1f} deg".format(self.position_angle.to(u.deg).value))
            if self.position_angle:
                JPEGcommand.append("-rotate")
                JPEGcommand.append(str(self.position_angle.to(u.deg).value))
                if self.image_flipped:
                    JPEGcommand.append("-flop")
            else:
                self.logger.warning("No position angle value found.  Not rotating JPEG.")
        ## Write label describing marking of pointing
        if markPointing and self.image_WCS and self.coordinate_from_header:
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            JPEGcommand.append("text 200,40 'Blue circle centered on target is {:.1f} arcmin diameter.'".format(self.tel.pointing_marker_size.to(u.arcmin).value))
        ## Use background subtracted image generated by SExtractor
        if not backgroundSubtracted:
            JPEGcommand.append(self.working_file)
        else:
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            JPEGcommand.append("text 200,120 'Background Subtracted Image'")
            JPEGcommand.append(self.check_image_file)
        if markDetectedStars and nStarsMarked > nStarsLimit:
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            JPEGcommand.append("text 200,80 'Marked {} brightest stars out of {}.'".format(nStarsLimit, self.n_stars_SExtracted))
        JPEGcommand.append(jpegFile)
        self.logger.debug("Issuing convert command to create jpeg.")
        try:
            ConvertSTDOUT = subprocess.check_output(JPEGcommand, stderr=subprocess.STDOUT, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            self.logger.error("Failed to create jpeg.")
            for line in e.output.split("\n"):
                self.logger.error(line)
        except OSError as e:
            self.logger.error("Failed to create jpeg.")
            for line in e.strerror.split("\n"):
                self.logger.error(line)
        except:
            self.logger.error("Convert process failed: {0} {1} {2}".format(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
        else:
            for line in ConvertSTDOUT.split("\n"):
                if len(line) > 0:
                    self.logger.debug(line)
            self.jpeg_file_names.append(jpegFileName)


    ##-------------------------------------------------------------------------
    ## Clean Up by Deleting Temporary Files
    ##-------------------------------------------------------------------------
    def CleanUp(self):
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
    def AddWebLogEntry(self, htmlImageList, fields=None):
        '''
        This function adds one line to the HTML table of images.  The line
        contains the image info extracted by IQMon.
        '''
        if not fields: fields=["Date and Time", "Filename", "Alt", "Az", "Airmass", "moon_sep", "MoonIllum", "FWHM", "ellipticity", "Background", "PErr", "PosAng", "ZeroPoint", "nStars", "total_process_time"]
        ## If HTML file does not yet exist, create it and insert header
        ## from template file.
        self.logger.info('Adding results to HTML table.')
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
                      '    <h2>IQMon Results for {}</h2>'.format(self.tel.long_name),
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
            if "moon_sep" in fields:
                header.append('        <th style="width:50px">Moon Sep (deg)</th>')
            if "MoonIllum" in fields:
                header.append('        <th style="width:50px">Moon Illum. (%)</th>')
            if "FWHM" in fields:
                header.append('        <th style="width:60px">FWHM ({})</th>'.format(str(self.tel.units_for_FWHM.unit)))
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
            if "total_process_time" in fields:
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
                IsEndTable = re.match("\s*</table>\s*", line)
                IsEndBody = re.match("\s*</body>\s*", line)
                IsEndHTML = re.match("\s*</html>\s*", line)
                if not IsEndTable and not IsEndBody and not IsEndHTML:
                    HTML.write(line+"\n")
        ## Write Lines for this Image to HTML File
        HTML = open(htmlImageList, 'a')
        HTML.write("    <tr>\n")
        ## Write Observation Date and Time
        if "Date and Time" in fields:
            HTML.write("      <td style='color:black;text-align:left'>{0}</td>\n".format(self.dateObs))
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
            if self.PSF_plotfile:
                PSFplot_html = " (<a href='{}'>PSF</a>)".format(os.path.join("..", "..", "Plots", self.PSF_plot_filename))
            else:
                PSFplot_html = ""
            if self.zeroPoint_plotfile:
                ZPplot_html = " (<a href='{}'>ZP</a>)".format(os.path.join("..", "..", "Plots", self.zeroPoint_plotfile))
            else:
                ZPplot_html = ""
            htmlline = "      <td style='color:black;text-align:left'>" + JPEG1_html + "{}</a>".format(self.raw_file_basename) + JPEG2_html + JPEG3_html + PSFplot_html + ZPplot_html + "</td>\n"
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
        if "moon_sep" in fields:
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
                if self.FWHM > self.tel.threshold_FWHM.to(u.pix, equivalencies=self.tel.pixel_scale_equivalency):
                    colorFWHM = "#FF5C33"
                else:
                    colorFWHM = "#70DB70"
                ## Convert FWHM value to appropriate units for HTML output
                if self.tel.units_for_FWHM.unit == u.arcsec:
                    FWHM_for_HTML = (self.FWHM * u.radian.to(u.arcsec)*self.tel.pixel_size.to(u.mm)/self.tel.focal_length.to(u.mm)).value
                else:
                    FWHM_for_HTML = self.FWHM.value
                HTML.write("      <td style='background-color:{0}'>{1:.2f}</td>\n".format(colorFWHM, FWHM_for_HTML))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("#FF5C33", ""))
        if "ellipticity" in fields:
            if self.ellipticity:
                ## Decide whether to flag ellipticity value with red color
                if self.ellipticity > self.tel.threshold_ellipticity:
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
                if self.pointing_error.arcminute > self.tel.threshold_pointing_err.to(u.arcmin).value:
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
            if self.zeroPoint:
                HTML.write("      <td style='color:{}'>{:.2f}</td>\n".format("black", self.zeroPoint))
            else:
                HTML.write("      <td style='color:{}'>{}</td>\n".format("black", ""))
        ## Write number of stars detected by SExtractor
        if "nStars" in fields:
            if self.n_stars_SExtracted:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", self.n_stars_SExtracted))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Write process time
        if "total_process_time" in fields:
            if self.total_process_time:
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.total_process_time))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Complete Table
        HTML.write("    </tr>\n")
        HTML.write("  </table>\n")
        HTML.write("</body>\n")
        HTML.write("</html>\n")
        HTML.close()


    ##-------------------------------------------------------------------------
    ## Append Line With Image Info to Summary Text File
    ##-------------------------------------------------------------------------
    def AddSummaryEntry(self, summaryFile):
        self.logger.info("Writing Summary File Entry.")
        self.logger.debug("  Summary File: {0}".format(summaryFile))
        ## Read in previous data
        if not os.path.exists(summaryFile):
            self.logger.info("  Making new astropy table object")
            SummaryTable = table.Table(names=("ExpStart", "File", "FWHM (pix)", "Ellipticity", 
                                       "Alt (deg)", "Az (deg)", "Airmass", "pointing_error (arcmin)", 
                                       "ZeroPoint", "nStars", "Background", "Background RMS"),
                                 dtype=('S22', 'S100', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4', 'f4', 'f4'),
                                 masked=True)
        else:
            self.logger.info("  Reading astropy table object from file: {0}".format(summaryFile))
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
                self.logger.critical("Failed to read summary file: {0} {1} {2}".format(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
        ## Astropy table writer can not write None to table initialized
        ## with type.  If any outputs are None, change to some value.
        tableMask = np.zeros(12)
        ## dateObs
        if self.dateObs: dateObs = self.dateObs
        else: 
            dateObs = ""
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
        if self.zeroPoint: zeroPoint = self.zeroPoint
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
        SummaryTable.add_row((dateObs, raw_file_name,
                              FWHM, ellipticity,
                              target_alt, target_az,
                              airmass, pointing_error,
                              zeroPoint, n_stars_SExtracted,
                              SExtractor_background, SExtractor_background_RMS),
                              mask=tableMask)
        ## Write Table to File
        self.logger.info("  Writing new summary file.")
        ascii.write(SummaryTable, summaryFile,
                    Writer=ascii.basic.Basic)


    ##-------------------------------------------------------------------------
    ## Calcualte Process Time
    ##-------------------------------------------------------------------------
    def Calculatetotal_process_time(self):
        self.end_process_time = datetime.datetime.now()
        self.total_process_time = self.end_process_time - self.start_process_time
        self.logger.info("IQMon processing time = {0:.1f} seconds".format(self.total_process_time))

    