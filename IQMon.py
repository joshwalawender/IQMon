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
import time
import subprocess
import logging
import math
import numpy as np

## Import Astronomy Specific Tools
import ephem
import astropy.units as u
from astropy.io import fits
import astropy.coordinates as coords
from astropy import table
from astropy import wcs
from astropy.io import ascii


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
#             IsCatalogPath = re.match("CATALOGPATH\s=\s([\w/\-\.]+)", line)
#             if IsCatalogPath:
#                 self.pathCatalog = os.path.abspath(IsCatalogPath.group(1))

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
      longName:      A long name for the telescope
      focalLength:   The focal length of the telescope
      pixelSize:     The size of the pixels in the CCD camera
      pixelScale:    The estimated pixel scale of the telescope in arcsec per
                     pixel.  This should be an astropy Unit object with a value
                     and units.
      fRatio:        The F/ratio of the telescope.
      aperture:      The aperture of the telescope.
      gain:          The gain (in electrons per ADU) of the CCD.
      nXPix:         The number of pixels in the X direction of the CCD.
      nYPix:         The number of pixels in the Y direction of the CCD.
      site:          A pyephem site object defining the site of the telescope
      unitsForFWHM:  A quantity with the units for the FWHM output.  The units
                     should be reducible to either pixels or arcseconds.
      ROI:           The "region of interest" to crop the image to.  Format
                     should be '[x1:x2,y1:y2]' or 'x1:x2,y1:y2'.
      thresholdFWHM: The FWHM above which the image should be flagged in the
                     HTML output log.
      thresholdPointingErr:   
      thresholdEllipticity:   
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
        self.longName = None
        self.focalLength = None
        self.pixelSize = None
        self.aperture = None
        self.gain = None
        self.nXPix = None
        self.nYPix = None
        self.unitsForFWHM = None
        self.ROI = None
        self.thresholdFWHM = None
        self.thresholdPointingErr = None
        self.thresholdEllipticity = None
        self.pixelScale = None
        self.fRatio = None
        self.SExtractorPhotAperture = None
        self.SExtractorSeeing = None
        self.SExtractorSaturation = None
        self.site = None
        self.pointingMarkerSize = 1*u.arcmin
        
    def CheckUnits(self):
        '''
        Checks whether the telescope properties have the right type.  If a unit
        is expected, checks whether the input has units and whether it is
        reducible to the expected unit.  If input has no units and they are
        expected, then adds the default unit.
        '''
        ## name is a string
        assert type(self.name) == str
        ## longName is a string
        assert type(self.longName) == str
        ## Default focalLength units to mm
        if type(self.focalLength) == u.quantity.Quantity:
            assert self.focalLength.to(u.mm)
        else:
            self.focalLength *= u.mm
        ## Default pixelSize units to microns
        if type(self.pixelSize) == u.quantity.Quantity:
            assert self.pixelSize.to(u.micron)
        else:
            self.pixelSize *= u.micron
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
        ## Default SExtractorSaturation to units of ADU
        if type(self.SExtractorSaturation) == u.quantity.Quantity:
            assert self.SExtractorSaturation.to(u.adu)
        else:
            self.SExtractorSaturation *= u.adu
        ## Default unitsForFWHM to units of arcsec
        if type(self.unitsForFWHM) == u.quantity.Quantity:
            assert self.unitsForFWHM.unit in [u.arcsec, u.pix]
        else:
            self.unitsForFWHM *= u.pix
        ## ROI is string
        assert type(self.ROI) == str
        ## Default thresholdFWHM to same units as unitsForFWHM
        if type(self.thresholdFWHM) == u.quantity.Quantity:
            assert self.thresholdFWHM.unit in [u.arcsec, u.pix]
        else:
            self.thresholdFWHM *= u.pix
        ## Default thresholdPointingErr to units of arcmin
        if type(self.thresholdPointingErr) == u.quantity.Quantity:
            assert self.thresholdPointingErr.to(u.arcmin)
        else:
            self.thresholdPointingErr *= u.arcmin
        ## Default thresholdEllipticity to dimensionless
        if type(self.thresholdEllipticity) == u.quantity.Quantity:
            assert self.thresholdEllipticity.to(u.dimensionless_unscaled)
        else:
            assert float(self.thresholdEllipticity) >=0 and float(self.thresholdEllipticity) <= 1.
        ## Default pixelScale to units of arcsec per pixel
        if type(self.pixelScale) == u.quantity.Quantity:
            assert self.pixelScale.to(u.arcsec / u.pix)
        else:
            self.pixelScale *= u.arcsec / u.pix
        ## Default fRatio to dimensionless
        if type(self.fRatio) == u.quantity.Quantity:
            assert self.fRatio.to(u.dimensionless_unscaled)
        else:
            assert float(self.fRatio)


    ##-------------------------------------------------------------------------
    ## Define astropy.units Equivalency for Arcseconds and Pixels
    ##-------------------------------------------------------------------------
    def DefinePixelScale(self):
        self.pixelScaleEquivalency = [(u.pix, u.arcsec,
                       lambda pix: (pix*u.radian.to(u.arcsec)*self.pixelSize/self.focalLength).decompose().value,
                       lambda arcsec: (arcsec/u.radian.to(u.arcsec)*self.focalLength/self.pixelSize).decompose().value
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
        self.startProcessTime = time.time()
        if os.path.exists(input):
            FitsFileDirectory, FitsFilename = os.path.split(input)
            self.rawFile = input
            self.rawFileName = FitsFilename
            self.rawFileDirectory = FitsFileDirectory
            self.rawFileBasename, self.fileExt = os.path.splitext(FitsFilename)
        else:
            self.rawFile = None
            self.rawFileName = None
            self.rawFileDirectory = None
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
        self.workingFile = None
        self.header = None
        self.exptime = None
        self.filter = None
        self.focusPos = None
        self.objectName = None
        self.astrometrySolved = None
        self.coordinate_WCS = None
        self.coordinate_header = None
        self.nSExtracted = None
        self.SExBackground = None
        self.SExBRMS = None
        self.tempFiles = []
        self.SExtractorResults = None
        self.nStarsSEx = None
        self.positionAngle = None
        self.zeroPoint = None
        self.processTime = None
        self.FWHM = None
        self.ellipticity = None
        self.pointingError = None
        self.imageFlipped = None
        self.jpegFileNames = []
        self.CheckImageFile = None
        self.cropped = False
        self.crop_x1 = None
        self.crop_x2 = None
        self.crop_y1 = None
        self.crop_y2 = None
        self.original_nXPix = None
        self.original_nXPix = None

    ##-------------------------------------------------------------------------
    ## Make Logger Object
    ##-------------------------------------------------------------------------
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
    ## Get Header
    ##-------------------------------------------------------------------------
    def GetHeader(self):
        '''
        Get information from the image fits header.
        '''
        hdulist = fits.open(self.workingFile, ignore_missing_end=True)
        self.header = hdulist[0].header
        self.image = hdulist[0].data
        hdulist.close()
        self.logger.info("Reading image header.")
        
        ## Get exposure time from header (assumes seconds)
        try:
            self.exptime = float(self.header['EXPTIME']) * u.s
        except:
            self.exptime = None
            self.logger.debug("No exposure time value found in header")
        else:
            self.logger.debug("Exposure time = {0:.1f} s".format(self.exptime.to(u.s).value))
        ## Get filter from header
        try:
            self.filter = self.header['FILTER']
        except:
            self.filter = None
            self.logger.debug("No filter keyword found in header")
        else:
            self.logger.debug("Filter = {0}".format(self.filter))
        ## Get focus position from header
        try:
            self.focusPos = self.header['FOCUSPOS']
        except:
            self.focusPos = None
            self.logger.debug("No focus position value found in header")
        else:
            self.logger.debug("Focus position = {0}".format(self.focusPos))
        ## Get object name from header
        try:
            self.objectName = self.header["OBJECT"]
        except:
            self.objectName = None
            self.logger.debug("No object value found in header")
        else:
            self.logger.debug("Header object name = {0}".format(self.objectName))
        ## Get airmass from header
        try:
            self.headerAirmass = self.header["AIRMASS"]
        except:
            self.headerAirmass = None
            self.logger.debug("No airmass value found in header")
        else:
            self.logger.debug("Header airmass = {0:.2f}".format(self.headerAirmass))
        ## Get Observation Date and Time from header
        ## (assumes YYYY-MM-DDTHH:MM:SS format)
        try:
            self.dateObs = self.header["DATE-OBS"]
        except:
            self.dateObs = None
            self.logger.debug("No date value found in header")
        else:
            self.logger.debug("Header date = {0}".format(self.dateObs))
        ## Get Site Latitude from header (assumes decimal degrees)
        try:
            self.latitude = self.header["LAT-OBS"] * u.deg
        except:
            self.latitude = None
            self.logger.debug("No latitude value found in header")
        else:
            self.logger.debug("Header latitude = {0:.4f} deg".format(self.latitude.to(u.deg).value))
        ## Get Site Longitude from header (assumes decimal degrees)
        try:
            self.longitude = self.header["LONG-OBS"] * u.deg
        except:
            self.longitude = None
            self.logger.debug("No longitiude value found in header")
        else:
            self.logger.debug("Header longitiude = {0:.4f} deg".format(self.longitude.to(u.deg).value))
        ## Get Site Altitude from header (assumes meters)
        try:
            self.altitude = self.header["ALT-OBS"] * u.meter
        except:
            self.altitude = None
            self.logger.debug("No altitude value found in header")
        else:
            self.logger.debug("Header altitude = {0:.0f} meters".format(self.altitude.to(u.meter).value))


        ## Determine Image Size in Pixels
        self.nYPix, self.nXPix = self.image.shape
        self.logger.debug('Image size is: {},{}'.format(self.nXPix, self.nYPix))

        ## Read Header Coordinates in to astropy coordinates object
        ImageRA  = self.header['RA']
        if len(ImageRA.split(":")) != 3:
            if len(ImageRA.split(" ")) == 3:
                ImageRA = ":".join(ImageRA.split(" "))
        ImageDEC = self.header['DEC']    
        if len(ImageDEC.split(":")) != 3:
            if len(ImageDEC.split(" ")) == 3:
                ImageDEC = ":".join(ImageDEC.split(" "))
        self.logger.debug("Read pointing info from header: "+ImageRA+" "+ImageDEC)
        try:
            self.coordinate_header = coords.ICRS(ImageRA+" "+ImageDEC,
                                                 unit=(u.hour, u.degree))
        except:
            self.logger.warning("Failed to read pointing info from header.")
            self.coordinate_header = None

        ## Read WCS
        try:
            self.imageWCS = wcs.WCS(self.header)
        except:
            self.imageWCS = None
            self.logger.info("No WCS found in image header")
        else:
            self.logger.debug("Found WCS in image header.")

        ## Determine PA of Image
        try:
            PC11 = float(self.imageWCS.to_header()['PC1_1'])
            PC12 = float(self.imageWCS.to_header()['PC1_2'])
            PC21 = float(self.imageWCS.to_header()['PC2_1'])
            PC22 = float(self.imageWCS.to_header()['PC2_2'])
        except:
            self.logger.debug("Could not find PCn_m values in WCS.")
            try:
                PC11 = float(self.imageWCS.to_header()['CD1_1'])
                PC12 = float(self.imageWCS.to_header()['CD1_2'])
                PC21 = float(self.imageWCS.to_header()['CD2_1'])
                PC22 = float(self.imageWCS.to_header()['CD2_2'])
            except:
                self.logger.debug("Could not find CDn_m values in WCS.")
                self.imageWCS = None
        if self.imageWCS:
            if (abs(PC21) > abs(PC22)) and (PC21 >= 0): 
                North = "Right"
                self.positionAngle = 270.*u.deg + math.degrees(math.atan(PC22/PC21))*u.deg
            elif (abs(PC21) > abs(PC22)) and (PC21 < 0):
                North = "Left"
                self.positionAngle = 90.*u.deg + math.degrees(math.atan(PC22/PC21))*u.deg
            elif (abs(PC21) < abs(PC22)) and (PC22 >= 0):
                North = "Up"
                self.positionAngle = 0.*u.deg + math.degrees(math.atan(PC21/PC22))*u.deg
            elif (abs(PC21) < abs(PC22)) and (PC22 < 0):
                North = "Down"
                self.positionAngle = 180.*u.deg + math.degrees(math.atan(PC21/PC22))*u.deg
            if (abs(PC11) > abs(PC12)) and (PC11 > 0): East = "Right"
            if (abs(PC11) > abs(PC12)) and (PC11 < 0): East = "Left"
            if (abs(PC11) < abs(PC12)) and (PC12 > 0): East = "Up"
            if (abs(PC11) < abs(PC12)) and (PC12 < 0): East = "Down"
            if North == "Up" and East == "Left": self.imageFlipped = False
            if North == "Up" and East == "Right": self.imageFlipped = True
            if North == "Down" and East == "Left": self.imageFlipped = True
            if North == "Down" and East == "Right": self.imageFlipped = False
            if North == "Right" and East == "Up": self.imageFlipped = False
            if North == "Right" and East == "Down": self.imageFlipped = True
            if North == "Left" and East == "Up": self.imageFlipped = True
            if North == "Left" and East == "Down": self.imageFlipped = False
            self.logger.debug("Position angle of WCS is {0:.1f} degrees.".format(self.positionAngle.to(u.deg).value))
            self.logger.debug("Image orientation is North {0}, East {1}.".format(North, East))
            if self.imageFlipped:
                self.logger.debug("Image is mirrored.")
        else:
            self.positionAngle = None
            self.imageFlipped = None


        ## Determine Alt, Az, Moon Sep, Moon Illum using ephem module
        if self.dateObs and self.latitude and self.longitude:
            ## Populate site object properties
            SiteDate = "/".join(self.dateObs[0:10].split("-"))
            SiteTime = self.dateObs[11:]        
            self.tel.site.date = ephem.Date(SiteDate+" "+SiteTime)
            self.tel.site.lat = str(self.latitude.to(u.deg).value)
            self.tel.site.lon = str(self.longitude.to(u.deg).value)
            if self.altitude: self.tel.site.elevation = self.altitude.to(u.meter).value
            ## Do calculations using ephem
            TargetObject = ephem.readdb("Target,f|M|F7,"+ImageRA+","+ImageDEC+",2.02,2000")
            TargetObject.compute(self.tel.site)
            self.targetAlt = TargetObject.alt * 180./ephem.pi * u.deg
            self.targetAz = TargetObject.az * 180./ephem.pi * u.deg
            self.logger.debug("Target Alt, Az = {0:.1f}, {1:.1f}".format(self.targetAlt.to(u.deg).value, self.targetAz.to(u.deg).value))
            self.zenithAngle = 90.*u.deg - self.targetAlt
            self.airmass = 1.0/math.cos(self.zenithAngle.to(u.radian).value)*(1.0 - 0.0012*(1.0/(math.cos(self.zenithAngle.to(u.radian).value)**2 - 1.0)))
            self.logger.debug("Target airmass (calculated) = {0:.2f}".format(self.airmass))
            ## Calculate Moon Position and Illumination
            TheMoon = ephem.Moon()
            TheMoon.compute(self.tel.site)
            self.moonPhase = TheMoon.phase
            self.moonSep = ephem.separation(TargetObject, TheMoon) * 180./ephem.pi * u.deg
            self.moonAlt = TheMoon.alt * 180./ephem.pi * u.deg
            self.logger.debug("A {0:.0f} percent illuminated Moon is {1:.0f} deg from target.".format(self.moonPhase, self.moonSep.to(u.deg).value))
        else:
            self.targetAlt = None
            self.targetAz = None
            self.moonPhase = None
            self.moonSep = None
            self.moonAlt = None
            self.zenithAngle = None
            self.airmass = None
            self.logger.warning("Object position and Moon position not calculated.")


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
        self.workingFile = os.path.join(self.config.pathTemp, self.rawFileName)
        shutil.copy2(self.rawFile, self.workingFile)
        self.tempFiles.append(self.workingFile)

    ##-------------------------------------------------------------------------
    ## Dark Subtract Image
    ##-------------------------------------------------------------------------
    def DarkSubtract(self, Darks):
        '''
        Create master dark and subtract from image.
        
        Input the filename of the appropriate master dark.  May want to write
        own function to make the master dark given input file data.
        '''
        self.logger.debug("Dark subtracting image.  Opening image data.")
        hdulist_image = fits.open(self.workingFile, mode='update')
        ## Load master dark if provided, but if multiple files input, combine
        ## them in to master dark, then load combined master dark.
        if len(Darks) == 1:
            self.logger.debug("Found master dark.  Opening master dark data.")
            hdulist_dark = fits.open(Darks[0])
            MasterDarkData = hdulist_dark[0].data
        elif len(Darks) > 1:
            self.logger.info("Multiple input darks detected.  Median combining {0} darks.".format(len(Darks)))
            ## Combine multiple darks frames
            DarkData = []
            for Dark in Darks:
                hdulist = fits.open(Dark)
                DarkData.append(hdulist[0].data)
            DarkData = np.array(DarkData)
            MasterDarkData = np.median(DarkData, axis=0)
            ## Save Master Dark to Fits File
            DataPath = os.path.split(self.rawFile)[0]
            DataNightString = os.path.split(DataPath)[1]
            MasterDarkFilename = "MasterDark_"+self.tel.name+"_"+DataNightString+"_"+str(int(math.floor(self.exptime.to(u.s).value)))+".fits"
            MasterDarkFile  = os.path.join(self.config.pathTemp, MasterDarkFilename)    
            hdu_MasterDark = fits.PrimaryHDU(MasterDarkData)
            hdulist_MasterDark = fits.HDUList([hdu_MasterDark])
            hdulist_MasterDark.header = hdulist[0].header
            hdulist_MasterDark.header['history'] = "Combined {0} images to make this master dark.".format(len(Darks))
            self.logger.info("Writing master dark file: {0}".format(MasterDarkFile))
            hdulist_MasterDark.writeto(MasterDarkFile)
        else:
            self.logger.error("No input dark files detected.")
        ## Now Subtract MasterDark from Image
        self.logger.info("Subtracting dark from image.")
        ImageData = hdulist_image[0].data
        DifferenceImage = ImageData - MasterDarkData
        hdulist_image[0].data = DifferenceImage
        hdulist_image.flush()
#         self.logger.debug("Median level of image = {0}".format(np.median(ImageData)))
#         self.logger.debug("Median level of dark = {0}".format(np.median(MasterDarkData)))
#         self.logger.debug("Median level of dark subtracted = {0}".format(np.median(DifferenceImage)))


    ##-------------------------------------------------------------------------
    ## Crop Image
    ##-------------------------------------------------------------------------
    def Crop(self):
        '''
        Crop working image to region of interest.
        '''
        if self.tel.ROI:
            ## Parse ROI String
            try:
                MatchROI = re.match("\[?(\d{1,5}):(\d{1,5}),(\d{1,5}):(\d{1,5})\]?", self.tel.ROI)
            except:
                self.logger.warning("Could not parse ROI string in telescope object.")
            else:
                self.crop_x1 = int(MatchROI.group(1))
                self.crop_x2 = int(MatchROI.group(2))
                self.crop_y1 = int(MatchROI.group(3))
                self.crop_y2 = int(MatchROI.group(4))
                self.logger.info("Cropping Image To [{0}:{1},{2}:{3}]".format(self.crop_x1, self.crop_x2, self.crop_y1, self.crop_y2))
                hdulist = fits.open(self.workingFile, mode="update")
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
        self.logger.info("Attempting to create WCS using Astrometry.net solver.")
        AstrometryCommand = ["solve-field", "-l", "5", "-O", "-p",
                             "-L", str(self.tel.pixelScale.value*0.90),
                             "-H", str(self.tel.pixelScale.value*1.10),
                             "-u", "arcsecperpix", "-z", "4", self.workingFile]
        AstrometrySTDOUT = ""

        try:
            StartTime = time.time()
            AstrometrySTDOUT = subprocess.check_output(AstrometryCommand, 
                               stderr=subprocess.STDOUT)
            EndTime = time.time()
        except subprocess.CalledProcessError as e:
            self.logger.warning("Astrometry.net failed.")
            for line in e.output.split("\n"):
                self.logger.error(line)
            self.astrometrySolved = False
        except:
            self.logger.error("solve-field process failed: {0} {1} {2}".format(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
        else:
            ProcessTime = EndTime - StartTime
            self.logger.debug("Astrometry.net Processing Time: %.1f s", ProcessTime)
            pos = AstrometrySTDOUT.find("Field center: (RA H:M:S, Dec D:M:S) = ")
            if pos != -1:
                IsFieldCenter = re.match("\s*(\d{1,2}:\d{2}:\d{2}\.\d+,\s-?\d{1,2}:\d{2}:\d{2}\.\d+).*", 
                                         AstrometrySTDOUT[pos+40:pos+75])
                if IsFieldCenter:
                    self.logger.info("Astrometry.net field center is: %s", IsFieldCenter.group(1))
                else:
                    self.logger.warning("Could not parse field center from astrometry.net output.")
                    for line in AstrometrySTDOUT.split("\n"):
                        self.logger.warning("  %s" % line)
            else:
                for line in AstrometrySTDOUT.split("\n"):
                    self.logger.warning("  %s" % line)
            NewFile = self.workingFile.replace(self.fileExt, ".new")
            NewFitsFile = self.workingFile.replace(self.fileExt, ".new.fits")
            if not os.path.exists(NewFile):
                self.logger.warning("No new file created by astrometry.net")
                self.astrometrySolved = False
            else:
                self.logger.debug("Astrometry.net succeeded")
                if os.path.exists(NewFitsFile): os.remove(NewFitsFile)
                os.rename(NewFile, NewFitsFile)
                self.astrometrySolved = True
                self.workingFile = NewFitsFile
                ## Update header history
#                 hdulist = fits.open(self.workingFile, mode="update", ignore_missing_end=True)
#                 now = time.gmtime()
#                 hdulist[0].header['history'] = "Solved by Astrometry.net at {0}".format(time.strftime("%Y-%m-%dT%H:%M:%S UTC"))
#                 hdulist.flush()
#                 hdulist.close()
            ## Add files created by astrometry.net to tempFiles list
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+".axy"))
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+".wcs"))
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+".solved"))
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+".rdls"))
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+".match"))
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+".corr"))
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+".new.fits"))
            self.tempFiles.append(os.path.join(self.config.pathTemp, self.rawFileBasename+"-indx.xyls"))

    ##-------------------------------------------------------------------------
    ## Refine WCS
    ##-------------------------------------------------------------------------
    def RefineWCS(self):
        '''
        Refine the WCS of the image to have accurate distortions.
        '''
        pass

    ##-------------------------------------------------------------------------
    ## Determine Pointing Error
    ##-------------------------------------------------------------------------
    def DeterminePointingError(self):
        '''
        Determine pointing error (difference between objects coordinates and
        solved WCS).
        '''
        self.logger.info("Detemining pointing error based on WCS solution")
        if self.imageWCS and self.coordinate_header:
            centerWCS = self.imageWCS.wcs_pix2world([[self.nXPix/2, self.nYPix/2]], 1)
            self.logger.debug("Using coordinates of center point: {0} {1}".format(centerWCS[0][0], centerWCS[0][1]))
            self.coordinate_WCS = coords.ICRS(ra=centerWCS[0][0],
                                                   dec=centerWCS[0][1],
                                                   unit=(u.degree, u.degree))
            self.pointingError = self.coordinate_WCS.separation(self.coordinate_header)
            self.logger.debug("Target Coordinates are:  %s %s",
                         self.coordinate_header.ra.format(u.hour, sep=":", precision=1),
                         self.coordinate_header.dec.format(u.degree, sep=":", precision=1, alwayssign=True))
            self.logger.debug("WCS of Central Pixel is: %s %s",
                         self.coordinate_WCS.ra.format(u.hour, sep=":", precision=1),
                         self.coordinate_WCS.dec.format(u.degree, sep=":", precision=1, alwayssign=True))
            self.logger.info("Pointing Error is %.2f arcmin", self.pointingError.arcminute)
        else:
            self.logger.warning("Pointing error not calculated.")

    ##-------------------------------------------------------------------------
    ## Run SExtractor
    ##-------------------------------------------------------------------------
    def RunSExtractor(self):
        '''
        Run SExtractor on image.
        '''
        assert type(self.tel.gain) == u.quantity.Quantity
        assert type(self.tel.pixelScale) == u.quantity.Quantity
        assert type(self.tel.SExtractorSeeing) == u.quantity.Quantity
        assert type(self.tel.SExtractorPhotAperture) == u.quantity.Quantity
        if self.tel.gain and self.tel.pixelScale and self.tel.SExtractorSeeing and self.tel.SExtractorPhotAperture:
            ## Set up file names
            SExtractorConfigFile = os.path.join(self.config.pathTemp, self.rawFileBasename+".sex")
            self.tempFiles.append(SExtractorConfigFile)
            SExtractorCatalog = os.path.join(self.config.pathTemp, self.rawFileBasename+".cat")
            self.tempFiles.append(SExtractorCatalog)
            PhotometryCatalogFile_xy = os.path.join(self.config.pathTemp, self.rawFileBasename+"PhotCat_xy.txt")
            self.tempFiles.append(PhotometryCatalogFile_xy)
            CheckImageType = "-BACKGROUND"
            self.CheckImageFile = os.path.join(self.config.pathPlots, self.rawFileBasename+"_bksub.fits")
            self.tempFiles.append(self.CheckImageFile)

            ## Create PhotometryCatalogFile_xy file for SExtractor Association
            if os.path.exists(PhotometryCatalogFile_xy): os.remove(PhotometryCatalogFile_xy)
            PhotCatFileObject = open(PhotometryCatalogFile_xy, 'w')
            PhotCatFileObject.write("# No Existing WCS Found for this image\n")
            PhotCatFileObject.write("# This is a dummy file to keep SExtractor happy\n")
            PhotCatFileObject.write("0.0  0.0  0.0  0.0\n")
            PhotCatFileObject.close()
            
            ## Make edits To default.sex based on telescope:
            ## Read in default config file
            DefaultConfig = subprocess.check_output(["sex", "-dd"]).split("\n")
            NewConfig     = open(SExtractorConfigFile, 'w')
            backgroundFilterSize = max(5.*self.tel.SExtractorSeeing.to(u.arcsec).value / self.tel.pixelScale.value, 5.)
            self.logger.debug("Using background filter size of 5x seeing = {0:.1f} pixels.".format(backgroundFilterSize))
            for line in DefaultConfig:
                newline = line
                if re.match("CATALOG_NAME\s+", line):
                    newline = "CATALOG_NAME     "+SExtractorCatalog+"\n"
                if re.match("CATALOG_TYPE\s+", line):
                    newline = "CATALOG_TYPE     "+"FITS_LDAC"+"\n"
                if re.match("PARAMETERS_NAME\s+", line):
                    newline = "PARAMETERS_NAME  "+os.path.join(self.config.pathIQMonExec, "default.param")+"\n"
                if re.match("DETECT_MINAREA\s+", line) and (2.*self.tel.pixelScale.value > self.tel.SExtractorSeeing.to(u.arcsec).value):
                    newline = "DETECT_MINAREA   "+"4"+"\n"
                if re.match("DETECT_THRESH\s+", line):
                    newline = "DETECT_THRESH    "+"5.0"+"\n"
                if re.match("ANALYSIS_THRESH\s+", line):
                    newline = "ANALYSIS_THRESH  "+"5.0"+"\n"
                if re.match("FILTER\s+", line):
                    newline = "FILTER           "+"N"+"\n"
                if re.match("BACK_SIZE\s+", line):
                    newline = "BACK_SIZE        {0:.1f}\n".format(backgroundFilterSize)
                if re.match("ASSOC_NAME\s+", line):
                    newline = "ASSOC_NAME       "+PhotometryCatalogFile_xy+"\n"
                if re.match("ASSOC_NAME\s+", line):
                    newline = "ASSOC_NAME       "+PhotometryCatalogFile_xy+"\n"
                if re.match("ASSOCSELEC_TYPE\s+", line):
                    newline = "ASSOCSELEC_TYPE  "+"ALL"+"\n"
                if re.match("CHECKIMAGE_TYPE\s+", line):
                    newline = "CHECKIMAGE_TYPE  "+CheckImageType+"\n"
                if re.match("CHECKIMAGE_NAME\s+", line):
                    newline = "CHECKIMAGE_NAME  "+self.CheckImageFile+"\n"
                if re.match("PHOT_APERTURES\s+", line):
                    newline = "PHOT_APERTURES   "+str(self.tel.SExtractorPhotAperture.to(u.pix).value)+"\n"
                if re.match("GAIN\s+", line):
                    newline = "GAIN             "+str(self.tel.gain.value)+"\n"
                if re.match("PIXEL_SCALE\s+", line):
                    newline = "PIXEL_SCALE      "+str(self.tel.pixelScale.value)+"\n"
                if self.tel.SExtractorSaturation:
                    if re.match("SATUR_LEVEL\s+", line):
                        newline = "SATUR_LEVEL      "+str(self.tel.SExtractorSaturation.to(u.adu).value)+"\n"
                if re.match("SEEING_FWHM\s+", line):
                    newline = "SEEING_FWHM      "+str(self.tel.SExtractorSeeing.to(u.arcsec).value)+"\n"
                NewConfig.write(newline+"\n")
            NewConfig.close()

            ## Run SExtractor
            SExtractorCommand = ["sex", self.workingFile, "-c", SExtractorConfigFile]
            self.logger.info("Invoking SExtractor")
            self.logger.debug("SExtractor command: {}".format(repr(SExtractorCommand)))
            try:
                SExSTDOUT = subprocess.check_output(SExtractorCommand, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                self.logger.error("SExtractor failed.")
                for line in e.output.split("\n"):
                    self.logger.error(line)
            except:
                self.logger.error("SExtractor process failed: {0} {1} {2}".format(sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2]))
            else:
                for line in SExSTDOUT.split("\n"):
                    line.replace("[1A", "")
                    line.replace("[1M>", "")
                    if not re.match(".*Setting up background map.*", line) and not re.match(".*Line:\s[0-9]*.*", line):
                        self.logger.debug("  "+line)
                ## Extract Number of Stars from SExtractor Output
                pos = SExSTDOUT.find("sextracted ")
                IsSExCount = re.match("\s*([0-9]+)\s+", SExSTDOUT[pos+11:pos+21])
                if IsSExCount:
                    self.nSExtracted = int(IsSExCount.group(1))
                    self.logger.info("SExtractor found {0} sources.".format(self.nSExtracted))
                else:
                    self.nSExtracted = None
                ## Extract Background Level from SExtractor Output
                pos = SExSTDOUT.find("Background: ")
                IsSExBkgnd = re.match("\s*([0-9\.]+)\s*", SExSTDOUT[pos+11:pos+21])
                if IsSExBkgnd:
                    self.SExBackground = float(IsSExBkgnd.group(1))
                    self.logger.info("SExtractor background is {0:.1f}".format(self.SExBackground))
                else:
                    self.SExBackground = None
                ## Extract Background RMS from SExtractor Output
                IsSExBRMS = re.match("\s*RMS:\s([0-9\.]+)\s*", SExSTDOUT[pos+21:pos+37])
                if IsSExBRMS:
                    self.SExBRMS = float(IsSExBRMS.group(1))
                    self.logger.info("SExtractor background RMS is {0:.1f}".format(self.SExBRMS))
                else:
                    self.SExBRMS = None

                ## If No Output Catalog Created ...
                if not os.path.exists(SExtractorCatalog):
                    self.logger.warning("SExtractor failed to create catalog.")
                    self.SExtractorCatalog = None
                else:
                    self.SExtractorCatalog = SExtractorCatalog

                ## Read FITS_LDAC SExtractor Catalog
                self.logger.debug("Reading SExtractor output catalog.")
                hdu = fits.open(self.SExtractorCatalog)
                self.SExtractorResults = table.Table(hdu[2].data)
#                 self.SExtractorResults = ascii.read(self.SExtractorCatalog, Reader=ascii.sextractor.SExtractor)
                SExImageRadius = []
                SExAngleInImage = []
                for star in self.SExtractorResults:
                    SExImageRadius.append(math.sqrt((self.nXPix/2-star['X_IMAGE'])**2 + (self.nYPix/2-star['Y_IMAGE'])**2))
                    SExAngleInImage.append(math.atan((star['X_IMAGE']-self.nXPix/2)/(self.nYPix/2-star['Y_IMAGE']))*180.0/math.pi)
                self.SExtractorResults.add_column(table.Column(data=SExImageRadius, name='ImageRadius'))
                self.SExtractorResults.add_column(table.Column(data=SExAngleInImage, name='AngleInImage'))
                self.nStarsSEx = len(self.SExtractorResults)
                self.logger.info("Read in {0} stars from SExtractor catalog.".format(self.nStarsSEx))

        else:
            self.logger.warning("Telescope proerties not set.")


    ##-------------------------------------------------------------------------
    ## Determine Image FWHM from SExtractor Catalog
    ##-------------------------------------------------------------------------
    def DetermineFWHM(self):
        '''
        Determine typical FWHM of image from SExtractor results.
        '''
        if self.nStarsSEx > 1:
            IQRadiusFactor = 1.0
            DiagonalRadius = math.sqrt((self.nXPix/2)**2+(self.nYPix/2)**2)
            IQRadius = DiagonalRadius*IQRadiusFactor
            CentralFWHMs = [star['FWHM_IMAGE'] for star in self.SExtractorResults if star['ImageRadius'] <= IQRadius]
            CentralEllipticities = [star['ELLIPTICITY'] for star in self.SExtractorResults if star['ImageRadius'] <= IQRadius]
#             CentralFWHMs = []
#             CentralEllipticities = []
#             for star in self.SExtractorResults:
#                 if star['ImageRadius'] <= IQRadius:
#                     CentralFWHMs.append(star['FWHM_IMAGE'])
#                     CentralEllipticities.append(star['ELLIPTICITY'])   
            if len(CentralFWHMs) > 3:
                self.FWHM = np.median(CentralFWHMs) * u.pix
                self.ellipticity = np.median(CentralEllipticities)
            else:
                self.logger.warning("Not enough stars detected in central region of image to form median FWHM.")
            self.logger.debug("Using {0} stars in central region to determine FWHM and ellipticity.".format(len(CentralFWHMs)))
            self.logger.info("Median FWHM in inner region is {0:.2f} pixels".format(self.FWHM.to(u.pix).value))
            self.logger.info("Median Ellipticity in inner region is {0:.2f}".format(self.ellipticity))
        else:
            self.FWHM = None
            self.ellipticity = None


    ##-------------------------------------------------------------------------
    ## Determine Zero Point from SExtractor Catalog
    ##-------------------------------------------------------------------------
    def DetermineZeroPoint(self):
        '''
        Determine zero point by comparing measured magnitudes with catalog
        magnitudes.
        '''
        pass


    ##-------------------------------------------------------------------------
    ## Make JPEG of Image
    ##-------------------------------------------------------------------------
    def MakeJPEG(self, jpegFileName, markStars=False, markPointing=False, rotate=False, binning=1, backgroundSubtracted=False):
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
        if markPointing and self.imageWCS and self.coordinate_header:
            self.logger.debug("Marking target pointing in jpeg.")
            ## Make markSize 180 arcseconds (3 arcmin)
#             markSize = (180*u.arcsec/self.tel.pixelScale).value/binning
            markSize = (self.tel.pointingMarkerSize.to(u.arcsec)/self.tel.pixelScale).value/binning
            ## Mark Central Pixel with a White Cross
            JPEGcommand.append("-stroke")
            JPEGcommand.append("white")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("3")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            pixelCenter = [self.nXPix/2/binning, self.nYPix/2/binning]
            self.logger.debug("Marking central pixel of JPEG: {:.1f},{:.1f}".format(pixelCenter[0], pixelCenter[1]))
#             JPEGcommand.append('-draw')
#             JPEGcommand.append("circle %d,%d %d,%d" % (pixelCenter[0], pixelCenter[1],
#                                pixelCenter[0]+markSize, pixelCenter[1]))
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
            ## Mark WCS of Target with a blue Circle
            JPEGcommand.append("-stroke")
            JPEGcommand.append("blue")
            JPEGcommand.append("-strokewidth")
            JPEGcommand.append("3")
            JPEGcommand.append("-fill")
            JPEGcommand.append("none")
            ## This next block of code seems to make the call to wcs_world2pix
            ## happy, but I'm not sure I understand why.
            foo = np.array([[self.coordinate_header.ra.hour*15., self.coordinate_header.dec.radian*180./math.pi], 
                            [self.coordinate_header.ra.hour*15., self.coordinate_header.dec.radian*180./math.pi]])
            targetPixel = (self.imageWCS.wcs_world2pix(foo, 1)[0])
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
        if markStars and self.SExtractorResults:
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
            sortedSExtractorResults = np.sort(self.SExtractorResults, order=['MAG_AUTO'])
            for star in sortedSExtractorResults:
                nStarsMarked += 1
                if nStarsMarked <= nStarsLimit:
                    MarkXPos = star['X_IMAGE']
                    MarkYPos = self.nYPix - star['Y_IMAGE']
                    JPEGcommand.append('-draw')
                    JPEGcommand.append("circle %d,%d %d,%d" % (MarkXPos, MarkYPos, MarkXPos+MarkRadius, MarkYPos))
                else:
                    self.logger.warning("Only marked brightest {} stars found in image.".format(nStarsLimit))
                    break
        if rotate and self.positionAngle:
            self.logger.debug("Rotating jpeg by {0:.1f} deg".format(self.positionAngle.to(u.deg).value))
            if self.positionAngle:
                JPEGcommand.append("-rotate")
                JPEGcommand.append(str(self.positionAngle.to(u.deg).value))
                if self.imageFlipped:
                    JPEGcommand.append("-flop")
            else:
                self.logger.warning("No position angle value found.  Not rotating JPEG.")
        if markPointing and self.imageWCS and self.coordinate_header:
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            JPEGcommand.append("text 200,40 'Blue circle centered on target is {:.1f} arcmin diameter.'".format(self.tel.pointingMarkerSize.to(u.arcmin).value))
        if not backgroundSubtracted:
            JPEGcommand.append(self.workingFile)
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
            JPEGcommand.append(self.CheckImageFile)
        if markStars and nStarsMarked > nStarsLimit:
            JPEGcommand.append("-stroke")
            JPEGcommand.append("none")
            JPEGcommand.append("-fill")
            JPEGcommand.append("white")
            JPEGcommand.append("-pointsize")
            JPEGcommand.append("28")
            JPEGcommand.append('-font')
            JPEGcommand.append('fixed')
            JPEGcommand.append('-draw')
            JPEGcommand.append("text 200,80 'Marked {} brightest stars out of {}.'".format(nStarsLimit, self.nSExtracted))
        JPEGcommand.append(jpegFile)
        self.logger.debug("Issuing convert command to create jpeg.")
#         self.logger.debug("Command: {}".format(repr(JPEGcommand)))
        try:
            ConvertSTDOUT = subprocess.check_output(JPEGcommand, stderr=subprocess.STDOUT)
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
            self.jpegFileNames.append(jpegFileName)


    ##-------------------------------------------------------------------------
    ## Clean Up by Deleting Temporary Files
    ##-------------------------------------------------------------------------
    def CleanUp(self):
        '''
        Clean up by deleting temporary files.
        '''
        self.logger.info("Cleaning Up Temporary Files.")
        for item in self.tempFiles:
            if os.path.exists(item):
                self.logger.debug("Deleting {0}".format(item))
                os.remove(item)


    ##-------------------------------------------------------------------------
    ## Append Line With Image Info to HTML File List
    ##-------------------------------------------------------------------------
    def AddWebLogEntry(self, htmlImageList, fields=None):
        '''
        This function adds one line to the HTML table of images.  The line
        contains the image info extracted by IQMon.
        '''
        if not fields: fields=["Date and Time", "Filename", "Alt", "Az", "Airmass", "MoonSep", "MoonIllum", "FWHM", "ellipticity", "Background", "PErr", "PosAng", "ZeroPoint", "nStars", "ProcessTime"]
        ## If HTML file does not yet exist, create it and insert header
        ## from template file.
        if not os.path.exists(htmlImageList):
            self.logger.debug("HTML file does not exist.  Creating it.")
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
                      '    <h2>IQMon Results for {}</h2>'.format(self.tel.longName),
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
                header.append('        <th style="width:60px">FWHM ({})</th>'.format(str(self.tel.unitsForFWHM.unit)))
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
                IsEndTable = re.match("\s*</table>\s*", line)
                IsEndBody = re.match("\s*</body>\s*", line)
                IsEndHTML = re.match("\s*</html>\s*", line)
                if not IsEndTable and not IsEndBody and not IsEndHTML:
                    HTML.write(line+"\n")
        ## Write Lines for this Image to HTML File
        self.logger.info("Adding image data to HTML log file.")
        HTML = open(htmlImageList, 'a')
        HTML.write("    <tr>\n")
        ## Write Observation Date and Time
        if "Date and Time" in fields:
            HTML.write("      <td style='color:black;text-align:left'>{0}</td>\n".format(self.dateObs))
        ## Write Filename (and links to jpegs)
        if "Filename" in fields:
            if len(self.jpegFileNames) == 0:
                HTML.write("      <td style='color:black;text-align:left'>{0}</td>\n".format(self.rawFileBasename))
            elif len(self.jpegFileNames) == 1:
                HTML.write("      <td style='color:black;text-align:left'><a href='{0}'>{1}</a></td>\n".format(os.path.join("..", "..", "Plots", self.jpegFileNames[0]), self.rawFileBasename))
            elif len(self.jpegFileNames) == 2:
                HTML.write("      <td style='color:black;text-align:left'><a href='{0}'>{1}</a> (<a href='{2}'>JPEG2</a>)</td>\n".format(os.path.join("..", "..", "Plots", self.jpegFileNames[0]), self.rawFileBasename, os.path.join("..", "..", "Plots", self.jpegFileNames[1])))                
            elif len(self.jpegFileNames) >= 3:
                HTML.write("      <td style='color:black;text-align:left'><a href='{0}'>{1}</a> (<a href='{2}'>JPEG2</a>)(<a href='{3}'>JPEG3</a>)</td>\n".format(os.path.join("..", "..", "Plots", self.jpegFileNames[0]), self.rawFileBasename, os.path.join("..", "..", "Plots", self.jpegFileNames[1]), os.path.join("..", "..", "Plots", self.jpegFileNames[2])))
        ## Write Target Name
        if "Target" in fields:
            if self.objectName:
                HTML.write("      <td style='color:black'>{0:}</td>\n".format(self.objectName))
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
            if self.targetAlt:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.targetAlt.to(u.deg).value))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        if "Az" in fields:
            if self.targetAz:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.targetAz.to(u.deg).value))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        if "Airmass" in fields:
            if self.targetAz:
                HTML.write("      <td style='color:{0}'>{1:.2f}</td>\n".format("black", self.airmass))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        if "MoonSep" in fields:
            if self.moonSep:
                if self.moonAlt > 0:
                    HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.moonSep.to(u.deg).value))
                else:
                    HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", "down"))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        if "MoonIllum" in fields:
            if self.moonPhase:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.moonPhase))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
        ## Write FWHM and ellipticity
        if "FWHM" in fields:
            if self.FWHM:
                ## Decide whether to flag FWHM value with red color
                if self.FWHM > self.tel.thresholdFWHM.to(u.pix, equivalencies=self.tel.pixelScaleEquivalency):
                    colorFWHM = "#FF5C33"
                else:
                    colorFWHM = "#70DB70"
                ## Convert FWHM value to appropriate units for HTML output
                if self.tel.unitsForFWHM.unit == u.arcsec:
                    FWHM_for_HTML = (self.FWHM * u.radian.to(u.arcsec)*self.tel.pixelSize.to(u.mm)/self.tel.focalLength.to(u.mm)).value
                else:
                    FWHM_for_HTML = self.FWHM.value
                HTML.write("      <td style='background-color:{0}'>{1:.2f}</td>\n".format(colorFWHM, FWHM_for_HTML))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("#FF5C33", ""))
        if "ellipticity" in fields:
            if self.ellipticity:
                ## Decide whether to flag ellipticity value with red color
                if self.ellipticity > self.tel.thresholdEllipticity:
                    colorEllipticity = "#FF5C33"
                else:
                    colorEllipticity = "#70DB70"
                HTML.write("      <td style='background-color:{0}'>{1:.2f}</td>\n".format(colorEllipticity, self.ellipticity))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("#FF5C33", ""))
        ## Write SExtractor background and background RMS
        if "Background" in fields:
            if self.SExBackground and self.SExBRMS:
                HTML.write("      <td style='color:{0}'>{1:.0f} [{2:.0f}]</td>\n".format("black", self.SExBackground, self.SExBRMS))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Write pointing error
        if "PErr" in fields:
            if self.pointingError:
                ## Decide whether to flag pointing error value with red color
                if self.pointingError.arcminute > self.tel.thresholdPointingErr.to(u.arcmin).value:
                    colorPointingError = "#FF5C33"
                else:
                    colorPointingError = "#70DB70"
                ## Write HTML
                HTML.write("      <td style='background-color:{0}'>{1:.1f}</td>\n".format(colorPointingError, self.pointingError.arcminute))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("#FF5C33", ""))
        ## Write WCS position angle
        if "PosAng" in fields:
            if self.positionAngle:
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.positionAngle.to(u.deg).value))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Write zero point
        if "ZeroPoint" in fields:
            if self.zeroPoint:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", self.zeroPoint))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Write number of stars detected by SExtractor
        if "nStars" in fields:
            if self.nStarsSEx:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", self.nStarsSEx))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
        ## Write process time
        if "ProcessTime" in fields:
            if self.processTime:
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.processTime))
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
        self.logger.debug("Summary File: {0}".format(summaryFile))
        ## Read in previous data
        if not os.path.exists(summaryFile):
            self.logger.info("Making new astropy table object")
            SummaryTable = table.Table(names=("ExpStart", "File", "FWHM (pix)", "Ellipticity", 
                                       "Alt (deg)", "Az (deg)", "Airmass", "PointingError (arcmin)", 
                                       "ZeroPoint", "nStars", "Background", "Background RMS"),
                                 dtypes=('S22', 'S100', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4', 'f4', 'f4'),
                                 masked=True)
        else:
            self.logger.info("Reading astropy table object from file: {0}".format(summaryFile))
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
                                          'PointingError (arcmin)': [ascii.convert_numpy('f4')],
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
        if self.rawFileName: rawFileName = self.rawFileName
        else: 
            rawFileName = ""
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
        if self.targetAlt: targetAlt = self.targetAlt.to(u.deg).value
        else:
            targetAlt = 0.
            tableMask[4] = True
        ## Target Az
        if self.targetAz: targetAz = self.targetAz.to(u.deg).value
        else:
            targetAz = 0.
            tableMask[5] = True
        ## Airmass
        if self.airmass: airmass = self.airmass
        else:
            airmass = 0.
            tableMask[6] = True
        ## Pointing Error
        if self.pointingError: pointingError = self.pointingError.arcminute
        else:
            pointingError = 0.
            tableMask[7] = True
        ## Zero Point
        if self.zeroPoint: zeroPoint = self.zeroPoint
        else:
            zeroPoint = 0.
            tableMask[8] = True
        ## nStarsSEx
        if self.nStarsSEx: nStarsSEx = self.nStarsSEx
        else: 
            nStarsSEx = 0.
            tableMask[9] = True
        ## SExtractor Background
        if self.SExBackground: SExBackground = self.SExBackground
        else:
            SExBackground = 0.
            tableMask[10] = True
        ## SExtractor Background RMS
        if self.SExBRMS: SExBRMS = self.SExBRMS
        else:
            SExBRMS = 0.
            tableMask[11] = True
        ## Add row to table
        self.logger.debug("Writing new row to log table.  Filename: {0}".format(rawFileName))
        SummaryTable.add_row((dateObs, rawFileName,
                              FWHM, ellipticity,
                              targetAlt, targetAz,
                              airmass, pointingError,
                              zeroPoint, nStarsSEx,
                              SExBackground, SExBRMS),
                              mask=tableMask)
        ## Write Table to File
        self.logger.info("Writing new summary file.")
        ascii.write(SummaryTable, summaryFile,
                    Writer=ascii.basic.Basic)


    ##-------------------------------------------------------------------------
    ## Calcualte Process Time
    ##-------------------------------------------------------------------------
    def CalculateProcessTime(self):
        self.endProcessTime = time.time()
        self.processTime = self.endProcessTime - self.startProcessTime
        self.logger.info("IQMon processing time = {0:.1f} seconds".format(self.processTime))

    