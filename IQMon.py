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
import subprocess32
import math
import numpy as np

## Import Astronomy Specific Tools
import ephem
import astropy.units as u
import astropy.io.fits as fits
import astropy.coordinates as coords
import astropy.table as table
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
    - location of configuration file for IQMon
    '''
    _singletons = dict()

    def __new__(cls):
        if not cls in cls._singletons:
            cls._singletons[cls] = object.__new__(cls)
        return cls._singletons[cls]

    def __init__(self):
        '''
        Read and parse configuration file.
        - Currently assumes that file is .IQMonConfig in the user's home
        directory.
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
            IsPythonPath = re.match("IQMONPYTHON\s=\s([\w/\-\.]+)", line)
            if IsPythonPath:
                PythonPath = os.path.abspath(IsPythonPath.group(1))
            IsCatalogPath = re.match("CATALOGPATH\s=\s([\w/\-\.]+)", line)
            if IsCatalogPath:
                self.pathCatalog = os.path.abspath(IsCatalogPath.group(1))


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
    - name
    - longName
    - focalLength
    - pixelSize
    - aperture
    - gain
    - nXPix
    - nYPix
    - unitsForFWHM
    - ROI
    - thresholdFWHM
    - thresholdPointingErr
    - thresholdEllipticity
    - pixelScale
    - fRatio
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
        self.site = None
        
    def CheckUnits(self, logger):
        '''
        Method to check whether the properties of the object have units and to
        add units if the value has no units.  Does not yet check whether the
        unit is reasonable (i.e. does focalLength have length units).
        '''
        ## name is a string
        ## longName is a string
        ## Default focalLength to units of mm
        if self.focalLength and not hasattr(self.focalLength, 'unit'):
            if logger:
                logger.debug(
                    "focalLength is unitless value.  Adding mm.")
            self.focalLength *= u.mm
        ## Default pixelSize to units of microns
        if self.pixelSize and not hasattr(self.pixelSize, 'unit'):
            if logger:
                logger.debug(
                    "pixelSize is unitless value.  Adding microns")
            self.pixelSize *= u.micron
        ## Default aperture to units of mm
        if self.aperture and not hasattr(self.aperture, 'unit'):
            if logger:
                logger.debug(
                    "aperture is unitless value.  Adding mm")
            self.aperture *= u.mm
        ## Default gain to units of 1/ADU
        if self.gain and not hasattr(self.gain, 'unit'):
            if logger:
                logger.debug(
                    "gain is unitless value.  Adding 1/ADU")
            self.gain /= u.adu
        ## Default nXPix to units of pixels
        if self.nXPix and not hasattr(self.nXPix, 'unit'):
            if logger:
                logger.debug(
                    "nXPix is unitless value.  Adding pixels")
            self.nXPix *= u.pix
        ## Default nYPix to units of pixels
        if self.nYPix and not hasattr(self.nYPix, 'unit'):
            if logger:
                logger.debug(
                    "nYPix is unitless value.  Adding pixels")
            self.nYPix *= u.pix
        ## Default unitsForFWHM to units of arcsec
        if self.unitsForFWHM and not hasattr(self.unitsForFWHM, 'unit'):
            if logger:
                logger.debug(
                    "unitsForFWHM is unitless value.  Adding arcsec")
            self.unitsForFWHM *= u.arcsec
        ## ROI is string
        ## Default thresholdFWHM to units of arcsec
        if self.thresholdFWHM and not hasattr(self.thresholdFWHM, 'unit'):
            if logger:
                logger.debug(
                    "thresholdFWHM is unitless value.  Adding arcsec")
            self.thresholdFWHM *= u.arcsec
        ## Default thresholdPointingErr to units of arcmin
        if self.thresholdPointingErr and not hasattr(self.thresholdPointingErr, 'unit'):
            if logger:
                logger.debug(
                    "thresholdPointingErr is unitless value.  Adding arcmin")
            self.thresholdPointingErr *= u.arcmin
        ## Default thresholdEllipticity to dimensionless
        if self.thresholdEllipticity and not hasattr(self.thresholdEllipticity, 'unit'):
            if logger:
                logger.debug(
                    "thresholdEllipticity is unitless value.  Adding dimensionless")
            self.thresholdEllipticity *= u.dimensionless_unscaled
        ## Default pixelScale to units of arcsec per pixel
        if self.pixelScale and not hasattr(self.pixelScale, 'unit'):
            if logger:
                logger.debug(
                    "pixelScale is unitless value.  Adding arcsec / pixel")
            self.pixelScale *= u.arcsec / u.pix
        ## Default fRatio to dimensionless
        if self.fRatio and not hasattr(self.fRatio, 'unit'):
            if logger:
                logger.debug(
                    "fRatio is unitless value.  Adding dimensionless")
            self.fRatio *= u.dimensionless_unscaled


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
    def __init__(self, input):
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
        self.htmlImageList = None
        self.positionAngle = None
        self.zeroPoint = None
        self.processTime = None
        self.FWHM = None
        self.ellipticity = None
        self.pointingError = None
        self.imageFlipped = None
        self.jpegFileNames = []
        self.summaryFile = None
        self.htmlImageList = None


    ##-------------------------------------------------------------------------
    ## Get Header
    ##-------------------------------------------------------------------------
    def GetHeader(self, tel, logger):
        '''
        Get information from the image fits header.
        '''
        hdulist = fits.open(self.workingFile)
        self.header = hdulist[0].header
        self.image = hdulist[0].data
        hdulist.close()
        logger.info("Reading image header.")
        
        ## Get exposure time from header (assumes seconds)
        try:
            self.exptime = float(self.header['EXPTIME']) * u.s
        except:
            self.exptime = None
            logger.debug("No exposure time value found in header")
        else:
            logger.debug("Exposure time = {0:.1f} s".format(self.exptime.to(u.s).value))
        ## Get filter from header
        try:
            self.filter = self.header['FILTER']
        except:
            self.filter = None
            logger.debug("No filter keyword found in header")
        else:
            logger.debug("Filter = {0}".format(self.filter))
        ## Get focus position from header
        try:
            self.focusPos = self.header['FOCUSPOS']
        except:
            self.focusPos = None
            logger.debug("No focus position value found in header")
        else:
            logger.debug("Focus position = {0}".format(self.focusPos))
        ## Get object name from header
        try:
            self.objectName = self.header["OBJECT"]
        except:
            self.objectName = None
            logger.debug("No object value found in header")
        else:
            logger.debug("Header object name = {0}".format(self.objectName))
        ## Get airmass from header
        try:
            self.headerAirmass = self.header["AIRMASS"]
        except:
            self.headerAirmass = None
            logger.debug("No airmass value found in header")
        else:
            logger.debug("Header airmass = {0:.2f}".format(self.headerAirmass))
        ## Get Observation Date and Time from header
        ## (assumes YYYY-MM-DDTHH:MM:SS format)
        try:
            self.dateObs = self.header["DATE-OBS"]
        except:
            self.dateObs = None
            logger.debug("No date value found in header")
        else:
            logger.debug("Header date = {0}".format(self.dateObs))
        ## Get Site Latitude from header (assumes decimal degrees)
        try:
            self.latitude = self.header["LAT-OBS"] * u.deg
        except:
            self.latitude = None
            logger.debug("No latitude value found in header")
        else:
            logger.debug("Header latitude = {0:.4f} deg".format(self.latitude.to(u.deg).value))
        ## Get Site Longitude from header (assumes decimal degrees)
        try:
            self.longitude = self.header["LONG-OBS"] * u.deg
        except:
            self.longitude = None
            logger.debug("No longitiude value found in header")
        else:
            logger.debug("Header longitiude = {0:.4f} deg".format(self.longitude.to(u.deg).value))
        ## Get Site Altitude from header (assumes meters)
        try:
            self.altitude = self.header["ALT-OBS"] * u.meter
        except:
            self.altitude = None
            logger.debug("No altitude value found in header")
        else:
            logger.debug("Header altitude = {0:.0f} meters".format(self.altitude.to(u.meter).value))


        ## Determine Image Size in Pixels
        self.nYPix, self.nXPix = self.image.shape

        ## Read Header Coordinates in to astropy coordinates object
        ImageRA  = self.header['RA']
        if len(ImageRA.split(":")) != 3:
            if len(ImageRA.split(" ")) == 3:
                ImageRA = ":".join(ImageRA.split(" "))
        ImageDEC = self.header['DEC']    
        if len(ImageDEC.split(":")) != 3:
            if len(ImageDEC.split(" ")) == 3:
                ImageDEC = ":".join(ImageDEC.split(" "))
        logger.debug("Read pointing info from header: "+ImageRA+" "+ImageDEC)
        try:
            self.coordinate_header = coords.ICRSCoordinates(
                                                      ImageRA+" "+ImageDEC,
                                                   unit=(u.hour, u.degree))
        except:
            logger.warning("Failed to read pointing info from header.")
            self.coordinate_header = None

        ## Read WCS
        try:
            self.imageWCS = wcs.WCS(self.header)
        except:
            self.imageWCS = None
            logger.info("No WCS found in image header")
        else:
            logger.debug("Found WCS in image header.")

        ## Determine PA of Image
        try:
            PC11 = float(self.imageWCS.to_header()['PC1_1'])
            PC12 = float(self.imageWCS.to_header()['PC1_2'])
            PC21 = float(self.imageWCS.to_header()['PC2_1'])
            PC22 = float(self.imageWCS.to_header()['PC2_2'])
        except:
            logger.debug("Could not find PCn_m values in WCS.")
            try:
                PC11 = float(self.imageWCS.to_header()['CD1_1'])
                PC12 = float(self.imageWCS.to_header()['CD1_2'])
                PC21 = float(self.imageWCS.to_header()['CD2_1'])
                PC22 = float(self.imageWCS.to_header()['CD2_2'])
            except:
                logger.debug("Could not find CDn_m values in WCS.")
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
            logger.debug("Position angle of WCS is {0:.1f} degrees.".format(self.positionAngle.to(u.deg).value))
            logger.debug("Image orientation is North {0}, East {1}.".format(North, East))
            if self.imageFlipped:
                logger.debug("Image is mirrored.")
        else:
            self.positionAngle = None
            self.imageFlipped = None


        ## Determine Alt, Az, Moon Sep, Moon Illum using ephem module
        if self.dateObs and self.latitude and self.longitude:
            ## Populate site object properties
            SiteDate = "/".join(self.dateObs[0:10].split("-"))
            SiteTime = self.dateObs[11:]        
            tel.site.date = ephem.Date(SiteDate+" "+SiteTime)
            tel.site.lat = str(self.latitude.to(u.deg).value)
            tel.site.lon = str(self.longitude.to(u.deg).value)
            if self.altitude: tel.site.elevation = self.altitude.to(u.meter).value
            ## Do calculations using ephem
            TargetObject = ephem.readdb("Target,f|M|F7,"+ImageRA+","+ImageDEC+",2.02,2000")
            TargetObject.compute(tel.site)
            self.targetAlt = TargetObject.alt * 180./ephem.pi * u.deg
            self.targetAz = TargetObject.az * 180./ephem.pi * u.deg
            logger.debug("Target Alt, Az = {0:.1f}, {1:.1f}".format(self.targetAlt.to(u.deg).value, self.targetAz.to(u.deg).value))
            self.zenithAngle = 90.*u.deg - self.targetAlt
            self.airmass = 1.0/math.cos(self.zenithAngle.to(u.radian).value)*(1.0 - 0.0012*(1.0/(math.cos(self.zenithAngle.to(u.radian).value)**2 - 1.0)))
            logger.debug("Target airmass (calculated) = {0:.2f}".format(self.airmass))
            ## Calculate Moon Position and Illumination
            TheMoon = ephem.Moon()
            TheMoon.compute(tel.site)
            self.moonPhase = TheMoon.phase
            self.moonSep = ephem.separation(TargetObject, TheMoon) * 180./ephem.pi * u.deg
            self.moonAlt = TheMoon.alt * 180./ephem.pi * u.deg
            logger.debug("A {0:.0f} percent illuminated Moon is {1:.0f} deg from target.".format(self.moonPhase, self.moonSep.to(u.deg).value))
        else:
            self.targetAlt = None
            self.targetAz = None
            self.moonPhase = None
            self.moonSep = None
            self.moonAlt = None
            self.zenithAngle = None
            self.airmass = None
            logger.warning("Object position and Moon position not calculated.")


    ##-------------------------------------------------------------------------
    ## Read Image
    ##-------------------------------------------------------------------------
    def ReadImage(self, config):
        '''
        Read the raw image and write out a working image in the IQMon temporary
        directory.
        
        - For the moment, this only copies a fits file from the original
          location to the IQMon tmp directory.
        - Later implement file format conversion from CRW, CR2, DNG, etc to
          fits using dcraw.
        '''
        self.workingFile = os.path.join(config.pathTemp, self.rawFileName)
        shutil.copy2(self.rawFile, self.workingFile)
        self.tempFiles.append(self.workingFile)

    ##-------------------------------------------------------------------------
    ## Dark Subtract Image
    ##-------------------------------------------------------------------------
    def DarkSubtract(self, Darks, tel, config, logger):
        '''
        Create master dark and subtract from image.
        
        Input the filename of the appropriate master dark.  May want to write
        own function to make the master dark given input file data.
        '''
        logger.debug("Dark subtracting image.  Opening image data.")
        hdulist_image = fits.open(self.workingFile, mode='update')
        ## Load master dark if provided, but if multiple files input, combine
        ## them in to master dark, then load combined master dark.
        if len(Darks) == 1:
            logger.debug("Found master dark.  Opening master dark data.")
            hdulist_dark = fits.open(Darks[0])
            MasterDarkData = hdulist_dark[0].data
        elif len(Darks) > 1:
            logger.info("Multiple input darks detected.  Median combining {0} darks.".format(len(Darks)))
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
            MasterDarkFilename = "MasterDark_"+tel.name+"_"+DataNightString+"_"+str(int(math.floor(self.exptime.to(u.s).value)))+".fits"
            MasterDarkFile  = os.path.join(config.pathTemp, MasterDarkFilename)    
            hdu_MasterDark = fits.PrimaryHDU(MasterDarkData)
            hdulist_MasterDark = fits.HDUList([hdu_MasterDark])
            hdulist_MasterDark.header = hdulist[0].header
            hdulist_MasterDark.header['history'] = "Combined {0} images to make this master dark.".format(len(Darks))
            logger.info("Writing master dark file: {0}".format(MasterDarkFile))
            hdulist_MasterDark.writeto(MasterDarkFile)
        else:
            logger.error("No input dark files detected.")
        ## Now Subtract MasterDark from Image
        logger.info("Subtracting dark from image.")
        ImageData = hdulist_image[0].data
        DifferenceImage = ImageData - MasterDarkData
        hdulist_image[0].data = DifferenceImage
        hdulist_image.flush()
#         logger.debug("Median level of image = {0}".format(np.median(ImageData)))
#         logger.debug("Median level of dark = {0}".format(np.median(MasterDarkData)))
#         logger.debug("Median level of dark subtracted = {0}".format(np.median(DifferenceImage)))


    ##-------------------------------------------------------------------------
    ## Crop Image
    ##-------------------------------------------------------------------------
    def Crop(self, tel, logger):
        '''
        Crop working image to region of interest.
        '''
        if tel.ROI:
            ## Parse ROI String
            try:
                MatchROI = re.match("\[?(\d{1,5}):(\d{1,5}),(\d{1,5}):(\d{1,5})\]?", tel.ROI)
            except:
                logger.warning("Could not parse ROI string in telescope object.")
            else:
                x1 = int(MatchROI.group(1))
                x2 = int(MatchROI.group(2))
                y1 = int(MatchROI.group(3))
                y2 = int(MatchROI.group(4))
                logger.info("Cropping Image To [{0}:{1},{2}:{3}]".format(x1, x2, y1, y2))
                hdulist = fits.open(self.workingFile, mode="update")
                hdulist[0].data = hdulist[0].data[x1:x2,y1:y2]
                hdulist.flush()
                hdulist.close()


    ##-------------------------------------------------------------------------
    ## Solve Astrometry Using astrometry.net
    ##-------------------------------------------------------------------------
    def SolveAstrometry(self, tel, config, logger):
        '''
        Solve astrometry in the working image using the astrometry.net solver.
        '''
        logger.info("Attempting to create WCS using Astrometry.net solver.")
        AstrometryCommand = ["solve-field", "-l", "5", "-O", "-p",
                             "-L", str(tel.pixelScale*0.90),
                             "-H", str(tel.pixelScale*1.10),
                             "-u", "arcsecperpix", "-z", "4", self.workingFile]
        AstrometrySTDOUT = ""

        try:
            StartTime = time.time()
            AstrometrySTDOUT = subprocess32.check_output(AstrometryCommand, 
                               stderr=subprocess32.STDOUT, timeout=20)
            EndTime = time.time()
#         except TimeoutExpired:
#             logger.warning("Astrometry.net timed out")
#             self.astrometrySolved = False
        except:
            logger.warning("Astrometry.net failed.")
            self.astrometrySolved = False
        else:
            ProcessTime = EndTime - StartTime
            logger.debug("Astrometry.net Processing Time: %.1f s", ProcessTime)
            pos = AstrometrySTDOUT.find("Field center: (RA H:M:S, Dec D:M:S) = ")
            if pos != -1:
                IsFieldCenter = re.match("\s*(\d{1,2}:\d{2}:\d{2}\.\d+,\s-?\d{1,2}:\d{2}:\d{2}\.\d+).*", 
                                         AstrometrySTDOUT[pos+40:pos+75])
                if IsFieldCenter:
                    logger.info("Astrometry.net field center is: %s", IsFieldCenter.group(1))
            else:
                for line in AstrometrySTDOUT:
                    logger.warning("  %s" % line)
            NewFile = self.workingFile.replace(self.fileExt, ".new")
            NewFitsFile = self.workingFile.replace(self.fileExt, ".new.fits")
            if not os.path.exists(NewFile):
                logger.warning("No new file created by astrometry.net")
                self.astrometrySolved = False
            else:
                logger.debug("Astrometry.net succeeded")
                if os.path.exists(NewFitsFile): os.remove(NewFitsFile)
                os.rename(NewFile, NewFitsFile)
                self.astrometrySolved = True
                ## Update header history
                hdulist = fits.open(self.workingFile, mode="update")
                now = time.gmtime()
                hdulist[0].header['history'] = "Solved by Astrometry.net at {0}".format(time.strftime("%Y-%m-%dT%H:%M:%S UTC"))
                hdulist.close()
            ## Add files created by astrometry.net to tempFiles list
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+".axy"))
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+".wcs"))
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+".solved"))
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+".rdls"))
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+".match"))
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+".corr"))
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+".new.fits"))
            self.tempFiles.append(os.path.join(config.pathTemp, self.rawFileBasename+"-indx.xyls"))

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
    def DeterminePointingError(self, logger):
        '''
        Determine pointing error (difference between objects coordinates and
        solved WCS).
        '''
        logger.info("Detemining pointing error based on WCS solution")
        if self.imageWCS and self.coordinate_header:
            centerWCS = self.imageWCS.wcs_pix2world([[self.nXPix/2, self.nYPix/2]], 1)
            logger.debug("Using coordinates of center point: {0} {1}".format(centerWCS[0][0], centerWCS[0][1]))
            self.coordinate_WCS = coords.ICRSCoordinates(ra=centerWCS[0][0],
                                                   dec=centerWCS[0][1],
                                                   unit=(u.degree, u.degree))
            self.pointingError = self.coordinate_WCS.separation(self.coordinate_header)
            logger.debug("Target Coordinates are:  %s %s",
                         self.coordinate_header.ra.format(u.hour, sep=":", precision=1),
                         self.coordinate_header.dec.format(u.degree, sep=":", precision=1, alwayssign=True))
            logger.debug("WCS of Central Pixel is: %s %s",
                         self.coordinate_WCS.ra.format(u.hour, sep=":", precision=1),
                         self.coordinate_WCS.dec.format(u.degree, sep=":", precision=1, alwayssign=True))
            logger.info("Pointing Error is %.2f arcmin", self.pointingError.arcmins)
        else:
            logger.warning("Pointing error not calculated.")

    ##-------------------------------------------------------------------------
    ## Run SExtractor
    ##-------------------------------------------------------------------------
    def RunSExtractor(self, tel, config, logger):
        '''
        Run SExtractor on image.
        
        - need to check that tel.SExtractorPhotAperture is set
        - need to check that tel.gain is set
        - need to check that tel.pixelScale is set
        - need to check that tel.SExtractorSeeing is set
        '''
        ## Set up file names
        SExtractorDefaultFile = os.path.join(config.pathIQMonExec, "default.sex")
        SExtractorConfigFile = os.path.join(config.pathTemp, self.rawFileBasename+".sex")
        self.tempFiles.append(SExtractorConfigFile)
        SExtractorCatalog = os.path.join(config.pathTemp, self.rawFileBasename+".cat")
        self.tempFiles.append(SExtractorCatalog)
        PhotometryCatalogFile_xy = os.path.join(config.pathTemp, self.rawFileBasename+"PhotCat_xy.txt")
        self.tempFiles.append(PhotometryCatalogFile_xy)

        ## Create PhotometryCatalogFile_xy file for SExtractor Association
        if os.path.exists(PhotometryCatalogFile_xy): os.remove(PhotometryCatalogFile_xy)
        PhotCatFileObject = open(PhotometryCatalogFile_xy, 'w')
        PhotCatFileObject.write("# No Existing WCS Found for this image\n")
        PhotCatFileObject.write("# This is a dummy file to keep SExtractor happy\n")
        PhotCatFileObject.write("0.0  0.0  0.0  0.0\n")
        PhotCatFileObject.close()

        ## Make edits To default.sex based on telescope:
        ## Read in default config file        
        DefaultConfig = open(SExtractorDefaultFile, 'r')
        NewConfig     = open(SExtractorConfigFile, 'w')
        for line in DefaultConfig:
            newline = line
            if re.match("CATALOG_NAME\s+", line):
                newline = "CATALOG_NAME     "+SExtractorCatalog+"\n"
            if re.match("PARAMETERS_NAME\s+", line):
                newline = "PARAMETERS_NAME  "+os.path.join(config.pathIQMonExec, "default.param")+"\n"
            if re.match("PHOT_APERTURES\s+", line):
                newline = "PHOT_APERTURES   "+str(tel.SExtractorPhotAperture.to(u.pix).value)+"\n"
            if re.match("GAIN\s+", line):
                newline = "GAIN             "+str(tel.gain.value)+"\n"
            if re.match("PIXEL_SCALE\s+", line):
                newline = "PIXEL_SCALE      "+str(tel.pixelScale.value)+"\n"
            if re.match("SEEING_FWHM\s+", line):
                newline = "SEEING_FWHM      "+str(tel.SExtractorSeeing.to(u.arcsec).value)+"\n"
            if re.match("ASSOC_NAME\s+", line):
                newline = "ASSOC_NAME       "+PhotometryCatalogFile_xy+"\n"
            NewConfig.write(newline)
        DefaultConfig.close()
        NewConfig.close()

        ## Run SExtractor
        logger.info("Invoking SExtractor.")
        SExtractorCommand = ["sex", self.workingFile, "-c", SExtractorConfigFile]
        try:
            SExSTDOUT = subprocess32.check_output(SExtractorCommand, stderr=subprocess32.STDOUT, timeout=30)
        except:
            logger.warning("SExtractor error.")
        else:
            for line in SExSTDOUT.split("\n"):
                line.replace("[1A", "")
                line.replace("[1M>", "")
                if not re.match(".*Setting up background map.*", line) and not re.match(".*Line:\s[0-9]*.*", line):
                    logger.debug("  "+line)
            ## Extract Number of Stars from SExtractor Output
            pos = SExSTDOUT.find("sextracted ")
            IsSExCount = re.match("\s*([0-9]+)\s+", SExSTDOUT[pos+11:pos+21])
            if IsSExCount:
                self.nSExtracted = int(IsSExCount.group(1))
                logger.info("SExtractor found {0} sources.".format(self.nSExtracted))
            else:
                self.nSExtracted = None
            ## Extract Background Level from SExtractor Output
            pos = SExSTDOUT.find("Background: ")
            IsSExBkgnd = re.match("\s*([0-9\.]+)\s*", SExSTDOUT[pos+11:pos+21])
            if IsSExBkgnd:
                self.SExBackground = float(IsSExBkgnd.group(1))
                logger.info("SExtractor background is {0:.1f}".format(self.SExBackground))
            else:
                self.SExBackground = None
            ## Extract Background RMS from SExtractor Output
            IsSExBRMS = re.match("\s*RMS:\s([0-9\.]+)\s*", SExSTDOUT[pos+21:pos+37])
            if IsSExBRMS:
                self.SExBRMS = float(IsSExBRMS.group(1))
                logger.info("SExtractor background RMS is {0:.1f}".format(self.SExBRMS))
            else:
                self.SExBRMS = None

            ## If No Output Catalog Created ...
            if not os.path.exists(SExtractorCatalog):
                logger.warning("SExtractor failed to create catalog.")
                self.SExtractorCatalog = None
            else:
                self.SExtractorCatalog = SExtractorCatalog

            ## Read Catalog
            logger.debug("Reading SExtractor output catalog.")
            self.SExtractorResults = ascii.read(self.SExtractorCatalog, Reader=ascii.sextractor.SExtractor)
            SExImageRadius = []
            SExAngleInImage = []
            for star in self.SExtractorResults:
                SExImageRadius.append(math.sqrt((self.nXPix/2-star['X_IMAGE'])**2 + (self.nYPix/2-star['Y_IMAGE'])**2))
                SExAngleInImage.append(math.atan((star['X_IMAGE']-self.nXPix/2)/(self.nYPix/2-star['Y_IMAGE']))*180.0/math.pi)
            self.SExtractorResults.add_column(table.Column(data=SExImageRadius, name='ImageRadius'))
            self.SExtractorResults.add_column(table.Column(data=SExAngleInImage, name='AngleInImage'))
            self.nStarsSEx = len(self.SExtractorResults)
            logger.info("Read in {0} stars from SExtractor catalog.".format(self.nStarsSEx))


    ##-------------------------------------------------------------------------
    ## Determine Image FWHM from SExtractor Catalog
    ##-------------------------------------------------------------------------
    def DetermineFWHM(self, logger):
        '''
        Determine typical FWHM of image from SExtractor results.
        '''
        if self.nStarsSEx > 1:
            IQRadiusFactor = 1.0
            DiagonalRadius = math.sqrt((self.nXPix/2)**2+(self.nYPix/2)**2)
            IQRadius = DiagonalRadius*IQRadiusFactor
            CentralFWHMs = []
            CentralEllipticities = []
            for star in self.SExtractorResults:
                if star['ImageRadius'] <= IQRadius:
                    CentralFWHMs.append(star['FWHM_IMAGE'])
                    CentralEllipticities.append(star['ELLIPTICITY'])   
            if len(CentralFWHMs) > 3:
                self.FWHM = np.median(CentralFWHMs) * u.pix
                self.ellipticity = np.median(CentralEllipticities)
            else:
                logger.warning("Not enough stars detected in central region of image to form median FWHM.")
            logger.info("Median FWHM in inner region is {0:.2f} pixels".format(self.FWHM.to(u.pix).value))
            logger.info("Median Ellipticity in inner region is {0:.2f}".format(self.ellipticity))
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
    def MakeJPEG(self, jpegFileName, tel, config, logger, marked=False, rotate=False, binning=1):
        '''
        Make jpegs of image.
        '''
        jpegFile = os.path.join(config.pathPlots, jpegFileName)
        if marked:
            logger.info("Making marked jpeg with binning factor of {0}.".format(binning))
        else:
            logger.info("Making jpeg with binning factor of {0}.".format(binning))
        if os.path.exists(jpegFile): os.remove(jpegFile)
        binningString = str(1./binning*100)+"%"
        JPEGcommand = ["convert", "-contrast-stretch", "0.9%,1%", "-compress", "JPEG", "-quality", "70", "-stroke", "red", "-fill", "none", "-resize", binningString]
        if marked:
            if self.FWHM:
                MarkRadius=max([4, 2*math.ceil(self.FWHM.value)])
            else:
                MarkRadius = 4
            for star in self.SExtractorResults:
                MarkXPos = star['X_IMAGE']
                MarkYPos = self.nXPix - star['Y_IMAGE']
                JPEGcommand.append('-draw')
                JPEGcommand.append("circle %d,%d %d,%d" % (MarkXPos, MarkYPos, MarkXPos+MarkRadius, MarkYPos))
        if rotate:
            if self.positionAngle:
                JPEGcommand.append("-rotate")
                JPEGcommand.append(str(self.positionAngle.to(u.deg).value))
                if self.imageFlipped:
                    JPEGcommand.append("-flop")
            else:
                logger.warning("No position angle value found.  Not rotating JPEG.")
        JPEGcommand.append(self.workingFile)
        JPEGcommand.append(jpegFile)
        try:        
            ConvertSTDOUT = subprocess32.check_output(JPEGcommand, stderr=subprocess32.STDOUT, timeout=30)
        except:
            logger.warning("Failed to create jpeg.")
        else:
            self.jpegFileNames.append(jpegFileName)


    ##-------------------------------------------------------------------------
    ## Clean Up by Deleting Temporary Files
    ##-------------------------------------------------------------------------
    def CleanUp(self, logger):
        '''
        Clean up by deleting temporary files.
        '''
        logger.info("Cleaning Up Temporary Files.")
        for item in self.tempFiles:
            if os.path.exists(item):
                logger.debug("Deleting {0}".format(item))
                os.remove(item)

    ##-------------------------------------------------------------------------
    ## Append Line With Image Info to HTML File List
    ##-------------------------------------------------------------------------
    def AddWebLogEntry(self, tel, config, logger):
        '''
        This function adds one line to the HTML table of images.  The line
        contains the image info extracted by IQMon.
        '''
        if self.htmlImageList:
            ## If HTML file does not yet exist, create it and insert header
            ## from template file.
            if not os.path.exists(self.htmlImageList):
                logger.debug("HTML files does not exist.  Creating it.")
                HTML = open(self.htmlImageList, 'w')
                HTMLheader = open(os.path.join(config.pathIQMonExec, "ImageListHeader.html"), 'r')
                header = HTMLheader.read()
                header = header.replace("telescopename", tel.longName)
                header = header.replace("FWHMunits", str(tel.unitsForFWHM.unit))
                HTMLheader.close()
                HTML.write(header)
                HTML.close()
            ## If HTML file does exist, we need to strip off the lines which
            ## end the file, so we can append more data to the table.
            else:
                logger.debug("HTML file exists.  Copying contents.")
                HTML = open(self.htmlImageList, 'r')
                existingContent = HTML.read().split("\n")
                HTML.close()
                HTML = open(self.htmlImageList, 'w')
                for line in existingContent:
                    IsEndTable = re.match("\s*</table>\s*", line)
                    IsEndBody = re.match("\s*</body>\s*", line)
                    IsEndHTML = re.match("\s*</html>\s*", line)
                    if not IsEndTable and not IsEndBody and not IsEndHTML:
                        HTML.write(line)
            ## Write Lines for this Image to HTML File
            logger.info("Adding image data to HTML log file.")
            HTML = open(self.htmlImageList, 'a')
            HTML.write("    <tr>\n")
            HTML.write("      <td style='color:black;text-align:left'>{0}</td>\n".format(self.dateObs))
            if len(self.jpegFileNames) == 0:
                HTML.write("      <td style='color:black;text-align:left'>{0}</td>\n".format(self.rawFileBasename))
            elif len(self.jpegFileNames) == 1:
                HTML.write("      <td style='color:black;text-align:left'><a href='{0}'>{1}</a></td>\n".format(os.path.join("..", "..", "Plots", self.jpegFileNames[0]), self.rawFileBasename))
            elif len(self.jpegFileNames) >= 2:
                HTML.write("      <td style='color:black;text-align:left'><a href='{0}'>{1}</a> (<a href='{2}'>JPEG2</a>)</td>\n".format(os.path.join("..", "..", "Plots", self.jpegFileNames[0]), self.rawFileBasename, os.path.join("..", "..", "Plots", self.jpegFileNames[1])))                
            if self.targetAlt and self.targetAz and self.airmass and self.moonSep and self.moonPhase:
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.targetAlt.to(u.deg).value))
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.targetAz.to(u.deg).value))
                HTML.write("      <td style='color:{0}'>{1:.2f}</td>\n".format("black", self.airmass))
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.moonSep.to(u.deg).value))
                HTML.write("      <td style='color:black'>{0:.1f}</td>\n".format(self.moonPhase))
            else:
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
                HTML.write("      <td style='color:black'>{0}</td>\n".format(""))
            if self.FWHM and self.ellipticity:
                if tel.unitsForFWHM.unit == u.arcsec:
                    FWHM_for_HTML = (self.FWHM * u.radian.to(u.arcsec)*tel.pixelSize.to(u.mm)/tel.focalLength.to(u.mm)).value
                else:
                    FWHM_for_HTML = self.FWHM.value
                HTML.write("      <td style='color:{0}'>{1:.2f}</td>\n".format("black", FWHM_for_HTML))
                HTML.write("      <td style='color:{0}'>{1:.2f}</td>\n".format("black", self.ellipticity))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
            if self.SExBackground and self.SExBRMS:
                HTML.write("      <td style='color:{0}'>{1:.1f} [{2:.1f}]</td>\n".format("black", self.SExBackground, self.SExBRMS))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
            if self.pointingError:
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.pointingError.arcmins))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
            if self.positionAngle:
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.positionAngle.to(u.deg).value))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
            if self.zeroPoint:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", self.zeroPoint))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
            if self.nStarsSEx:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", self.nStarsSEx))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
            if self.processTime:
                HTML.write("      <td style='color:{0}'>{1:.1f}</td>\n".format("black", self.processTime))
            else:
                HTML.write("      <td style='color:{0}'>{1}</td>\n".format("black", ""))
            HTML.write("    </tr>\n")
            HTML.write("  </table>\n")
            HTML.write("</body>\n")
            HTML.write("</html>\n")
            HTML.close()


    ##-------------------------------------------------------------------------
    ## Append Line With Image Info to Summary Text File
    ##-------------------------------------------------------------------------
    def AddSummaryEntry(self, logger):
        if self.summaryFile:
            logger.info("Writing Summary File Entry.")
            logger.debug("Summary File: {0}".format(self.summaryFile))
            ## Read in previous data
            if not os.path.exists(self.summaryFile):
                logger.info("Making new astropy table object")
                SummaryTable = table.Table(names=("ExpStart", "File", "FWHM (pix)", "Ellipticity", 
                                           "Alt (deg)", "Az (deg)", "Airmass", "PointingError (arcmin)", 
                                           "ZeroPoint", "nStars", "Background", "Background RMS"),
                                     dtypes=('a22', 'a120', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4', 'f4', 'f4'),
                                     masked=True)
            else:
                logger.info("Reading astropy table object from file: {0}".format(self.summaryFile))
                SummaryTable = ascii.read(self.summaryFile)
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
            if self.pointingError: pointingError = self.pointingError.arcmins
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
            logger.debug("Writing new row to log table.")
            SummaryTable.add_row((dateObs, rawFileName,
                                  FWHM, ellipticity,
                                  targetAlt, targetAz,
                                  airmass, pointingError,
                                  zeroPoint, nStarsSEx,
                                  SExBackground, SExBRMS),
                                  mask=tableMask)
            ## Write Table to File
            logger.info("Writing new summary file.")
            ascii.write(SummaryTable, self.summaryFile,
                        Writer=ascii.basic.Basic)


    ##-------------------------------------------------------------------------
    ## Calcualte Process Time
    ##-------------------------------------------------------------------------
    def CalculateProcessTime(self, logger):
        self.endProcessTime = time.time()
        self.processTime = self.endProcessTime - self.startProcessTime
        logger.info("IQMon processing time = {0:.1f} seconds".format(self.processTime))


if __name__ == '__main__':
    unittest.main()
