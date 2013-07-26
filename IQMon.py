#!/usr/bin/env python
# encoding: utf-8
"""
IQMon.py

Created by Josh Walawender on 2013-07-22.
Copyright (c) 2013 . All rights reserved.
"""
from __future__ import division, print_function

import sys
import os
import re
import unittest

import astropy.units as u

##--------------------------------------------------------------------------------------------------
## Define Config object to hold IQMon configuration information
##--------------------------------------------------------------------------------------------------
class Config(object):
	'''
	Contains configuration information for IQMon package that can be passed to methods and functions.
	
	Properties:
	- location of configuration file for IQMon
	'''
	def __init__(self):
		'''
		Read and parse configuration file.
		- Currently assumes that file is .IQMonConfig in the user's home directory.
		- No defaults set, if config file not read, then values default to None.
		'''
		IQMonExecPath = None
		LogPath = None
		PlotsPath = None
		tmpPath = None
		PythonPath = None
		DataPath = None
		CatalogPath = None
		
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
			IsIQMonExecPath = re.match("IQMONPATH\s=\s([a-zA-z0-9/\-_\.]+)", line)
			if IsIQMonExecPath: IQMonExecPath = os.path.abspath(IsIQMonExecPath.group(1))
			IsLogPath = re.match("IQMONLOGS\s=\s([a-zA-z0-9/\-_\.]+)", line)
			if IsLogPath: LogPath = os.path.abspath(IsLogPath.group(1))
			IsPlotsPath = re.match("IQMONPLOTS\s=\s([a-zA-z0-9/\-_\.]+)", line)
			if IsPlotsPath: PlotsPath = os.path.abspath(IsPlotsPath.group(1))
			IstmpPath = re.match("IQMONTMP\s=\s([a-zA-z0-9/\-_\.]+)", line)
			if IstmpPath: tmpPath = os.path.abspath(IstmpPath.group(1))
			IsPythonPath = re.match("IQMONPYTHON\s=\s([a-zA-z0-9/\-_\.]+)", line)
			if IsPythonPath: PythonPath = os.path.abspath(IsPythonPath.group(1))
			IsDataPath = re.match("VYSOS20DATAPATH\s=\s([a-zA-z0-9/\-_\.]+)", line)
			if IsDataPath: DataPath = os.path.abspath(IsDataPath.group(1))
			IsCatalogPath = re.match("CATALOGPATH\s=\s([a-zA-z0-9/\-_\.]+)", line)
			if IsCatalogPath: CatalogPath = os.path.abspath(IsCatalogPath.group(1))
		
		self.IQMonExecPath = IQMonExecPath
		self.LogPath       = LogPath
		self.PlotsPath     = PlotsPath
		self.tmpPath       = tmpPath
		self.PythonPath    = PythonPath
		self.DataPath      = DataPath
		self.CatalogPath   = CatalogPath 

##--------------------------------------------------------------------------------------------------
## Define Telescope object to hold telescope information
##--------------------------------------------------------------------------------------------------
class Telescope(object):
	'''
	Contains information about the telescope that can be passed to methods and functions.  The concept
	for operation is that the user will write a simple script which creates a telescope object and
	assigned values to all it's properties (or sets them to None).  The processing is then just a 
	sequence of calls to the IQMon.Image object and it's methods.

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
		if not cls._singletons.has_key(cls):
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

	def CheckUnits(self):
		'''
		Method to check whether the properties of the object have units and to add units if the
		value has no units.  Does not yet check whether the unit is reasonable (i.e. does
		focalLength have length units).
		'''
		## name is a string
		## longName is a string
		## Default focalLength to units of mm
		if self.focalLength and not hasattr(self.focalLength, 'unit'):
			print("focalLength is unitless value.  Adding units of mm.")
			self.focalLength = self.focalLength * u.mm
		## Default pixelSize to units of microns
		if self.pixelSize and not hasattr(self.pixelSize, 'unit'):
			print("pixelSize is unitless value.  Adding units of microns")
			self.pixelSize = self.pixelSize * u.micron
		## Default aperture to units of mm
		if self.aperture and not not hasattr(self.aperture, 'unit'):
			print("aperture is unitless value.  Adding units of mm")
			self.aperture = self.aperture * u.mm
		## Default gain to units of 1/ADU
		if self.gain and not not hasattr(self.gain, 'unit'):
			print("gain is unitless value.  Adding units of 1/ADU")
			self.gain = self.gain / u.adu
		## Default nXPix to units of pixels
		if self.nXPix and not not hasattr(self.nXPix, 'unit'):
			print("nXPix is unitless value.  Adding units of pixels")
			self.nXPix = self.nXPix * u.pix
		## Default nYPix to units of pixels
		if self.nYPix and not not hasattr(self.nYPix, 'unit'):
			print("nYPix is unitless value.  Adding units of pixels")
			self.nYPix = self.nYPix * u.pix
		## Default unitsForFWHM to units of arcsec
		if self.unitsForFWHM and not not hasattr(self.unitsForFWHM, 'unit'):
			print("unitsForFWHM is unitless value.  Adding units of arcsec")
			self.unitsForFWHM = self.unitsForFWHM * u.arcsec
		## ROI is string
		## Default thresholdFWHM to units of arcsec
		if self.thresholdFWHM and not not hasattr(self.thresholdFWHM, 'unit'):
			print("thresholdFWHM is unitless value.  Adding units of arcsec")
			self.thresholdFWHM = self.thresholdFWHM * u.arcsec
		## Default thresholdPointingErr to units of arcmin
		if self.thresholdPointingErr and not not hasattr(self.thresholdPointingErr, 'unit'):
			print("thresholdPointingErr is unitless value.  Adding units of arcmin")
			self.thresholdPointingErr = self.thresholdPointingErr * u.arcmin
		## Default thresholdEllipticity to dimensionless
		if self.thresholdEllipticity and not not hasattr(self.thresholdEllipticity, 'unit'):
			print("thresholdEllipticity is unitless value.  Adding units of dimensionless")
			self.thresholdEllipticity = self.thresholdEllipticity * u.dimensionless_unscaled
		## Default pixelScale to units of arcsec per pixel
		if self.pixelScale and not not hasattr(self.pixelScale, 'unit'):
			print("pixelScale is unitless value.  Adding units of arcsec / pixel")
			self.pixelScale = self.pixelScale * u.arcsec / u.pix
		## Default fRatio to dimensionless
		if self.fRatio and not not hasattr(self.fRatio, 'unit'):
			print("fRatio is unitless value.  Adding units of dimensionless")
			self.fRatio = self.fRatio * u.dimensionless_unscaled

			
		

##--------------------------------------------------------------------------------------------------
## Define Image object which holds information on the image and methods for analysis
##--------------------------------------------------------------------------------------------------
class Image(object):
	'''
	The IQMon.Image object represents a single input image to the IQMon process.

	When defined, the image objects requires both a filename to a valid fits file and an IQMon.Config object.

	Properties:
	- original file name
	- working file name
	- target object name
	- target RA
	- target Dec
	- target alt
	- target az
	- image WCS
	- image exposure time
	- moon alt
	- moon separation from target
	- moon illumination

	Methods:
	- extract header info
	- read input file (rfits, dcraw, etc.)
	- dark subtract image
	- solve image using astrometry.net
	- crop image
	- filter cosmic rays
	- refine WCS
	- find stars (run sextractor) (detemine FWHM, ellipticity)
	- determine pointing error
	- determine zero point
	- make image plots
	- make jpeg (full frame or cropped)
	'''
	def __init__(self, filename):
		if os.path.exists(filename):
			self.filename = filename
		else:
			self.filename = None
			raise IOError("File {0} does not exist".format(filename))

	def GetHeader(self):
		'''
		Get information from the image fits header.
		'''
		pass

	def ReadImage(self):
		'''
		Read the raw image and write out a working image in the IQMon temporary directory.
		'''
		pass

	def DarkSubtract(self):
		'''
		Create master dark and subtract from image.
		'''
		pass

	def FilterCosmicRays(self):
		'''
		Use IRAF to filter cosmic rays from image (not supported)
		'''
		pass

	def Crop(self):
		'''
		Crop working image to region of interest
		'''
		pass

	def SolveAstrometry(self):
		'''
		Solve astrometry in the working image using the astrometry.net solver.
		'''
		pass

	def RefineWCS(self):
		'''
		Refine the WCS of the image to have accurate distortions.
		'''
		pass

	def DeterminePointingError(self):
		'''
		Determine pointing error (difference between objects coordinates and solved WCS).
		'''
		pass

	def RunSExtractor(self):
		'''
		Run SExtractor on image.
		'''
		pass

	def DetermineZeroPoint(self):
		'''
		Determine zero point by comparing measured magnitudes with catalog magnitudes.
		'''
		pass

	def MakePlots(self):
		'''
		Make plots for image.
		'''
		pass

	def MakeJpegs(self):
		'''
		Make jpegs of cropped version and full frame of image.
		'''
		pass
		

class ConfigTests(unittest.TestCase):
	def setUp(self):
		pass



class ImageTests(unittest.TestCase):
	def setUp(self):
		pass



if __name__ == '__main__':
	unittest.main()