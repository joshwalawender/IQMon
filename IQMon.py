#!/usr/bin/env python
# encoding: utf-8
"""
IQMon.py

Created by Josh Walawender on 2013-07-22.
Copyright (c) 2013 . All rights reserved.
"""

import sys
import os
import re
import unittest


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