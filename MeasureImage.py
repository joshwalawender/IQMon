#!/usr/bin/env python
# encoding: utf-8
"""
MeasureImage.py

Created by Josh Walawender on 2012-10-19.
Copyright (c) 2012 . All rights reserved.
"""

import sys
import os
import getopt
import time
import datetime
import re
import StringIO
import subprocess32
import shutil
import math
import numpy
import logging

import pylab
import matplotlib.pyplot as pyplot

import astropy
import astropy.coordinates
from astropy import units as u
import astropy.io
import astropy.io.ascii
import astropy.table

from pyraf import iraf
import pyfits
import ephem

import IQMonTools
import SiteCustomization
import Focus

help_message = '''
This program is intended to be executed any time a new fits file appears in a watched directory.
The program will open the fits file, find stars, measure the average FWHM of the stars and record 
the image quality in a log file.
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg	

##############################################################
##
##  Main Program
##
##############################################################
def main(argv=None):
	StartTime = time.time()
	
	## Initialize Some Variables
	infile = ""
	MakeImagePlots = False
	CleanUp        = True
	CosmicRayClean = False
	BadThings      = False
	ExistingWCS    = True
	WCSFail        = False
	Clobber        = False
	NoCoords       = False
	ForceWCSSolve  = False
	NoJPEG         = False
	AstrometryNetSolved = False
	RefineWCS      = False
	MeasureZeroPoint = False
	DarkSubtract   = True
	ImageIsCropped = False
	telescope      = ""
	pyplot.ioff()
	
	
	
	##############################################################
	## Parse Command Line Arguments and Options
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hpi:jwzd ", ["help", "MakePlots", "input=", "no-jpeg", "wcs", "zeropoint", "no-dark", "clobber", "cosmicrays", "no-cleanup", "refineWCS"])
		except getopt.error, msg:
			raise Usage(msg)
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-p", "--MakePlots"):
				MakeImagePlots = True
			if option in ("-i", "--input"):
				infile = value
			if option in ("-j", "--no-jpeg"):
				NoJPEG = True
			if option in ("-w", "--wcs"):
				ForceWCSSolve = True
			if option in ("-z", "--zeropoint"):
				MeasureZeroPoint = True
			if option in ("-d", "--no-dark"):
				DarkSubtract = False
			if option in  ("--clobber"):
				Clobber = True
			if option in ("--cosmicrays"):
				CosmicRayClean = True
			if option in ("--no-cleanup"):
				CleanUp = False
			if option in ("--refineWCS"):
				RefineWCS = True
			
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	if infile == "":
		infile = args[0]
		
	if len(args) > 1:
		print "Only One Input Image Allowed, %d found" % len(args)

	## FITS File
	FitsFile  = os.path.abspath(infile)
	if not os.path.exists(FitsFile):
		print "Unable to find input file: %s" % FitsFile
		BadThings = True
	FitsFileDirectory, FitsFilename = os.path.split(FitsFile)
	FitsBasename, FitsExt   = os.path.splitext(FitsFilename)


	
	##############################################################
	## 1)  Read Configuration File to get locations of key IQMon files and folders
	IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, LogBuffer = IQMonTools.ReadConfigFile()
	


	##############################################################
	## 2)  Get Properties of the Telescope
	## - uses information in SiteCustomization.py
	## - user should edit SiteCustomization.py to set values
	tel = SiteCustomization.GetTelescopeProperties(FitsFileDirectory, FitsFilename, IQMonExecPath)
	
	
	
	##############################################################
	## 3)  Form File Names for Summary and Log Files and various temporary files
	## - determine the UT date of the image by looking at the path of the input file
	IQMonLogFileName  = os.path.join(LogPath, tel['LongName'], tel['DataNight']+"_"+tel['name']+"_IQMonLog.txt")
	SummaryFileName   = os.path.join(LogPath, tel['LongName'], tel['DataNight']+"_"+tel['name']+"_Summary.txt")
	SummaryHTMLFile   = os.path.join(LogPath, tel['LongName'], tel['DataNight']+"_"+tel['name']+".html")
	## Individual Image Plot Filenames
	Plot_filename        = os.path.join(PlotsPath, tel['name']+"_"+tel['DataNight'], tel['FileID']+".png")
	FWHM_filename        = os.path.join(PlotsPath, tel['name']+"_"+tel['DataNight'], tel['FileID']+"_FWHM_Histogram.png")
	ZeroPoint_filename   = os.path.join(PlotsPath, tel['name']+"_"+tel['DataNight'], tel['FileID']+"_ZeroPoint.png")
	
	TemporaryHTMLFile = os.path.join(tmpPath, "tmp_"+tel['DataNight']+".html")
	ImageLink  = 'http://72.235.176.178/IQMon_Plots/'+tel['name']+"_"+tel['DataNight']+"/"+tel['FileID']+".jpg"
	CropImageLink  = 'http://72.235.176.178/IQMon_Plots/'+tel['name']+"_"+tel['DataNight']+"/"+tel['FileID']+"_crop.jpg"
	GraphsLink = 'http://72.235.176.178/IQMon_Plots/'+tel['name']+"_"+tel['DataNight']+"/"+tel['FileID']+".png"
	HistLink   = 'http://72.235.176.178/IQMon_Plots/'+tel['name']+"_"+tel['DataNight']+"/"+tel['FileID']+"_FWHM_Histogram.png"
	
	# if DoFocus:
	# 	IQMonLogFileName = os.path.join(LogPath, tel['LongName'], "FocusLog_"+tel['name']+"_"+tel['DataNight']+".txt")
	# 	SummaryFileName  = os.path.join("/Volumes", "Data_"+tel['name'], "Logs", tel['DataNight'], "FocusSummary_"+tel['name']+"_"+tel['DataNight']+".txt")
	# 	DOSFileName      = os.path.join("/Volumes", "Data_"+tel['name'], "Logs", tel['DataNight'], "FocusSummary_"+tel['name']+"_"+tel['DataNight']+"_dos.txt")
	if Clobber:
		if os.path.exists(IQMonLogFileName): os.remove(IQMonLogFileName)
		if os.path.exists(SummaryFileName): os.remove(SummaryFileName)
		if os.path.exists(SummaryHTMLFile): os.remove(SummaryHTMLFile)
	## form filenames for temporary fits files
	if not BadThings:
		WorkingFilename = tel['FileID']+"_working.fits"
		WorkingFile     = os.path.join(tmpPath, WorkingFilename)
		CroppedFilename = tel['FileID']+"_crop.fits"
		CroppedFile     = os.path.join(tmpPath, CroppedFilename)
		RFITSFilename   = tel['FileID']+"_rfits.fits"
		RFITSFile       = os.path.join(tmpPath, RFITSFilename)
		DarkSubFilename = tel['FileID']+"_darksub.fits"
		DarkSubFile     = os.path.join(tmpPath, DarkSubFilename)
		CRFilename      = tel['FileID']+"_cleaned.fits"
		CRFile          = os.path.join(tmpPath, CRFilename)
	## Set temporary file names
	if not BadThings:
		SExtractorDefaultFile    = os.path.join(IQMonExecPath, "default.sex")
		PhotometryCatalogFile_xy = os.path.join(tmpPath, tel['FileID']+"PhotCat_xy.txt")
		centerpix_file           = os.path.join(tmpPath, tel['FileID']+"_centerpix.txt")
		phot_file                = os.path.join(tmpPath, tel['FileID']+"_phot.txt")
		psfmeasurein_file        = os.path.join(tmpPath, tel['FileID']+"_psfmeasurein.txt")
		psfmeasureout_file       = os.path.join(tmpPath, tel['FileID']+"_psfmeasureout.txt")
		graphcur_file            = os.path.join(tmpPath, tel['FileID']+"_graphcur.txt")
		SExtractorConfigFile     = os.path.join(tmpPath, tel['FileID']+".sex")
		SExtractorCatalog        = os.path.join(tmpPath, tel['FileID']+".cat")
	## Set up IO strings for errors
	IRAFerrors = StringIO.StringIO()

	## Make sure paths for files exist
	if not BadThings:
		if not os.path.exists(LogPath): 
			os.mkdir(LogPath)
		if not os.path.exists(os.path.join(LogPath, tel['LongName'])): 
			os.mkdir(os.path.join(LogPath, tel['LongName']))
		if not os.path.exists(tmpPath): 
			os.mkdir(tmpPath)
		if not os.path.exists(PlotsPath): 
			os.mkdir(PlotsPath)
		if not os.path.exists(os.path.join(PlotsPath, tel['name']+"_"+tel['DataNight'])): 
			os.mkdir(os.path.join(PlotsPath, tel['name']+"_"+tel['DataNight']))



	##############################################################
	## 4)  Configure Log and Write Staring Info to Log
	##     - set logging to go to file based on telescope and date
	##     - also set up console to print log messages to screen
	logger = logging.getLogger('IQMonLogger')
	logger.setLevel(logging.DEBUG)
	LogFileHandler = logging.FileHandler(IQMonLogFileName)
	LogFileHandler.setLevel(logging.DEBUG)
	LogConsoleHandler = logging.StreamHandler()
	LogConsoleHandler.setLevel(logging.INFO)
	LogFormat = logging.Formatter('%(asctime)23s %(levelname)8s: %(message)s')
	LogFileHandler.setFormatter(LogFormat)
	LogConsoleHandler.setFormatter(LogFormat)
	logger.addHandler(LogConsoleHandler)
	logger.addHandler(LogFileHandler)
		
	logger.info("###### Processing Image:  %s ######", FitsFilename)
	## Check to Make sure input Fits File Exists
	if os.path.exists(FitsFile):
		logger.info("  Found file %s ", FitsFile)
	else:
		logger.warning("  Could not find file %s", FitsFile)
		BadThings = True
	logger.info("  Setting telescope variable to %s", tel['name'])
	for line in LogBuffer:
		logger.info(line)



	##############################################################
	##  5)  Determine Basic Properties of Image from Header
	##  - exptime
	##  - telescope
	##  - observatory
	##  - object
	##  - image size in pixels
	##  - Whether it has a WCS (crpix and crval)
	##  - Pointing RA, Dec
	##  - Exposure Start Date and Time
	##  - Observatory Location (lat, long, alt)
	if not BadThings:
		## Get image data and header
		hdulist = pyfits.open(FitsFile)
		header = hdulist[0].header
		image = hdulist[0].data
		hdulist.close()
		logger.info("  Reading image header")

		## Get exposure time from header
		try:
			exptime = header['EXPTIME']
			logger.info("  Exposure Time = %.1f" % exptime)
		except:
			logger.info("  No Exposure Time Value Found in Header")
			exptime = 0.0		
		## Get filter from header
		try:
			CCDfilter = header['FILTER']
			logger.info("  Filter = %s" % CCDfilter)
		except:
			logger.info("  No Filter Keyword Found in Header")
			CCDfilter = 'PSr'
		## Get focus position from header
		try:
			focuspos = header['FOCUSPOS']
			logger.info("  Focus Position = %d" % focuspos)
		except:
			logger.info("  No Focus Position Value Found in Header")
			focuspos = 0		
		## Get telescope from header
		try:
			headertel = header['TELESCOP']
			logger.info("  Header Telescope = "+headertel)
		except:
			headertel = ""
		## Get observatory from header
		try:
			headerobs = header['OBSERVAT']
			logger.info("  Header Observatory = "+headerobs)
		except:
			headerobs = ""
		## Get object name from header
		try:
			headerobject = "OBJECT"
			logger.info("  Header Object Name = "+headerobject)
		except:
			headerobject = ""
		
		## Determine Image Size in Pixels
		NYPix, NXPix = image.shape
		XCenterPix = NXPix/2
		YCenterPix = NYPix/2
		## Is the Image cropped or binned in camera?
		if NXPix < tel['FullXPix'] or NYPix < tel['FullYPix']:
			ImageIsCropped = True
		
				
		## Check for WCS
		WCS = IQMonTools.ExtractWCS(header)
		if WCS['Extracted'] == True:
			logger.info("  Found WCS in Header.")
			logger.info("  Image Orientation is %s and the image flip state is %s" % (WCS['Orientation'], WCS['Flipped']))
			logger.info("  Image Position Angle is %.1f" % WCS['PA'])
			ExistingWCS = True
		else:
			logger.info("  No WCS in Header.")
			ExistingWCS = False
				
		## Create ephem style Object
		try:
			ImageRA  = header['RA']
			NoCoords = False
			if len(ImageRA.split(":")) != 3:
				if len(ImageRA.split(" ")) == 3:
					ImageRA = ":".join(ImageRA.split(" "))
			ImageDEC = header['DEC']	
			if len(ImageDEC.split(":")) != 3:
				if len(ImageDEC.split(" ")) == 3:
					ImageDEC = ":".join(ImageDEC.split(" "))
			logger.info("  Reading pointing info from header: "+ImageRA+" "+ImageDEC)
			
			try:
				HeaderCoordinate = astropy.coordinates.ICRSCoordinates(ImageRA+" "+ImageDEC, unit=(u.hour, u.degree))
			except:
				BadCoords = True
		except:
			ImageRA  = "0:00:00.0"
			ImageDEC = "0:00:00.0"

		## Create ephem style Site
		site = ephem.Observer()
		## Grab UT Date and Time of Exposure Start
		ExpStartTime = header['DATE-OBS']
		SiteDate = "/".join(ExpStartTime[0:10].split("-"))
		SiteTime = ExpStartTime[11:]		
		site.date = ephem.Date(SiteDate+" "+SiteTime)
		## Latitude
		try:
			Latitude = header['LAT-OBS']
			site.lat = str(Latitude)
		except KeyError:
			Latitude = float('NaN')
		##   Longitude
		try:
			Longitude = header['LONG-OBS']
			site.lon = str(Longitude)
		except KeyError:
			Longitude = float('NaN')
		##   Elevation
		try:
			Elevation = header['ALT-OBS']
			site.elevation = Elevation
		except KeyError:
			Elevation = float('NaN')
		
		## Calculate Elevation and Azimuth of Objects at Time of Observation
		TargetObject = ephem.readdb("Target,f|M|F7,"+ImageRA+","+ImageDEC+",2.02,2000")
		TheMoon = ephem.Moon()
		TargetObject.compute(site)
		TheMoon.compute(site)
		TargetAlt = TargetObject.alt * 180. / ephem.pi
		TargetAz = TargetObject.az * 180. / ephem.pi
		MoonPhase = TheMoon.phase
		MoonSep = ephem.separation(TargetObject, TheMoon) * 180. / ephem.pi
		MoonAlt = TheMoon.alt * 180. / ephem.pi
		logger.info("  A %0.f percent illuminated Moon is %.1f degrees from the target." % (MoonPhase, MoonSep))
		
		## Calculate Airmass
		## - take value from header if available
		## - if not, calculate from Target Altitude
		if not NoCoords:
			try:
				airmass = header['AIRMASS']
			except KeyError:
				ZenithAngle = math.radians(90. - TargetAlt)
				airmass = 1.0/math.cos(ZenithAngle)*(1.0 - 0.0012*(1.0/(math.cos(ZenithAngle))**2 - 1.0))

		## Write Info About Image to HTML Log
		logger.info("  Writing Preliminary HTML Line")
		if os.path.exists(SummaryHTMLFile):
			shutil.copy(SummaryHTMLFile, TemporaryHTMLFile)
			ExistingHTML = open(TemporaryHTMLFile, 'r')
			HTML = open(SummaryHTMLFile, 'w')
			MatchEndTable = re.compile("\s*</table>")
			MatchEndBody  = re.compile("</body>")
			MatchEndHTML  = re.compile("</html>")
			for line in ExistingHTML:
				IsEndTable = MatchEndTable.match(line)
				IsEndBody  = MatchEndBody.match(line)
				IsEndHTML  = MatchEndHTML.match(line)
				if not IsEndTable and not IsEndBody and not IsEndHTML:
					HTML.write(line)
			ExistingHTML.close()
			os.remove(TemporaryHTMLFile)
		else:
			HTML = open(SummaryHTMLFile, 'w')
			HTMLheader = open(os.path.join(IQMonExecPath, "header.html"), 'r')
			header = HTMLheader.read()
			header = header.replace("telescopename", tel['LongName'])
			header = header.replace("datanight", tel['DataNight'])
			header = header.replace("FWHMunits", tel['FWHMunits'])
			HTMLheader.close()
			HTML.write(header)
		
		## Add Info for This Image
		HTML.write("    <tr>\n")
		HTML.write("      <td style='color:black;text-align:left'>%-21s</td>\n" % ExpStartTime.strip())
		## Write Image Name Field
		if MakeImagePlots and not NoJPEG:
			HTML.write("      <td style='color:black;text-align:left'><a href='%s'>%-50s</a> (<a href='%s'>C</a>) (<a href='%s'>G</a>) (<a href='%s'>H</a>)</td>\n" % (ImageLink, FitsFilename.strip(), CropImageLink, GraphsLink, HistLink))
		elif not MakeImagePlots and not NoJPEG:
			HTML.write("      <td style='color:black;text-align:left'><a href='%s'>%-50s</a> (<a href='%s'>C</a>)</td>\n" % (ImageLink, FitsFilename.strip(), CropImageLink))
		elif MakeImagePlots and NoJPEG:
			HTML.write("      <td style='color:black;text-align:left'>%-50s (<a href='%s'>G</a>) (<a href='%s'>H</a>)</td>\n" % (FitsFilename.strip(), GraphsLink, HistLink))
		elif not MakeImagePlots and NoJPEG:
			HTML.write("      <td style='color:black;text-align:left'>%-50s</td>\n" % (FitsFilename.strip()))
		## Write Target Alt, Az, airmass fields
		HTML.write("      <td style='color:black'>%.1f</td>\n" % TargetAlt)
		HTML.write("      <td style='color:black'>%.1f</td>\n" % TargetAz)
		HTML.write("      <td style='color:black'>%.2f</td>\n" % airmass)
		## Write Moon Separation and Illumination Fields
		HTML.write("      <td style='color:black'>%.1f</td>\n" % MoonSep)
		HTML.write("      <td style='color:black'>%.1f</td>\n" % MoonPhase)
	
		HTML.write("    </table>\n")
		HTML.write("</body>\n")
		HTML.write("</html>\n")
		HTML.close()
	
	
	##############################################################
	## 6)  Convert Fits File using RFITS in IRAF
	if not BadThings:
		if os.path.exists(WorkingFile): os.remove(WorkingFile)
		try:
			logger.info("  Converting FITS file using IRAF.RFITS ...")
			RFITSoutput = iraf.rfits(fits_file=FitsFile, file_list="0", iraf_file=WorkingFile, datatype="u", Stdout=1)
		except:
			for line in RFITSoutput:
				logger.warning("  IRAF.RFITS failed.")
		## Check to see if new file created
		if os.path.exists(WorkingFile):
			logger.info("  Converted FITS file written to %s", WorkingFilename)
		else:
			logger.warning("  Problem Running IRAF.RFITS")
			BadThings = True



	##############################################################
	## 7)  Dark Subtract Image
	if DarkSubtract and not BadThings:
		## Create Master Dark
		MasterDarkFile = SiteCustomization.ObtainMasterDark(tel['name'], tel['DataNight'], exptime, tmpPath, FitsFileDirectory, logger)
		if MasterDarkFile != "":
			logger.info("  Master Dark File is %s" % MasterDarkFile)
			logger.info("  Subtracting Master Dark from Image.")
			if os.path.exists(DarkSubFile): os.remove(DarkSubFile)
			iraf.imarith(WorkingFile, "-", MasterDarkFile, DarkSubFile)
			if os.path.exists(DarkSubFile):
				if os.path.exists(WorkingFile): os.remove(WorkingFile)
				shutil.copy2(DarkSubFile, WorkingFile)
		else:
			logger.warning("  Problem creating master dark file.  Continuing without dark subtraction.")
			


	##############################################################
	## 8)  Add WCS from astrometry.net solver if needed
	if not ExistingWCS:
		logger.info("  No exisiting WCS found in FITS file.")
	if ForceWCSSolve or not ExistingWCS:
		try:
			AstrometryNetSolved = IQMonTools.RunAstrometryDotNet(WorkingFile, tel, logger)
		except:
			logger.warning("  Call to Astrometry.net wrapper failed.")
		if AstrometryNetSolved:
			NewFitsFile = os.path.join(tmpPath, WorkingFile.replace(".fits", ".new.fits"))
			if os.path.exists(WorkingFile): os.remove(WorkingFile)			
			shutil.copy2(NewFitsFile, WorkingFile)
			hdulist = pyfits.open(WorkingFile)
			header = hdulist[0].header
			hdulist.close()
			## Check for WCS
			WCS = IQMonTools.ExtractWCS(header)
			if WCS['Extracted'] == True:
				logger.info("  Found WCS in Header after Astrometry.net solve.")
				logger.info("  Image Orientation is %s" % WCS['Orientation'])
				ExistingWCS = True
			else:
				logger.warning("  No WCS in Header after Astrometry.net solve.")
				ExistingWCS = False
		else:
			logger.warning("  Astrometry.net solve failed.")
			ExistingWCS = False


	##############################################################
	## 9)  If ROI is Called for, crop file using IRAF.IMCOPY
	if tel['UseROI'] and not ImageIsCropped:
		## Rename PWMvalue to CroppedFile
		os.rename(WorkingFile, CroppedFile)
		logger.info("  Cropping FITS file using IRAF.IMCOPY with ROI %s", tel['ROI'])
		IMCOPYoutput = iraf.imcopy(input=CroppedFile+tel['ROI'], output=WorkingFile, Stdout=1)
		## Get image data and header for cropped image
		hdulist = pyfits.open(WorkingFile)
		header = hdulist[0].header
		image = hdulist[0].data
		hdulist.close()
		## Determine Image Size in Pixels
		NYPix, NXPix = image.shape
		XCenterPix = NXPix/2
		YCenterPix = NYPix/2
		## Check for WCS
		WCS = IQMonTools.ExtractWCS(header)
		if WCS['Extracted'] == True:
			logger.info("  Found WCS in Header after cropping.")
			ExistingWCS = True
		else:
			logger.info("  No WCS in Header after cropping.")
			ExistingWCS = False		
		
		
		
	##############################################################
	## 10)  Use IRAF to filter out cosmic rays from image
	if CosmicRayClean and not BadThings:
		if os.path.exists(CRFile): os.remove(CRFile)
		Imstat = iraf.imstat(images=WorkingFile, fields="stddev", 
		                     nclip=4, lsigma=5, usigma=5, format="no", Stdout=1)
		ImageSTDDEV = float(Imstat[0])
		iraf.imred(_doprint=0)
		iraf.crutil(_doprint=0)
		logger.info("  Cleaning Cosmic Rays ...")
		iraf.cosmicrays.threshold=5.0*ImageSTDDEV
		iraf.cosmicrays.fluxratio=6.0
		iraf.cosmicrays.interactive="no"
		CRout = iraf.cosmicrays(input=WorkingFile, output=CRFile, answer="yes", Stdout=1, Stderr=IRAFerrors)
		if os.path.exists(CRFile):
			if os.path.exists(WorkingFile): os.remove(WorkingFile)
			shutil.copy2(CRFile, WorkingFile)
			logger.info("  Cleaned FITS image written to %s", CRFilename)
		else:
			logger.warning("  Problem Running IRAF.COSMICRAYS")
			BadThings = True
		


	##############################################################
	## 11) Refine Existing WCS
	if ExistingWCS and RefineWCS and not BadThings:
		IQMonTools.RefineWCS(WorkingFile, tel, WCS, NXPix, NYPix, logger)
	
	

	##############################################################
	## 12)  Determine Pointing Error
	## - compare WCS coordinates for center pixel, to target coordinates in image header
	if ExistingWCS and not BadThings:
		logger.info("  Detemining pointing error based on WCS solution.")
		if os.path.exists(centerpix_file): os.remove(centerpix_file)
		CenterPixFile = open(centerpix_file, 'w')
		CenterPixFile.write("%d  %d\n" % (XCenterPix, YCenterPix))
		CenterPixFile.close()
		CenterCoords = iraf.wcsctran(input=centerpix_file, output="STDOUT", image=WorkingFile, 
		                             inwcs="logical", outwcs="world", format="%H %h", 
		                             Stdout=1, Stderr=IRAFerrors)
		MatchRADEC = re.compile("\s?([0-9]{1,2}:[0-9]{2}:[0-9\.]{2,9})\s+([0-9\-\+]{1,3}:[0-9]{2}:[0-9\.]{2,9})")
		IsRADEC = MatchRADEC.match(CenterCoords[3])
		if IsRADEC:
			WCS_RA = IsRADEC.group(1)
			WCS_DEC = IsRADEC.group(2)
			WCSCoordinate = astropy.coordinates.ICRSCoordinates(ra=WCS_RA, dec=WCS_DEC, unit=(u.hour, u.degree))
			delta = WCSCoordinate.separation(HeaderCoordinate)
			PointingError = delta.arcmins
			logger.info("  Target Coordinates are:  %s %s", HeaderCoordinate.ra.format(u.hour, sep=":"), HeaderCoordinate.dec.format(u.degree, sep=":", alwayssign=True))
			logger.info("  WCS of Central Pixel is: %s %s", WCSCoordinate.ra.format(u.hour, sep=":"), WCSCoordinate.dec.format(u.degree, sep=":", alwayssign=True))
			logger.info("  ===> Pointing Error is %.2f arcmin", PointingError)
		else:
			logger.warning("  ERROR:  Could not determine center coords of WCS.")
			WCSFail = True
			PointingError = None
	else:
		PointingError = None
	
	
	##############################################################
	## 13)  Find Stars In Image Using SExtractor
	## - sectractor results give FWHM of stars
	## - sextractor results give instrumental magnitude which
	##   can be compared to catalog magnitude to determine zero point (only if RefineWCS selected)
	if not BadThings:
		## 13a)  Create PhotometryCatalogFile_xy file for SExtractor Association
		if MeasureZeroPoint:
			pass
			## Get magnitudes of catalog stars and give to SExtractor
		else:
			if os.path.exists(PhotometryCatalogFile_xy): os.remove(PhotometryCatalogFile_xy)
			PhotCatFileObject = open(PhotometryCatalogFile_xy, 'w')
			PhotCatFileObject.write("# No Existing WCS Found for this image\n")
			PhotCatFileObject.write("# This is a false file to keep SExtractor happy\n")
			PhotCatFileObject.write("0.0  0.0  0.0  0.0\n")
			PhotCatFileObject.close()
		
		
		## 13b) Make edits To default.sex based on telescope:
		## Read in default config file		
		DefaultConfig = open(SExtractorDefaultFile, 'r')
		NewConfig     = open(SExtractorConfigFile, 'w')
		for line in DefaultConfig:
			newline = line
			if re.match("CATALOG_NAME\s+", line):
				newline = "CATALOG_NAME     "+SExtractorCatalog+"\n"
			if re.match("PARAMETERS_NAME\s+", line):
				newline = "PARAMETERS_NAME  "+os.path.join(IQMonExecPath, "default.param")+"\n"
			if re.match("PHOT_APERTURES\s+", line):
				newline = "PHOT_APERTURES   "+str(tel['PhotAperture'])+"\n"
			if re.match("GAIN\s+", line):
				newline = "GAIN             "+str(tel['Gain'])+"\n"
			if re.match("PIXEL_SCALE\s+", line):
				newline = "PIXEL_SCALE      "+str(tel['PixelScale'])+"\n"
			if re.match("SEEING_FWHM\s+", line):
				newline = "SEEING_FWHM      "+str(tel['Seeing'])+"\n"
			if re.match("ASSOC_NAME\s+", line):
				newline = "ASSOC_NAME       "+PhotometryCatalogFile_xy+"\n"
			NewConfig.write(newline)
		DefaultConfig.close()
		NewConfig.close()


		## 13c) Run SExtractor
		logger.info("  Invoking SExtractor.")
		SExtractorCommand = ["sex", WorkingFile, "-c", SExtractorConfigFile]
		try:
			SExSTDOUT = subprocess32.check_output(SExtractorCommand, stderr=subprocess32.STDOUT, timeout=30)
			for line in SExSTDOUT.split("\n"):
				line.replace("[1A", "")
				line.replace("[1M>", "")
				if not re.match(".*Setting up background map.*", line) and not re.match(".*Line:\s[0-9]*.*", line):
					logger.info("  "+line)
		except:
			logger.warning("  SExtractor error.")
		## Extract Number of Stars from SExtractor Output
		pos = SExSTDOUT.find("sextracted ")
		IsSExCount = re.match("\s*([0-9]+)\s+", SExSTDOUT[pos+11:pos+21])
		if IsSExCount:
			nSExtracted = int(IsSExCount.group(1))
			logger.info("  SExtractor found %d sources.", nSExtracted)
		else:
			nSExtracted = None
		## Extract Background Level from SExtractor Output
		pos = SExSTDOUT.find("Background: ")
		IsSExBkgnd = re.match("\s*([0-9\.]+)\s*", SExSTDOUT[pos+11:pos+21])
		if IsSExBkgnd:
			SExBackground = float(IsSExBkgnd.group(1))
			logger.info("  SExtractor background is %f", SExBackground)
		else:
			SExBackground = None
		## Extract Background RMS from SExtractor Output
		IsSExBRMS = re.match("\s*RMS:\s([0-9\.]+)\s*", SExSTDOUT[pos+21:pos+37])
		if IsSExBRMS:
			SExBRMS = float(IsSExBRMS.group(1))
			logger.info("  SExtractor background RMS is %f", SExBRMS)
		else:
			SExBRMS = None

		## If No Output Catalog Created ...
		if not os.path.exists(SExtractorCatalog):
			logger.warning("  SExtractor failed to create catalog.")
			BadThings = True
	
	
	
	##############################################################
	## 14) Read in SExtractor Results
	if not BadThings:
		SExtractorResults = astropy.io.ascii.read(SExtractorCatalog, Reader=astropy.io.ascii.sextractor.SExtractor)
		SExImageRadius = []
		SExAngleInImage = []
		for star in SExtractorResults:
			SExImageRadius.append(math.sqrt((XCenterPix-star['X_IMAGE'])**2 + (YCenterPix-star['Y_IMAGE'])**2))
			SExAngleInImage.append(math.atan((star['X_IMAGE']-XCenterPix)/(YCenterPix-star['Y_IMAGE']))*180.0/math.pi)
		SExtractorResults.add_column(astropy.table.Column(data=SExImageRadius, name='ImageRadius'))
		SExtractorResults.add_column(astropy.table.Column(data=SExAngleInImage, name='AngleInImage'))
		nStarsSEx = len(SExtractorResults)		
		
		
		
	##############################################################
	## 15)  Determine Zero Point
	## - Read in USNO catalog (which was converted to pixel coordinates)
	## - The x pixel value is the matched value listed by sextractor ASSOC result
	if MeasureZeroPoint and not BadThings:
		MatchedData = []
		countAssoc = 0
		## Read USNO pixel coordinate file in to local variables
		USNO_x = []
		USNO_y = []
		USNO_magr = []
		USNOfile = open(USNOpixfile, 'r')
		MatchUSNO = re.compile("([0-9\.\-]{5,13})\s+([0-9\.\-]{5,13})\s+([0-9\.\-]{5,7})\s+([0-9\.\-]{5,7})")
		for line in USNOfile:
			IsMatch = MatchUSNO.match(line)
			if IsMatch:
				USNO_x.append(float(IsMatch.group(1)))
				USNO_magr.append(float(IsMatch.group(3)))
		USNOfile.close()
		for i in range(0, nStarsSEx):
			if Associated[i] != 0:
				countAssoc = countAssoc + 1
				x_distance = []
				for j in range(0, len(USNO_x)):
					x_distance.append([math.fabs(USNO_x[j]-Associated[i]), USNO_magr[j]])
				USNO_magr_associated = min(x_distance)[1]
				zp1 = USNO_magr_associated - Mag_Aper[i]
				zp2 = USNO_magr_associated - Mag_Auto[i]
				MatchedData.append([Associated[i], min(x_distance)[0], min(x_distance)[1], 
				                   Mag_Aper[i], Mag_Auto[i], zp1, zp2, USNO_magr_associated])
		nMatched = len(MatchedData)
		logger.info("  SExtractor associated %d stars to USNO catalog.", nMatched)
		
		## Calulate Zero Point for each star as difference between catalog and instrumental magnitude
		## - get two zero points.  one for MAG_APER, one for MAG_AUTO from SExtractor
		ZPbinsize = 0.1
		InstrumentalMag = [row[4] for row in MatchedData]
		CatalogMag      = [row[7] for row in MatchedData]
		ZeroPoint1_data = [row[5] for row in MatchedData]
		ZeroPoint2_data = [row[6] for row in MatchedData]
		ZeroPoint2_med = numpy.median(ZeroPoint2_data)
		ZeroPoint2_avg = numpy.average(ZeroPoint2_data)
		ZeroPoint2  = ZeroPoint2_med
		ZPErr2      = numpy.std(ZeroPoint2_data)
		ZPfit, ZPcov, nRejected = IQMonTools.FitZPSigmaReject(InstrumentalMag, CatalogMag, 2.0, 0.0003, 5, plot=False)
		logger.info("  ===> Median Zero Point SExtractor Automatic Photometry:  %.2f (with Std. Dev. of %.1f for %d stars)", ZeroPoint2_med, ZPErr2, nStarsSEx)
		logger.info("  ===> Zero Point Fit from SExtractor Photometry: %.2f (with covariance of %.1f and %d stars rejected)", ZPfit, ZPcov, nRejected)
	else:
		ZeroPoint2_med = 0.0
		ZPfit = 0.0
	
	
	
	##############################################################
	## 16)  Determine FWHM & Ellipticity from SExtractor
	if not BadThings:
		if nStarsSEx > 1:
			DiagonalRadius = math.sqrt(XCenterPix**2+YCenterPix**2)
			SquareRadius   = (XCenterPix+YCenterPix)/2.
			IQRadius = DiagonalRadius/tel['IQRadiusFactor']
			CentralFWHMs = []
			CentralEllipticities = []
			for star in SExtractorResults:
				if star['ImageRadius'] <= IQRadius:
					CentralFWHMs.append(star['FWHM_IMAGE'])
					CentralEllipticities.append(star['ELLIPTICITY'])		
			IQ_SEx_median = numpy.median(CentralFWHMs)
			IQ_SEx_modes = IQMonTools.modes(CentralFWHMs, 0.1)
			if len(IQ_SEx_modes) > 1:
				logger.warning("  Multiple Modes for FWHM found.  %d", len(IQ_SEx_modes))
				IQ_SEx_mode = numpy.mean(IQ_SEx_modes)
			else:
				IQ_SEx_mode = IQ_SEx_modes[0]
			IQ_SEx = IQ_SEx_median
			IQ_Ellipticity = numpy.median(CentralEllipticities)
			logger.info("  ===> Median FWHM in inner region is %-4.2f pixels", IQ_SEx_median)
			logger.info("  ===> Mode FWHM in inner region is %-4.2f pixels", IQ_SEx_mode)
			logger.info("  ===> Median Ellipticity in inner region is %-.2f", IQ_Ellipticity)
		else:
			IQ_SEx = None
			IQ_Ellipticity = None
			
	
	
	
	##############################################################
	## 17)  Make plots of various values
	if not BadThings and MakeImagePlots:		
		## Generate radial averages for FWHM, Elongation, Ellipticity
		logger.info("  Making FWHM, Ellipticity, and Theta Plots")
		minR = 0
		maxR = math.sqrt(NXPix**2 + NYPix**2)/2.0
		nbins = 25
		binsize = maxR/nbins
		FWHM_plot = []
		FWHMmed_plot = []
		FWHMerr_plot = []
		Radius_plot = []
		nInRadius_plot = []
		Ellipticity_plot = []
		EllipticityErr_plot = []
		FWHMvalues         = numpy.array(SExtractorResults['FWHM_IMAGE'])
		RadiusValues       = numpy.array(SExtractorResults['ImageRadius'])
		EllipticityValues  = numpy.array(SExtractorResults['ELLIPTICITY'])
		for i in range(0,nbins):
			## FWHM & Radius
			FWHMvalues_InBin   = FWHMvalues[(RadiusValues>i*binsize) & (RadiusValues<(i+1)*binsize)]
			RadiusValues_InBin = RadiusValues[(RadiusValues>i*binsize) & (RadiusValues<(i+1)*binsize)]
			FWHM_plot.append(numpy.average(FWHMvalues_InBin))
			FWHMmed_plot.append(numpy.median(FWHMvalues_InBin))
			FWHMerr_plot.append(3.0*numpy.std(FWHMvalues_InBin)/math.sqrt(len(FWHMvalues_InBin)))
			Radius_plot.append((i+0.5)*binsize)
			nInRadius_plot.append(len(FWHMvalues_InBin))
			## Ellipticity
			Ellipticity_InBin  = EllipticityValues[(RadiusValues>i*binsize) & (RadiusValues<(i+1)*binsize)]
			Ellipticity_plot.append(numpy.average(Ellipticity_InBin))
			EllipticityErr_plot.append(3.0*numpy.std(Ellipticity_InBin)/math.sqrt(len(Ellipticity_InBin)))

		## Generate averages for Theta and Position Angle plot
		minPA = -90
		maxPA = +90
		nbinsPA = 18
		binsizePA = (maxPA - minPA)/nbinsPA
		PAvalues    = numpy.array(SExtractorResults['AngleInImage'])
		ThetaValues = numpy.array(SExtractorResults['THETA_IMAGE'])
		PA_plot = []
		Theta_plot = []
		ThetaErr_plot = []
		for j in range(nbinsPA):
			PAvalues_InBin = PAvalues[(PAvalues > minPA+j*binsizePA) & (PAvalues < minPA+(j+1)*binsizePA)]
			Theta_InBin    = ThetaValues[(PAvalues > minPA+j*binsizePA) & (PAvalues < minPA+(j+1)*binsizePA)]
			PA_plot.append(minPA+(j+0.5)*binsizePA)
			Theta_plot.append(numpy.average(Theta_InBin))
			ThetaErr_plot.append(3.0*numpy.std(Theta_InBin)/math.sqrt(len(Theta_InBin)))


		## 17a)  FWHM vs. Radius Plot
		pyplot.figure(figsize=(6,8))
		pyplot.subplot(311)
		pyplot.xlabel("Radius (pixels)", size="small")
		pyplot.ylabel("FWHM (pixels)", size="small")
		pyplot.errorbar(Radius_plot, FWHM_plot, yerr=FWHMerr_plot, 
		                ecolor="blue", elinewidth=1,
		                label="Avg. FWHM (3 sigma errors)")
		pyplot.errorbar(Radius_plot, FWHMmed_plot, yerr=FWHMerr_plot, 
		                ecolor="green", elinewidth=1,
		                label="Median FWHM (3 sigma errors)")
		pyplot.plot([0,IQRadius], [IQ_SEx, IQ_SEx], 
		                color="black", linewidth=1,
		                label="Median FWHM for inner region")
		# pyplot.legend(loc="upper left", fontsize="small")
		pyplot.ylim(0,max(FWHM_plot)*1.5)
		pyplot.xlim(minR,maxR*1.02)
		pyplot.yticks(size="small")
		pyplot.xticks(size="small")
		pyplot.grid()
	
		## 17b)  Ellipticity vs. Radius Plot
		pyplot.subplot(312)
		pyplot.xlabel("Radius (pixels)")
		pyplot.ylabel("Ellipticity")
		pyplot.errorbar(Radius_plot, Ellipticity_plot, yerr=EllipticityErr_plot, ecolor="blue", elinewidth=1)
		pyplot.plot([0,IQRadius], [IQ_Ellipticity, IQ_Ellipticity], 
		                color="black", linewidth=1,
		                label="Median Ellipticity for inner region")
		pyplot.ylim(0,max(Ellipticity_plot)*1.5)
		pyplot.xlim(minR,maxR*1.02)
		pyplot.grid()
	
		## 17c)  Theta vs. Position Angle Plot
		pyplot.subplot(313)
		pyplot.xlabel("Position Angle (from Center of CCD)")
		pyplot.ylabel("Theta (PSF Angle)")
		pyplot.errorbar(PA_plot, Theta_plot, yerr=ThetaErr_plot, ecolor="blue", elinewidth=1)
		pyplot.plot([-90,90], [-90,90], color="black", linewidth=1)
		pyplot.xlim(-90,90)
		pyplot.ylim(-90,90)
		pyplot.grid()
		
		pyplot.savefig(Plot_filename, dpi=100, bbox_inches='tight', pad_inches=0.10)
				
		## 17d)  Zero Point Plots
		if MeasureZeroPoint:			
			logger.info("  Making Zero Point Histograms")
			xfit = [min(InstrumentalMag), max(InstrumentalMag)]
			yfit = [ZPfit + 1.0*min(InstrumentalMag),ZPfit + 1.0*max(InstrumentalMag)]
			
			pyplot.figure(figsize=(10,12))
			pyplot.subplot(211)
			# ZPmin = math.floor(min([min(ZeroPoint1_data), min(ZeroPoint2_data)]))
			# ZPmax = math.ceil(max([max(ZeroPoint1_data), max(ZeroPoint2_data)]))
			ZPmin = math.floor(ZeroPoint2 - 3.0*ZPErr2)
			ZPmax = math.ceil(ZeroPoint2 + 3.0*ZPErr2)
			ZPrange = ZPmax - ZPmin
			ZPnbins = math.ceil(ZPrange/ZPbinsize)
			ZPbins = pylab.linspace(ZPmin, ZPmax, ZPnbins)
			pyplot.title("Zero Points from Aperture Photometry")
			pyplot.xlabel("Zero Point (magnitude)")
			pyplot.ylabel("Number of Stars")
			# pyplot.hist(ZeroPoint1_data, bins=ZPbins, 
			#             alpha=1.0, color="red", label="Aperture Photometry")
			pyplot.hist(ZeroPoint2_data, bins=ZPbins, 
			            alpha=1.0, color="blue", label="Automatic Photometry")
			pyplot.plot([ZeroPoint2_med, ZeroPoint2_med], [0,200], 'r-')
			pyplot.plot([ZPfit, ZPfit], [0,200], 'k-')
			pyplot.xlim([ZeroPoint2_med-2.0*ZPErr2,ZeroPoint2_med+2.0*ZPErr2])
			# pyplot.legend(loc="upper left")
			pyplot.grid()
			
			pyplot.subplot(212)
			pyplot.plot(InstrumentalMag, CatalogMag, 'bo')			
			pyplot.plot(xfit, yfit, 'k-')
			pyplot.xlim(min(InstrumentalMag)-0.2, max(InstrumentalMag)+0.2)
			pyplot.ylim(min(CatalogMag)-0.2, max(CatalogMag)+0.2)
			xticklow = int(math.floor(min(InstrumentalMag)))
			xtickhigh = int(math.ceil(max(InstrumentalMag)))
			pyplot.xticks(range(xticklow, xtickhigh, 1))
			yticklow = int(math.floor(min(CatalogMag)))
			ytickhigh = int(math.ceil(max(CatalogMag)))
			pyplot.yticks(range(yticklow, ytickhigh, 1))
			pyplot.title("Photometry Correlation")
			pyplot.xlabel("Instrumental Magnitude")
			pyplot.ylabel("Catalog Magnitude")
			pyplot.grid()
			
			pyplot.savefig(ZeroPoint_filename, dpi=100, bbox_inches='tight', pad_inches=0.10)
			
			
		## 17e)  FWHM Histogram
		logger.info("  Making FWHM Histograms")
		pyplot.figure(figsize=(6,8))
		pyplot.subplot(211)
		pyplot.title("Histogram of FWHM Values", size="small")
		pyplot.xlabel("FWHM (pixels)", size="small")
		pyplot.ylabel("Number of Stars", size="small")
		FWHMbinsize = 0.10
		FWHMmin = min(SExtractorResults['FWHM_IMAGE'])
		FWHMmax = max(SExtractorResults['FWHM_IMAGE'])
		FWHMrange = FWHMmax - FWHMmin
		FWHMnbins = math.ceil(FWHMrange/FWHMbinsize)
		FWHMbins = pylab.linspace(FWHMmin, FWHMmax, FWHMnbins)
		pyplot.hist(SExtractorResults['FWHM_IMAGE'], bins=FWHMbins, log=True)
		pyplot.plot([IQ_SEx_median, IQ_SEx_median], [.1, 10000], 'r-')
		pyplot.plot([IQ_SEx_mode, IQ_SEx_mode], [.1, 10000], 'k-')
		pyplot.yticks(size="small")
		pyplot.xticks(size="small")
		pyplot.grid()

		pyplot.subplot(212)
		pyplot.title("Histogram of FWHM Values", size="small")
		pyplot.xlabel("FWHM (pixels)", size="small")
		pyplot.ylabel("Number of Stars", size="small")
		pyplot.hist(SExtractorResults['FWHM_IMAGE'], bins=FWHMbins, log=True)
		pyplot.plot([IQ_SEx_median, IQ_SEx_median], [.1, 10000], 'r-')
		pyplot.plot([IQ_SEx_mode, IQ_SEx_mode], [.1, 10000], 'k-')
		pyplot.yticks(size="small")
		pyplot.xticks(size="small")
		pyplot.grid()
		MaxPlot = math.ceil(max([IQ_SEx_median, IQ_SEx_mode])*1.5)
		pyplot.xlim(0,MaxPlot)
		
		pyplot.savefig(FWHM_filename, dpi=100, bbox_inches='tight', pad_inches=0.10)



	##############################################################
	## 18)  Make JPEGs of Image
	if not NoJPEG:
		## 18a)  Make Cropped and Marked Image
		logger.info("  Creating cropped and marked jpeg version of image for web")
		JPEGFileName2 = os.path.join(PlotsPath, tel['name']+"_"+tel['DataNight'], tel['FileID']+"_crop.jpg")
		if os.path.exists(JPEGFileName2): os.remove(JPEGFileName2)
		JPEGcommand = ["convert", "-contrast-stretch", "0.9%,1%", "-compress", "JPEG", "-quality", "70", "-stroke", "red", "-fill", "none"]
		MarkRadius=max([4, 2*math.ceil(IQ_SEx)])   ## radius of circle on jpeg is 2x seeing (in pix), but no smaller than 4 pixels
		for star in SExtractorResults:
			MarkXPos = star['X_IMAGE']
			MarkYPos = NYPix - star['Y_IMAGE']
			JPEGcommand.append('-draw')
			JPEGcommand.append("circle %d,%d %d,%d" % (MarkXPos, MarkYPos, MarkXPos+MarkRadius, MarkYPos))
		JPEGcommand.append(WorkingFile)
		JPEGcommand.append(JPEGFileName2)
		try:		
			ConvertSTDOUT = subprocess32.check_output(JPEGcommand, stderr=subprocess32.STDOUT, timeout=30)
		except:
			logger.warning("  Failed to create jpeg.")
	
		## 18b)  Make Full Field Image
		JRF = 4  ## Jpeg Reduction factor
		logger.info("  Creating full frame jpeg version of image for web")
		JPEGFileName = os.path.join(PlotsPath, tel['name']+"_"+tel['DataNight'], tel['FileID']+".jpg")
		if os.path.exists(JPEGFileName): os.remove(JPEGFileName)
		JRFstring = str(1./JRF*100)
		if ExistingWCS:
			JPEGcommand = ["convert", "-contrast-stretch", "0.9%,1%", "-compress", "JPEG", "-quality", "70", "-rotate", str(WCS['PA']) ,"-resize", JRFstring+"%"]
			if WCS['Flipped']:
				JPEGcommand.append("-flop")
		else:
			JPEGcommand = ["convert", "-contrast-stretch", "0.9%,1%", "-compress", "JPEG", "-quality", "70" ,"-resize", JRFstring+"%"]
		## Use Dark Subtracted file if it is available, otherwise use raw
		if os.path.exists(DarkSubFile):
			JPEGcommand.append(DarkSubFile)
		else:
			JPEGcommand.append(FitsFile)
		JPEGcommand.append(JPEGFileName)
		try:		
			ConvertSTDOUT = subprocess32.check_output(JPEGcommand, stderr=subprocess32.STDOUT, timeout=30)
		except:
			logger.warning("  Failed to create jpeg.")
	
	

	##############################################################
	## 19)  Cleanup Files in Temporary Directory
	if CleanUp:
		logger.info("  Cleaning Up.  Deleting temporary files.")
		if os.path.exists(WorkingFile):
			os.remove(WorkingFile)
		if os.path.exists(RFITSFile):
			os.remove(RFITSFile)
		if os.path.exists(DarkSubFile):
			os.remove(DarkSubFile)
		if os.path.exists(CroppedFile):
			os.remove(CroppedFile)
		if os.path.exists(CRFile):
			os.remove(CRFile)
		if os.path.exists(centerpix_file):
			os.remove(centerpix_file)
		if os.path.exists(phot_file):
			os.remove(phot_file)
		if os.path.exists(psfmeasurein_file):
			os.remove(psfmeasurein_file)
		if os.path.exists(psfmeasureout_file):
			os.remove(psfmeasureout_file)
		if os.path.exists(graphcur_file):
			os.remove(graphcur_file)
		if os.path.exists(SExtractorConfigFile):
			os.remove(SExtractorConfigFile)
		if os.path.exists(SExtractorCatalog):
			os.remove(SExtractorCatalog)
		if os.path.exists(PhotometryCatalogFile_xy):
			os.remove(PhotometryCatalogFile_xy)
	## Close out open streams and files
	IRAFerrors.close()
	

	
	##############################################################
	## 20)  Focus
	## Calculate Move
	logger.info("  Calculating Focus Move")
	Move = Focus.Reactive(SummaryFileName, tel)
	logger.info("  Focus Move = %d" % Move)
		
	
			
	##############################################################
	## 21)  Make Summary File Entry for This Image
	if not BadThings:
		logger.info("  Writing Summary File Entry")		
		## Read in previous data
		if not os.path.exists(SummaryFileName):
			SummaryTable = astropy.table.Table(names=("ExpStart", "File", "FWHM", "Ellipticity", 
			                            "Alt", "Az", "Airmass", "PointingError", 
			                            "ZeroPointMedian", "ZeroPointFit", "nStars", "Move", "Background", "Background RMS"),
			                     dtypes=('a22', 'a120', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'i4', 'i4', 'f4', 'f4'),
			                     masked=False)
		else:
			SummaryTable = astropy.io.ascii.read(SummaryFileName)
		## Write new line to table
		if tel['FWHMunits'] == "arcsec":
			SummaryTable.add_row((ExpStartTime, FitsFilename, IQ_SEx*tel['PixelScale'], IQ_Ellipticity, 
			                      TargetAlt, TargetAz, airmass, PointingError, 
			                      ZeroPoint2_med, ZPfit, nStarsSEx, Move, SExBackground, SExBRMS))
		else:
			SummaryTable.add_row((ExpStartTime, FitsFilename, IQ_SEx, IQ_Ellipticity, 
			                      TargetAlt, TargetAz, airmass, PointingError, 
			                      ZeroPoint2_med, ZPfit, nStarsSEx, Move, SExBackground, SExBRMS))
		## Write Table to File
		astropy.io.ascii.write(SummaryTable, SummaryFileName, Writer=astropy.io.ascii.basic.Basic)



	##############################################################
	## 22)  Determine Run Time for IQMon
	EndTime = time.time()
	ProcessTime = EndTime - StartTime
	logger.info("  IQMon Processing Time: %.1f", ProcessTime)



	##############################################################
	## 23)  Make Summary HTML File Entry for this Image
	if not BadThings:
		## Test to See if the Image Should be Flagged
		if IQ_SEx >= tel['IQLimit']:
			IQFlag = True
		else:
			IQFlag = False
		if ExistingWCS and not WCSFail:
			if PointingError >= tel['PELimit']:
				PEFlag = True
			else:
				PEFlag = False
		if IQ_Ellipticity > tel['EllipLimit']:
			EllipticityFlag = True
		else:
			EllipticityFlag = False
		
		## Write HTML
		logger.info("  Writing HTML File Entry")
		if os.path.exists(SummaryHTMLFile):
			shutil.copy(SummaryHTMLFile, TemporaryHTMLFile)
			ExistingHTML = open(TemporaryHTMLFile, 'r')
			HTML = open(SummaryHTMLFile, 'w')
			MatchEndTable = re.compile("\s*</table>")
			MatchEndBody  = re.compile("</body>")
			MatchEndHTML  = re.compile("</html>")
			for line in ExistingHTML:
				IsEndTable = MatchEndTable.match(line)
				IsEndBody  = MatchEndBody.match(line)
				IsEndHTML  = MatchEndHTML.match(line)
				if not IsEndTable and not IsEndBody and not IsEndHTML:
					HTML.write(line)
			ExistingHTML.close()
			os.remove(TemporaryHTMLFile)
		else:
			HTML = open(SummaryHTMLFile, 'w')
			HTMLheader = open(os.path.join(IQMonExecPath, "header.html"), 'r')
			header = HTMLheader.read()
			header = header.replace("telescopename", tel['LongName'])
			header = header.replace("datanight", tel['DataNight'])
			header = header.replace("FWHMunits", tel['FWHMunits'])
			HTMLheader.close()
			HTML.write(header)


		## Write FWHM Field			
		try:
			if tel['FWHMunits'] == "pixels":
				if IQFlag:
					HTML.write("      <td style='color:red'>%.2f</td>\n" % IQ_SEx)
				else:
					HTML.write("      <td style='color:black'>%.2f</td>\n" % IQ_SEx)
			elif tel['FWHMunits'] == "arcsec":
				if IQFlag:
					HTML.write("      <td style='color:red'>%.2f</td>\n" % FWHMarcsec)
				else:
					HTML.write("      <td style='color:black'>%.2f</td>\n" % FWHMarcsec)
			else:
				if IQFlag:
					HTML.write("      <td style='color:red'>%.2f</td>\n" % IQ_SEx)
				else:
					HTML.write("      <td style='color:black'>%.2f</td>\n" % IQ_SEx)
		except TypeError:
			HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		## Write Ellipticity Field			
		try:
			if EllipticityFlag:
				HTML.write("      <td style='color:red'>%.2f</td>\n" % IQ_Ellipticity)
			else:
				HTML.write("      <td style='color:black'>%.2f</td>\n" % IQ_Ellipticity)
		except TypeError:
			HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		## Write Background Field
		try:
			HTML.write("      <td style='color:black'>%.1f [%.1f]</td>\n" % (SExBackground, SExBRMS))
		except TypeError:
			HTML.write("      <td style='color:black'>%s</td>\n" % ("n/a"))
		## Write PointingError Field
		if ExistingWCS and not WCSFail:
			if PEFlag:
				HTML.write("      <td style='color:red'>%.1f</td>\n" % PointingError)
				HTML.write("      <td style='color:black'>%.1f</td>\n" % WCS['PA'])
			else:
				HTML.write("      <td style='color:black'>%.1f</td>\n" % PointingError)
				HTML.write("      <td style='color:black'>%.1f</td>\n" % WCS['PA'])
		else:
			HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
			HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		## Write Number of Stars and Process Time Fields
		HTML.write("      <td style='color:black'>%d</td>\n" % nStarsSEx)
		HTML.write("      <td style='color:black'>%.1f</td>\n" % ProcessTime)
		HTML.write("    </tr>\n")
		HTML.write("    </table>\n")
		HTML.write("</body>\n")
		HTML.write("</html>\n")
		HTML.close()
	else:
		## Write Dummy Info if BadThings flag is set
		logger.info("  Writing Dummy HTML File Entry")
		if os.path.exists(SummaryHTMLFile):
			shutil.copy(SummaryHTMLFile, TemporaryHTMLFile)
			ExistingHTML = open(TemporaryHTMLFile, 'r')
			HTML = open(SummaryHTMLFile, 'w')
			MatchEndTable = re.compile("\s*</table>")
			MatchEndBody  = re.compile("</body>")
			MatchEndHTML  = re.compile("</html>")
			for line in ExistingHTML:
				IsEndTable = MatchEndTable.match(line)
				IsEndBody  = MatchEndBody.match(line)
				IsEndHTML  = MatchEndHTML.match(line)
				if not IsEndTable and not IsEndBody and not IsEndHTML:
					HTML.write(line)
			ExistingHTML.close()
			os.remove(TemporaryHTMLFile)
		else:
			HTML = open(SummaryHTMLFile, 'w')
			HTMLheader = open(os.path.join(IQMonExecPath, "header.html"), 'r')
			header = HTMLheader.read()
			header = header.replace("telescopename", tel['LongName'])
			header = header.replace("datanight", tel['DataNight'])
			header = header.replace("FWHMunits", tel['FWHMunits'])
			HTMLheader.close()
			HTML.write(header)
		## Write FWHM Field
		HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		## Write Ellipticity Field
		HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		## Write Background Field
		HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		## Write PointingError Field
		HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		## Write nStars and ProcessTime
		HTML.write("      <td style='color:black'>%s</td>\n" % "n/a")
		HTML.write("      <td style='color:black'>%.1f</td>\n" % ProcessTime)
		HTML.write("    </tr>\n")
		HTML.write("    </table>\n")
		HTML.write("</body>\n")
		HTML.write("</html>\n")
		HTML.close()
		
	
	return None
			
	

if __name__ == "__main__":
	sys.exit(main())
