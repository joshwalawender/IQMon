#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2012-10-11.
Copyright (c) 2012 . All rights reserved.
"""

import sys
import getopt
import os
import re
import string
import fnmatch
import numpy
import math
import astropy
from astropy import io
from astropy.io import ascii

import asciitable
import matplotlib.pyplot as pyplot
import ephem
import datetime

import IQMonTools
import SiteCustomization

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def TimeStringToDecimal(TimeString):
	hms = string.split(TimeString, ":")
	DecimalTime = float(hms[0])+ float(hms[1])/60.0 + float(hms[2])/60.0/60.0
	return DecimalTime

def ConvertHSTtoUTString(TimeString):
	hmsHST = string.split(TimeString, ":")
	if int(hmsHST[0]) >= 14:
		UTString = str(int(hmsHST[0])+10-24)+":"+hmsHST[1]+":"+hmsHST[2]
	else:
		UTString = str(int(hmsHST[0])+10)+":"+hmsHST[1]+":"+hmsHST[2]
	return UTString

def CompareIQ(DateString, telescope):
	print "#### Making Nightly Plots for "+telescope+" on the Night of "+DateString+" ####"
	
	FoundACPLog       = False
	FoundCCDIFile     = False
	FoundIQMonFile    = False
	FoundFocusMaxFile = False
	FoundEnvLogFile   = False
	pyplot.ioff()
		
	##############################################################
	## Read Configuration File to get the following items
	## - IQMONEXECPATH
	## - IQMONLOGS
	## - IQMONPLOTS
	## - IQMONTMP
	IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, LogBuffer = IQMonTools.ReadConfigFile()
	
	##############################################################
	## Set up pathnames and filenames
	if telescope == "V5":
		VYSOSDATAPath = V5DataPath
		PixelScale = 2.53
		telname = "VYSOS-5"
	if telescope == "V20":
		VYSOSDATAPath = V20DataPath
		PixelScale = 0.44
		telname = "VYSOS-20"
	
	## Set File Name
	PlotFileName = DateString+"_"+telescope+".png"
	PlotFile = os.path.join(LogPath, telname, PlotFileName)
	EnvPlotFileName = DateString+"_"+telescope+"_Env.png"
	EnvPlotFile = os.path.join(LogPath, telname, EnvPlotFileName)
	
	## Compile Various Regular Expressions for File Name Matching and ACP Log Parsing
	MatchDir           = re.compile(DateString)
	MatchStartOfImage  = re.compile("(\d{2}:\d{2}:\d{2})\s*Imaging\sto\s([a-zA-Z0-9@\-_]*)")
	MatchPointingError = re.compile("(\d{2}:\d{2}:\d{2})\s*Pointing\serror\sis\s([0-9\.]+)\sarcmin.*")
	MatchImageFWHM     = re.compile("(\d{2}:\d{2}:\d{2})\s*Image\sFWHM\sis\s[0-9\.]{2,5}\sarcsec\s\(([0-9\.]{2,5})\spixels\)")
	MatchAvgFWHM       = re.compile("(\d{2}:\d{2}:\d{2})\s*\(avg\sFWHM\s=\s([0-9\.]{2,5})\sarcsec\)")
	MatchRunComplete   = re.compile("(\d{2}:\d{2}:\d{2})\s*Run\scomplete")
	MatchCCDILine      = re.compile("\d{2}/\d{2}/\d{2}\s*(\d{2}:\d{2}:\d{2})\s*([a-zA-Z0-9@_\-]{5,50}\.fts)\s+([0-9\.]{2,6})px\s+.*")
	
	
	##############################################################
	## Use pyephem determine sunrise and sunset times
	Observatory = ephem.Observer()
	Observatory.lon = "-155:34:33.9"
	Observatory.lat = "+19:32:09.66"
	Observatory.elevation = 3400.0
	Observatory.temp = 10.0
	Observatory.pressure = 680.0
	Observatory.date = DateString[0:4]+"/"+DateString[4:6]+"/"+DateString[6:8]+" 10:00:00.0"

	Observatory.horizon = '0.0'
	SunsetTime  = Observatory.previous_setting(ephem.Sun()).datetime()
	SunriseTime = Observatory.next_rising(ephem.Sun()).datetime()
	SunsetDecimal = float(datetime.datetime.strftime(SunsetTime, "%H"))+float(datetime.datetime.strftime(SunsetTime, "%M"))/60.+float(datetime.datetime.strftime(SunsetTime, "%S"))/3600.
	SunriseDecimal = float(datetime.datetime.strftime(SunriseTime, "%H"))+float(datetime.datetime.strftime(SunriseTime, "%M"))/60.+float(datetime.datetime.strftime(SunriseTime, "%S"))/3600.
	Observatory.horizon = '-6.0'
	EveningCivilTwilightTime = Observatory.previous_setting(ephem.Sun(), use_center=True).datetime()
	MorningCivilTwilightTime = Observatory.next_rising(ephem.Sun(), use_center=True).datetime()
	EveningCivilTwilightDecimal = float(datetime.datetime.strftime(EveningCivilTwilightTime, "%H"))+float(datetime.datetime.strftime(EveningCivilTwilightTime, "%M"))/60.+float(datetime.datetime.strftime(EveningCivilTwilightTime, "%S"))/3600.
	MorningCivilTwilightDecimal = float(datetime.datetime.strftime(MorningCivilTwilightTime, "%H"))+float(datetime.datetime.strftime(MorningCivilTwilightTime, "%M"))/60.+float(datetime.datetime.strftime(MorningCivilTwilightTime, "%S"))/3600.
	Observatory.horizon = '-12.0'
	EveningNauticalTwilightTime = Observatory.previous_setting(ephem.Sun(), use_center=True).datetime()
	MorningNauticalTwilightTime = Observatory.next_rising(ephem.Sun(), use_center=True).datetime()
	EveningNauticalTwilightDecimal = float(datetime.datetime.strftime(EveningNauticalTwilightTime, "%H"))+float(datetime.datetime.strftime(EveningNauticalTwilightTime, "%M"))/60.+float(datetime.datetime.strftime(EveningNauticalTwilightTime, "%S"))/3600.
	MorningNauticalTwilightDecimal = float(datetime.datetime.strftime(MorningNauticalTwilightTime, "%H"))+float(datetime.datetime.strftime(MorningNauticalTwilightTime, "%M"))/60.+float(datetime.datetime.strftime(MorningNauticalTwilightTime, "%S"))/3600.
	Observatory.horizon = '-18.0'
	EveningAstronomicalTwilightTime = Observatory.previous_setting(ephem.Sun(), use_center=True).datetime()
	MorningAstronomicalTwilightTime = Observatory.next_rising(ephem.Sun(), use_center=True).datetime()
	EveningAstronomicalTwilightDecimal = float(datetime.datetime.strftime(EveningAstronomicalTwilightTime, "%H"))+float(datetime.datetime.strftime(EveningAstronomicalTwilightTime, "%M"))/60.+float(datetime.datetime.strftime(EveningAstronomicalTwilightTime, "%S"))/3600.
	MorningAstronomicalTwilightDecimal = float(datetime.datetime.strftime(MorningAstronomicalTwilightTime, "%H"))+float(datetime.datetime.strftime(MorningAstronomicalTwilightTime, "%M"))/60.+float(datetime.datetime.strftime(MorningAstronomicalTwilightTime, "%S"))/3600.
	
	PlotStartUT = math.floor(SunsetDecimal)
	PlotEndUT = math.ceil(SunriseDecimal)
	# PlotStartUT = int(4)
	# PlotEndUT = int(16)
	nUTHours = PlotEndUT-PlotStartUT+1
	
	now = datetime.datetime.utcnow()
	DecimalTime = now.hour+now.minute/60.+now.second/3600.
		
	###########################################################
	## Read ACP Logs
	## - extract ACP FWHM and Pointing Error
	ACPLogDirectory = os.path.join(VYSOSDATAPath, "Logs", DateString)
	ACPdata = []
	if os.path.exists(ACPLogDirectory):
		print "  Found ACP Log Directory: "+ACPLogDirectory
		FoundACPLog = True
		ACPLogFiles = os.listdir(ACPLogDirectory)
		## Loop through all log files
		for LogFile in ACPLogFiles:
			if fnmatch.fnmatch(LogFile,"*.log"):
				input = open(os.path.join(ACPLogDirectory, LogFile), 'r')
				ImageFile = ""
				ImageFWHM = float("NaN")
				AvgFWHM = float("NaN")
				TimeString = ""
				TimeDecimal = float("NaN")
				ACPPointingError = float("NaN")
				for line in input:
					IsStartOfImage = MatchStartOfImage.match(line)
					if IsStartOfImage:
						if ImageFile != "":
							ACPdata.append([TimeDecimal, ImageFile, TimeString, ImageFWHM, AvgFWHM, ACPPointingError])
						ImageFile = IsStartOfImage.group(2)
						TimeString = IsStartOfImage.group(1)
						TimeDecimal = TimeStringToDecimal(TimeString)
					IsPointingError = MatchPointingError.match(line)
					if IsPointingError:
						ACPPointingError = float(IsPointingError.group(2))
					IsImageFWHM = MatchImageFWHM.match(line)
					if IsImageFWHM:
						ImageFWHM = float(IsImageFWHM.group(2))
					IsAvgFWHM = MatchAvgFWHM.match(line)
					if IsAvgFWHM:
						AvgFWHM = float(IsAvgFWHM.group(2))/PixelScale
					IsRunComplete = MatchRunComplete.match(line)
					if IsRunComplete:
						ACPdata.append([TimeDecimal, ImageFile, TimeString, ImageFWHM, AvgFWHM, ACPPointingError])
				input.close()
		ACPImages = []
		for item in ACPdata:
			if not re.match(".*Empty.*", item[1]):
				ACPImages.append(item[1])
		nACPImages = len(ACPImages)
		print "    Data for %d images (Empty filter images excluded) extracted from ACP Logs." % nACPImages
	else:
		print "  Failed to Find ACP Log Directory: "+ACPLogDirectory
			
	###########################################################
	## Read IQMon Logs
	## - extract IQMon FWHM, ellipticity, pointing error
	if os.path.exists(os.path.join(LogPath, telname)):
		print "  Found IQMon Logs "+os.path.join(LogPath, telname)
		Files = os.listdir(os.path.join(LogPath, telname))
		if telescope == "V5":
			MatchIQMonFile = re.compile("([0-9]{8}UT)_V5_Summary\.txt")
		if telescope == "V20":
			MatchIQMonFile = re.compile("([0-9]{8}UT)_V20_Summary\.txt")
		for File in Files:
			IsIQMonFile = MatchIQMonFile.match(File)
			if IsIQMonFile:
				FullIQMonFile = os.path.join(LogPath, telname, File)
				IQMonFileDate = IsIQMonFile.group(1)
				if IQMonFileDate == DateString:
					print "    IQMon Summary file for "+DateString+" is "+File
					FoundIQMonFile = True
					ColStarts = [ 0, 11, 21, 71, 83,  95, 107, 119, 131, 146, 158, 170]
					ColEnds   = [ 9, 20, 70, 82, 94, 106, 118, 130, 145, 157, 169, 181]
					IQMonSummary = astropy.io.ascii.read(FullIQMonFile, data_start=2, Reader=ascii.FixedWidth, col_starts=ColStarts, col_ends=ColEnds, guess=False, comment=";", header_start=0)

					IQMonDate        = IQMonSummary['# Exposure']
					IQMonTime        = IQMonSummary['Start']
					IQMonFile        = IQMonSummary['Image File']
					IQMonFWHM        = IQMonSummary['FWHM']
					IQMonEllipticity = IQMonSummary['Ellipticity']
					IQMonAlt         = IQMonSummary['Alt']
					IQMonAz          = IQMonSummary['Az']
					IQMonAirmass     = IQMonSummary['Airmass']
					IQMonPErr        = IQMonSummary['PointingError']
					IQMonZP1         = IQMonSummary['ZeroPoint1']
					IQMonZP2         = IQMonSummary['ZeroPoint2']
					
					IQMonData = zip(IQMonFile, IQMonTime, IQMonFWHM, IQMonPErr, IQMonZP1, IQMonZP2, IQMonEllipticity, IQMonAlt, IQMonAz, IQMonAirmass)										
					print "    Data for %d images extracted from IQMon summary file." % len(IQMonData)
	else:
		print "  Failed to Find IQMon Logs: "+os.path.join(LogPath, telname)
	
			
	###########################################################
	## Read FocusMax Logs
	## - extract times of focus operations
	UTDate = int(DateString[0:8])
	LocalNight = UTDate-1
	LocalNightString = str(LocalNight)
	FocusMaxLogPath = os.path.join(VYSOSDATAPath, "Logs", "FocusMax")
	if os.path.exists(FocusMaxLogPath):
		FMLogFiles = os.listdir(FocusMaxLogPath)
		MatchLogFile = re.compile(LocalNightString+"_([0-9]{5,6})\.log")
		RightLogFiles = []
		for file in FMLogFiles:
			IsRightLogFile = MatchLogFile.match(file)
			if IsRightLogFile:
				FoundFocusMaxFile = True
				RightLogFiles.append(file)
		for FocusMaxLogFileName in RightLogFiles:
			print "  Found FocusMax Log File: "+FocusMaxLogFileName
			FocusMaxLogFile = open(os.path.join(FocusMaxLogPath, FocusMaxLogFileName), 'r')
			MatchStartingAS  = re.compile("\s?([0-9]{1,2}:[0-9]{2}:[0-9]{2})\s\s\s\*\*\sStarting AcquireStar Sequence")
			MatchBestFocus   = re.compile("\s?([0-9]{1,2}:[0-9]{2}:[0-9]{2})\s\s\sBest Focus is:\s([0-9]{1,5})")
			MatchTemperature = re.compile("\s?([0-9]{1,2}:[0-9]{2}:[0-9]{2})\s\s\sTemperature = ([0-9]{1,3}\.[0-9]{1,3})")
			MatchASCompleted = re.compile("\s?([0-9]{1,2}:[0-9]{2}:[0-9]{2})\s\s\sAcquireStar completed")
			MatchASnotCompleted = re.compile("\s?([0-9]{1,2}:[0-9]{2}:[0-9]{2})\s\s\sAcquireStar not completed")
			Sequence = 0
			FocusRuns = []
			for line in FocusMaxLogFile:
				IsStartingAS = MatchStartingAS.match(line)
				if IsStartingAS:
					Sequence = 1
					StartTime = IsStartingAS.group(1)
					StartUT = TimeStringToDecimal(StartTime)+10.0
					if StartUT >= 24.0: StartUT = StartUT - 24.0
				if Sequence == 1:
					IsBestFocus = MatchBestFocus.match(line)
					if IsBestFocus:
						BestFocus = int(IsBestFocus.group(2))
						Sequence = 2
					IsASnotCompleted = MatchASnotCompleted.match(line)
					if IsASnotCompleted:
						EndTime = IsASnotCompleted.group(1)
						EndUT = TimeStringToDecimal(EndTime)+10.0
						if EndUT >= 24.0: EndUT = EndUT - 24.0					
						Sequence = 0
				if Sequence == 2:
					IsTemperature = MatchTemperature.match(line)
					if IsTemperature:
						FocusMaxTemperature = float(IsTemperature.group(2))
						Sequence = 3
					IsASnotCompleted = MatchASnotCompleted.match(line)
					if IsASnotCompleted:
						EndTime = IsASnotCompleted.group(1)
						EndUT = TimeStringToDecimal(EndTime)+10.0
						if EndUT >= 24.0: EndUT = EndUT - 24.0
						Sequence = 0
				if Sequence == 3:
					IsASCompleted = MatchASCompleted.match(line)
					if IsASCompleted:
						EndTime = IsASCompleted.group(1)
						EndUT = TimeStringToDecimal(EndTime)+10.0
						if EndUT >= 24.0: EndUT = EndUT - 24.0
						Sequence = 4
					IsASnotCompleted = MatchASnotCompleted.match(line)
					if IsASnotCompleted:
						Sequence = 0
				if Sequence == 4:
					FocusRuns.append([StartTime, EndTime, StartUT, EndUT, BestFocus, FocusMaxTemperature])
					print "      Focus: %8s %8s %5d %8.2f" % (StartTime, EndTime, BestFocus, FocusMaxTemperature)
					Sequence = 0
	else:
		print "  Failed to Find FocusMax Log File"
	
	
	
	###########################################################			
	## Get Environmental Data
	print "  Reading Environmental Logs"
	V20EnvLogFileName = os.path.join(V20DataPath, "Logs", DateString, "EnvironmentalLog.txt")
	V5EnvLogFileName  = os.path.join(V5DataPath,  "Logs", DateString, "EnvironmentalLog.txt")
	FoundV20EnvLogFile = False
	FoundV5EnvLogFile  = False
	FoundOtherLogFile = False
	if os.path.exists(V20EnvLogFileName):
		print "  Found VYSOS-20 Environmental Logs"
		FoundV20EnvLogFile = True
		ColStarts = [ 0, 11, 22, 32, 42, 52, 62, 72, 82,  92, 102, 112, 122, 132, 142]
		ColEnds   = [ 9, 18, 31, 41, 51, 61, 71, 81, 91, 101, 111, 121, 131, 141, 151]
		V20Night = astropy.io.ascii.read(V20EnvLogFileName, data_start=2, Reader=ascii.FixedWidth, col_starts=ColStarts, col_ends=ColEnds, guess=False, comment=";", header_start=0)
		
		V20Time      = V20Night['col1']
		V20TubeT     = V20Night['Tube']
		V20PrimaryT  = V20Night['Primary']
		V20FanPower  = V20Night['FanPwr']
		V20FocusPos  = V20Night['Focus']
		V20SkyT      = V20Night['Sky']
		V20OutsideT  = V20Night['Outside']
		V20WindSpeed = V20Night['WindSpd']
		V20Humidity  = V20Night['Humid']
		V20Wetness   = V20Night['Wetness']
		V20SkyDiff   = V20Night['Sky'] - V20Night['Outside']
		
		V20TimeDecimal = []
		for atime in V20Time:
			V20TimeDecimal.append(TimeStringToDecimal(atime[0:8]))
		
		V20WetnessFilt = []
		V20WetnessTime = []
		V20HumidityFilt = []
		for values in V20Night:
			try:
				V20WetnessFilt.append(int(values['Wetness']))
				V20WetnessTime.append(values['col1'])
				V20HumidityFilt.append(values['Humid'])
			except:
				pass
						
	if os.path.exists(V5EnvLogFileName):
		print "  Found VYSOS-5 Environmental Logs"
		FoundV5EnvLogFile = True
		ColStarts = [ 0, 11, 22, 32, 42, 52, 62, 72, 82,  92, 102, 112]
		ColEnds   = [ 9, 18, 31, 41, 51, 61, 71, 81, 91, 101, 111, 121]
		V5Night = astropy.io.ascii.read(V5EnvLogFileName, data_start=2, Reader=ascii.FixedWidth, col_starts=ColStarts, col_ends=ColEnds, guess=False, comment=";", header_start=0)
		
		V5Time      = V5Night['col1']
		V5TubeT     = V5Night['Tube']
		V5FocusPos  = V5Night['Focus']
		V5SkyT      = V5Night['Sky']
		V5OutsideT  = V5Night['Outside']
		V5WindSpeed = V5Night['WindSpd']
		V5Humidity  = V5Night['Humid']
		V5Wetness   = V5Night['Wetness']
		V5SkyDiff   = V5Night['Sky'] - V5Night['Outside']
		
		V5TimeDecimal = []
		for time in V5Time:
			V5TimeDecimal.append(TimeStringToDecimal(time[0:8]))
		
		V5WetnessFilt = []
		V5WetnessTime = []
		V5HumidityFilt = []
		for values in V20Night:
			try:
				V5WetnessFilt.append(int(values['Wetness']))
				V5WetnessTime.append(values['col1'])
				V5HumidityFilt.append(values['Humid'])
			except:
				pass
							
	if telescope == "V20" and FoundV20EnvLogFile:
		if FoundV20EnvLogFile: FoundEnvLogFile = True
		if FoundV5EnvLogFile: FoundOtherLogFile = True
		OtherTelescope = "V5"
		EnvLogTime      = numpy.array(V20TimeDecimal)
		EnvLogTubeT     = numpy.array(V20TubeT)
		EnvLogPrimaryT  = numpy.array(V20PrimaryT)
		EnvLogOutsideT  = numpy.array(V20OutsideT)
		EnvLogFanPower  = numpy.array(V20FanPower)
		EnvLogFocusPos  = numpy.array(V20FocusPos)
		EnvLogWindSpeed = numpy.array(V20WindSpeed)
		EnvLogSkyDiff   = numpy.array(V20SkyDiff)
		EnvLogHumidity  = numpy.array(V20Humidity)
		EnvLogWetness   = numpy.array(V20Wetness)
		EnvLogWetnessFilt  = numpy.array(V20WetnessFilt)
		EnvLogWetnessTime  = numpy.array(V20WetnessTime)
		EnvLogHumidityFilt = numpy.array(V20HumidityFilt)
		if FoundV5EnvLogFile:
			OtherTime      = numpy.array(V5TimeDecimal)
			OtherOutsideT  = numpy.array(V5OutsideT)
			OtherWindSpeed = numpy.array(V5WindSpeed)
			OtherSkyDiff   = numpy.array(V5SkyDiff)
			OtherHumidity  = numpy.array(V5Humidity)
	if telescope == "V5" and FoundV5EnvLogFile:
		if FoundV5EnvLogFile: FoundEnvLogFile = True
		if FoundV20EnvLogFile: FoundOtherLogFile = True
		OtherTelescope = "V20"
		EnvLogTime      = numpy.array(V5TimeDecimal)
		EnvLogTubeT     = numpy.array(V5TubeT)
		EnvLogOutsideT  = numpy.array(V5OutsideT)
		EnvLogFanPower  = numpy.array(V5FocusPos)
		EnvLogWindSpeed = numpy.array(V5WindSpeed)
		EnvLogSkyDiff   = numpy.array(V5SkyDiff)
		EnvLogFocusPos  = numpy.array(V5FocusPos)
		EnvLogHumidity  = numpy.array(V5Humidity)
		EnvLogWetness   = numpy.array(V5Wetness)
		EnvLogWetnessFilt  = numpy.array(V5WetnessFilt)
		EnvLogWetnessTime  = numpy.array(V5WetnessTime)
		EnvLogHumidityFilt = numpy.array(V5HumidityFilt)
		if FoundV20EnvLogFile:
			OtherTime      = numpy.array(V20TimeDecimal)
			OtherOutsideT  = numpy.array(V20OutsideT)
			OtherWindSpeed = numpy.array(V20WindSpeed)
			OtherSkyDiff   = numpy.array(V20SkyDiff)
			OtherHumidity  = numpy.array(V20Humidity)

	if FoundEnvLogFile:
		EnvLogTimePlot  = EnvLogTime[(EnvLogTime>PlotStartUT) & (EnvLogTime<PlotEndUT)]
		## Filter Out SkyTemp of -998 and Wind of -1.2
		EnvLogTimeFilt      = EnvLogTime[(EnvLogWindSpeed > 0) & (EnvLogSkyDiff > -140)]
		EnvLogSkyDiffFilt   = EnvLogSkyDiff[(EnvLogWindSpeed > 0) & (EnvLogSkyDiff > -140)]
		EnvLogWindSpeedFilt = EnvLogWindSpeed[(EnvLogWindSpeed > 0) & (EnvLogSkyDiff > -140)]

		if FoundOtherLogFile:
			## Filter Out SkyTemp of -998 and Wind of -1.2
			OtherTimeFilt      = OtherTime[(OtherWindSpeed > 0) & (OtherSkyDiff > -140)]	
			OtherWindSpeedFilt = OtherWindSpeed[(OtherWindSpeed > 0) & (OtherSkyDiff > -140)]
			OtherSkyDiffFilt   = OtherSkyDiff[(OtherWindSpeed > 0) & (OtherSkyDiff > -140)]						
			
	###########################################################
	## Match up ACP Log and IQMon Results Based on filename
	if len(ACPdata) >= 5:					
		if FoundIQMonFile and FoundACPLog:
			data = []
			for ACPentry in ACPdata:
				FoundIQMonMatch = False
				IsEmptyFilter = re.match(".*Empty.*", ACPentry[1])
				if IsEmptyFilter:
					FoundIQMonMatch = True
				else:
					for IQMonEntry in IQMonData:
						if IQMonEntry[0] == ACPentry[1]+".fts":
							FoundIQMonMatch = True
							if not FoundCCDIFile:
								CCDIinfo = [0,0,0]
							else:
								for CCDIentry in CCDIdata:
									if CCDIentry[0] == ACPentry[1]+".fts":
										CCDIinfo = CCDIentry
							data.append([ACPentry[0], ACPentry[1], ACPentry[2], ACPentry[3], ACPentry[4], ACPentry[5],
							             CCDIinfo[0], CCDIinfo[2], 
							             IQMonEntry[0], IQMonEntry[1], IQMonEntry[2], IQMonEntry[3], 
							             IQMonEntry[4], IQMonEntry[5], IQMonEntry[6], IQMonEntry[7], IQMonEntry[8], IQMonEntry[9]
							             ])	
				if not FoundIQMonMatch:
					print "  - Could not find IQMon results for "+ACPentry[1]+".fts."
	


		###########################################################			
		## Gather ACP and IQMon FWHM Data
		if FoundACPLog and FoundIQMonFile:
			Time = []
			Time_PErr = []
			Time_ZP = []
			Time_ACPPErr = []
			FWHM_ACPimage = []
			FWHM_ACPavg = []
			FWHM_CCDI = []
			FWHM_IQMon = []
			Ellipticity = []
			ZP1_IQMon = []
			ZP2_IQMon = []
			PErr_IQMon = []
			PErr_ACP = []
			Alt = []
			Az = []
			Airmass = []
			for entry in sorted(data):
				Time.append(entry[0])
				FWHM_ACPimage.append(entry[3])
				FWHM_ACPavg.append(entry[4])
				if not math.isnan(entry[5]):
					PErr_ACP.append(entry[5])
					Time_ACPPErr.append(entry[0])
				Ellipticity.append(entry[14])
				if FoundCCDIFile:
					FWHM_CCDI.append(entry[7])				
				FWHM_IQMon.append(entry[10])
				if not math.isnan(entry[11]):
					PErr_IQMon.append(entry[11])
					Time_PErr.append(entry[0])
				if not math.isnan(entry[12]) or math.isnan(entry[13]):
					ZP1_IQMon.append(entry[12])
					ZP2_IQMon.append(entry[13])
					Time_ZP.append(entry[0])
					Alt.append(entry[15])
					Az.append(entry[16])
					Airmass.append(entry[17])
			FWHM_IQMon    = numpy.array(FWHM_IQMon)
			FWHM_ACPimage = numpy.array(FWHM_ACPimage)
			FWHM_ACPavg   = numpy.array(FWHM_ACPavg)
			FWHM_CCDI     = numpy.array(FWHM_CCDI)	
			Ellipticity   = numpy.array(Ellipticity)

	###########################################################			
	## Make Nightly Sumamry Plot (show only night time)
	###########################################################			
	print "  Writing Output File: "+PlotFileName
	dpi=100
	Figure = pyplot.figure(figsize=(15,12), dpi=dpi)
	
	###########################################################			
	## Temperatures
	if FoundEnvLogFile:
		Figure.add_axes([0.0, 0.765, 0.47, 0.235])
		pyplot.title("Environmental Data for "+telescope + " on the Night of " + DateString)
		pyplot.plot(EnvLogTime, EnvLogTubeT, 'g-', label="Tube Temp")
		pyplot.plot(EnvLogTime, EnvLogOutsideT, 'k-', label="Outside Temp ("+telescope+")")
		if FoundOtherLogFile:
			pyplot.plot(OtherTime, OtherOutsideT, 'k-', alpha=0.7, label="Outside Temp ("+OtherTelescope+")")
		if telescope == "V20":
			pyplot.plot(EnvLogTime, EnvLogPrimaryT, 'r-', label="Mirror Temp")
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Temperature (F)")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.grid()
		
		## Overplot Twilights
		pyplot.axvspan(SunsetDecimal, EveningCivilTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.1)
		pyplot.axvspan(EveningCivilTwilightDecimal, EveningNauticalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.2)
		pyplot.axvspan(EveningNauticalTwilightDecimal, EveningAstronomicalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.3)
		pyplot.axvspan(EveningAstronomicalTwilightDecimal, MorningAstronomicalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.5)
		pyplot.axvspan(MorningAstronomicalTwilightDecimal, MorningNauticalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.3)
		pyplot.axvspan(MorningNauticalTwilightDecimal, MorningCivilTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.2)
		pyplot.axvspan(MorningCivilTwilightDecimal, SunriseDecimal, ymin=0, ymax=1, color='blue', alpha=0.1)
		
		## Add Fan Power (if VYSOS-20)
		if telescope == "V20":
			Figure.add_axes([0.0, 0.695, 0.47, 0.05], xticklabels=[])
			pyplot.plot(EnvLogTime, EnvLogFanPower, 'b-', label="Fan Power (%)")
			# pyplot.xlabel("Time in Hours UT")
			pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
			pyplot.xlim(PlotStartUT,PlotEndUT)
			pyplot.ylim(-10,110)
			pyplot.yticks(numpy.linspace(0,100,3,endpoint=True))
			# pyplot.ylabel("Power")
			pyplot.legend(loc='best', prop={'size':10})
			pyplot.grid()
	
	###########################################################			
	## Humidity
	if FoundEnvLogFile:
		if telescope == "V5":
			Figure.add_axes([0.0, 0.510, 0.47, 0.235], xticklabels=[])
		if telescope == "V20":
			Figure.add_axes([0.0, 0.510, 0.47, 0.175])
		if FoundOtherLogFile:
			pyplot.plot(OtherTime, OtherHumidity, 'k-', alpha=0.6, label="Humidity ("+OtherTelescope+")")
		pyplot.plot(EnvLogTime, EnvLogHumidity, 'b-', label="Humidity ("+telescope+")")
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Humidity (%)")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylim(-5,105)
		pyplot.grid()
				
		pyplot.fill_between(EnvLogTime, 0, EnvLogHumidity, where=EnvLogWetness==2, color='blue', alpha=0.3)
		pyplot.fill_between(EnvLogTime, 0, EnvLogHumidity, where=EnvLogWetness==3, color='blue', alpha=0.6)
	
		# pyplot.fill_between(EnvLogWetnessTime, 0, EnvLogHumidityFilt, where=EnvLogWetnessFilt==2, color='blue', alpha=0.3)
		# pyplot.fill_between(EnvLogWetnessTime, 0, EnvLogHumidityFilt, where=EnvLogWetnessFilt==3, color='blue', alpha=0.6)
	

	###########################################################			
	## Sky Condition (Cloudiness)
	if FoundEnvLogFile:
		Figure.add_axes([0.0, 0.255, 0.47, 0.235])
		if FoundOtherLogFile:
			pyplot.plot(OtherTimeFilt, OtherSkyDiffFilt, 'k-', alpha=0.6, label="Sky Condition ("+OtherTelescope+")")
		pyplot.plot(EnvLogTimeFilt, EnvLogSkyDiffFilt, 'b-', label="Sky Condition ("+telescope+")")
		if telescope == "V20":
			VeryCloudyThreshold = -53.6
			CloudinessThreshold = -70.0
		if telescope == "V5":
			VeryCloudyThreshold = -47.3
			CloudinessThreshold = -65.2
		pyplot.fill_between(EnvLogTimeFilt, -140, EnvLogSkyDiffFilt, where=EnvLogSkyDiffFilt<CloudinessThreshold, color='green', alpha=0.5)
		pyplot.fill_between(EnvLogTimeFilt, -140, EnvLogSkyDiffFilt, where=(EnvLogSkyDiffFilt<VeryCloudyThreshold)&(EnvLogSkyDiffFilt>CloudinessThreshold), color='yellow', alpha=0.7)
		pyplot.fill_between(EnvLogTimeFilt, -140, EnvLogSkyDiffFilt, where=EnvLogSkyDiffFilt>=VeryCloudyThreshold, color='red', alpha=0.5)

		pyplot.plot([PlotStartUT,PlotEndUT], [VeryCloudyThreshold, VeryCloudyThreshold], 'r-')
		pyplot.plot([PlotStartUT,PlotEndUT], [CloudinessThreshold, CloudinessThreshold], 'y-')
		pyplot.text(PlotStartUT+0.25, VeryCloudyThreshold+5, 'Very Cloudy')
		pyplot.text(PlotStartUT+0.25, CloudinessThreshold+5, 'Cloudy')
		pyplot.text(PlotStartUT+0.25, CloudinessThreshold-10, 'Clear')

		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Temperature Difference (F)")
		# pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylim(-140,-20)
		pyplot.grid()

	###########################################################			
	## Wind Speed
	if FoundEnvLogFile:
		Figure.add_axes([0.0, 0.000, 0.47, 0.235])
		if FoundOtherLogFile:
			pyplot.plot(OtherTimeFilt, OtherWindSpeedFilt, 'k-', alpha=0.6, label="Wind Speed ("+OtherTelescope+")")
		pyplot.plot(EnvLogTimeFilt, EnvLogWindSpeedFilt, 'b-', label="Wind Speed ("+telescope+")")
		if telescope == "V20":
			WindinessThreshold = 22.0*0.621371
			VeryWindyThreshold = 40.0*0.621371
		if telescope == "V5":
			WindinessThreshold = 30.0*0.621371
			VeryWindyThreshold = 50.0*0.621371

		pyplot.fill_between(EnvLogTimeFilt, 0, EnvLogWindSpeedFilt, where=EnvLogWindSpeedFilt<WindinessThreshold, color='green', alpha=0.5)
		pyplot.fill_between(EnvLogTimeFilt, 0, EnvLogWindSpeedFilt, where=(EnvLogWindSpeedFilt>WindinessThreshold)&(EnvLogWindSpeedFilt<VeryWindyThreshold), color='yellow', alpha=0.7)
		pyplot.fill_between(EnvLogTimeFilt, 0, EnvLogWindSpeedFilt, where=EnvLogWindSpeedFilt>=VeryWindyThreshold, color='red', alpha=0.5)

		pyplot.plot([PlotStartUT,PlotEndUT], [VeryWindyThreshold, VeryWindyThreshold], 'r-')
		pyplot.plot([PlotStartUT,PlotEndUT], [WindinessThreshold, WindinessThreshold], 'y-')
		pyplot.text(PlotStartUT+0.25, VeryWindyThreshold+1, 'Very Windy')
		pyplot.text(PlotStartUT+0.25, WindinessThreshold+3, 'Windy')
		pyplot.text(PlotStartUT+0.25, WindinessThreshold-3, 'Calm')
		
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Wind Speed (mph)")
		pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylim(0,max([max(EnvLogWindSpeedFilt)*1.1,35.]))
		pyplot.grid()
	
	
	
	###########################################################			
	## FWHM vs. Time
	if FoundIQMonFile:
		Figure.add_axes([0.53, 0.765, 0.47, 0.235])			
		pyplot.title("IQ Mon Results for "+telescope + " on the Night of " + DateString)
		if FoundACPLog:
			if telescope == "V5":
				pyplot.plot(Time, FWHM_ACPimage, 'g.-', drawstyle="steps-post", label="FWHM (ACP Image)", alpha=0.5)
			if telescope == "V20":
				pyplot.plot(Time, FWHM_ACPimage*PixelScale, 'g.-', drawstyle="steps-post", label="FWHM (ACP Image)", alpha=0.5)
		pyplot.plot(Time, FWHM_IQMon, 'k.-', drawstyle="steps-post", label="FWHM (IQMon)")
		pyplot.ylabel("FWHM (arcsec)")
		pyplot.yticks(numpy.linspace(0,15,16,endpoint=True))
		pyplot.ylim(0,5)
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.grid()
		if FoundFocusMaxFile:
			for FocusRun in FocusRuns:
				pyplot.axvspan(FocusRun[2], FocusRun[3], color="gray", alpha=0.7)
		pyplot.legend(loc='best', prop={'size':10})
		# pyplot.xlabel("Time in Hours UT")

		## Overplot Twilights
		pyplot.axvspan(SunsetDecimal, EveningCivilTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.1)
		pyplot.axvspan(EveningCivilTwilightDecimal, EveningNauticalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.2)
		pyplot.axvspan(EveningNauticalTwilightDecimal, EveningAstronomicalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.3)
		pyplot.axvspan(EveningAstronomicalTwilightDecimal, MorningAstronomicalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.5)
		pyplot.axvspan(MorningAstronomicalTwilightDecimal, MorningNauticalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.3)
		pyplot.axvspan(MorningNauticalTwilightDecimal, MorningCivilTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.2)
		pyplot.axvspan(MorningCivilTwilightDecimal, SunriseDecimal, ymin=0, ymax=1, color='blue', alpha=0.1)


	###########################################################			
	## Focus Position
	if FoundEnvLogFile and FoundIQMonFile:
		Figure.add_axes([0.53, 0.510, 0.47, 0.235])
		pyplot.plot(EnvLogTime, EnvLogFocusPos, 'b.-', label="Focus Position")
		# pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		## Clip the extreme values off of the focus position list (bad values?)
		nFocusPos = len(EnvLogFocusPos)
		ClippingFactor = 0.02
		nClipped = int(nFocusPos*ClippingFactor)
		ylimlower = (sorted(EnvLogFocusPos))[nClipped]
		ylimupper = (sorted(EnvLogFocusPos))[nFocusPos-nClipped]
		ylimrange = max([100 , ylimupper - ylimlower])
		pyplot.ylim(ylimlower-0.15*ylimrange, ylimupper+0.15*ylimrange)
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.grid()
	

	###########################################################			
	## Ellipticity vs. Time
	if FoundIQMonFile:
		Figure.add_axes([0.53, 0.255, 0.47, 0.235])
		pyplot.plot(Time_ZP, Ellipticity, 'b.-', drawstyle="steps-post", label="Ellipticity")
		pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylabel("Ellipticity")
		pyplot.ylim(0,1)
		pyplot.grid()


	###########################################################			
	## Pointing Error vs. Time
	if FoundIQMonFile:
		Figure.add_axes([0.53, 0.000, 0.47, 0.235])
		if FoundACPLog:
			pyplot.plot(Time_ACPPErr, PErr_ACP, 'g.-', drawstyle="steps-post", label="ACP Log", alpha=0.6)
		pyplot.plot(Time_PErr, PErr_IQMon, 'b.-', drawstyle="steps-post", label="IQMon")
		# pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylabel("Pointing Error (arcmin)")
		if telescope == "V5":
			PErrPlotMax = 10
		if telescope == "V20":
			PErrPlotMax = 10
		pyplot.ylim(0,PErrPlotMax)
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.grid()
	
	pyplot.savefig(PlotFile, dpi=dpi, bbox_inches='tight', pad_inches=0.01)
	


	###########################################################			
	## Make Environmental Plot (show entire day)
	###########################################################			
	if FoundEnvLogFile:
		print "  Writing Output File: "+EnvPlotFileName
		dpi=100
		Figure = pyplot.figure(figsize=(12,12), dpi=dpi)

		###########################################################			
		## Temperatures
		if FoundEnvLogFile:
			Figure.add_axes([0.0, 0.765, 1.0, 0.235])
			pyplot.title("Environmental Data for "+telescope + " on the Night of " + DateString)
			pyplot.plot(EnvLogTime, EnvLogTubeT, 'g-', label="Tube Temp")
			pyplot.plot(EnvLogTime, EnvLogOutsideT, 'k-', label="Outside Temp ("+telescope+")")
			if FoundOtherLogFile:
				pyplot.plot(OtherTime, OtherOutsideT, 'k-', alpha=0.7, label="Outside Temp ("+OtherTelescope+")")
			if telescope == "V20":
				pyplot.plot(EnvLogTime, EnvLogPrimaryT, 'r-', label="Mirror Temp")
			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Temperature (F)")
			# if telescope == "V5":
			# 	pyplot.xlabel("Time in Hours UT")
			pyplot.xticks(numpy.linspace(0,24,25,endpoint=True))
			pyplot.xlim(0,24)
			pyplot.grid()

			## Overplot Twilights
			pyplot.axvspan(SunsetDecimal, EveningCivilTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.1)
			pyplot.axvspan(EveningCivilTwilightDecimal, EveningNauticalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.2)
			pyplot.axvspan(EveningNauticalTwilightDecimal, EveningAstronomicalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.3)
			pyplot.axvspan(EveningAstronomicalTwilightDecimal, MorningAstronomicalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.5)
			pyplot.axvspan(MorningAstronomicalTwilightDecimal, MorningNauticalTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.3)
			pyplot.axvspan(MorningNauticalTwilightDecimal, MorningCivilTwilightDecimal, ymin=0, ymax=1, color='blue', alpha=0.2)
			pyplot.axvspan(MorningCivilTwilightDecimal, SunriseDecimal, ymin=0, ymax=1, color='blue', alpha=0.1)

			## Add Fan Power (if VYSOS-20)
			if telescope == "V20":
				Figure.add_axes([0.0, 0.695, 1.0, 0.05], xticklabels=[])
				pyplot.plot(EnvLogTime, EnvLogFanPower, 'b-', label="Fan Power (%)")
				# pyplot.xlabel("Time in Hours UT")
				pyplot.xticks(numpy.linspace(0,24,25,endpoint=True))
				pyplot.xlim(0,24)
				pyplot.ylim(-10,110)
				pyplot.yticks(numpy.linspace(0,100,3,endpoint=True))
				# pyplot.ylabel("Power")
				pyplot.legend(loc='best', prop={'size':10})
				pyplot.grid()

		###########################################################			
		## Humidity
		if FoundEnvLogFile:
			if telescope == "V5":
				Figure.add_axes([0.0, 0.510, 1.0, 0.235], xticklabels=[])
			if telescope == "V20":
				Figure.add_axes([0.0, 0.510, 1.0, 0.175])
			if FoundOtherLogFile:
				pyplot.plot(OtherTime, OtherHumidity, 'k-', alpha=0.6, label="Humidity ("+OtherTelescope+")")
			pyplot.plot(EnvLogTime, EnvLogHumidity, 'b-', label="Humidity ("+telescope+")")
			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Humidity (%)")
			pyplot.xticks(numpy.linspace(0,24,25,endpoint=True))
			pyplot.xlim(0,24)
			pyplot.ylim(-5,105)
			pyplot.grid()
			
			pyplot.fill_between(EnvLogTime, 0, EnvLogHumidity, where=EnvLogWetness==2, color='blue', alpha=0.3)
			pyplot.fill_between(EnvLogTime, 0, EnvLogHumidity, where=EnvLogWetness==3, color='blue', alpha=0.6)

			# pyplot.fill_between(EnvLogWetnessTime, 0, EnvLogHumidityFilt, where=EnvLogWetnessFilt==2, color='blue', alpha=0.3)
			# pyplot.fill_between(EnvLogWetnessTime, 0, EnvLogHumidityFilt, where=EnvLogWetnessFilt==3, color='blue', alpha=0.6)

		###########################################################			
		## Sky Condition (Cloudiness)
		if FoundEnvLogFile:
			Figure.add_axes([0.0, 0.255, 1.0, 0.235])
			if FoundOtherLogFile:
				pyplot.plot(OtherTimeFilt, OtherSkyDiffFilt, 'k-', alpha=0.6, label="Sky Condition ("+OtherTelescope+")")
			pyplot.plot(EnvLogTimeFilt, EnvLogSkyDiffFilt, 'b-', label="Sky Condition ("+telescope+")")
			if telescope == "V20":
				VeryCloudyThreshold = -53.6
				CloudinessThreshold = -70.0
			if telescope == "V5":
				VeryCloudyThreshold = -47.3
				CloudinessThreshold = -65.2
			pyplot.fill_between(EnvLogTimeFilt, -140, EnvLogSkyDiffFilt, where=EnvLogSkyDiffFilt<CloudinessThreshold, color='green', alpha=0.5)
			pyplot.fill_between(EnvLogTimeFilt, -140, EnvLogSkyDiffFilt, where=(EnvLogSkyDiffFilt<VeryCloudyThreshold)&(EnvLogSkyDiffFilt>CloudinessThreshold), color='yellow', alpha=0.7)
			pyplot.fill_between(EnvLogTimeFilt, -140, EnvLogSkyDiffFilt, where=EnvLogSkyDiffFilt>=VeryCloudyThreshold, color='red', alpha=0.5)

			pyplot.plot([0,24], [VeryCloudyThreshold, VeryCloudyThreshold], 'r-')
			pyplot.plot([0,24], [CloudinessThreshold, CloudinessThreshold], 'y-')
			pyplot.text(0+0.25, VeryCloudyThreshold+5, 'Very Cloudy')
			pyplot.text(0+0.25, CloudinessThreshold+5, 'Cloudy')
			pyplot.text(0+0.25, CloudinessThreshold-10, 'Clear')

			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Temperature Difference (F)")
			# pyplot.xlabel("Time in Hours UT")
			pyplot.xticks(numpy.linspace(0,24,25,endpoint=True))
			pyplot.xlim(0,24)
			pyplot.ylim(-140,-20)
			pyplot.grid()

		###########################################################			
		## Wind Speed
		if FoundEnvLogFile:
			Figure.add_axes([0.0, 0.000, 1.0, 0.235])
			if FoundOtherLogFile:
				pyplot.plot(OtherTimeFilt, OtherWindSpeedFilt, 'k-', alpha=0.6, label="Wind Speed ("+OtherTelescope+")")
			pyplot.plot(EnvLogTimeFilt, EnvLogWindSpeedFilt, 'b-', label="Wind Speed ("+telescope+")")
			if telescope == "V20":
				WindinessThreshold = 22.0*0.621371
				VeryWindyThreshold = 40.0*0.621371
			if telescope == "V5":
				WindinessThreshold = 30.0*0.621371
				VeryWindyThreshold = 50.0*0.621371

			pyplot.fill_between(EnvLogTimeFilt, 0, EnvLogWindSpeedFilt, where=EnvLogWindSpeedFilt<WindinessThreshold, color='green', alpha=0.5)
			pyplot.fill_between(EnvLogTimeFilt, 0, EnvLogWindSpeedFilt, where=(EnvLogWindSpeedFilt>WindinessThreshold)&(EnvLogWindSpeedFilt<VeryWindyThreshold), color='yellow', alpha=0.7)
			pyplot.fill_between(EnvLogTimeFilt, 0, EnvLogWindSpeedFilt, where=EnvLogWindSpeedFilt>=VeryWindyThreshold, color='red', alpha=0.5)

			pyplot.plot([0,24], [VeryWindyThreshold, VeryWindyThreshold], 'r-')
			pyplot.plot([0,24], [WindinessThreshold, WindinessThreshold], 'y-')
			pyplot.text(0+0.25, VeryWindyThreshold+1, 'Very Windy')
			pyplot.text(0+0.25, WindinessThreshold+3, 'Windy')
			pyplot.text(0+0.25, WindinessThreshold-3, 'Calm')

			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Wind Speed (mph)")
			pyplot.xlabel("Time in Hours UT")
			pyplot.xticks(numpy.linspace(0,24,25,endpoint=True))
			pyplot.xlim(0,24)
			pyplot.ylim(0,max([max(EnvLogWindSpeedFilt)*1.1,35.]))
			pyplot.grid()
	
		pyplot.savefig(EnvPlotFile, dpi=dpi, bbox_inches='tight', pad_inches=0.01)
	
	if os.path.exists(PlotFile):
		return True
	else:
		return False




	
def main(argv=None):
	DateString = ""	
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hd:t:", ["help", "date=", "telescope="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-d", "--date"):
				DateString = value
			if option in ("-t", "--telescope"):
				telescope = value
			
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2
	
	## Set date to tonight if not specified
	now = datetime.datetime.utcnow()
	
	if (DateString == "tonight") or (DateString == "Tonight") or (DateString == "lastnight") or (DateString == "LastNight") or (DateString ==""):
		DateString = now.strftime("%Y%m%dUT")
	
	## Input Should be Folder Name of ACP Logs in format yyyymmddUT
	## which indicates the night on which the data was taken.
	MatchACPLogFolderName = re.compile("[0-9]{8}UT")
	IsACPLogFolderName = MatchACPLogFolderName.match(DateString)
	if telescope == "V5" or telescope == "VYSOS5" or telescope == "VYSOS-5":
		telescope = "V5"
	if telescope == "V20" or telescope == "VYSOS20" or telescope == "VYSOS-20":
		telescope = "V20"
	Success = CompareIQ(DateString, telescope)
	
	# if Success:
		


if __name__ == '__main__':
	main()

