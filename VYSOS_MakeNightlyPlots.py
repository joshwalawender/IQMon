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
import astropy.io
import astropy.io.ascii
import astropy.table

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
	
	
###########################################################
## Read ACP Logs
## - extract ACP FWHM and Pointing Error
##
## VYSOSDATAPath
## DateString
## FoundACPLog (should be output)
## ACPdata (output)
def ReadACPLog(DateString, VYSOSDATAPath, PixelScale):
	MatchStartOfImage  = re.compile("(\d{2}:\d{2}:\d{2})\s*Imaging\sto\s([a-zA-Z0-9@\-_\+]*)")
	MatchPointingError = re.compile("(\d{2}:\d{2}:\d{2})\s*Pointing\serror\sis\s([0-9\.]+)\sarcmin.*")
	MatchImageFWHM     = re.compile("(\d{2}:\d{2}:\d{2})\s*Image\sFWHM\sis\s[0-9\.]{2,5}\sarcsec\s\(([0-9\.]{2,5})\spixels\)")
	MatchAvgFWHM       = re.compile("(\d{2}:\d{2}:\d{2})\s*\(avg\sFWHM\s=\s([0-9\.]{2,5})\sarcsec\)")
	MatchRunComplete   = re.compile("(\d{2}:\d{2}:\d{2})\s*Run\scomplete")
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
						if (ImageFile != "") and (not re.match(".*Empty.*", ImageFile)):
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
		nACPImages = len(ACPdata)
		print "    Data for %d images (Empty filter images excluded) extracted from ACP Logs." % nACPImages
	else:
		print "  Failed to Find ACP Log Directory: "+ACPLogDirectory
		ACPdata = []
	
	return ACPdata



###########################################################
## Read IQMon Logs
## - extract IQMon FWHM, ellipticity, pointing error
def ReadIQMonLog(LogPath, telescope, DateString):
	FoundIQMonFile = False
	# ColStarts = [ 0, 11, 21, 71, 83,  95, 107, 119, 131, 146, 158, 170]
	# ColEnds   = [ 9, 20, 70, 82, 94, 106, 118, 130, 145, 157, 169, 181]
	# ColNames  = ['Date', 'TimeString', 'File', 'FWHM', 'Ellipticity', 'Alt', 'Az', 'Airmass', 'PointingError', 'ZeroPoint1', 'ZeroPoint2', 'nStars']
	if telescope == "V5":
		telname = "VYSOS-5"
	if telescope == "V20":
		telname = "VYSOS-20"
	if os.path.exists(os.path.join(LogPath, telname)):
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
					try:
						IQMonTable = astropy.io.ascii.read(FullIQMonFile)
						# IQMonTable = astropy.io.ascii.read(FullIQMonFile, data_start=2, Reader=astropy.io.ascii.FixedWidth, 
						#              col_starts=ColStarts, col_ends=ColEnds, names=ColNames, 
						#              guess=False, comment=";", header_start=0)
						IQMonTimeDecimals = []
						for i in range(0,len(IQMonTable),1):
							IQMonTimeDecimals.append(TimeStringToDecimal(IQMonTable[i]['ExpStart'][11:19]))
						IQMonTable.add_column(astropy.table.Column(data=IQMonTimeDecimals, name='Time'))
						print "    Data for %d images extracted from IQMon summary file." % len(IQMonTable)
					except:
						print "    Failed to Read IQMon Log File"
						IQMonTable = astropy.table.Table(names=ColNames)
	if not FoundIQMonFile:
		print "  Failed to Find IQMon Logs: "+os.path.join(LogPath, telname)
		IQMonTable = astropy.table.Table()

	return IQMonTable
	
	
	
###########################################################
## Read FocusMax Logs
## - extract times of focus operations
def ReadFocusMaxLog(VYSOSDATAPath, DateString):
	FocusRuns = []
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
		FocusRuns = []

	return FocusRuns
	
	
###########################################################			
## Read Environmental Logs
def ReadEnvironmentalLogs(DateString, telescope, V5DataPath, V20DataPath):
	print "  Reading Environmental Logs"
	V20EnvLogFileName = os.path.join(V20DataPath, "Logs", DateString, "EnvironmentalLog.txt")
	V5EnvLogFileName  = os.path.join(V5DataPath,  "Logs", DateString, "EnvironmentalLog.txt")
	FoundV20EnvLogFile = False
	FoundV5EnvLogFile  = False
	FoundOtherLogFile = False
	if os.path.exists(V20EnvLogFileName):
		print "  Found VYSOS-20 Environmental Logs"
		FoundV20EnvLogFile = True
		ColStarts = [ 0, 11, 22, 32, 42, 52, 62, 72, 82,  92, 102, 112, 122, 132, 142, 152, 162]
		ColEnds   = [ 9, 18, 31, 41, 51, 61, 71, 81, 91, 101, 111, 121, 131, 141, 151, 161, 171]
		ColNames  = ['Date', 'TimeString', 'TubeTemp', 'PrimaryTemp', 'SecTemp', 'FanPower', 'FocusPos', 
		             'SkyTemp', 'OutsideTemp', 'WindSpeed', 'Humidity', 'DewPoint', 'Alt', 'Az', 'Condition', 'DomeTemp', 'DomeFanState']
		V20EnvTable = astropy.io.ascii.read(V20EnvLogFileName, data_start=2, Reader=astropy.io.ascii.FixedWidth,
		              col_starts=ColStarts, col_ends=ColEnds, names=ColNames, 
		              guess=False, comment=";", header_start=0)
		V20SkyDiff   = V20EnvTable['SkyTemp'] #- V20EnvTable['OutsideTemp']
		V20EnvTable.add_column(astropy.table.Column(data=V20SkyDiff, name='SkyDiff'))
		V20DomeFan = []
		V20TimeDecimal = []
		V20Wetness = []
		V20Cloudiness = []
		V20Windiness = []
		for i in range(0,len(V20EnvTable),1):
			## Make Time Decimal
			V20TimeDecimal.append(TimeStringToDecimal(V20EnvTable[i]['TimeString'][0:8]))
			## Parse Condition String
			ConditionMatch = re.match("([0-3])([0-3])([0-3])", str(V20EnvTable[i]['Condition']))
			if ConditionMatch:
				V20Wetness.append(ConditionMatch.group(1))
				V20Cloudiness.append(ConditionMatch.group(2))
				V20Windiness.append(ConditionMatch.group(3))
			else:
				# print "no match to |"+str(V20EnvTable[i]['Condition'])+"| in V20 Env Log"
				V20Wetness.append("-1")
				V20Cloudiness.append("-1")
				V20Windiness.append("-1")
			## Filter Out Bad Sky Diff Values
			if V20EnvTable[i]['SkyTemp'] < -100.: V20EnvTable[i]['SkyTemp'] = float("nan")
			## Parse Dome Fan State
			FanMatch = re.match(".*([0-1])([0-1]).*", str(V20EnvTable[i]['DomeFanState']))
			if FanMatch:
				V20DomeFan.append(float(FanMatch.group(1))*100)
			else:
				V20DomeFan.append(float(0))
		V20EnvTable.add_column(astropy.table.Column(data=V20TimeDecimal, name='Time'))
		V20EnvTable.add_column(astropy.table.Column(data=V20Wetness, name='WetCondition'))
		V20EnvTable.add_column(astropy.table.Column(data=V20Cloudiness, name='CloudCondition'))
		V20EnvTable.add_column(astropy.table.Column(data=V20Windiness, name='WindCondition'))
		V20EnvTable.add_column(astropy.table.Column(data=V20DomeFan, name='DomeFan'))
	else:
		ColNames  = ['Date', 'TimeString', 'TubeTemp', 'PrimaryTemp', 'SecTemp', 'FanPower', 'FocusPos', 
		             'SkyTemp', 'OutsideTemp', 'WindSpeed', 'Humidity', 'DewPoint', 'Alt', 'Az', 'Condition', 'DomeTemp', 'DomeFanState']
		V20EnvTable = astropy.table.Table(names=ColNames)
		
	if os.path.exists(V5EnvLogFileName):
		print "  Found VYSOS-5 Environmental Logs"
		FoundV5EnvLogFile = True
		ColStarts = [ 0, 11, 22, 32, 42, 52, 62, 72, 82,  92, 102, 112]
		ColEnds   = [ 9, 18, 31, 41, 51, 61, 71, 81, 91, 101, 111, 121]
		ColNames  = ['Date', 'TimeString', 'TubeTemp', 'FocusPos', 
		             'SkyTemp', 'OutsideTemp', 'WindSpeed', 'Humidity', 'DewPoint', 'Alt', 'Az', 'Condition']
		V5EnvTable = astropy.io.ascii.read(V5EnvLogFileName, data_start=2, Reader=astropy.io.ascii.FixedWidth, 
		             col_starts=ColStarts, col_ends=ColEnds, names=ColNames, 
		             guess=False, comment=";", header_start=0)
		V5SkyDiff   = V5EnvTable['SkyTemp'] - V5EnvTable['OutsideTemp']
		V5EnvTable.add_column(astropy.table.Column(data=V5SkyDiff, name='SkyDiff'))
		V5TimeDecimal = []
		V5Wetness = []
		V5Cloudiness = []
		V5Windiness = []
		for i in range(0,len(V5EnvTable),1):
			## Make Time Decimal
			V5TimeDecimal.append(TimeStringToDecimal(V5EnvTable[i]['TimeString'][0:8]))
			## Parse Condition String
			ConditionMatch = re.match("\s*([\-0-9])([\-0-9])([\-0-9])\s*", str(V5EnvTable[i]['Condition']))
			if ConditionMatch:
				V5Wetness.append(ConditionMatch.group(1))
				V5Cloudiness.append(ConditionMatch.group(2))
				V5Windiness.append(ConditionMatch.group(3))
			else:
				# print "no match to |"+str(V5EnvTable[i]['Condition'])+"| in V5 Env Log"
				V5Wetness.append("-1")
				V5Cloudiness.append("-1")
				V5Windiness.append("-1")
			## Filter Out Bad Sky Diff Values
			if V5EnvTable[i]['SkyTemp'] < -100.: V5EnvTable[i]['SkyTemp'] = float("nan")
		V5EnvTable.add_column(astropy.table.Column(data=V5TimeDecimal, name='Time'))
		V5EnvTable.add_column(astropy.table.Column(data=V5Wetness, name='WetCondition'))
		V5EnvTable.add_column(astropy.table.Column(data=V5Cloudiness, name='CloudCondition'))
		V5EnvTable.add_column(astropy.table.Column(data=V5Windiness, name='WindCondition'))
	else:
		ColNames  = ['Date', 'TimeString', 'TubeTemp', 'FocusPos', 
		             'SkyTemp', 'OutsideTemp', 'WindSpeed', 'Humidity', 'DewPoint', 'Alt', 'Az', 'Condition']
		V5EnvTable = astropy.table.Table(names=ColNames)
	
		
	return V20EnvTable, V5EnvTable
		
							




###########################################################			
## Make Plots
def MakePlots(DateString, telescope):
	print "#### Making Nightly Plots for "+telescope+" on the Night of "+DateString+" ####"
	
	FoundACPLog       = False
	FoundIQMonFile    = False
	FoundFocusMaxFile = False
	FoundV20Env       = False
	FoundV5Env        = False
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
		OtherTelescope = "V20"
	if telescope == "V20":
		VYSOSDATAPath = V20DataPath
		PixelScale = 0.44
		telname = "VYSOS-20"
		OtherTelescope = "V5"
	
	## Set File Name
	PlotFileName = DateString+"_"+telescope+".png"
	PlotFile = os.path.join(LogPath, telname, PlotFileName)
	EnvPlotFileName = DateString+"_"+telescope+"_Env.png"
	EnvPlotFile = os.path.join(LogPath, telname, EnvPlotFileName)
	RecentPlotFileName = "Recent_"+telname+"_Conditions.png"
	RecentPlotFile = os.path.join(LogPath, telname, RecentPlotFileName)
	
	## Compile Various Regular Expressions for File Name Matching and ACP Log Parsing
	MatchDir           = re.compile(DateString)
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
	nUTHours = PlotEndUT-PlotStartUT+1
	
	Observatory.date = DateString[0:4]+"/"+DateString[4:6]+"/"+DateString[6:8]+" 0:00:01.0"
	TheMoon = ephem.Moon()
	TheMoon.compute(Observatory)
	MoonsetTime  = Observatory.next_setting(ephem.Moon()).datetime()
	MoonriseTime = Observatory.next_rising(ephem.Moon()).datetime()
	MoonsetDecimal = float(datetime.datetime.strftime(MoonsetTime, "%H"))+float(datetime.datetime.strftime(MoonsetTime, "%M"))/60.+float(datetime.datetime.strftime(MoonsetTime, "%S"))/3600.
	MoonriseDecimal = float(datetime.datetime.strftime(MoonriseTime, "%H"))+float(datetime.datetime.strftime(MoonriseTime, "%M"))/60.+float(datetime.datetime.strftime(MoonriseTime, "%S"))/3600.		
	
	MoonTimes = numpy.arange(0,24,0.1)
	MoonAlts = []
	for MoonTime in MoonTimes:
		TimeString = "%02d:%02d:%04.1f" % (math.floor(MoonTime), math.floor((MoonTime % 1)*60), ((MoonTime % 1 * 60) % 1)*60.0)
		Observatory.date = DateString[0:4]+"/"+DateString[4:6]+"/"+DateString[6:8]+" "+TimeString
		TheMoon.compute(Observatory)
		MoonAlts.append(TheMoon.alt * 180. / ephem.pi)
	MoonAlts = numpy.array(MoonAlts)
	
	MoonPeakAlt = max(MoonAlts)
	MoonPeakTime = (MoonTimes[(MoonAlts == MoonPeakAlt)])[0]
	MoonPeakTimeString = "%02d:%02d:%04.1f" % (math.floor(MoonPeakTime), math.floor((MoonPeakTime % 1)*60), ((MoonPeakTime % 1 * 60) % 1)*60.0)
	Observatory.date = DateString[0:4]+"/"+DateString[4:6]+"/"+DateString[6:8]+" "+MoonPeakTimeString
	TheMoon.compute(Observatory)
	MoonPhase = TheMoon.phase
		
	now = datetime.datetime.utcnow()
	DecimalTime = now.hour+now.minute/60.+now.second/3600.
		
	###########################################################
	## Read ACP Logs
	ACPdata = ReadACPLog(DateString, VYSOSDATAPath, PixelScale)
	if len(ACPdata) > 1: FoundACPLog = True

	###########################################################
	## Read IQMon Logs
	## - extract IQMon FWHM, ellipticity, pointing error
	IQMonTable = ReadIQMonLog(LogPath, telescope, DateString)
	if len(IQMonTable) > 1: FoundIQMonFile = True
	# IQMonTable = IQMonTable.sort('Time')

	###########################################################
	## Read FocusMax Logs
	## - extract times of focus operations
	FocusRuns = ReadFocusMaxLog(VYSOSDATAPath, DateString)
	if len(FocusRuns) > 1: FoundFocusMaxFile = True

	###########################################################			
	## Get Environmental Data
	V20EnvTable, V5EnvTable = ReadEnvironmentalLogs(DateString, telescope, V5DataPath, V20DataPath)
	if len(V20EnvTable) > 1: FoundV20Env = True
	if len(V5EnvTable) > 1:  FoundV5Env = True
	# V20EnvTable = V20EnvTable.sort('Time')
	# V5EnvTable = V5EnvTable.sort('Time')
		
	###########################################################
	## Match up ACP Log and IQMon Results Based on filename
	print "  Matching IQMon and ACP data"
	if FoundIQMonFile and FoundACPLog:
		MatchedData = astropy.table.Table(
		              names=('ACP Time', 'ACP File', 'ACP FWHM', 'ACP PErr', 'IQMon Time', 'IQMon File', 'IQMon FWHM', 'IQMon PErr'),
		              dtypes=('f',       'S',        'f',        'f',        'S',           'S',         'f',          'f')
		              )
		for ACPentry in ACPdata:
			FoundIQMonMatch = False
			for IQMonEntry in IQMonTable:
				if re.match(IQMonEntry['File']+".*", ACPentry[1]+".fts"):
					FoundIQMonMatch = True							
					MatchedData.add_row([ACPentry[0], ACPentry[1], ACPentry[3], ACPentry[5], 
					                     IQMonEntry['ExpStart'], IQMonEntry['File'], IQMonEntry['FWHM'], IQMonEntry['PointingError']])	
			if not FoundIQMonMatch:
				ACPFile = os.path.join(VYSOSDATAPath, "Images", DateString, ACPentry[1]+".fts")
				if os.path.exists(ACPFile):
					print "  - Could not find IQMon results for "+ACPentry[1]+".fts."
				else:
					print "  - Could not find file for ACP log entry "+ACPentry[1]+".fts."
	# MatchedData = MatchedData.sort('ACP Time')
		

	###########################################################			
	## Make Nightly Sumamry Plot (show only night time)
	###########################################################			
	print "  Writing Output File: "+PlotFileName
	dpi=100
	Figure = pyplot.figure(figsize=(11,11), dpi=dpi)
	
	###########################################################			
	## Temperatures
	if FoundV20Env or FoundV5Env:
		TemperatureAxes = pyplot.axes([0.0, 0.765, 0.46, 0.235])
		pyplot.title("Environmental Data for "+telescope + " on the Night of " + DateString)
	
		if telescope == "V20" and FoundV20Env:
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['OutsideTemp'], 'k-', alpha=0.5, drawstyle="steps-post", label="Outside Temp ("+OtherTelescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['TubeTemp'], 'g-', drawstyle="steps-post", label="Tube Temp")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['OutsideTemp'], 'k-', drawstyle="steps-post", label="Outside Temp ("+telescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['PrimaryTemp'], 'r-', drawstyle="steps-post", label="Mirror Temp")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['DomeTemp'], 'c-', drawstyle="steps-post", label="Dome Temp")
			
		
		if telescope == "V5" and FoundV5Env:
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['OutsideTemp'], 'k-', alpha=0.5, drawstyle="steps-post", label="Outside Temp ("+OtherTelescope+")")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['TubeTemp'], 'g-', drawstyle="steps-post", label="Tube Temp")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['OutsideTemp'], 'k-', drawstyle="steps-post", label="Outside Temp ("+telescope+")")
		
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

		## Overplot Moon Up Time
		MoonAxes = TemperatureAxes.twinx()
		MoonAxes.set_ylabel('Moon Alt (%.0f%% full)' % MoonPhase, color='y')
		pyplot.plot(MoonTimes, MoonAlts, 'y-')
		pyplot.ylim(0,100)
		pyplot.yticks([10,30,50,70,90], color='y')
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		MoonFill = MoonPhase/100.*0.5+0.05
		pyplot.fill_between(MoonTimes, 0, MoonAlts, where=MoonAlts>0, color='yellow', alpha=MoonFill)

		## Add Fan Power (if VYSOS-20)
		if telescope == "V20":
			Figure.add_axes([0.0, 0.675, 0.46, 0.07], xticklabels=[])
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['DomeFan'], 'c-', drawstyle="steps-post", label="Dome Fan State")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['FanPower'], 'b-', drawstyle="steps-post", label="Mirror Fans (%)")
			pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
			pyplot.xlim(PlotStartUT,PlotEndUT)
			pyplot.ylim(-10,110)
			pyplot.yticks(numpy.linspace(0,100,3,endpoint=True))
			pyplot.legend(loc='best', prop={'size':10})
			pyplot.grid()
	
	###########################################################			
	## Sky Condition (Cloudiness)
	if FoundV20Env or FoundV5Env:
		# Figure.add_axes([0.0, 0.255, 0.46, 0.235])
		if telescope == "V20" and FoundV20Env:
			Figure.add_axes([0.0, 0.430, 0.46, 0.235])
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['SkyTemp'], 'k-', alpha=0.5, drawstyle="steps-post", label="Cloudiness ("+OtherTelescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['SkyTemp'], 'b-', drawstyle="steps-post", label="Cloudiness ("+telescope+")")
			pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="3"), color='red', alpha=0.8)
		if telescope == "V5" and FoundV5Env:
			Figure.add_axes([0.0, 0.510, 0.46, 0.235])
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['SkyTemp'], 'k-', alpha=0.5, drawstyle="steps-post", label="Cloudiness ("+OtherTelescope+")")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['SkyTemp'], 'b-', drawstyle="steps-post", label="Cloudiness ("+telescope+")")
			pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="3"), color='red', alpha=0.8)
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Temperature Difference (F)")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylim(-100,-20)
		pyplot.grid()

	###########################################################			
	## Humidity
	if FoundV20Env or FoundV5Env:
		if telescope == "V5" and FoundV5Env:
			# Figure.add_axes([0.0, 0.510, 0.46, 0.235])
			Figure.add_axes([0.0, 0.255, 0.46, 0.235])
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['Humidity'], 'b-', drawstyle="steps-post", label="Humidity ("+telescope+")")
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['Humidity'], 'k-', alpha=0.5, drawstyle="steps-post", label="Humidity ("+OtherTelescope+")")
			pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="2"), color='red', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="3"), color='red', alpha=0.8)			
		if telescope == "V20" and FoundV20Env:
			# Figure.add_axes([0.0, 0.510, 0.46, 0.155])
			Figure.add_axes([0.0, 0.255, 0.46, 0.155])
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['Humidity'], 'b-', drawstyle="steps-post", label="Humidity ("+telescope+")")
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['Humidity'], 'k-', drawstyle="steps-post", alpha=0.5, label="Humidity ("+OtherTelescope+")")
			pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="2"), color='red', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="3"), color='red', alpha=0.8)
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Humidity (%)")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylim(-5,105)
		pyplot.grid()

	###########################################################			
	## Wind Speed
	if FoundV20Env or FoundV5Env:
		Figure.add_axes([0.0, 0.000, 0.46, 0.235])
		if telescope == "V20" and FoundV20Env:
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['WindSpeed'], 'k-', alpha=0.5, drawstyle="steps-post", label="Wind Speed ("+OtherTelescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['WindSpeed'], 'b-', drawstyle="steps-post", label="Wind Speed ("+telescope+")")
			pyplot.ylim(0,max([max(V20EnvTable['WindSpeed'])*1.1,35.]))
			pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="3"), color='red', alpha=0.8)
		if telescope == "V5" and FoundV5Env:
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['WindSpeed'], 'k-', alpha=0.5, label="Wind Speed ("+OtherTelescope+")")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['WindSpeed'], 'b-', label="Wind Speed ("+telescope+")")
			pyplot.ylim(0,max([max(V5EnvTable['WindSpeed'])*1.1,35.]))
			pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="3"), color='red', alpha=0.8)
					
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Wind Speed (mph)")
		pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.grid()
	
	
	
	###########################################################			
	## FWHM vs. Time
	if FoundIQMonFile:
		Figure.add_axes([0.54, 0.765, 0.46, 0.235])			
		pyplot.title("IQ Mon Results for "+telescope + " on the Night of " + DateString)
		if FoundACPLog:
			if telescope == "V5":
				pyplot.plot(MatchedData['ACP Time'], MatchedData['ACP FWHM'], 'g.', drawstyle="steps-post", label="FWHM (ACP Image)", alpha=0.5)
				pyplot.ylabel("FWHM (pixels)")
			if telescope == "V20":
				pyplot.plot(MatchedData['ACP Time'], MatchedData['ACP FWHM']*PixelScale, 'g.', drawstyle="steps-post", label="FWHM (ACP Image)", alpha=0.5)
				pyplot.ylabel("FWHM (arcsec)")
		pyplot.plot(MatchedData['ACP Time'], MatchedData['IQMon FWHM'], 'k.', drawstyle="steps-post", label="FWHM (IQMon)")
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
	if (telescope == "V5" and FoundV5Env and FoundIQMonFile):
		Figure.add_axes([0.54, 0.510, 0.46, 0.235])
		pyplot.plot(V5EnvTable['Time'], V5EnvTable['FocusPos'], 'b-', drawstyle="steps-post", label="Focus Position")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		## Clip the extreme values off of the focus position list (bad values?)
		nFocusPos = len(V5EnvTable['FocusPos'])
		ClippingFactor = 0.02
		nClipped = int(nFocusPos*ClippingFactor)
		ylimlower = (sorted(V5EnvTable['FocusPos']))[nClipped]
		ylimupper = (sorted(V5EnvTable['FocusPos']))[nFocusPos-nClipped]
		ylimrange = max([100 , ylimupper - ylimlower])
		pyplot.ylim(ylimlower-0.15*ylimrange, ylimupper+0.15*ylimrange)
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.grid()
	if (telescope == "V20" and FoundV20Env and FoundIQMonFile):
		Figure.add_axes([0.54, 0.510, 0.46, 0.235])
		pyplot.plot(V20EnvTable['Time'], V20EnvTable['FocusPos'], 'b-', drawstyle="steps-post", label="Focus Position")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		## Clip the extreme values off of the focus position list (bad values?)
		nFocusPos = len(V20EnvTable['FocusPos'])
		ClippingFactor = 0.02
		nClipped = int(nFocusPos*ClippingFactor)
		ylimlower = (sorted(V20EnvTable['FocusPos']))[nClipped]
		ylimupper = (sorted(V20EnvTable['FocusPos']))[nFocusPos-nClipped]
		ylimrange = max([100 , ylimupper - ylimlower])
		pyplot.ylim(ylimlower-0.15*ylimrange, ylimupper+0.15*ylimrange)
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.grid()
	

	###########################################################			
	## Ellipticity vs. Time
	if FoundIQMonFile:
		Figure.add_axes([0.54, 0.255, 0.46, 0.235])
		pyplot.plot(IQMonTable['Time'], IQMonTable['Ellipticity'], 'b.-', drawstyle="steps-post", label="Ellipticity")
		pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylabel("Ellipticity")
		pyplot.ylim(0,1)
		pyplot.grid()


	###########################################################			
	## Pointing Error vs. Time
	if FoundIQMonFile:
		Figure.add_axes([0.54, 0.000, 0.46, 0.235])
		if FoundACPLog:
			pyplot.plot(MatchedData['ACP Time'], MatchedData['ACP PErr'], 'g.', drawstyle="steps-post", label="ACP Log", alpha=0.6)
		pyplot.plot(MatchedData['ACP Time'], MatchedData['IQMon PErr'], 'b.', drawstyle="steps-post", label="IQMon")
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
	
	pyplot.savefig(PlotFile, dpi=dpi, bbox_inches='tight', pad_inches=0.10)
	


	###########################################################			
	## Make Environmental Plot (show entire day)
	###########################################################			
	print "  Writing Output File: "+EnvPlotFileName
	dpi=100
	Figure = pyplot.figure(figsize=(11,11), dpi=dpi)
	PlotStartUT = 0
	PlotEndUT = 24
	nUTHours = 25
	
	###########################################################			
	## Temperatures
	if FoundV20Env or FoundV5Env:
		TemperatureAxes = pyplot.axes([0.0, 0.765, 1.0, 0.235])		
		pyplot.title("Environmental Data for "+telescope + " on the Night of " + DateString)
	
		if telescope == "V20" and FoundV20Env:
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['OutsideTemp'], 'k-', alpha=0.5, drawstyle="steps-post", label="Outside Temp ("+OtherTelescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['TubeTemp'], 'g-', drawstyle="steps-post", label="Tube Temp")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['OutsideTemp'], 'k-', drawstyle="steps-post", label="Outside Temp ("+telescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['PrimaryTemp'], 'r-', drawstyle="steps-post", label="Mirror Temp")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['DomeTemp'], 'c-', drawstyle="steps-post", label="Dome Temp")
		
		if telescope == "V5" and FoundV5Env:
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['OutsideTemp'], 'k-', alpha=0.5, drawstyle="steps-post", label="Outside Temp ("+OtherTelescope+")")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['TubeTemp'], 'g-', drawstyle="steps-post", label="Tube Temp")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['OutsideTemp'], 'k-', drawstyle="steps-post", label="Outside Temp ("+telescope+")")
		
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
		
		## Overplot Moon Up Time
		MoonAxes = TemperatureAxes.twinx()
		MoonAxes.set_ylabel('Moon Alt (%.0f%% full)' % MoonPhase, color='y')
		pyplot.plot(MoonTimes, MoonAlts, 'y-')
		pyplot.ylim(0,100)
		pyplot.yticks([10,30,50,70,90], color='y')
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.fill_between(MoonTimes, 0, MoonAlts, where=MoonAlts>0, color='yellow', alpha=MoonFill)		
				
		## Add Fan Power (if VYSOS-20)
		if telescope == "V20":
			Figure.add_axes([0.0, 0.675, 1.0, 0.07], xticklabels=[])
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['DomeFan'], 'c-', drawstyle="steps-post", label="Dome Fan State")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['FanPower'], 'b-', drawstyle="steps-post", label="Mirror Fans (%)")
			
			pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
			pyplot.xlim(PlotStartUT,PlotEndUT)
			pyplot.ylim(-10,110)
			pyplot.yticks(numpy.linspace(0,100,3,endpoint=True))
			pyplot.legend(loc='best', prop={'size':10})
			pyplot.grid()
	

	###########################################################			
	## Sky Condition (Cloudiness)
	if FoundV20Env or FoundV5Env:
		# Figure.add_axes([0.0, 0.255, 1.0, 0.235])
		if telescope == "V20" and FoundV20Env:
			Figure.add_axes([0.0, 0.430, 1.0, 0.235])
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['SkyTemp'], 'k-', drawstyle="steps-post", alpha=0.5, label="Cloudiness ("+OtherTelescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['SkyTemp'], 'b-', drawstyle="steps-post", label="Cloudiness ("+telescope+")")
			pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="3"), color='red', alpha=0.8)
		if telescope == "V5" and FoundV5Env:
			Figure.add_axes([0.0, 0.510, 1.0, 0.235])
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['SkyTemp'], 'k-', drawstyle="steps-post", alpha=0.5, label="Cloudiness ("+OtherTelescope+")")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['SkyTemp'], 'b-', drawstyle="steps-post", label="Cloudiness ("+telescope+")")
			pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="3"), color='red', alpha=0.8)
			
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Temperature Difference (F)")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylim(-100,-20)
		pyplot.grid()

	###########################################################			
	## Humidity
	if FoundV20Env or FoundV5Env:
		if telescope == "V5" and FoundV5Env:
			# Figure.add_axes([0.0, 0.510, 1.0, 0.235])
			Figure.add_axes([0.0, 0.255, 1.0, 0.235])
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['Humidity'], 'k-', drawstyle="steps-post", alpha=0.5, label="Humidity ("+OtherTelescope+")")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['Humidity'], 'b-', drawstyle="steps-post", label="Humidity ("+telescope+")")
			pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="2"), color='red', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="3"), color='red', alpha=0.8)			
		if telescope == "V20" and FoundV20Env:
			# Figure.add_axes([0.0, 0.510, 1.0, 0.155])
			Figure.add_axes([0.0, 0.255, 1.0, 0.155])
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['Humidity'], 'k-', drawstyle="steps-post", alpha=0.5, label="Humidity ("+OtherTelescope+")")			
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['Humidity'], 'b-', drawstyle="steps-post", label="Humidity ("+telescope+")")
			pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="2"), color='red', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="3"), color='red', alpha=0.8)			
		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Humidity (%)")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.ylim(-5,105)
		pyplot.grid()

	###########################################################			
	## Wind Speed
	if FoundV20Env or FoundV5Env:
		Figure.add_axes([0.0, 0.000, 1.0, 0.235])
		if telescope == "V20" and FoundV20Env:
			if FoundV5Env:
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['WindSpeed'], 'k-', drawstyle="steps-post", alpha=0.5, label="Wind Speed ("+OtherTelescope+")")
			pyplot.plot(V20EnvTable['Time'], V20EnvTable['WindSpeed'], 'b-', drawstyle="steps-post", label="Wind Speed ("+telescope+")")
			pyplot.ylim(0,max([max(V20EnvTable['WindSpeed'])*1.1,35.]))
			pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="3"), color='red', alpha=0.8)
		if telescope == "V5" and FoundV5Env:
			if FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['WindSpeed'], 'k-', drawstyle="steps-post", alpha=0.5, label="Wind Speed ("+OtherTelescope+")")
			pyplot.plot(V5EnvTable['Time'], V5EnvTable['WindSpeed'], 'b-', drawstyle="steps-post", label="Wind Speed ("+telescope+")")
			pyplot.ylim(0,max([max(V5EnvTable['WindSpeed'])*1.1,35.]))
			pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="1"), color='green', alpha=0.5)
			pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="2"), color='yellow', alpha=0.8)
			pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="3"), color='red', alpha=0.8)

		pyplot.legend(loc='best', prop={'size':10})
		pyplot.ylabel("Wind Speed (mph)")
		pyplot.xlabel("Time in Hours UT")
		pyplot.xticks(numpy.linspace(PlotStartUT,PlotEndUT,nUTHours,endpoint=True))
		pyplot.xlim(PlotStartUT,PlotEndUT)
		pyplot.grid()

	pyplot.savefig(EnvPlotFile, dpi=dpi, bbox_inches='tight', pad_inches=0.10)



	###########################################################			
	## Make Recent Conditions Plot (Last 2 hours)
	###########################################################			
	print "  Writing Output File: "+RecentPlotFileName
	dpi=100
	Figure = pyplot.figure(figsize=(11.2,5.8), dpi=dpi)
	now = datetime.datetime.utcnow()
	nowDateString = "%04d%02d%02dUT" % (now.year, now.month, now.day)
	nowDecimal = now.hour + now.minute/60. + now.second/3600.
	
	if nowDateString == DateString:
		if nowDecimal < 2:
			PlotStartUT = 0
			PlotEndUT = 2
		else:
			PlotStartUT = nowDecimal-2.
			PlotEndUT = nowDecimal
		nUTHours = 3
	
		###########################################################			
		## Temperatures
		if FoundV20Env or FoundV5Env:
			TemperatureAxes = pyplot.axes([0.0, 0.53, 0.45, 0.47])
			pyplot.title("Recent Environmental Data for "+telescope)
	
			if telescope == "V20" and FoundV20Env:
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['TubeTemp'], 'g-', drawstyle="steps-post", label="Tube Temp")
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['OutsideTemp'], 'k-', drawstyle="steps-post", label="Outside Temp ("+telescope+")")
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['PrimaryTemp'], 'r-', drawstyle="steps-post", label="Mirror Temp")
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['DomeTemp'], 'c-', drawstyle="steps-post", label="Dome Temp")
		
			if telescope == "V5" and FoundV5Env:
				if FoundV20Env:
					pyplot.plot(V20EnvTable['Time'], V20EnvTable['OutsideTemp'], 'k-', alpha=0.5, drawstyle="steps-post", label="Outside Temp ("+OtherTelescope+")")
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['TubeTemp'], 'g-', drawstyle="steps-post", label="Tube Temp")
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['OutsideTemp'], 'k-', drawstyle="steps-post", label="Outside Temp ("+telescope+")")
		
			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Temperature (F)")
			pyplot.xticks(numpy.arange(math.floor(PlotStartUT),math.ceil(PlotEndUT),0.25))
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
		
			# ## Overplot Moon Up Time
			# MoonAxes = TemperatureAxes.twinx()
			# MoonAxes.set_yticklabels([])
			# # MoonAxes.set_ylabel('Moon Alt (%.0f%% full)' % MoonPhase, color='y')
			# pyplot.plot(MoonTimes, MoonAlts, 'y-')
			# pyplot.ylim(0,100)
			# pyplot.yticks([10,30,50,70,90], color='y')
			# pyplot.xticks(numpy.arange(math.floor(PlotStartUT),math.ceil(PlotEndUT),0.25))
			# pyplot.xlim(PlotStartUT,PlotEndUT)
			# pyplot.fill_between(MoonTimes, 0, MoonAlts, where=MoonAlts>0, color='yellow', alpha=MoonFill)		
		
			## Add Fan Power (if VYSOS-20)
			if telescope == "V20":
				Figure.add_axes([0.0, 0.34, 0.45, 0.13], xticklabels=[])
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['DomeFan'], 'c-', drawstyle="steps-post", label="Dome Fan State")
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['FanPower'], 'b-', drawstyle="steps-post", label="Mirror Fans (%)")

				pyplot.xticks(numpy.arange(math.floor(PlotStartUT),math.ceil(PlotEndUT),0.25))
				pyplot.xlim(PlotStartUT,PlotEndUT)
				pyplot.ylim(-10,110)
				pyplot.yticks(numpy.linspace(0,100,3,endpoint=True))
				pyplot.legend(loc='center left', prop={'size':10})
				pyplot.grid()	


		###########################################################			
		## Humidity
		if FoundV20Env or FoundV5Env:
			if telescope == "V5" and FoundV5Env:
				Figure.add_axes([0.0, 0.0, 0.45, 0.47])
				if FoundV20Env:
					pyplot.plot(V20EnvTable['Time'], V20EnvTable['Humidity'], 'k-', drawstyle="steps-post", alpha=0.5, label="Humidity ("+OtherTelescope+")")
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['Humidity'], 'b-', drawstyle="steps-post", label="Humidity ("+telescope+")")
				pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="1"), color='green', alpha=0.5)
				pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="2"), color='red', alpha=0.5)
				pyplot.fill_between(V5EnvTable['Time'], -5, V5EnvTable['Humidity'], where=(V5EnvTable['WetCondition']=="3"), color='red', alpha=0.8)			
			if telescope == "V20" and FoundV20Env:
				Figure.add_axes([0.0, 0.0, 0.45, 0.33])
				if FoundV5Env:
					pyplot.plot(V5EnvTable['Time'], V5EnvTable['Humidity'], 'k-', drawstyle="steps-post", alpha=0.5, label="Humidity ("+OtherTelescope+")")			
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['Humidity'], 'b-', drawstyle="steps-post", label="Humidity ("+telescope+")")
				pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="1"), color='green', alpha=0.5)
				pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="2"), color='red', alpha=0.5)
				pyplot.fill_between(V20EnvTable['Time'], -5, V20EnvTable['Humidity'], where=(V20EnvTable['WetCondition']=="3"), color='red', alpha=0.8)			
			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Humidity (%)")
			pyplot.xticks(numpy.arange(math.floor(PlotStartUT),math.ceil(PlotEndUT),0.25))
			pyplot.xlim(PlotStartUT,PlotEndUT)
			pyplot.ylim(-5,105)
			pyplot.grid()

		###########################################################			
		## Sky Condition (Cloudiness)
		if FoundV20Env or FoundV5Env:
			Figure.add_axes([0.53, 0.53, 0.45, 0.47])
			if telescope == "V20" and FoundV20Env:
				if FoundV5Env:
					pyplot.plot(V5EnvTable['Time'], V5EnvTable['SkyTemp'], 'k-', drawstyle="steps-post", alpha=0.5, label="Cloudiness ("+OtherTelescope+")")
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['SkyTemp'], 'b-', drawstyle="steps-post", label="Cloudiness ("+telescope+")")
				pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="1"), color='green', alpha=0.5)
				pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="2"), color='yellow', alpha=0.8)
				pyplot.fill_between(V20EnvTable['Time'], -140, V20EnvTable['SkyTemp'], where=(V20EnvTable['CloudCondition']=="3"), color='red', alpha=0.8)
			if telescope == "V5" and FoundV5Env:
				if FoundV20Env:
					pyplot.plot(V20EnvTable['Time'], V20EnvTable['SkyTemp'], 'k-', drawstyle="steps-post", alpha=0.5, label="Cloudiness ("+OtherTelescope+")")
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['SkyTemp'], 'b-', drawstyle="steps-post", label="Cloudiness ("+telescope+")")
				pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="1"), color='green', alpha=0.5)
				pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="2"), color='yellow', alpha=0.8)
				pyplot.fill_between(V5EnvTable['Time'], -140, V5EnvTable['SkyTemp'], where=(V5EnvTable['CloudCondition']=="3"), color='red', alpha=0.8)

			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Temperature Difference (F)")
			pyplot.xticks(numpy.arange(math.floor(PlotStartUT),math.ceil(PlotEndUT),0.25))
			pyplot.xlim(PlotStartUT,PlotEndUT)
			pyplot.ylim(-100,-20)
			pyplot.grid()

		###########################################################			
		## Wind Speed
		if FoundV20Env or FoundV5Env:
			Figure.add_axes([0.53, 0.0, 0.45, 0.47])
			if telescope == "V20" and FoundV20Env:
				if FoundV5Env:
					pyplot.plot(V5EnvTable['Time'], V5EnvTable['WindSpeed'], 'k-', alpha=0.5, drawstyle="steps-post", label="Wind Speed ("+OtherTelescope+")")
				pyplot.plot(V20EnvTable['Time'], V20EnvTable['WindSpeed'], 'b-', drawstyle="steps-post", label="Wind Speed ("+telescope+")")
				pyplot.ylim(0,max([max(V20EnvTable['WindSpeed'])*1.1,35.]))
				pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="1"), color='green', alpha=0.5)
				pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="2"), color='yellow', alpha=0.8)
				pyplot.fill_between(V20EnvTable['Time'], 0, V20EnvTable['WindSpeed'], where=(V20EnvTable['WindCondition']=="3"), color='red', alpha=0.8)
			if telescope == "V5" and FoundV5Env:
				if FoundV20Env:
					pyplot.plot(V20EnvTable['Time'], V20EnvTable['WindSpeed'], 'k-', alpha=0.5, drawstyle="steps-post", label="Wind Speed ("+OtherTelescope+")")
				pyplot.plot(V5EnvTable['Time'], V5EnvTable['WindSpeed'], 'b-', drawstyle="steps-post", label="Wind Speed ("+telescope+")")
				pyplot.ylim(0,max([max(V5EnvTable['WindSpeed'])*1.1,35.]))
				pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="1"), color='green', alpha=0.5)
				pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="2"), color='yellow', alpha=0.8)
				pyplot.fill_between(V5EnvTable['Time'], 0, V5EnvTable['WindSpeed'], where=(V5EnvTable['WindCondition']=="3"), color='red', alpha=0.8)

			pyplot.legend(loc='best', prop={'size':10})
			pyplot.ylabel("Wind Speed (mph)")
			pyplot.xlabel("Time in Hours UT")
			pyplot.xticks(numpy.arange(math.floor(PlotStartUT),math.ceil(PlotEndUT),0.25))
			pyplot.xlim(PlotStartUT,PlotEndUT)
			pyplot.grid()

		pyplot.savefig(RecentPlotFile, dpi=dpi, bbox_inches='tight', pad_inches=0.10)
	
	
	
	
	
	
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
	
	if telescope == "V5" or telescope == "VYSOS5" or telescope == "VYSOS-5":
		telescope = "V5"
	if telescope == "V20" or telescope == "VYSOS20" or telescope == "VYSOS-20":
		telescope = "V20"

	Success = MakePlots(DateString, telescope)
	
		


if __name__ == '__main__':
	main()

