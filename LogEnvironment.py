#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2012-11-16.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import getopt
import os
import re
import sys
import time
import string
import logging
import shutil
import win32com.client
import math
import numpy
import matplotlib.pyplot as pyplot
import asciitable
import datetime
# import ephem
import urllib
from xml.dom import minidom

class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

################################################################
## TimeStringToDecimal
def TimeStringToDecimal(TimeString):
	hms = string.split(TimeString, ":")
	DecimalTime = float(hms[0])+ float(hms[1])/60.0 + float(hms[2])/60.0/60.0
	return DecimalTime


################################################################
## ConvertHSTtoUTString
def ConvertHSTtoUTString(TimeString):
	hmsHST = string.split(TimeString, ":")
	if int(hmsHST[0]) >= 14:
		UTString = str(int(hmsHST[0])+10-24)+":"+hmsHST[1]+":"+hmsHST[2]
	else:
		UTString = str(int(hmsHST[0])+10)+":"+hmsHST[1]+":"+hmsHST[2]
	return UTString


################################################################
## GetFMTemp
def GetFMTemp(FocusMax):
	if not FocusMax.Link:
		try:
			FocusMax.Link = True
		except:
			return -999.0
		
	FocuserTemps = []
	FocuserTemp = float("NaN")
	for i in range(0,4,1):
		try:
			newTemp = FocusMax.Temperature*9./5. + 32.
		except:
			newTemp = -999.
		if not (newTemp >= 90.0 or newTemp <= 10.0):
			FocuserTemps.append(newTemp)
		else:
			pass
			# print "## FocusMax Returned Temperature of %f" % newTemp
	if len(FocuserTemps) > 5:
		FocuserTemp = numpy.median(FocuserTemps)
	
	return FocuserTemp


################################################################
## GetFMPos
def GetFMPos(FocusMax):
	if not FocusMax.Link:
		try:
			FocusMax.Link = True
		except:
			return -1
		
	FocuserPositions = []
	for i in range(0,4,1):
		newPos = FocusMax.Position
		FocuserPositions.append(newPos)
	FocuserPos = numpy.median(FocuserPositions)
					
	return FocuserPos

################################################################
## GetRCOSData
def GetRCOSData(RCOST, RCOSF):
	
	## Get Ambient Temperature
	TrussTemps = []
	for i in range(0,4,1):
		try:
			newTrussTemp = RCOST.AmbientTemp
		except:
			newTrussTemp = -999.
		if not (newTrussTemp >= 99.0 or newTrussTemp <= 10.0):
			TrussTemps.append(newTrussTemp)
		else:
			print "## Focuser Returned Ambient Temperature of %f" % newTrussTemp
	if len(TrussTemps) > 2:
		TrussTemp = numpy.median(TrussTemps)
	else:
		TrussTemp = float("NaN")
		
	## Get Primary Temperature
	PriTemps = []
	for i in range(0,4,1):
		try:
			newPriTemp = RCOST.PrimaryTemp
		except:
			newPriTemp = -999.
		if not (newPriTemp >= 99.0 or newPriTemp <= 10.0):
			PriTemps.append(newPriTemp)
		else:
			print "## Focuser Returned Primary Temperature of %f" % newPriTemp
	if len(PriTemps) > 2:
		PriTemp = numpy.median(PriTemps)
	else:
		PriTemp = float("NaN")

	## Get Secondary Temperature
	SecTemps = []
	for i in range(0,4,1):
		try:
			newSecTemp = RCOST.SecondaryTemp
		except:
			newSecTemp = -999.
		if not (newSecTemp >= 99.0 or newSecTemp <= 10.0):
			SecTemps.append(newSecTemp)
		else:
			print "## Focuser Returned Secondary Temperature of %f" % newSecTemp
	if len(SecTemps) > 2:
		SecTemp = numpy.median(SecTemps)
	else:
		SecTemp = float("NaN")
	
	## Get Fan Speed
	FanSpeeds = []
	for i in range(0,4,1):
		try:
			newFanSpeed = RCOST.FanSpeed
		except:
			newFanSpeed = -999.
		if (newFanSpeed > 100.0 or newFanSpeed < 0.0):
			print "## Focuser Returned Fan Speed of %f" % newFanSpeed
		else:
			FanSpeeds.append(newFanSpeed)
	if len(FanSpeeds) > 2:
		FanSpeed = numpy.median(FanSpeeds)
	else:
		FanSpeed = -1

	## Get FocusPosition
	Positions = []
	for i in range(0,4,1):
		try:
			newPos = RCOSF.Position
		except:
			newPos = -999.
		if (newPos > 100000.0 or newPos < 0.0):
			print "## Focuser Returned Position of %f" % newPos
		else:
			Positions.append(newPos)
	if len(Positions) > 2:
		FocusPos = numpy.median(Positions)
	else:
		FocusPos = -1
	
	
	return 	TrussTemp, PriTemp, SecTemp, FanSpeed, FocusPos



################################################################
## GetClarity
def GetClarity(DecimalTime):
	ClarityDataFile = os.path.join("C:\\", "Users", "vysosuser", "Documents", "ClarityII", "ClarityData.txt")
	
	SkyTempValues = []
	AmbTempValues = []
	WindSpeedValues = []
	HumidityValues = []
	DewPointValues = []
	CloudConditionValues = []
	WindConditionValues = []
	RainConditionValues = []
	DayConditionValues = []
	RoofCloseValues = []
	
	
	nAverage = 6
	for i in range(0,nAverage,1):
		ClarityGood = False
		ClarityData = asciitable.read(ClarityDataFile, guess=False, delimiter=" ", comment="#", data_start=0, Reader=asciitable.NoHeader)
		Date             = ClarityData[0][0]
		Time             = ClarityData[0][1]
		TempUnits        = ClarityData[0][2]
		WindUnits        = ClarityData[0][3]
		SkyTemp          = ClarityData[0][4]
		AmbTemp          = ClarityData[0][5]
		SensorTemp       = ClarityData[0][6]
		WindSpeed        = ClarityData[0][7]
		Humidity         = ClarityData[0][8]
		DewPoint         = ClarityData[0][9]
		Heater           = ClarityData[0][10]
		RainFlag         = ClarityData[0][11]
		WetFlag          = ClarityData[0][12]
		Since            = ClarityData[0][13]
		VB6Now           = ClarityData[0][14]
		CloudCondition   = ClarityData[0][15]
		WindCondition    = ClarityData[0][16]
		RainCondition    = ClarityData[0][17]
		DayCondition     = ClarityData[0][18]
		RoofClose        = ClarityData[0][19]
		## Cloud Condition   | Wind Condition   | Rain Condition | Day Condition    | RoofClose
		##   0 = unknown     |   0 = unknown    |   0 = unknown  |   0 = unknown    |  0 = not requested
		##   1 = clear       |   1 = calm       |   1 = dry      |   1 = dark       |  1 = requested
		##   2 = cloudy      |   2 = windy      |   2 = wet      |   2 = light      |
		##   3 = very cloudy |   3 = very windy |   3 = rain     |   3 = very light | 

		## Check that Data File is not Stale
		ClarityTime       = TimeStringToDecimal(ConvertHSTtoUTString(Time))
		if abs(ClarityTime - DecimalTime) <= 1.0/60.0:  ## Time Difference of Less Than 1.0 Minute
			ClarityGood = True

		## Standardize Temperature Units to Farenheit (F)
		if TempUnits == "C":
			## Change Units to F
			SkyTemp    = SkyTemp*9./5. + 32.
			AmbTemp    = AmbTemp*9./5. + 32.
			DewPoint   = DewPoint*9./5. + 32.
		elif TempUnits == "F":
			pass
		else:
			ClarityGood=False
		## Standardize Wind Speed Units to Miles Per Hour (MPH)
		if WindUnits == "M":
			## units are MPH
			pass
		elif WindUnits == "K":
			## units are KPH
			WindSpeed = WindSpeed * 0.621
		elif WindUnits == "m":
			WindSpeed = WindSpeed * 2.237
		else:
			ClarityGood == False
		## Append Values to Arrays for averaging
		if ClarityGood:
			SkyTempValues.append(SkyTemp)
			AmbTempValues.append(AmbTemp)
			WindSpeedValues.append(WindSpeed)
			HumidityValues.append(Humidity)
			DewPointValues.append(DewPoint)
			CloudConditionValues.append(CloudCondition)
			WindConditionValues.append(WindCondition)
			RainConditionValues.append(RainCondition)
			DayConditionValues.append(DayCondition)
			RoofCloseValues.append(RoofClose)
		## If Problems, set values to either NaN or 0
		# if ClarityGood == False:
		# 	SkyTemp    = float("NaN")
		# 	AmbTemp    = float("NaN")
		# 	WindSpeed  = float("NaN")
		# 	Humidity   = float("NaN")
		# 	DewPoint   = float("NaN")
		# 	CloudCondition = 0
		# 	WindCondition = 0
		# 	RainCondition = 0
		# 	DayCondition = 0
		# 	RoofClose = 0
		## Wait For Next Reading
		# print "%6.1f %6.1f %6.1f %6.1f %6.1f %1d %1d %1d %1d %1d" % (SkyTemp, AmbTemp, WindSpeed, Humidity, DewPoint, CloudCondition, WindCondition, RainCondition, DayCondition, RoofClose)
		time.sleep(3)
		
	## Average Values
	SkyTemp   = numpy.median(SkyTempValues)
	AmbTemp   = numpy.median(AmbTempValues)
	WindSpeed = numpy.median(WindSpeedValues)
	Humidity  = numpy.median(HumidityValues)
	DewPoint  = numpy.median(DewPointValues)

	CloudCondition = max(CloudConditionValues)
	WindCondition  = max(WindConditionValues)
	RainCondition  = max(RainConditionValues)
	DayCondition   = max(DayConditionValues)
	RoofClose      = max(RoofCloseValues)
	return [SkyTemp, AmbTemp, WindSpeed, Humidity, DewPoint, CloudCondition, WindCondition, RainCondition, DayCondition, RoofClose]


################################################################
## GetACPStatus
def GetACPStatus(ACP):
	Connected = ACP.Connected
	try:
		Altitude = ACP.Altitude
		Azimuth  = ACP.Azimuth
		AtPark   = ACP.AtPark
		Slewing  = ACP.Slewing
		Tracking = ACP.Tracking
	except:
		Altitude = float("NaN")
		Azimuth  = float("NaN")
		AtPark   = "Unknown"
		Slewing  = "Unknown"
		Tracking = "Unknown"
		
	return [Altitude, Azimuth, AtPark, Slewing, Tracking, Connected]
		
		
		
################################################################
## GetTemperatureModuleInfo
def GetTemperatureModuleInfo():
	IPaddress = "192.168.1.115"

	try:
		page = urllib.urlopen("http://"+IPaddress+"/state.xml")
		contents = page.read()
		ContentLines = contents.split("\n")	
	
		xmldoc = minidom.parseString(contents)
		units = str(xmldoc.getElementsByTagName('units')[0].firstChild.nodeValue)
		temp1 = float(xmldoc.getElementsByTagName('sensor1temp')[0].firstChild.nodeValue)
		temp2 = float(xmldoc.getElementsByTagName('sensor2temp')[0].firstChild.nodeValue)
		r1state = int(xmldoc.getElementsByTagName('relay1state')[0].firstChild.nodeValue)
		r2state = int(xmldoc.getElementsByTagName('relay2state')[0].firstChild.nodeValue)
		if units == "C":
			temp1 = temp1*9./5. + 32.
			temp2 = temp2*9./5. + 32.
		return [units, temp1, temp2, r1state, r2state]
	except:
		return ["", float("nan"), float("nan"), -1, -1]

################################################################
## SetTemperatureModuleState
## - if UT < sunset or UT > sunrise
## - if InsideTemp > OutsideTemp + DeadbandHigh:  Change State to On
## - if InsideTemp < OutsideTemp + DeadbandLow:  Change State to Off
def SetTemperatureModuleState(InsideTemp, OutsideTemp, RelayState, Enable, logger):
	IPaddress = "192.168.1.115"
	DeadbandHigh = 3.0
	DeadbandLow = 0.25
	Day = False
	
	## Use pyephem determine sunrise and sunset times
	# Observatory = ephem.Observer()
	# Observatory.lon = "-155:34:33.9"
	# Observatory.lat = "+19:32:09.66"
	# Observatory.elevation = 3400.0
	# Observatory.temp = 10.0
	# Observatory.pressure = 680.0
	# Observatory.date = DateString[0:4]+"/"+DateString[4:6]+"/"+DateString[6:8]+" 10:00:00.0"
	# 
	# Observatory.horizon = '0.0'
	# SunsetTime  = Observatory.previous_setting(ephem.Sun()).datetime()
	# SunriseTime = Observatory.next_rising(ephem.Sun()).datetime()
	# SunsetDecimal = float(datetime.datetime.strftime(SunsetTime, "%H"))+float(datetime.datetime.strftime(SunsetTime, "%M"))/60.+float(datetime.datetime.strftime(SunsetTime, "%S"))/3600.
	# SunriseDecimal = float(datetime.datetime.strftime(SunriseTime, "%H"))+float(datetime.datetime.strftime(SunriseTime, "%M"))/60.+float(datetime.datetime.strftime(SunriseTime, "%S"))/3600.
	SunsetDecimal = 5.0
	SunriseDecimal = 16.0
	
	## Is it Daytime
	now = datetime.datetime.utcnow()
	DecimalTime = now.hour+now.minute/60.+now.second/3600.
	if (DecimalTime < SunsetDecimal) or (DecimalTime > SunriseDecimal):	Day = True

	## Temperature Difference
	DeltaT = abs(InsideTemp-OutsideTemp)

	## if Override is Set Turn Fans Off
	if Enable == 0:
		if (RelayState == 1):
			logger.info("Remote Control Set to Off.  Turning Dome Fan Off.  DeltaT = %.1f" % DeltaT)
			page = urllib.urlopen("http://"+IPaddress+"/state.xml?relay1State=0")
		elif (RelayState == 0):
			logger.info("Remote Control Set to Off.  Leaving Dome Fan Off.  DeltaT = %.1f" % DeltaT)
	else:
		## Turn on Fans if Inside Temperature is High
		if Day and (RelayState == 0) and (InsideTemp > OutsideTemp + DeadbandHigh):
			logger.info("Turning Dome Fan On.  DeltaT = %.1f" % DeltaT)
			page = urllib.urlopen("http://"+IPaddress+"/state.xml?relay1State=1")
		## Turn off Fans if Night
		elif not Day and (RelayState == 1):
			logger.info("Turning Dome Fan Off for Night.  DeltaT = %.1f" % DeltaT)
			page = urllib.urlopen("http://"+IPaddress+"/state.xml?relay1State=0")
		## Turn off Fans if Inside Temperature is Low
		elif (RelayState == 1) and (InsideTemp < OutsideTemp + DeadbandLow):
			logger.info("Turning Dome Fan Off.  DeltaT = %.1f" % DeltaT)
			page = urllib.urlopen("http://"+IPaddress+"/state.xml?relay1State=0")
		else:
			if RelayState == 1:
				logger.info("Leaving Dome Fan On.  DeltaT = %.1f", DeltaT)
			if RelayState == 0:
				logger.info("Leaving Dome Fan Off.  DeltaT = %.1f", DeltaT)
		
		
	


################################################################
## Main
def main(argv=None):
	MakePlot = False
	telescope = ""
	
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ht:", ["help", "telescope="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-t", "--telescope"):
				telescope = value
	
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2
	
	###########
	
	## Set up logger
	logger = logging.getLogger('EnvLogger')
	logger.setLevel(logging.DEBUG)
	LogConsoleHandler = logging.StreamHandler()
	LogConsoleHandler.setLevel(logging.INFO)
	LogFormat = logging.Formatter('%(asctime)23s: %(message)s')
	LogConsoleHandler.setFormatter(LogFormat)
	logger.addHandler(LogConsoleHandler)
	
	now = time.gmtime()
	DecimalTime = now[3]+now[4]/60.+now[5]/3600.
	TimeString = time.strftime("%Y/%m/%d %H:%M:%SUT", now)
	DateString = time.strftime("%Y%m%dUT", now)

	ACP = win32com.client.Dispatch("ACP.Telescope")
	if telescope == "V5":
		FocusMax = win32com.client.Dispatch("FocusMax.Focuser")
		# FeatherTouch = win32com.client.Dispatch("FeatherTouch.Focuser")
	elif telescope == "V20":
		RCOST = win32com.client.Dispatch("RCOS_AE.Temperature")
		RCOSF = win32com.client.Dispatch("RCOS_AE.Focuser")
		

	# if telescope == "V5":
	# 	print "%-21s %7s %6s %7s %7s %5s %5s %6s %6s %6s" % ("## Time", "FM_Temp", "FM_Pos", "SkyTemp", "AmbTemp", "Wind", "Humid", "DewPt", "Alt", "Az")
	# if telescope == "V20":
	# 	print "%-21s %7s %7s %7s %6s %8s %7s %7s %5s %5s %6s %6s %6s" % ("## Time", "AirTemp", "PriTemp", "SecTemp", "FanSpd", "FocusPos", "SkyTemp", "AmbTemp", "Wind", "Humid", "DewPt", "Alt", "Az")
	while telescope != "":
		## Get Time
		now = time.gmtime()
		DecimalTime = now[3]+now[4]/60.+now[5]/3600.
		TimeString = time.strftime("%Y/%m/%d %H:%M:%SUT", now)
		DateString = time.strftime("%Y%m%dUT", now)

		## Set Log File name and create directory if needed
		if telescope == "V5":
			LogFilePath = os.path.join("C:\\", "Data_V5", "Logs", DateString)
		if telescope == "V20":
			LogFilePath = os.path.join("C:\\", "Data_V20", "Logs", DateString)
		
		if not os.path.exists(LogFilePath): 
			os.mkdir(LogFilePath)
		LogFileName = "EnvironmentalLog.txt"
		LogFile = os.path.join(LogFilePath, LogFileName)

		##
		## Get Focus Temperature & Position
		##
		if telescope == "V5":
			logger.info("Getting Data from FocusMax")
			FocuserTemp = GetFMTemp(FocusMax)
			FocuserPos  = GetFMPos(FocusMax)
			logger.info("  Temp = %.1f, Position = %d" % (FocuserTemp, FocuserPos))
		if telescope == "V20":
			logger.info("Getting Data from RCOS TCC")
			TrussTemp, PriTemp, SecTemp, FanSpeed, FocusPos = GetRCOSData(RCOST, RCOSF)
			logger.info("  TrussTemp=%.1f, PriTemp=%.1f, Position=%.0f, Fan=%.0f" % (TrussTemp, PriTemp, FocusPos, FanSpeed))
		
		##
		## Get ACP Status
		##
		logger.info("Getting ACP Status")
		Altitude, Azimuth, AtPark, Slewing, Tracking, Connected = GetACPStatus(ACP)	
		if Connected:
			logger.info("  Alt=%.1f, Az=%.1f, Park=%s, Slew=%s, Track=%s" % (Altitude, Azimuth, AtPark, Slewing, Tracking))
		else:
			logger.info("  ACP Not Connected.")

		##
		## Get Clarity Data
		##
		logger.info("Getting Clarity/Boltwood Data")
		ClaritySkyTemp, ClarityAmbTemp, ClarityWindSpeed, ClarityHumidity, ClarityDewPoint, CloudCondition, WindCondition, RainCondition, DayCondition, RoofClose = GetClarity(DecimalTime)
		TempDiff = ClaritySkyTemp - ClarityAmbTemp
		logger.info("  Sky=%.1f, Ambient=%.1f" % (ClaritySkyTemp, ClarityAmbTemp))
		logger.info("  Safety Conditions (CWRDR): %d%d%d%d%d" % (CloudCondition, WindCondition, RainCondition, DayCondition, RoofClose))

		# print CloudCondition, WindCondition, RainCondition, DayCondition, RoofClose

		##
		## Get Temperature Module Info
		##
		if telescope == "V20":
			logger.info("Getting Temperature Module Info")
			TMInfo = GetTemperatureModuleInfo()
			logger.info("  Temperature = %.1f %1s" % (TMInfo[1], TMInfo[0]))
			logger.info("  Relay Status = %d %d" % (TMInfo[3], TMInfo[4]))
			SetTemperatureModuleState(TrussTemp, ClarityAmbTemp, TMInfo[3], TMInfo[4], logger)

		##
		## Write HTML snippet for inclusion in status web page
		##
		logger.info("Writing HTML snippet for status page.")
		HTMLfilename = os.path.join("C:\\", "Data_"+telescope, "Boltwood_"+telescope+".html")
		HTMLsnippet = open(HTMLfilename, 'w')
		## Telescope Name With Roof Close Condition
		HTMLsnippet.write("<tr>\n")
		if RoofClose == 0:
			if telescope == "V5":
				HTMLsnippet.write("  <td colspan=1 style='background-color:green;align:center'><a style='background-color:green;align:center' href='VYSOS5.html'>VYSOS-5</a></td>\n")
			if telescope == "V20":
				HTMLsnippet.write("  <td colspan=1 style='background-color:green;align:center'><a style='background-color:green;align:center' href='VYSOS20.html'>VYSOS-20</a></td>\n")
			if Connected:
				HTMLsnippet.write(  "<td>Alt = %.1f, Az = %.1f</td>\n" % (Altitude, Azimuth))
				StatusString = ""
				if AtPark:
					StatusString = "Parked"
				elif Slewing:
					StatusString = "Slewing"
				elif Tracking:
					StatusString = "Tracking"
				else:
					StatusString = "Not Tracking"
				HTMLsnippet.write(  "<td>%s</td>\n" % StatusString)
			else:
				HTMLsnippet.write(  "<td>Telescope not connected</td>\n")
				HTMLsnippet.write(  "<td>Telescope not connected</td>\n")
		if RoofClose == 1:
			AlertString = ""
			if CloudCondition == 3:
				AlertString += "Clouds"
			if WindCondition == 3:
				AlertString += "Wind"
			if RainCondition >= 2:
				AlertString += "Rain"
			if DayCondition == 3:
				AlertString += "Day"
			if AlertString == "":
				AlertString = "None"
			if telescope == "V5":
				HTMLsnippet.write("  <td colspan=1 style='background-color:red;align:center'><a style='background-color:red;align:center' href='VYSOS5.html'>VYSOS-5</a> (Alerts = %s)</td>\n" % AlertString)
			if telescope == "V20":
				HTMLsnippet.write("  <td colspan=1 style='background-color:red;align:center'><a style='background-color:red;align:center' href='VYSOS20.html'>VYSOS-20</a> (Alerts = %s)</td>\n" % AlertString)
			if Connected:
				HTMLsnippet.write(  "<td>Alt = %.1f, Az = %.1f</td>\n" % (Altitude, Azimuth))
				StatusString = ""
				if AtPark:
					StatusString = "Parked"
				elif Slewing:
					StatusString = "Slewing"
				elif Tracking:
					StatusString = "Tracking"
				else:
					StatusString = "Not Tracking"				
				HTMLsnippet.write(  "<td>%s</td>\n" % StatusString)
			else:
				HTMLsnippet.write(  "<td>Telescope not connected</td>\n")
				HTMLsnippet.write(  "<td>Telescope not connected</td>\n")
		HTMLsnippet.write("</tr>\n")
		HTMLsnippet.write("<tr>\n")
		## Cloudiness
		if CloudCondition == 0:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>Cloudiness Unknown (%.1f F)</td>\n" % TempDiff)
		elif CloudCondition == 1:
			HTMLsnippet.write("  <td style='background-color:green;align:center'>Clear (%.1f F)</td>\n" % TempDiff)
		elif CloudCondition == 2:
			HTMLsnippet.write("  <td style='background-color:yellow;align:center'>Cloudy (%.1f F)</td>\n" % TempDiff)
		elif CloudCondition == 3:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>Very Cloudy (%.1f F)</td>\n" % TempDiff)
		else:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>State Reporting Error (%.1f F)</td>\n" % TempDiff)
		## Wind
		if RainCondition == 0:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>Wind Speed Unknown (%.1f mph)</td>\n" % ClarityWindSpeed)
		elif RainCondition == 1:
			HTMLsnippet.write("  <td style='background-color:green;align:center'>Calm (%.1f mph)</td>\n" % ClarityWindSpeed)
		elif RainCondition == 2:
			HTMLsnippet.write("  <td style='background-color:yellow;align:center'>Windy (%.1f mph)</td>\n" % ClarityWindSpeed)
		elif RainCondition == 3:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>Very Windy (%.1f mph)</td>\n" % ClarityWindSpeed)
		else:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>State Reporting Error (%.1f mph)</td>\n" % ClarityWindSpeed)
		## Rain
		if RainCondition == 0:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>Rain: Unknown</td>\n")
		elif RainCondition == 1:
			HTMLsnippet.write("  <td style='background-color:green;align:center'>Rain: Dry</td>\n")
		elif RainCondition == 2:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>Rain: Wet</td>\n")
		elif RainCondition == 3:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>Rain: Rain</td>\n")
		else:
			HTMLsnippet.write("  <td style='background-color:red;align:center'>error</td>\n")
		
		HTMLsnippet.write("</tr>\n")
		HTMLsnippet.close()
		
		##
		## Update Log File
		##
		logger.info("Writing data to log file.")
		if telescope != "":
			## Read in Previous Logs (to avoid clobbering them when new ones written)
			PreviousLogs = ""
			if not os.path.exists(LogFile):
				output = open(LogFile, 'w')
				if telescope == "V5":
					output.write("%-22s" % "#")
					output.write("%10s" % "Tube")
					output.write("%10s" % "Focus")
					output.write("%10s" % "Sky")
					output.write("%10s" % "Outside")
					output.write("%10s" % "WindSpd")
					output.write("%10s" % "Humid")
					output.write("%10s" % "DewPt")
					output.write("%10s" % "Alt")
					output.write("%10s" % "Az")
					output.write("%10s\n" % "WetCldWnd")
					output.write("%-22s" % "# Date & Time UT")
					output.write("%10s" % "(F)")
					output.write("%10s" % "Pos.")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(km/h)")
					output.write("%10s" % "(%)")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(deg)")
					output.write("%10s" % "(deg)")
					output.write("%10s\n" % "()")
				if telescope == "V20":
					output.write("%-22s" % "#")
					output.write("%10s" % "Tube")
					output.write("%10s" % "Primary")
					output.write("%10s" % "Sec.")
					output.write("%10s" % "FanPwr")
					output.write("%10s" % "Focus")
					output.write("%10s" % "Sky")
					output.write("%10s" % "Outside")
					output.write("%10s" % "WindSpd")
					output.write("%10s" % "Humid")
					output.write("%10s" % "DewPt")
					output.write("%10s" % "Alt")
					output.write("%10s" % "Az")
					output.write("%10s" % "Wetness")
					output.write("%10s" % "Dome")
					output.write("%10s\n" % "DomeFan")
					output.write("%-22s" % "# Date & Time UT")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(%)")
					output.write("%10s" % "Pos.")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(km/h)")
					output.write("%10s" % "(%)")
					output.write("%10s" % "(F)")
					output.write("%10s" % "(deg)")
					output.write("%10s" % "(deg)")
					output.write("%10s" % "()")
					output.write("%10s" % "(F)")
					output.write("%10s\n" % "()")
					
				output.close()

			input = open(LogFile, 'r')
			PreviousLogs = input.read()
			input.close()
			
			## Open Output and Write Previous Data
			output = open(LogFile, 'w')
			output.write(PreviousLogs)

			## Write Line of New Data to Log
			WetCldWnd = str(RainCondition)+str(CloudCondition)+str(WindCondition)
			if telescope == "V5":
				output.write("%-22s" % TimeString)
				output.write("%10.2f" % FocuserTemp)
				output.write("%10d" % FocuserPos)		
				output.write("%10.2f" % ClaritySkyTemp)
				output.write("%10.2f" % ClarityAmbTemp)
				output.write("%10.1f" % ClarityWindSpeed)
				output.write("%10.0f" % ClarityHumidity)
				output.write("%10.2f" % ClarityDewPoint)
				output.write("%10.2f" % Altitude)
				output.write("%10.2f" % Azimuth)
				output.write("%10s\n" % WetCldWnd)
			elif telescope == "V20":
				output.write("%-22s" % TimeString)
				output.write("%10.2f" % TrussTemp)
				output.write("%10.2f" % PriTemp)
				output.write("%10.2f" % SecTemp)
				output.write("%10d" % FanSpeed)
				output.write("%10d" % FocusPos)
				output.write("%10.2f" % ClaritySkyTemp)
				output.write("%10.2f" % ClarityAmbTemp)
				output.write("%10.1f" % ClarityWindSpeed)
				output.write("%10.0f" % ClarityHumidity)
				output.write("%10.2f" % ClarityDewPoint)
				output.write("%10.2f" % Altitude)
				output.write("%10.2f" % Azimuth)
				output.write("%10s" % WetCldWnd)
				output.write("%10.1f" % TMInfo[1])
				output.write("%9d%1d\n" % (TMInfo[3], TMInfo[4]))
			output.close()

			## Print Same Line to Screen
			# if telescope == "V5":
			# 	print "%-21s %7.2f %6d %7.2f %7.2f %5.1f %5.0f %6.2f %6.2f %6.2f" % (TimeString, FocuserTemp, FocuserPos, 
			# 	                        ClaritySkyTemp, ClarityAmbTemp, ClarityWindSpeed, ClarityHumidity, ClarityDewPoint, 
			# 	                        Altitude, Azimuth)
			# if telescope == "V20":
			# 	print "%-21s %7.2f %7.2f %7.2f %6d %8d %7.2f %7.2f %5.1f %5.0f %6.2f %6.2f %6.2f" % (TimeString, 
			# 	             TrussTemp, PriTemp, SecTemp, FanSpeed, FocusPos,
			# 	             ClaritySkyTemp, ClarityAmbTemp, ClarityWindSpeed, ClarityHumidity, ClarityDewPoint, 
			# 	             Altitude, Azimuth)
			
		##
		## Wait
		##
		logger.info("#### Waiting 5 seconds ... ####")
		time.sleep(5)

if __name__ == '__main__':
	main()

