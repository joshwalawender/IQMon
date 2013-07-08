import sys
import os
import re
import math
import numpy
import logging
import subprocess32
import StringIO
import time

import astropy
import astropy.io
import astropy.io.ascii
import pyfits
from scipy.optimize import curve_fit
import matplotlib.pyplot as pyplot
from collections import defaultdict
from pyraf import iraf




##############################################################
## Read Configuration File to get the following items
def ReadConfigFile():
	log = []
	## Look for configuration file
	HomePath = os.path.expandvars("$HOME")
	ConfigFilePath = os.path.join(HomePath, ".IQMonConfig")
	if os.path.exists(ConfigFilePath):
		ConfigFile = open(ConfigFilePath, 'r')
		ConfigFileLines = ConfigFile.readlines()
		ConfigFile.close()
	else:
		ConfigFileLines = []

	## read configuration file
	FoundIQMonExecPath = False
	FoundLogPath = False
	FoundPlotsPath = False
	FoundtmpPath = False
	FoundPythonPath = False
	FoundV5DataPath = False
	FoundV20DataPath = False
	FoundCatalogPath = False
	for line in ConfigFileLines:
		IsIQMonExecPath = re.match("IQMONPATH\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IsIQMonExecPath:
			FoundIQMonExecPath = True
			IQMonExecPath = os.path.abspath(IsIQMonExecPath.group(1))
		IsLogPath = re.match("IQMONLOGS\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IsLogPath:
			FoundLogPath = True
			LogPath = os.path.abspath(IsLogPath.group(1))
		IsPlotsPath = re.match("IQMONPLOTS\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IsPlotsPath:
			FoundPlotsPath = True
			PlotsPath = os.path.abspath(IsPlotsPath.group(1))
		IstmpPath = re.match("IQMONTMP\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IstmpPath:
			FoundtmpPath = True
			tmpPath = os.path.abspath(IstmpPath.group(1))
		IsPythonPath = re.match("IQMONPYTHON\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IsPythonPath:
			FoundPythonPath = True
			PythonPath = os.path.abspath(IsPythonPath.group(1))
		IsV5DataPath = re.match("VYSOS5DATAPATH\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IsV5DataPath:
			FoundV5DataPath = True
			V5DataPath = os.path.abspath(IsV5DataPath.group(1))
		IsV20DataPath = re.match("VYSOS20DATAPATH\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IsV20DataPath:
			FoundV20DataPath = True
			V20DataPath = os.path.abspath(IsV20DataPath.group(1))
		IsCatalogPath = re.match("CATALOGPATH\s=\s([a-zA-z0-9/\-_\.]+)", line)
		if IsCatalogPath:
			FoundCatalogPath = True
			CatalogPath = os.path.abspath(IsCatalogPath.group(1))


	## if configuration file not found or not read, then use defaults
	if not FoundIQMonExecPath:
		log.append("  IQMonExecPath being set to default")
		IQMonExecPath = os.path.join(HomePath, "bin", "IQMon")
	if not FoundLogPath:
		log.append("  LogPath being set to default")
		LogPath = os.path.join(HomePath, "IQMon", "Logs")
	if not FoundPlotsPath:
		log.append("  PlotsPath being set to default")
		PlotsPath = os.path.join(HomePath, "IQMon", "Plots")
	if not FoundtmpPath:
		log.append("  tmpPath being set to default")
		tmpPath = os.path.join(HomePath, "IQMon", "tmp")
	if not FoundPythonPath:
		log.append("  Using Default Python Path")
		PythonPath = ""
	if not FoundV5DataPath:
		log.append("  V5DataPath being set to default")
		V5DataPath = os.path.join("/Volumes", "Data_V5")
	if not FoundV20DataPath:
		log.append("  V20DataPath being set to default")
		V20DataPath = os.path.join("/Volumes", "Data_V20")
	if not FoundCatalogPath:
		log.append("  CatalogPath being set to default")
		CatalogPath = os.path.join(HomePath, "UCAC4", "u4b")
		
	return IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, log


##############################################################
## Extract Basic WCS Info from a Header Object
##############################################################
def ExtractWCS(header):
	## Check for WCS
	WCS = {}
	try:
		## Grab initial WCS from image (only zero point of WCS)
		WCS['CRPIX1'] = header['CRPIX1']
		WCS['CRPIX2'] = header['CRPIX2']
		WCS['CRVAL1'] = header['CRVAL1']
		WCS['CRVAL2'] = header['CRVAL2']
		WCS['CD1_1']  = header['CD1_1']
		WCS['CD1_2']  = header['CD1_2']
		WCS['CD2_1']  = header['CD2_1']
		WCS['CD2_2']  = header['CD2_2']
		WCS['Extracted'] = True
		if (abs(WCS['CD2_1']) > abs(WCS['CD2_2'])) and (WCS['CD2_1'] > 0): 
			North = "Right"
			WCS['PA'] = 270. + math.degrees(math.atan(WCS['CD2_2']/WCS['CD2_1']))
		elif (abs(WCS['CD2_1']) > abs(WCS['CD2_2'])) and (WCS['CD2_1'] < 0):
			North = "Left"
			WCS['PA'] = 90. + math.degrees(math.atan(WCS['CD2_2']/WCS['CD2_1']))
		elif (abs(WCS['CD2_1']) < abs(WCS['CD2_2'])) and (WCS['CD2_2'] > 0):
			North = "Up"
			WCS['PA'] = 0. + math.degrees(math.atan(WCS['CD2_1']/WCS['CD2_2']))
		elif (abs(WCS['CD2_1']) < abs(WCS['CD2_2'])) and (WCS['CD2_2'] < 0):
			North = "Down"
			WCS['PA'] = 180. + math.degrees(math.atan(WCS['CD2_1']/WCS['CD2_2']))
		else:
			print WCS['CD1_1']
			print WCS['CD1_2']
			print WCS['CD2_1']
			print WCS['CD2_2']
		if (abs(WCS['CD1_1']) > abs(WCS['CD1_2'])) and (WCS['CD1_1'] > 0): East = "Right"
		if (abs(WCS['CD1_1']) > abs(WCS['CD1_2'])) and (WCS['CD1_1'] < 0): East = "Left"
		if (abs(WCS['CD1_1']) < abs(WCS['CD1_2'])) and (WCS['CD1_2'] > 0): East = "Up"
		if (abs(WCS['CD1_1']) < abs(WCS['CD1_2'])) and (WCS['CD1_2'] < 0): East = "Down"
		WCS['Orientation'] = "North %s, East %s" % (North, East)
		if North == "Up" and East == "Left": WCS['Flipped'] = False
		if North == "Up" and East == "Right": WCS['Flipped'] = True
		if North == "Down" and East == "Left": WCS['Flipped'] = True
		if North == "Down" and East == "Right": WCS['Flipped'] = False
		if North == "Right" and East == "Up": WCS['Flipped'] = False
		if North == "Right" and East == "Down": WCS['Flipped'] = True
		if North == "Left" and East == "Up": WCS['Flipped'] = True
		if North == "Left" and East == "Down": WCS['Flipped'] = False
	except:
		WCS['Extracted'] = False

	return WCS



##############################################################
## Run Astrometry.net Solver
##############################################################
def RunAstrometryDotNet(FitsFile, tel, logger):
	logger.info("  Attempting to create WCS using Astrometry.net solver ...")
	FitsFileDirectory, FitsFilename = os.path.split(FitsFile)
	FitsBasename, FitsExt   = os.path.splitext(FitsFilename)
	## Set pixel scale estimates to 10% above and below nominal pixel scale
	PxScL = str(tel['PixelScale']*0.95)
	PxScH = str(tel['PixelScale']*1.05)
	AstrometryCommand = ["solve-field", "-l", "5", "-O", "-p", "-t", "3", "-L", PxScL, "-H", PxScH, "-u", "arcsecperpix", "-z", "4", FitsFile]
	AstrometrySTDOUT = ""
	
	try:
		StartTime = time.time()
		AstrometrySTDOUT = subprocess32.check_output(AstrometryCommand, stderr=subprocess32.STDOUT, timeout=30)
		EndTime = time.time()
		ProcessTime = EndTime - StartTime
		logger.info("  Astrometry.net Processing Time: %.1f", ProcessTime)
	except TimeoutExpired:
		logger.warning("  Astrometry.net timed out")
		return False
	except:
		logger.warning("  Astrometry.net failed.")
		return False
		
	pos = AstrometrySTDOUT.find("Field center: (RA H:M:S, Dec D:M:S) = ")
	if pos != -1:
		IsFieldCenter = re.match("\s*([0-9]{1,2}:[0-9]{2}:[0-9]{2}\.[0-9]+,\s-?[0-9]{1,2}:[0-9]{2}:[0-9]{2}\.[0-9]+).*", AstrometrySTDOUT[pos+40:pos+75])
		if IsFieldCenter:
			logger.info("  Astrometry.net field center is: %s", IsFieldCenter.group(1))
	else:
		for line in AstrometrySTDOUT:
			logger.warning("  %s" % line)
		

	NewFile        = FitsFile.replace(FitsExt, ".new")
	NewFitsFile    = FitsFile.replace(FitsExt, ".new.fits")
	if not os.path.exists(NewFile):
		logger.warning("  Astrometry.net failed.  No new file created by astrometry.net")
		return False
	else:
		logger.info("  Astrometry.net succeeded")
		if os.path.exists(NewFitsFile): os.remove(NewFitsFile)
		os.rename(NewFile, NewFitsFile)
		return True	




##############################################################
## Get Stars from UCAC4 Catalog
##############################################################
def GetUCAC4(RAcenter, DECcenter, RAsize, DECsize):
	## Get catalog of stars in image from UCAC4 catalog
	## Read Configuration File
	IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, log = ReadConfigFile()
	if os.path.exists("ucac4.txt"): os.remove("ucac4.txt")
	UCAC4command = [os.path.join(IQMonExecPath, "u4test"), "-h", str(RAcenter), str(DECcenter), str(RAsize), str(DECsize), CatalogPath]
	try:
		UCAC4output = subprocess32.check_output(UCAC4command, stderr=subprocess32.STDOUT, timeout=20)
	except:
		pass
	if os.path.exists("ucac4.txt"):
		ColStarts = [ 0, 11, 24, 37, 44, 50, 54, 57, 61, 77, 81, 85, 88, 91, 94, 101, 108, 117, 127, 134, 141, 148, 160, 169, 176, 183, 190, 197, 204, 210, 216, 222, 228]
		ColEnds   = [ 9, 22, 35, 42, 49, 53, 56, 60, 75, 79, 83, 86, 89, 92, 99, 106, 115, 125, 132, 139, 146, 158, 167, 174, 181, 188, 195, 202, 208, 214, 220, 226, 232]
		UCAC4catalog = astropy.io.ascii.read("ucac4.txt", data_start=1, Reader=ascii.FixedWidth, col_starts=ColStarts, col_ends=ColEnds)
		return UCAC4catalog
	else:
		return []


##############################################################
## Get Stars from USNO Catalog
##############################################################
def GetUSNO(RAcenter, DECcenter, boxsize, OutputFile):
	if OutputFile != False and OutputFile != "":
		WriteOutputFile = True
	else:
		WriteOutputFile = False
	FaintLimit = 16.0
	USNOcommand = ["findusnob1", str(RAcenter), str(DECcenter), "-b", str(boxsize*60.), "-m", "10000", "-lmr1", "0,%.1f" % FaintLimit]
	try:
		USNOoutput = subprocess32.check_output(USNOcommand, stderr=subprocess32.STDOUT, timeout=60)
		USNOoutput = USNOoutput.split("\n")[4:-3]
	except:
		return []
	if WriteOutputFile:
		OutputFileObject = open(OutputFile, 'w')
	## Parse Text Output
	data = []
	for line in USNOoutput:
		## Get RA Value
		try:
			RA = float(line[26:36])
		except:
			RA = float("NaN")
		##Get DEC Value
		try:
			DEC = float(line[36:46])
		except:
			DEC = float("Nan")
		## Get Rmag1 Value
		try:
			Rmag1 = float(line[97:103])
		except:
			Rmag1 = float("NaN")
		## Get Rmag2 Value
		try:
			Rmag2 = float(line[190:195])
		except:
			Rmag2 = float("NaN")
		## Add to Data Array
		data.append([RA, DEC, Rmag1, Rmag2])
		if WriteOutputFile:
			OutputFileObject.write("%11.6f %11.6f %6.2f %6.2f\n" % (RA, DEC, Rmag1, Rmag2))
		
	if WriteOutputFile:
		OutputFileObject.close()	
			
	return data


##############################################################
## - if RefineWCS is set
## - copy standard database solution from file based on telescope (V5 or V20)
## - copy CRPIX and CRVAl back from original header
## - get UCAC4 catalog in to temp file
## - use CCFIND to ID stars in catalog and create list of coordinates (x, y, RA, Dec)
## - use CCMAP to determine astrometric solution
##############################################################
def RefineWCS(WorkingFile, tel, WCS, NXPix, NYPix, logger):
	IRAFerrors = StringIO.StringIO()
	logger.info("  Refining WCS")
	WorkingFileDir, WorkingFilename = os.path.split(WorkingFile)
	## Read Configuration File
	IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, log = ReadConfigFile()

	##############################################################
	## Copy Standard Distortion Solution for this telescope
		
	try:
		iraf.hedit(WorkingFile, 'WCSDIM'  , tel['WCSDIM']  , add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'CTYPE1'  , tel['CTYPE1']  , add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'CTYPE2'  , tel['CTYPE2']  , add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT0_001', tel['WAT0_001'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT1_001', tel['WAT1_001'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT1_002', tel['WAT1_002'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT1_003', tel['WAT1_003'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT1_004', tel['WAT1_004'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT1_005', tel['WAT1_005'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT2_001', tel['WAT2_001'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT2_002', tel['WAT2_002'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT2_003', tel['WAT2_003'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT2_004', tel['WAT2_004'], add="yes", update="yes", verify="no", show="no")
		iraf.hedit(WorkingFile, 'WAT2_005', tel['WAT2_005'], add="yes", update="yes", verify="no", show="no")
	except:
		logger.warning("  Error installing standard distortion keywords for %s" % tel['LongName'])
		
	
	##############################################################
	## Get UCAC Stars In Image and Feed in to IRAF.CCFIND
	
	## Determine RA and DEC of Image Center and Size of Image in RA and DEC
	CDarray = numpy.array([[WCS['CD1_1'], WCS['CD1_2']], [WCS['CD2_1'], WCS['CD2_2']]])
	RefPixel     = numpy.array([[WCS['CRPIX1']], [WCS['CRPIX2']]])
	RefCoords    = numpy.array([[WCS['CRVAL1']], [WCS['CRVAL2']]])
	CenterPixel  = numpy.array([[NXPix/2], [NYPix/2]])
	CenterDelta  = numpy.dot(CDarray, CenterPixel-RefPixel)
	CenterCoords = RefCoords + CenterDelta
	LLPixel      = numpy.array([[1], [1]])
	LLDelta      = numpy.dot(CDarray, LLPixel-RefPixel)
	URPixel      = numpy.array([[NXPix], [NYPix]])
	URDelta      = numpy.dot(CDarray, URPixel-RefPixel)
	RAcenter = (CenterCoords[0])[0]
	DECcenter = (CenterCoords[1])[0]
	RAsize  = (abs(LLDelta[0] - URDelta[0]))[0]
	DECsize = (abs(LLDelta[1] - URDelta[1]))[0]
	radius = numpy.sqrt(numpy.square(RAsize/2.) + numpy.square(DECsize/2.))
	boxsize = max(RAsize, DECsize)

	## Get USNO Coordinates of Stars in Image for Plate Solving
	# logger.info("  Getting USNO catalog.")
	# USNOfile = os.path.join(tmpPath, tel['FileID']+"_USNO.txt")
	# if os.path.exists(USNOfile): os.remove(USNOfile)
	# USNOtable = GetUSNO(RAcenter, DECcenter, radius, USNOfile)
	
	## Get UCAC stars
	logger.info("  Getting UCAC catalog stars")
	UCACfile = os.path.join(tmpPath, tel['FileID']+"_UCAC.txt")
	UCAC4Table = GetUCAC4(RAcenter, DECcenter, RAsize, DECsize)
	if len(UCAC4Table) <= 1:
		logger.warning("  Failed to get UCAC4 Catalog")

	UCAC4data = []
	if os.path.exists(UCACfile): os.remove(UCACfile)
	UCACFO = open(UCACfile, "w")
	for item in UCAC4Table:
		UCAC4data.append([item['RA'], item['dec'], item['mag1'], item['mag2'], item['rmag'], item['imag']])
		UCACFO.write("%8.4f %8.4f\n" % (item['RA'], item['dec']))
	nStars = len(UCAC4data)
	UCACFO.close()

	## Filter UCAC catalog for target number of stars
	nStarsMax = 5000
	maglimit = 17.5
	if nStars > nStarsMax:
		while nStars > nStarsMax:
			maglimit = maglimit - 0.25
			UCAC4data = []
			if os.path.exists(UCACfile): os.remove(UCACfile)
			UCACFO = open(UCACfile, "w")
			for item in UCAC4Table:
				if item['mag1'] <= maglimit:
					UCAC4data.append([item['RA'], item['dec'], item['mag1'], item['mag2'], item['rmag'], item['imag']])
					UCACFO.write("%8.4f %8.4f\n" % (item['RA'], item['dec']))
			nStars = len(UCAC4data)
		UCACFO.close()
	logger.info("  Got %d UCAC catalog stars brighter than magnitude %.2f (mag1)" % (nStars, maglimit))
		
	## Run CCFIND
	logger.info("  Running CCFIND ...")
	CCFINDoutputFile = os.path.join(tmpPath, tel['FileID']+"_CCFINDout.txt")
	if os.path.exists(CCFINDoutputFile): os.remove(CCFINDoutputFile)
	CCFINDout = iraf.ccfind(UCACfile, CCFINDoutputFile, WorkingFile, Stdout=1, Stderr=IRAFerrors,
                            lngcol=1, latcol=2, lngunits="degrees", latunits="degrees", 
                            insystem="j2000", usewcs="yes", projection="tnx", center="yes")
	MatchCCFINDout = re.compile("Located ([0-9]{1,8}) objects in image")
	if len(CCFINDout) >= 8:
		IsCCFINDout = MatchCCFINDout.match(CCFINDout[8])
		if IsCCFINDout:
			nFound = int(IsCCFINDout.group(1))
			logger.info("  CCFIND:  Located %d objects in image." % nFound)
		else:
			nFound = -1
			logger.warning("  CCFIND failed.")

	## Run CCMAP
	logger.info("  Running CCMAP")
	DatabaseFile = os.path.join(tmpPath, tel['FileID']+"_db.txt")
	if os.path.exists(DatabaseFile): os.remove(DatabaseFile)
	
	CCMAP_STDOUT = iraf.ccmap(input=CCFINDoutputFile, database=DatabaseFile, Stdout=1, Stderr=IRAFerrors,
	                          xcol=3, ycol=4, lngcol=1, latcol=2, lngunits="degrees", latunits="degrees", 
	                          projection="tnx", fitgeometry="general", function="polynomial", xxo=4, xyo=4, yxo=4, yyo=4, xxterms="half", yxterms="half", 
	                          maxiter=3, reject=4.0,
	                          interactive="no", update="yes")
	logger.info("  CCMAP Solution:")
	for line in CCMAP_STDOUT:
		if re.match(".*Ra/Dec or Long/Lat fit rms.*", line):
			logger.info("    %s" % line)
		if re.match(".*X and Y scale.*", line):
			logger.info("    %s" % line)
		if re.match(".*X and Y axis rotation.*", line):
			logger.info("    %s" % line)

	logger.info("  Running CCSETWCS")
	CCSETWCS_STDOUT = iraf.ccsetwcs(images=WorkingFile, database=DatabaseFile, projection="tnx", Stdout=1, Stderr=1)
	
	## Clean Up
	IRAFerrors.close()
	if os.path.exists(UCACfile):
		logger.info("  Cleaning Up.  Deleting: %s", UCACfile)
		os.remove(UCACfile)
	if os.path.exists(CCFINDoutputFile):
		logger.info("  Cleaning Up.  Deleting: %s", CCFINDoutputFile)
		os.remove(CCFINDoutputFile)
	if os.path.exists(DatabaseFile):
		logger.info("  Cleaning Up.  Deleting: %s", DatabaseFile)
		os.remove(DatabaseFile)
	
	

##############################################################
## Fitting Routine for Zero Point
## - fits line with slope one to instrumental magnitude data
## - does sigma rejection to ignore extraneous points
##############################################################
def DetermineZeroPoint():
	pass



##############################################################
## Fitting Routine for Zero Point
## - fits line with slope one to instrumental magnitude data
## - does sigma rejection to ignore extraneous points
##############################################################
def func(x, a):
	return a + 1.0*x
def FitZPSigmaReject(xdata, ydata, RejectLevel, targetcov, maxiter, plot=False):
	x = numpy.array(xdata)
	y = numpy.array(ydata)
	popt, pcov = curve_fit(func, x, y)
	change = 1
	count = 0
	while (pcov[0][0] > targetcov) and (change >= 1) and (count < maxiter):
		count += 1
		nPointsStart = len(x)
		yerrs = abs(y - func(x, popt[0]))
		sigma = numpy.std(yerrs)
		x = x[(yerrs<RejectLevel*sigma)]
		y = y[(yerrs<RejectLevel*sigma)]
		popt, pcov = curve_fit(func, x, y)
		nPointsEnd = len(x)
		change = nPointsStart - nPointsEnd
	nRejected = len(xdata) - len(x)
	if plot:
		pyplot.figure()
		pyplot.plot(xdata, ydata, 'ro')
		pyplot.plot(x,y, 'bo')
		pyplot.savefig("plot.png", dpi=100)
	return popt[0], pcov[0][0], nRejected


##############################################################
## Determines mode(s) of an array
## - returns mode(s) given a binsize
##############################################################
def modes(array, binsize):
	values = []
	for element in array:
		#print element, element/binsize, round(element/binsize,0), round(element/binsize,0)*binsize
		values.append(round(element/binsize,0))
	count = defaultdict(int)
	for v in values:
		count[v] +=1
	best = max(count.values())
	return [k*binsize for k,v in count.items() if v == best]
