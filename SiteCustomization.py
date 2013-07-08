import os
import re
import logging
import datetime
import math
import ephem
from pyraf import iraf

def GetTelescopeProperties(FitsFileDirectory, FitsFilename, IQMonPath):
	## Determine Data Night from Folder Structure
	## - first, check to see if containing folder matches DataNight Pattern
	DataNightPattern = re.compile("([0-9]{4}[0-9]{2}[0-9]{2})U?T?")
	ISDNP1 = DataNightPattern.match(os.path.split(os.path.abspath(FitsFileDirectory))[1])
	if ISDNP1:
		DataNight = ISDNP1.group(1)+"UT"
	else:
		ISDNP2 = DataNightPattern.match(os.path.split(os.path.split(os.path.abspath(FitsFileDirectory))[0])[1])
		if ISDNP2:
			DataNight = ISDNP2.group(1)+"UT"
		else:
			print "  Could not Determine DataNight"
		
	## Match Filename
	MatchVYSOSFilename = re.compile("(V[0-9]{1,2})_([0-9a-zA-Z_\-]+)\-([0-9]{8})at([0-9]{6})\.fts")
	IsVYSOSFilename = MatchVYSOSFilename.match(FitsFilename)
	if IsVYSOSFilename:
		telescope    = IsVYSOSFilename.group(1)
		FileNameDate = IsVYSOSFilename.group(3)
		FileNameTime = IsVYSOSFilename.group(4)
		FileID       = telescope+"_"+FileNameDate+"at"+FileNameTime
	else:
		MatchOldVYSOSFilename = re.compile("(V[0-9]{1,2})_(.*)\.fts")
		IsOldVYSOSFilename = MatchOldVYSOSFilename.match(FitsFilename)
		if IsOldVYSOSFilename:
			telescope = IsOldVYSOSFilename.group(1)
			FileID       = FitsFilename.split(".")[0]
		else:
			telescope = "other"

	## Define Properties
	tel = dict()
	if telescope == "V20":
		## Telescope Properties
		tel['name'] = 'V20'
		tel['LongName'] = "VYSOS-20"
		tel['FocalLength'] = 4175.  # (mm)
		tel['PixelSize'] = 9.0      # (um)
		tel['Aperture'] = 20.0*25.4 # (mm)
		tel['Gain'] = 1.6           # (electrons per ADU)
		tel['FullXPix'] = 4096      # Full size of the CCD in pixels
		tel['FullYPix'] = 4096      # Full size of the CCD in pixels
		## Focus Algorithm Properties
		tel['StepSize'] = 100.      # (um) Focus distance per step
		tel['TargetFWHM'] = 1.5     # in units of tel['FWHMunits'] below
		tel['Aggressiveness'] = 0.5 # fraction of calculated focus move to actually execute
		tel['MaxMove'] = 10         # Maximum move in steps to make regardless of what focus algorithm calculates
		## SExtractor Properties
		tel['PhotAperture'] = 16.0  # (pixels)
		tel['Seeing'] = 2.0         # (arcsec)
		## IQMon Options and Properties
		tel['FWHMunits'] = "arcsec" # choice for reported units
		tel['IQRadiusFactor'] = 1.0 # take all of CCD for FWHM?
		tel['UseROI'] = True        # use only portion of image for IQMon analysis
		tel['ROI'] = "[1025:3072,1025:3072]"
		## Calculated Properties
		tel['PixelScale'] = 206.265*tel['PixelSize']/tel['FocalLength']
		tel['FRatio'] = tel['FocalLength']/tel['Aperture']
		tel['FileID'] = FileID
		tel['DataNight'] = DataNight
		## User Chosen Limits for Image Flagging
		tel['IQLimit'] = 3.0/tel['PixelScale'] # (pix) FWHM Limit before image is flagged
		tel['PELimit'] = 3.0        # (arcmin) Pointing error limit before image is flagged
		tel['EllipLimit'] = 0.35    # Ellipticity limit before image is flagged
		
	if telescope == "V5":
		## Telescope Properties
		tel['name'] = 'V5'
		tel['LongName'] = "VYSOS-5"
		tel['FocalLength'] = 735.   # (mm)
		tel['PixelSize'] = 9.0      # (um)
		tel['Aperture'] = 135.      # (mm)
		tel['Gain'] = 1.6           # (electrons per ADU)
		tel['FullXPix'] = 4096      # Full size of the CCD in pixels
		tel['FullYPix'] = 4096      # Full size of the CCD in pixels
		## Focus Algorithm Properties
		tel['StepSize'] = 7.28      # (um) Focus distance per step
		tel['TargetFWHM'] = 1.25    # in units of tel['FWHMunits'] below
		tel['Aggressiveness'] = 0.5 # fraction of calculated focus move to actually execute
		tel['MaxMove'] = 10         # Maximum move in steps to make regardless of what focus algorithm calculates
		## SExtractor Properties
		tel['PhotAperture'] = 6.0   # (pixels)
		tel['Seeing'] = 2.0         # (arcsec)
		## IQMon Options and Properties
		tel['FWHMunits'] = "pixels" # choice for reported units
		tel['IQRadiusFactor'] = 1.0 # take all of CCD for FWHM?
		tel['UseROI'] = True        # use only portion of image for IQMon analysis
		tel['ROI'] = "[1025:3072,1025:3072]"
		## Calculated Properties
		tel['PixelScale'] = 206.265*tel['PixelSize']/tel['FocalLength']
		tel['FRatio'] = tel['FocalLength']/tel['Aperture']
		tel['FileID'] = FileID
		tel['DataNight'] = DataNight
		## User Chosen Limits for Image Flagging
		tel['IQLimit'] = 2.5        # (pix) FWHM Limit before image is flagged
		tel['PELimit'] = 5.0        # (arcmin) Pointing error limit before image is flagged
		tel['EllipLimit'] = 0.35    # Ellipticity limit before image is flagged
		## Distortion Coefficients for WCS Refinement
		## - Uses IRAF TNX format for distortions
		tel['WCSDIM']   = 2
		tel['CTYPE1']   = 'RA---TNX'
		tel['CTYPE2']   = 'DEC--TNX'
		tel['WAT0_001'] = 'system=image'
		tel['WAT1_001'] = 'wtype=tnx axtype=ra lngcor = "3. 4. 4. 2. -1.398949147760993 1.46139'
		tel['WAT1_002'] = '5418573523 -1.362189763246309 1.521748524766638 -6.232240670682045E-'
		tel['WAT1_003'] = '4 0.01330577077620412 0.001084712871350142 -0.007314776300857944 1.5'
		tel['WAT1_004'] = '28463887188006E-4 0.001066055285100885 9.863395632396422E-6 3.592239'
		tel['WAT1_005'] = '621696304E-4 -0.007349566338775161 9.188017571342054E-7 "'
		tel['WAT2_001'] = 'wtype=tnx axtype=dec latcor = "3. 4. 4. 2. -1.398949147760993 1.4613'
		tel['WAT2_002'] = '95418573523 -1.362189763246309 1.521748524766638 -7.964375522791034E'
		tel['WAT2_003'] = '-4 1.831471606780815E-4 5.096231979490976E-4 -1.222779983305424E-6 0'
		tel['WAT2_004'] = '.01341963536388697 7.049581472363932E-4 -0.007373876338113699 0.0015'
		tel['WAT2_005'] = '06907175257101 -1.397714361455436E-5 -0.007296656826212161 "'
		
		
	return tel



##############################################################
## Obtain Flat Field Frame
## - Look for flats in /Volumes/Data/VYSOS-5/calFrames/YYYYMMDD
## - Find flat files, determine which are for this filter
##############################################################
def ObtainMasterFlat(telescope, CCDfilter, tel, logger):
	MasterFlatFilename = "Flat_"+telescope+"_"+CCDfilter+"_"+DataNightString+".fits"
	FlatFileBaseDirectory = os.path.join("/Volumes", "Data", "VYSOS-5", "calFrames")
	FoundFlats = False
	DataNightPattern = re.compile("([0-9]{4})([0-9]{2})([0-9]){2}UT")
	DataNight = datetime.datetime.strptime(DataNightString, "%Y%m%dUT")
	OneDay = datetime.timedelta(days=1)
	NewDate = DataNight
	NewDateString = datetime.datetime.strftime(NewDate, "%Y%m%dUT")
	## Check to see if a flat file exists for this data night
	if os.path.exists(MasterFlatFilename):
		logger.info("  Found master flat file: %s" % ObtainMasterFlat)
	else:
		## Look for individual files which can be made in to a master flat
		pass
		
	return MasterFlatFilename
		
	

##############################################################
## Obtain Master Dark Frame
## - look for raw dark frames from previous night
## - if not found, look back in time for MaxNights
## - once a qualifying set of darks is found (same telescope, same exptime)
## - check that CCD temp is within tolerance for darks and image
## - if not same temp, then continue going back in time
## - once a set of dark frames which satify the criterion is found, combine all via a median combine
## - write that combined dark to a file named for the DataNight (not the night it was taken)
##############################################################
def ObtainMasterDark(telescope, DataNightString, exptime, DarkPath, DataPath, logger):
	logger.info("  Looking for master dark frame.")
	DataNightPattern = re.compile("([0-9]{4})([0-9]{2})([0-9]){2}UT")
	DarkFilePattern = re.compile("Dark\-([0-9]{3})\-([0-9]{8})at([0-9]{6})\.fts")
	FoundDarks = False
	DataNight = datetime.datetime.strptime(DataNightString, "%Y%m%dUT")
	OneDay = datetime.timedelta(days=1)
	NewDate = DataNight - OneDay
	NewDateString = datetime.datetime.strftime(NewDate, "%Y%m%dUT")
	DateLimit = DataNight - datetime.timedelta(days=20)
	## Check to see if MasterDark Exists for this Date
	MasterDarkFilename = "MasterDark_"+telescope+"_"+DataNightString+"_"+str(int(math.floor(exptime)))+".fits"
	MasterDarkFile  = os.path.join(DarkPath, MasterDarkFilename)	
	## Is that Master Dark File does not exist, see if the raw files exit to build one.
	if os.path.exists(MasterDarkFile):
		logger.info("  Found Master Dark: %s" % MasterDarkFilename)
	else:
		logger.info("  Could Not Find Master Dark.  Looking for raw frames.")
		BaseDirectory, NightDirectory = os.path.split(DataPath)
		while FoundDarks == False and NewDate > DateLimit:
			SearchPath = os.path.join(BaseDirectory, NewDateString, "Calibration")
			if os.path.exists(SearchPath):
				Files = os.listdir(SearchPath)
				Darks = []
				RawDarks = []
				for File in Files:
					IsDark = DarkFilePattern.match(File)
					if IsDark:
						DarkExp = float(IsDark.group(1))
						if DarkExp == exptime:
							RawDarks.append(os.path.join(SearchPath, File))
							Darks.append(os.path.join(SearchPath, File.replace("fts", "fits")))
				if len(Darks) >= 3:
					## Make Master Dark
					logger.info("  Found %d dark files in %s" % (len(Darks), SearchPath))
					logger.info("  Converting %d raw dark frames using RFITS" % len(Darks))
					for i in range(0,len(RawDarks),1):
						if os.path.exists(Darks[i]): os.remove(Darks[i])
						RFITSoutput = iraf.rfits(RawDarks[i], 0, Darks[i], datatype="u", Stdout=1)
					iraf.imcombine(input=",".join(Darks), output=MasterDarkFile, combine="median", Stdout=1, Stderr=1)
					# try:
					# 	IMCOMBINEout = iraf.imcombine(input=Darks, output=MasterDarkFile, combine="median", Stdout=1, Stderr=1)
					# 	print "IMCOMBINEout:"
					# 	print IMCOMBINEout
					# except:
					# 	logger.warning("  IMCOMBINE failed.")
					for Dark in Darks:
						os.remove(Dark)
					if os.path.exists(MasterDarkFile):
						logger.info("  Combined %d dark files to make %s" % ((len(Darks), MasterDarkFilename)))
					else:
						logger.warning("  Creation of Master Dark failed.")
					FoundDarks = True
			NewDate = NewDate - OneDay
			NewDateString = datetime.datetime.strftime(NewDate, "%Y%m%dUT")

	if os.path.exists(MasterDarkFile):
		return MasterDarkFile
	else:
		return ""

