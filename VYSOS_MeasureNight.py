#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2012-10-29.
Copyright (c) 2012 . All rights reserved.
"""

import sys
import getopt
import os
import subprocess32
import re
import fnmatch
import numpy
import time

import IQMonTools

help_message = '''
The help message goes here.
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):
	telescope = ""
	Focus = False
	Clobber = True
	cosmicrays = False
	
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hfd:n:t:", ["help", "focus", "no-clobber", "date=", "night=", "telescope=", "cosmicrays"])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-f", "--focus"):
				Focus = True
			if option in ("-d", "--date", "-n", "--night"):
				input = value
			else:
				input = ""
			if option in ("-t", "--telescope"):
				telescope = value
			if option in ("--no-clobber"):
				Clobber = False
			if option in ("--cosmicrays"):
				cosmicrays = True
			
	
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2
	
	##############################################################
	## Read Configuration File to get the following items
	IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, LogBuffer = IQMonTools.ReadConfigFile()
	
	## Set date to tonight if not specified
	now = time.gmtime()
	DateString = time.strftime("%Y%m%dUT", now)
	if (input == "tonight") or (input == "Tonight") or (input == "lastnight") or (input == "LastNight") or (input ==""):
		input = DateString
	
	## Set Default Telescope to V5 if not set
	if telescope == "V20" or telescope == "VYSOS20" or telescope == "VYSOS-20":
		telescope = "V20"
	if telescope == "V5" or telescope == "VYSOS5" or telescope == "VYSOS-5":
		telescope = "V5"
	if telescope != "V20" and telescope != "V5":
		telescope = "V5"
	
		
	###########################
	## 
	if telescope == "V5":
		VYSOSDATAPath = V5DataPath
	if telescope == "V20":
		VYSOSDATAPath = V20DataPath
	ImagesDirectory = os.path.join(VYSOSDATAPath, "Images", input)
	LogsDirectory = os.path.join(VYSOSDATAPath, "Logs", input)
	PythonString = os.path.join(PythonPath, "python")
	
	print "Analyzing data for night of "+input
	if os.path.exists(ImagesDirectory) and os.path.exists(LogsDirectory):
		print "  Found "+ImagesDirectory+" and "+LogsDirectory
		##
		## Loop Through All Images in Images Directory
		##
		Files = os.listdir(ImagesDirectory)
		print "Found %d files in images directory" % len(Files)
		if len(Files) >= 1:
			## Parse filename for date and time
			MatchFilename = re.compile("(.*)\-([0-9]{8})at([0-9]{6})\.fts")
			MatchEmpty = re.compile(".*\-Empty\-.*\.fts")
			Properties = []
			for File in Files:
				IsMatch = MatchFilename.match(File)
				IsEmpty = MatchEmpty.match(File)
				if IsMatch and not IsEmpty:
					target = IsMatch.group(1)
					FNdate = IsMatch.group(2)
					FNtime = IsMatch.group(3)
					Properties.append([FNtime, FNdate, target, File])
				else:
					print "  File Rejected: %s" % File
		
			SortedImageTimes   = numpy.array([row[0] for row in sorted(Properties)])
			SortedImageDates   = numpy.array([row[1] for row in sorted(Properties)])
			SortedImageTargets = numpy.array([row[2] for row in sorted(Properties)])
			SortedImageFiles   = numpy.array([row[3] for row in sorted(Properties)])
		
			print "%d out of %d files meet selection criteria." % (len(SortedImageFiles), len(Files))
			for Image in SortedImageFiles:
				if fnmatch.fnmatch(Image, "*.fts"):
					now = time.gmtime()
					TimeString = time.strftime("%Y/%m/%d %H:%M:%S UT -", now)
					DateString = time.strftime("%Y%m%dUT", now)

					# ProcessCall = [PythonString, IQMonExecPath+"/MeasureImage.py", os.path.join(ImagesDirectory, Image)]
					ProcessCall = ["MeasureImage.py"]
					if Clobber and Image == SortedImageFiles[0]:
						ProcessCall.append("--clobber")
					if cosmicrays:
						ProcessCall.append("--cosmicrays")
					ProcessCall.append(os.path.join(ImagesDirectory, Image))
					print "  %s Calling MeasureImage.py with %s" % (TimeString, ProcessCall)
					
					try:
						MIoutput = subprocess32.check_output(ProcessCall, stderr=subprocess32.STDOUT, timeout=150)
						for line in MIoutput.split("\n"):
							print line
					except:
						print "Call to MeasureImage.py Failed"
		else:
			print "No image files found in directory: "+ImagesDirectory
	else:
		print "No Images or Logs directory for this night"
	if not Focus:
		# print "  Making Nightly Plots for "+input
		# MNPcommand = [PythonString, IQMonExecPath+"/VYSOS_MakeNightlyPlots.py", "--date="+input, "--telescope="+telescope]
		# MNPoutput = subprocess32.check_output(MNPcommand, stderr=subprocess32.STDOUT, timeout=600)
		print "  Updating index.html web pages with new nights"
		MWPcommand = [PythonString, IQMonExecPath+"/VYSOS_MakeWebPage.py", "--telescope="+telescope]
		MWPoutput = subprocess32.check_output(MWPcommand, stderr=subprocess32.STDOUT, timeout=600)

if __name__ == "__main__":
	sys.exit(main())
