#!/usr/bin/env python
# encoding: utf-8
"""
IQMon.py

Created by Josh Walawender on 2013-05-31.
Copyright (c) 2013 . All rights reserved.
"""

import sys
import os
import re
import getopt
import time
import datetime
import subprocess32

import IQMonTools

help_message = '''
The help message goes here.
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):
	telescope = ""
	Operate = True
	
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

	##############################################################
	## Read Configuration File to get the following items
	IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, LogBuffer = IQMonTools.ReadConfigFile()
	PythonString = os.path.join(PythonPath, "python")

	## Set date to tonight
	now = time.gmtime()
	DateString = time.strftime("%Y%m%dUT", now)
	
	## Set data path
	if telescope == "":
		telescope = V5
	if telescope == "V5":
		DataPath = os.path.join(V5DataPath, "Images", DateString)
	if telescope == "V20":
		DataPath = os.path.join(V20DataPath, "Images", DateString)
	
	## Look for Pre-existing Files
	if not os.path.exists(DataPath): os.mkdir(DataPath)
	PreviousFiles = os.listdir(DataPath)
	PreviousFilesTime = time.gmtime()
	
	## Operations Loop
	while Operate:
		## Set date to tonight
		now = time.gmtime()
		nowDecimalHours = now.tm_hour + now.tm_min/60. + now.tm_sec/3600.
		DateString = time.strftime("%Y%m%dUT", now)
		TimeString = time.strftime("%Y/%m/%d %H:%M:%S UT -", now)
		
		Files = os.listdir(DataPath)
		FilesTime = now
		
		time.sleep(1)
				
		if len(Files) > len(PreviousFiles):
			for File in Files:
				FileFound = False
				for PreviousFile in PreviousFiles:
					if File == PreviousFile:
						FileFound = True
				if not FileFound:
					if re.match(".*\.fi?ts", File) and not re.match(".*\-Empty\-.*\.fts", File):
						print "New fits File Found:  %s" % File
						Focus = False
						ProcessCall = [PythonString, IQMonExecPath+"/MeasureImage.py", "--input="+os.path.join(DataPath, File)]
						print "  %s Calling MeasureImage.py with %s" % (TimeString, ProcessCall[2:])
						try:
							MIoutput = subprocess32.check_output(ProcessCall, stderr=subprocess32.STDOUT, timeout=150)
							print "Call to MeasureImage.py Succeeded"
						except:
							print "Call to MeasureImage.py Failed"
		PreviousFiles = Files
		PreviousFilesTime = now
		time.sleep(5)
		if nowDecimalHours > 18.0:
			Operate = False
		
		

if __name__ == "__main__":
	sys.exit(main())
