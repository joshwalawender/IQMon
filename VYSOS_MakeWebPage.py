#!/usr/bin/env python
# encoding: utf-8
"""
VYSOS_LongTermLog.py

Created by Josh Walawender on 2013-05-01.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import getopt
import re
import shutil
import IQMonTools
import subprocess32
import datetime

def main(argv=None):
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
	
	
	
	##############################################################
	## Read Configuration File to get the following items
	## - IQMONEXECPATH
	## - IQMONLOGS
	## - IQMONPLOTS
	## - IQMONTMP
	IQMonExecPath, LogPath, PlotsPath, tmpPath, PythonPath, V5DataPath, V20DataPath, CatalogPath, LogBuffer = IQMonTools.ReadConfigFile()
	if telescope == "V20": telname = "VYSOS-20"
	if telescope == "V5":  telname = "VYSOS-5"
	NightSummariesDirectory = os.path.join(LogPath, telname)
	SummaryHTMLFile = os.path.join(NightSummariesDirectory, "index.html")
	TemporaryHTMLFile = os.path.join(NightSummariesDirectory, "index_tmp.html")



	##############################################################
	## Read Contents of Night Summaries Directory
	Files = os.listdir(NightSummariesDirectory)
	MatchDateOnFile  = re.compile("([0-9]{8}UT)_"+telescope+".*")
	
	## Make List of Dates for Files in Directory
	## - loop through files, extract date
	## - compare against list of dates already recorded
	## - if date not already recorded, add to list
	Dates = []
	for File in Files:
 		HasDate = MatchDateOnFile.match(File)
		if HasDate:
			Date = HasDate.group(1)
			DateAlreadyListed = False
			for ListedDate in Dates:
				if ListedDate[0] == Date:
					DateAlreadyListed = True
			if not DateAlreadyListed:
				Dates.append([Date, "", "", "", "", "", 0])

			IsPNGFile     = re.match(Date+"_"+telescope+"\.png", File)
			if IsPNGFile:
				print "Found Summary Graphs File for "+Date
				Dates[-1][1] = File
			IsEnvFile     = re.match(Date+"_"+telescope+"_Env\.png", File)
			if IsEnvFile:
				print "Found Environmantal Graphs File for "+Date
				Dates[-1][2] = File
			IsHTMLFile    = re.match(Date+"_"+telescope+"\.html", File)
			if IsHTMLFile:
				print "Found HTML File for "+Date
				Dates[-1][3] = File
			IsIQMonFile   = re.match(Date+"_"+telescope+"_IQMonLog\.txt", File)
			if IsIQMonFile:
				print "Found IQMonLog File for "+Date
				Dates[-1][4] = File
			IsSummaryFile = re.match(Date+"_"+telescope+"_Summary\.txt", File)
			if IsSummaryFile:
				print "Found Summary File for "+Date
				Dates[-1][5] = File
				wcSTDOUT = subprocess32.check_output(["wc", "-l", os.path.join(NightSummariesDirectory, File)], stderr=subprocess32.STDOUT, timeout=5)

				try:
					nLines = int(wcSTDOUT.strip().split(" ")[0])
					nImages = nLines - 1
				except:
					nImages = 0
				Dates[-1][6] = nImages
				
			
	# for Date in Dates:
	# 	print Date

	SortedDates = sorted(Dates, reverse=True)
	# for item in SortedDates:
	# 	print item

	##############################################################
	## Make index.html file
	HTML = open(SummaryHTMLFile, 'w')
	HTMLheader = open(os.path.join(IQMonExecPath, "VYSOS_ListOfNights.html"), 'r')
	header = HTMLheader.read()
	header = header.replace("telescopename", telname)
	HTMLheader.close()
	HTML.write(header)
	

	for DateInfo in SortedDates:
		HTML.write("    <tr>\n")
		## Write UT Date
		DateObject = datetime.datetime.strptime(DateInfo[0], "%Y%m%dUT")
		NiceDateString = DateObject.strftime("%Y/%m/%d UT")
		HTML.write("      <td style='text-align:center'>%-21s</td>\n" % NiceDateString)
		## Write Link to Night Summary Graphs
		if DateInfo[1] != "":
			HTML.write("      <td style='text-align:center'><a href='%s'>%-50s</a></td>\n" % (DateInfo[1], "Night Summary Graphs"))
		else:
			HTML.write("      <td style='text-align:center'></td>\n")
		## Write Link to Environmental Plots
		if DateInfo[2] != "":
			HTML.write("      <td style='text-align:center'><a href='%s'>%-50s</a></td>\n" % (DateInfo[2], "Enviornmental Graphs"))
		else:
			HTML.write("      <td style='text-align:center'></td>\n")
		## Write Link to HTML Summary
		if DateInfo[3] != "":
			HTML.write("      <td style='text-align:center'><a href='%s'>%-50s</a></td>\n" % (DateInfo[3], "Image Summary"))
		else:
			HTML.write("      <td style='text-align:center'></td>\n")
		## Write Link to IQMon Log
		if DateInfo[4] != "":
			HTML.write("      <td style='text-align:center'><a href='%s'>%-50s</a></td>\n" % (DateInfo[4], "IQMon Log"))
		else:
			HTML.write("      <td style='text-align:center'></td>\n")
		## Write Link to Text Summary
		if DateInfo[5] != "":
			HTML.write("      <td style='text-align:center'><a href='%s'>%-50s</a></td>\n" % (DateInfo[5], "Text Summary"))
		else:
			HTML.write("      <td style='text-align:center'></td>\n")
		## Write Number of Images
		HTML.write("      <td style='text-align:center'>%-5d</td>\n" % (DateInfo[6]))

	HTML.write("    </tr>\n")
	HTML.write("    </table>\n")
	HTML.write("</body>\n")
	HTML.write("</html>\n")
	HTML.close()


if __name__ == '__main__':
	main()

