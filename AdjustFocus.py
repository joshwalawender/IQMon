#!/usr/bin/env python
# encoding: utf-8
"""
AdjustFocus.py

Created by Josh Walawender on 2012-11-26.
Copyright (c) 2012 . All rights reserved.
"""

import sys
import os
import re
import time
import subprocess
import shutil
import logging
import asciitable
import win32com.client

def main():
	Continue = True
	UseFocusMax = False
	telescope = "V5"
	
	now = time.gmtime()
	DecimalTime = now[3]+now[4]/60.+now[5]/3600.
	TimeString = time.strftime("%Y/%m/%d %H:%M:%SUT", now)
	DateString = time.strftime("%Y%m%dUT", now)
	
	LastMoveTime = TimeString
	
	if UseFocusMax:
		FocuserObject = win32com.client.Dispatch("FocusMax.Focuser")
	else:
		FocuserObject = win32com.client.Dispatch("FeatherTouch.Focuser")
	
	while Continue:
		time.sleep(5)
		if not FocuserObject.Link:
			FocuserObject.Link = True
		
		now = time.gmtime()
		DecimalTime = now[3]+now[4]/60.+now[5]/3600.
		TimeString = time.strftime("%Y/%m/%d %H:%M:%SUT", now)
		DateString = time.strftime("%Y%m%dUT", now)
		
		DateString = "20130107UT"
		
		IQMonFilePath   = os.path.join("C:\\", "Data_V5", "Logs", DateString, "FocusSummary_"+telescope+"_"+DateString+".txt")
		DOSFilePath     = os.path.join("C:\\", "Data_V5", "Logs", DateString, "DOS_"+telescope+"_"+DateString+".txt")
		SummaryFileName = os.path.join("C:\\", "Data_V5", "Logs", DateString, "FocusMoves_"+telescope+"_"+DateString+".txt")
		logging.basicConfig(filename=SummaryFileName, format='%(asctime)23s %(levelname)8s: %(message)s', level=logging.DEBUG)
		
		if os.path.exists(IQMonFilePath):
			shutil.copy(IQMonFilePath, DOSFilePath)
			dos2unix_result = subprocess.Popen(["C:\\bin\\unix2dos.exe", DOSFilePath], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			IQMonSummary = asciitable.read(DOSFilePath, comment="#")

			LastLineTime = IQMonSummary[-1][0]
			if LastLineTime != LastMoveTime:
				Move = int(IQMonSummary[-1][7])
				LastMoveTime = LastLineTime
				CurrentPosition = int(FocuserObject.Position)
				Destination = CurrentPosition + Move
				FocuserObject.Move(Destination)
				Moving = True
				while Moving:
					time.sleep(1)
					CurrentPos = int(FocuserObject.Position)
					if CurrentPos == Destination:
						Moving = False
						LogString = "Moved from %d to %d (%d steps) based on FITS image taken at %s\n" % (CurrentPosition, Destination, Move, LastLineTime)
						logging.info(LogString)
						

		## Stop Loop After 7am HST (17 UT)
		if DecimalTime >= 17.0:
			Continue = False
			


if __name__ == '__main__':
	main()

