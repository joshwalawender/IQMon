#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2012-12-05.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import subprocess
import time
import math
import re

def TestDevice(address, nPings):	
	MatchPingResult = re.compile(".*([0-9]+)\spackets\stransmitted,\s([0-9]+)\spackets received,\s([0-9\.]+).\spacket\sloss.*")
	MatchPingStats  = re.compile(".*round\-trip\smin/avg/max/stddev\s=\s([0-9\.]+)/([0-9\.]+)/([0-9\.]+)/([0-9\.]+)\sms.*")
	
	result = subprocess.check_output(["ping", "-c "+str(nPings), address])
	foo = result.find("statistics ---") + len("statistics ---")
	result = result[foo+1:-1]
	IsMatch = MatchPingResult.match(result)
	if IsMatch:
		PacketLoss = float(IsMatch.group(3))
		bar = result.find("packet loss") + len("packet loss")
		statstring = result[bar+1:]
		IsStats = MatchPingStats.match(statstring)
		if IsStats:
			AvgRT = float(IsStats.group(2))
		else:
			AvgRT = float("NaN")
	else:
		PacketLoss = float("NaN")
		AvgRT = float("NaN")
	if not math.isnan(PacketLoss):
		if PacketLoss <= 2./7.*100.:
			Status = "up"
		else:
			Status = "down"
	else:
		Status = "unknown"
	return Status, PacketLoss, AvgRT

def main():
	names     = ['Router', 'Switch', 'OldRouter', 'Olivier', 'Altair', 'CCTV', 'MLOAllSky']
	addresses = ['192.168.1.1', '192.168.1.2', '192.168.1.10', '192.168.1.50', '192.168.1.102', '192.168.1.103', '192.168.1.104']
	results   = ['', '', '', '', '', '', '']
	Addresses = dict(zip(names, addresses))
	PingResults = dict(zip(names, results))
	nPings = 7
	Continue = True

	HeaderString = "%-22s " % "## Time"
	for Device in names:
		HeaderString += "%-15s " % Device
	print HeaderString

	while Continue:
		## Start Only on Integer Multiples of 5 min		
		Waiting = True
		while Waiting:
			time.sleep(1)
			now = time.gmtime()
			DecimalTime = now[3]+now[4]/60.+now[5]/3600.
			TimeString = time.strftime("%Y/%m/%d %H:%M:%SUT", now)
			DateString = time.strftime("%Y%m%dUT", now)
			# print TimeString, now[4], now[5], now[4]%2
			if now[4]%5 == 0 and now[5] <= 1.0:
				Waiting = False

		## Form log file name (must be within loop to keep track of date)
		IQMonPath = os.path.expandvars("$IQMONPATH")
		LogFile = os.path.join(IQMonPath, "Logs", "Pings", "Pings_"+DateString+".txt")

		## Read in existing log data
		if os.path.exists(LogFile):
			Log = open(LogFile, 'r')
			PreviousLogs = Log.read()
			Log.close()
		## if previous log did not exist, make header string
		else:
			PreviousLogs = HeaderString+"\n"
			
		Log = open(LogFile, 'w')
		Log.write(PreviousLogs)
	
		## Loop through devices and get ping results
		for Device in names:
			Status, PacketLoss, AvgRT = TestDevice(Addresses[Device], nPings)
			ResultString = "%(status)s(%(percent).0f%(sym)s,%(time).2fms)" % {"status": Status, "percent": PacketLoss, "sym": "%", "time": AvgRT}
			PingResults[Device] = ResultString

		## Loop though devices and write those results to the log
		LogLine = "%-22s " % TimeString
		Log.write("%-22s " % TimeString)
		for Device in names:
			Log.write("%-15s " % PingResults[Device])
			LogLine += "%-15s " % PingResults[Device]
		Log.write("\n")
		Log.close()
		print LogLine
		
		time.sleep(2)		

if __name__ == '__main__':
	main()

