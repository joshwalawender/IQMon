#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2012-11-27.
Copyright (c) 2012 . All rights reserved.
"""

import sys
import os
import string
import asciitable
from matplotlib import pyplot
import numpy

def TimeStringToDecimal(TimeString):
	hms = string.split(TimeString, ":")
	DecimalTime = float(hms[0])+ float(hms[1])/60.0 + float(hms[2])/60.0/60.0
	return DecimalTime

def main():
	nights = ["20121122UT", "20121123UT", "20121124UT", "20121125UT", "20121126UT", 
	          "20121127UT", "20121128UT", "20121129UT", "20121130UT", "20121201UT",
	          "20121202UT", "20121203UT", "20121204UT", "20121205UT"]
	nights = ["20121204UT", "20121205UT"]

	DiffTimes = []
	Differences = []
	V20AvgSD = []
	V5AvgSD = []
	V20AvgH = []
	V5AvgH = []
	V20AvgOT = []
	V5AvgOT = []
	V20AvgW = []
	V5AvgW = []

	nightcount = 0
	for night in nights:
		nightcount += 1
		print "Analyzing "+night+" ("+str(nightcount)+")"
		V20TimeDecimal = []
		V5TimeDecimal = []
		
		## Read in VYSOS-20 Logs
		LogFileName = os.path.join("/Volumes", "Data", "VYSOS-20", "Logs", night, "EnvironmentalLog.txt")
		V20Night = asciitable.read(LogFileName)
		V20Time      = [item[1] for item in V20Night]
		V20SkyT      = [item[7] for item in V20Night]
		V20OutsideT  = [item[8] for item in V20Night]
		V20WindSpeed = [item[9] for item in V20Night]
		V20Humidity  = [item[10] for item in V20Night]
		V20SkyDiff   = [item[7]-item[8] for item in V20Night]
		for time in V20Time:
			V20TimeDecimal.append(TimeStringToDecimal(time[0:8]))
		
		## Read in VYSOS-5 Logs
		LogFileName = os.path.join("/Volumes", "Data", "VYSOS-5", "Logs", night, "EnvironmentalLog.txt")
		V5Night = asciitable.read(LogFileName)
		V5Time      = [item[1] for item in V5Night]
		V5SkyT      = [item[4] for item in V5Night]
		V5OutsideT  = [item[5] for item in V5Night]
		V5WindSpeed = [item[6] for item in V5Night]
		V5Humidity  = [item[7] for item in V5Night]
		V5SkyDiff   = [item[4]-item[5] for item in V5Night]
		for time in V5Time:
			V5TimeDecimal.append(TimeStringToDecimal(time[0:8]))

		V20TimeDecimal = numpy.array(V20TimeDecimal)
		V20SkyDiff = numpy.array(V20SkyDiff)
		V20OutsideT = numpy.array(V20OutsideT)
		V20Humidity = numpy.array(V20Humidity)
		V20WindSpeed = numpy.array(V20WindSpeed)
		V5TimeDecimal = numpy.array(V5TimeDecimal)
		V5SkyDiff = numpy.array(V5SkyDiff)
		V5OutsideT = numpy.array(V5OutsideT)
		V5Humidity = numpy.array(V5Humidity)
		V5WindSpeed = numpy.array(V5WindSpeed)
		
		for interval in range(1,145):
			hour = interval/6.
			AvgV20SD = numpy.median(V20SkyDiff[(V20TimeDecimal>hour-1./6.) & (V20TimeDecimal<=hour)])
			AvgV5SD  = numpy.median(V5SkyDiff[(V5TimeDecimal>hour-1./6.) & (V5TimeDecimal<=hour)])
			AvgV20OT = numpy.median(V20OutsideT[(V20TimeDecimal>hour-1./6.) & (V20TimeDecimal<=hour)])
			AvgV5OT  = numpy.median(V5OutsideT[(V5TimeDecimal>hour-1./6.) & (V5TimeDecimal<=hour)])
			AvgV20H  = numpy.median(V20Humidity[(V20TimeDecimal>hour-1./6.) & (V20TimeDecimal<=hour)])
			AvgV5H   = numpy.median(V5Humidity[(V5TimeDecimal>hour-1./6.) & (V5TimeDecimal<=hour)])
			AvgV20W  = numpy.median(V20WindSpeed[(V20TimeDecimal>hour-1./6.) & (V20TimeDecimal<=hour)])
			AvgV5W   = numpy.median(V5WindSpeed[(V5TimeDecimal>hour-1./6.) & (V5TimeDecimal<=hour)])

			V20AvgSD.append(AvgV20SD)
			V5AvgSD.append(AvgV5SD)
			V20AvgOT.append(AvgV20OT)
			V5AvgOT.append(AvgV5OT)
			V20AvgH.append(AvgV20H)
			V5AvgH.append(AvgV5H)
			V20AvgW.append(AvgV20W)
			V5AvgW.append(AvgV5W)
			DiffTimes.append(hour + 24.*nightcount)
			Differences.append(AvgV20SD - AvgV5SD)



		pyplot.ioff()
		
		
		## Make Plot
		pyplot.figure(night+"OT")
		pyplot.title("Outside Temperature on "+night)
		pyplot.plot(V20TimeDecimal, V20OutsideT, 'b-', label="V20 Outside Temp")
		pyplot.plot(V5TimeDecimal, V5OutsideT, 'g-', label="V5 Outside Temp")
		pyplot.legend(loc="best")
		pyplot.xlabel("Time in Hours UT")
		pyplot.ylabel("Temperature (F)")
		pyplot.xlim(0,24)
		pyplot.savefig("Plots/OutsideTemp_"+night+".png", dpi=100)

		pyplot.figure(night+"ST")
		pyplot.title("Sky Temperature on "+night)
		pyplot.plot(V20TimeDecimal, V20SkyT, 'b.', label="V20 Sky Temp")
		pyplot.plot(V5TimeDecimal, V5SkyT, 'g.', label="V5 Sky Temp")
		pyplot.legend(loc="best")
		pyplot.xlabel("Time in Hours UT")
		pyplot.ylabel("Temperature (F)")
		pyplot.ylim(-50,50)
		pyplot.xlim(0,24)
		pyplot.savefig("Plots/SkyTemp_"+night+".png", dpi=100)

		pyplot.figure(night+"SD")
		pyplot.title("SkyTemp - AmbientTemp on "+night)
		pyplot.plot(V20TimeDecimal, V20SkyDiff, 'b.', label="V20 Sky Diff")
		pyplot.plot(V5TimeDecimal, V5SkyDiff, 'g.', label="V5 Sky Diff")
		pyplot.legend(loc="best")
		pyplot.xlabel("Time in Hours UT")
		pyplot.ylabel("Temperature Difference (F)")
		pyplot.ylim(-140,10)
		pyplot.grid()
		pyplot.xlim(0,24)
		pyplot.savefig("Plots/SkyDiff_"+night+".png", dpi=100)

		pyplot.figure(night+"H")
		pyplot.title("Humidity on "+night)
		pyplot.plot(V20TimeDecimal, V20Humidity, 'b-', label="V20 Humidity")
		pyplot.plot(V5TimeDecimal, V5Humidity, 'g-', label="V5 Humidity")
		pyplot.legend(loc="best")
		pyplot.xlabel("Time in Hours UT")
		pyplot.ylabel("Humidity (%)")
		pyplot.ylim(0,100)
		pyplot.xlim(0,24)
		pyplot.savefig("Plots/Humidity_"+night+".png", dpi=100)

	pyplot.figure("avgs")
	pyplot.title("SkyDiff Difference V20 - V5")
	pyplot.plot(DiffTimes, Differences, 'bo')
	pyplot.xlabel("Time in Hours")
	pyplot.ylabel("Temperature Difference (F)")
	pyplot.ylim(-40,10)
	pyplot.grid()
	pyplot.savefig("Plots/Difference.png", dpi=100)

	pyplot.figure("corrSD")
	pyplot.title("Correlation of Sky Difference")
	pyplot.plot(V20AvgSD, V5AvgSD, 'bo')
	pyplot.xlabel("VYSOS-20 Sky Difference Temp. (F)")
	pyplot.ylabel("VYSOS-5 Sky Difference Temp. (F)")
	pyplot.ylim(-100,0)
	pyplot.xlim(-100,0)
	pyplot.grid()
	pyplot.savefig("Plots/Correlation_SD.png", dpi=100)

	pyplot.figure("corrOT")
	pyplot.title("Correlation of Outside Temp")
	pyplot.plot(V20AvgOT, V5AvgOT, 'bo')
	pyplot.xlabel("VYSOS-20 Outside Temp. (F)")
	pyplot.ylabel("VYSOS-5 Outside Temp. (F)")
	pyplot.ylim(30,70)
	pyplot.xlim(30,70)
	pyplot.grid()
	pyplot.savefig("Plots/Correlation_OT.png", dpi=100)

	pyplot.figure("corrH")
	pyplot.title("Correlation of Humidity")
	pyplot.plot(V20AvgH, V5AvgH, 'bo')
	pyplot.xlabel("VYSOS-20 Humidity (%)")
	pyplot.ylabel("VYSOS-5 Humidity (%)")
	pyplot.ylim(0,100)
	pyplot.xlim(0,100)
	pyplot.grid()
	pyplot.savefig("Plots/Correlation_H.png", dpi=100)

	pyplot.figure("corrW")
	pyplot.title("Correlation of Wind Speed")
	pyplot.plot(V20AvgW, V5AvgW, 'bo')
	pyplot.xlabel("VYSOS-20 Wind Speed (kph)")
	pyplot.ylabel("VYSOS-5 Wind Speed (kph)")
	# pyplot.ylim(0,100)
	# pyplot.xlim(0,100)
	pyplot.grid()
	pyplot.savefig("Plots/Correlation_W.png", dpi=100)

if __name__ == '__main__':
	main()

