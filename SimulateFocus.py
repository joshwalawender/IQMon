#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2012-11-17.
Copyright (c) 2012 . All rights reserved.
"""

import sys
import os
import math
import numpy
import time
import random
import matplotlib.pyplot as pyplot

import Focus

def SimulateImageResult(index, StartingError, TargetFWHM, Aggressiveness):
	# print "index = %d" % index
	telescope = "V5"
	if telescope == "V5":
		FRatio = 5.4
		StepSize = 7.28
		PixelSize = 9.0
	
	error_sigma = 0.10
	MaxMove = 1.0*PixelSize*FRatio/StepSize  #MaxMove is steps to change FWHM by 1 pix
	
	now = time.gmtime()
	TimeString = time.strftime("%Y/%m/%d %H:%M:%S", now)
	DateString = time.strftime("%Y%m%dUT", now)
	
	StartingTemp = 10.0
	RealIQ = 1.2
	DeltaT = 15.0
	Temp = StartingTemp - float(index**2)/200**2.*DeltaT
	# print "Temp = %.4f" % Temp
	
	##
	## Results Without Focusing
	##
	DeltaL = 735.*24e-6*(Temp-StartingTemp)*1000.0 + StartingError*StepSize
	GeoBlur = abs(DeltaL)/5.4/9.0
	IQ_SEx = math.sqrt(RealIQ**2 + GeoBlur**2) + random.gauss(0,error_sigma)

	SummaryFileName = os.path.join("tmp", "NoFocus.txt")
	if os.path.exists(SummaryFileName):
		file = open(SummaryFileName, 'r')
		PreviousLogs = file.read()
		file.close()
	else:
		PreviousLogs = ""
	file = open(SummaryFileName, 'w')
	file.write(PreviousLogs)
	file.write("%-20s " % TimeString)
	file.write("%-49s " % "FitsFilename")
	file.write("%-12.2f\n" % IQ_SEx)
	file.close()

	##
	## Reactive Focus Algorithm
	##
	SummaryFileName = os.path.join("tmp", "WithFocus.txt")
	
	Move = Focus.Reactive(SummaryFileName, TargetFWHM, Aggressiveness, MaxMove, PixelSize, FRatio, StepSize)

	## Determine Total Moves to This Point for Simulation Only
	History = []
	Moves = []
	MoveSum = 0.0
	if os.path.exists(SummaryFileName):
		SumFile = open(SummaryFileName, 'r')
		for line in SumFile:
			if line[0] != "#":
				FWHM  = float(line[71:82])
				Moved = float(line[84:95])
				History.append(FWHM)
				Moves.append(Moved)
	MoveSum = numpy.sum(Moves)+Move
	
	GeoBlur = abs(DeltaL+MoveSum*StepSize)/5.4/9.0
	IQ_SEx = math.sqrt(RealIQ**2 + GeoBlur**2) + random.gauss(0,error_sigma)
	
	# if index == 0:
	# 	print "%4s %10s %10s %10s %10s %10s %10s %10s" % ("i", "S2LastFWHM", "LastFWHM", "DeltaL", "MS*StepSz", "GeoBlur", "Move", "IQ")
	# if index < 2:
	# 	LastFWHM = 0.0
	# 	S2LastFWHM = 0.0
	# print "%4d %10.2f %10.2f %10.3f %10.2f %10.2f %10d %10.2f" % (index, S2LastFWHM, LastFWHM, DeltaL, MoveSum*StepSize, GeoBlur, Move, IQ_SEx)	

	if os.path.exists(SummaryFileName):
		file = open(SummaryFileName, 'r')
		PreviousLogs = file.read()
		file.close()
	else:
		PreviousLogs = ""
	file = open(SummaryFileName, 'w')
	file.write(PreviousLogs)
	file.write("%-20s " % TimeString)
	file.write("%-49s " % "FitsFilename")
	file.write("%-12.2f " % IQ_SEx)
	file.write("%-12.2f\n" % Move)
	file.close()


def SimulateNights(StartingError, TargetFWHM, Aggressiveness, nNights, Plot):
	MaxAchieved = []
	MaxWorstCase = []
	MeanAchieved = []
	MeanWorstCase = []
	Scatter = []

	for j in range(0,nNights):
		if os.path.exists(os.path.join("tmp", "NoFocus.txt")): 
			os.remove(os.path.join("tmp", "NoFocus.txt"))
		if os.path.exists(os.path.join("tmp", "WithFocus.txt")): 
			os.remove(os.path.join("tmp", "WithFocus.txt"))
		for i in range(0,201):
			SimulateImageResult(i, StartingError, TargetFWHM, Aggressiveness)
		## Plot Results
		NoFocusFile = open(os.path.join("tmp", "NoFocus.txt"), 'r')
		NoFocusFWHM = []
		for line in NoFocusFile:
			if line[0] != "#":
				DateTime = line[0:20]
				File     = line[21:70]
				# print line
				# print "=>"+line[78:89]+"<="
				# print "=>"+line[91:102]+"<="
				FWHM  = float(line[71:82])
				NoFocusFWHM.append(FWHM)
		WithFocusFile = open(os.path.join("tmp", "WithFocus.txt"), 'r')
		WithFocusFWHM = []
		for line in WithFocusFile:
			if line[0] != "#":
				DateTime = line[0:20]
				File     = line[21:70]
				# print line
				# print "=>"+line[78:89]+"<="
				# print "=>"+line[91:102]+"<="
				FWHM  = float(line[71:82])
				WithFocusFWHM.append(FWHM)
		
		## Make Scatter Plot
		History = numpy.array(WithFocusFWHM)
		RunningDiff = []
		for i in range(3,len(History)):
			PreviousAvg = numpy.mean(History[i-3:i-1])
			CurrentValue = History[i]
			RunningDiff.append(CurrentValue - PreviousAvg)
		## Make Plots
		if Plot:
			PlotFile = "FWHM_T%02d_A%02d_SE%02d.png" % (int(TargetFWHM*10.0), int(Aggressiveness*10.0), int(StartingError))
			pyplot.ioff()
			pyplot.figure()
			TitleString = "StartingError = %d, TargetFWHM = %.1f, Aggressiveness = %.1f" % (StartingError, TargetFWHM, Aggressiveness)
			pyplot.title(TitleString)
			pyplot.plot(NoFocusFWHM, color="red", label="No Focus", marker="o")
			pyplot.plot(WithFocusFWHM, color="blue", label="Reactive Focus", marker="o")
			pyplot.plot(RunningDiff, color="black", label="FWHM Scatter", marker="o")
			pyplot.legend(loc='lower left')
			pyplot.xlabel("Image Number")
			pyplot.ylabel("FWHM (pix)")
			pyplot.ylim(-2,3)
			# pyplot.ylim(min([min(NoFocusFWHM), min(WithFocusFWHM), min(RunningDiff)])-0.3,max([max(NoFocusFWHM), max(WithFocusFWHM), max(RunningDiff)])+0.3)
			pyplot.savefig(PlotFile, dpi=100)
		##
		MaxAchieved.append(max(WithFocusFWHM[20:-1]))
		MaxWorstCase.append(max(NoFocusFWHM[20:-1]))
		MeanAchieved.append(numpy.mean(WithFocusFWHM[20:-1]))
		MeanWorstCase.append(numpy.mean(NoFocusFWHM[20:-1]))
		Scatter.append(numpy.std(RunningDiff[20:-1]))
	
	result = [TargetFWHM, Aggressiveness, 
	          numpy.mean(MaxAchieved), numpy.mean(MaxWorstCase), 
	          numpy.mean(MeanAchieved), numpy.mean(MeanWorstCase),
	          numpy.mean(Scatter)]
	return result
		

def main():
	nSims = 1
	StartingError = 13
	SummaryFileName = os.path.join("tmp", "WithFocus.txt")
	print "%-9s %9s %9s %9s %9s" % ("# Target", "Agg.", "Max Err", "Mean Err", "Scatter")
	print "%-9s %9s %9s %9s %9s" % ("#   FWHM",  "",    "Focus",   "Focus",    "")
	for j in range(0,7):
		for i in range(1, 7):
			Aggressiveness = i/10.0
			TargetFWHM = j/4.0
			## Run Simulations
			result = SimulateNights(StartingError, TargetFWHM, Aggressiveness, nSims, Plot=False)
			print "%9.2f %9.1f %9.2f %9.2f %9.2f" % (result[0], result[1], result[2]-1.2, result[4]-1.2, result[6])
			## One More Sim for Plot Only
			# foo = SimulateNights(StartingError, TargetFWHM, Aggressiveness, 1, Plot=True)


	
		
if __name__ == '__main__':
	main()

