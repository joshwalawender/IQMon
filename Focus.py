import os
import numpy
import astropy
import astropy.io
import astropy.io.ascii
import astropy.table

#############################################################################
## Focus by Simple Reactive Rules:
## - if FWHM is withing Target, then do nothing
## - if improving and FWHM large, move again in same direction
## - if not improving and FWHM large, move in opposite direction
## Move Amount is Aggressiveness times distace to move to target FWHM
#############################################################################
def Reactive(SummaryFileName, tel):
	History = []
	Moves = []
	MoveSum = 0.0
	if os.path.exists(SummaryFileName):
		IQMonTable = astropy.io.ascii.read(SummaryFileName)
		FWHM = IQMonTable['FWHM']
		Moved = IQMonTable['Move']
		History.append(FWHM)
		Moves.append(Moved)
		
		# SumFile = open(SummaryFileName, 'r')
		# for line in SumFile:
		# 	if line[0] != "#":
		# 		FWHM  = float(line[71:82])
		# 		Moved = float(line[84:95])
		# 		History.append(FWHM)
		# 		Moves.append(Moved)
		nMeasures = len(History)
		if nMeasures >= 2:
			LastFWHM = History[-1]
			S2LastFWHM = History[-2]
			LastMove = Moves[-1]
			## Is FWHM Improving?
			if History[-1] <= History[-2]:
				Improving = True
			else:
				Improving = False
			## Adjust Focus
			## - if improving and FWHM small, then do nothing
			if Improving and (LastFWHM <= tel['TargetFWHM']):
				Direction = 0.0
			## - if not improving and FWHM small, do nothing
			if not Improving and (LastFWHM <= tel['TargetFWHM']):
				Direction = 0.0
			## - if improving and FWHM large, move again in same direction
			if Improving and (LastFWHM >= tel['TargetFWHM']):
				if Moves[-1] < 0:
					Direction = -1
				if Moves[-1] >= 0:
					Direction = 1
			## - if not improving and FWHM large, move in opposite direction
			if not Improving and (LastFWHM >= tel['TargetFWHM']):
				if Moves[-1] < 0:
					Direction = 1
				if Moves[-1] >= 0:
					Direction = -1
			DeltaFWHM = abs(LastFWHM - tel['TargetFWHM'])
			Move = Direction*int(round(tel['Aggressiveness']*DeltaFWHM*tel['PixelSize']*tel['FRatio']/tel['StepSize']))
			MoveSum = numpy.sum(Moves)+Move
		else:
			Move = 0.0
	else:
		Move = 0.0

	if Move >= tel['MaxMove']:
		Move = tel['MaxMove']
	elif Move <= -1.*tel['MaxMove']:
		Move = -1.*tel['MaxMove']	
	
	return Move
	
	
#############################################################################
## Focus by Adaptive Reaction Rules:
##   Use same basic rules as the Reactive algorithm, except look to see if
##   recent moves have all been in same direction or not.  Change the 
##   Aggressiveness term based on that.  If similar amount of positive
##   and neg moves, then decrease aggressivness.  If most moves are in 
##   same direction, then increase Aggressiveness.
#############################################################################
# def AdaptiveReactive(SummaryFileName, tel):
	
#############################################################################
## Focus by Temperature Rules:
##   Works similar to the SimpleRules algorithm above except that in the case
##   where a move is desired AND dTdt is above a set threshold, move based on
##   dTdt first, then if that does not improve 
## - if FWHM is within target, do nothing
## - if FWHM is larger than target and dT/dt is large, move in direction indicated by dT/dt
##   - if last move was from dT/dt and FWHM increases, then reverse direction
## - if FWHM is larger then target and dT/dt is small, do nothing
##   - if FWHM is larger than target again and dT/dt is small, move
##     - if FWHM continues to increase, reverse direction
## - Move Amount is Aggressiveness times distace to move to target FWHM
#############################################################################
# def TemperatureReactive(SummaryFileName, EnvironmentalLogFile, tel):

	
	
	