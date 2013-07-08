#!/usr/bin/env python
# encoding: utf-8
"""
MeasureScatter.py

Created by Josh Walawender on 2012-11-21.
Copyright (c) 2012 . All rights reserved.
"""

import sys
import os
import re

import matplotlib.pyplot as pyplot

def ReadCCDInspector(CCDIFile):
	MatchLine = re.compile("[0-9]{2}/[0-9]{2}/[0-9]{2}\s[0-9]{2}:[0-9]{2}:[0-9]{2}\s+([a-zA-Z0-9\-_]{1,50}\.[a-zA-Z0-9]{3,4})\s+([0-9\.]{1,5})\s+([0-9]{1,3})\s+([0-9]{1,5})\s+([0-9\.]{1,7})")
	FWHMs = []
	file = open(CCDIFile, 'r')
	for line in file:
		IsMatch = MatchLine.match(line)
		if IsMatch:
			FWHMs.append(IsMatch.group(2))
	return FWHMs
	
def main():
	FWHMs = ReadCCDInspector("Logs/CCDInspector_20121010_slope12stepsOut.txt")

	pyplot.ioff()
	pyplot.figure()
	pyplot.plot(FWHMs, 'b,-')
	pyplot.savefig("Scatter.png", dpi=100)

if __name__ == '__main__':
	main()

