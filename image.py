#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2013-07-13.
Copyright (c) 2013 . All rights reserved.
"""

import sys
import os
import unittest


class Image(object):
	'''
	The IQMon.Image object represents a single input image to the IQMon process.
	
	When defined, the image objects requires both a filename to a valid fits file and an IQMon.Config object.
	
	Properties:
	- original file name
	- working file name
	- target object name
	- target RA
	- target Dec
	- target alt
	- target az
	- image WCS
	- image exposure time
	- moon alt
	- moon separation from target
	- moon illumination
	
	Methods:
	- extract header info
	- read input file (rfits, dcraw, etc.)
	- dark subtract image
	- solve image using astrometry.net
	- crop image
	- filter cosmic rays
	- refine WCS
	- find stars (run sextractor) (detemine FWHM, ellipticity)
	- determine pointing error
	- determine zero point
	- make image plots
	- make jpeg (full frame or cropped)
	- calculate focus move
	'''
	def __init__(self, filename, config):
		pass


class untitledTests(unittest.TestCase):
	def setUp(self):
		pass


if __name__ == '__main__':
	unittest.main()