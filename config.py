#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Josh Walawender on 2013-07-16.
Copyright (c) 2013 . All rights reserved.
"""

import sys
import os
import unittest


class Config(object):
	'''
	The IQMon.Config object contains configuration information for processing a particular image.
	
	Properties:
	- telescope name
	- telescope focal length
	- telescope aperture
	- telescope pixel scale
	- site name
	- site latitude
	- site longitude
	- site elevation
	
	Methods:
	- read config file
	'''
	def __init__(self):
		pass


class untitledTests(unittest.TestCase):
	def setUp(self):
		pass


if __name__ == '__main__':
	unittest.main()