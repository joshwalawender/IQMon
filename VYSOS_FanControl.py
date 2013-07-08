#!/usr/bin/env python
# encoding: utf-8
"""
VYSOS_FanControl.py

Created by Josh Walawender on 2013-06-17.
Copyright (c) 2013 . All rights reserved.
"""

import sys
import getopt


help_message = '''
The help message goes here.
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):
	## Set telescope to default value of VYSOS-20
	telescope = "V20"
	
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ht:", ["help", "telescope="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-t", "--telescope"):
				telescope = value
	
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	## VYSOS-20
	if telescope == "V20":
		
	## VYSOS-5
	if telescope == "V5":
		pass

if __name__ == "__main__":
	sys.exit(main())
