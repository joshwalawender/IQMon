# !usr/bin/env python

import sys
import os
import glob

from setuptools import setup

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])
#metadata = dict(conf.items('metadata'))


# MODIFY THE NAME OF THE PACKAGE to be the one chosen
NAME = 'IQMon'
VERSION = '2.0dev'
RELEASE = 'dev' not in VERSION

scripts = []

entry_points = {
    'console_scripts': [
#         "analyzeone = vysosdrp.script:analyze_one",
#         "watchdirectory = vysosdrp.script:watch_directory",
#         "qlcd = vysosdrp.script:change_directory",
#         "qlclear = vysosdrp.script:clear_queue",
#         "listqueue = vysosdrp.script:list_queue",
    ]
}


# modify the list of packages, to make sure that your package is defined correctly
setup(name=NAME,
      provides=NAME,
      version=VERSION,
      license='BSD2',
      description='Quick look imaging DRP.',
      long_description=open('README.txt').read(),
      author='Josh Walawender',
      author_email='jmwalawender@gmail.com',
      packages=['iqmon',],
      scripts=scripts,
      entry_points=entry_points,
      )
