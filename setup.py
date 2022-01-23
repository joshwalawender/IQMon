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
        # Ingest Pipeline
        "iqmon_ingest_monitor = iqmon.scripts.scripts:ingest_start_queue",
        "iqmon_ingest_one = iqmon.scripts.scripts:ingest_one",
        "iqmon_ingest_all = iqmon.scripts.scripts:ingest_all",
        "iqmon_ingest_clear = iqmon.scripts.scripts:ingest_clear",
        "iqmon_ingest_list = iqmon.scripts.scripts:ingest_list",
        # Analysis Pipeline
        "iqmon_analyze_monitor = iqmon.scripts.scripts:analyze_start_queue",
        "iqmon_analyze_one = iqmon.scripts.scripts:analyze_one",
        "iqmon_analyze_all = iqmon.scripts.scripts:analyze_all",
        "iqmon_analyze_cd = iqmon.scripts.scripts:analyze_cd",
        "iqmon_analyze_clear = iqmon.scripts.scripts:analyze_clear",
        "iqmon_analyze_list = iqmon.scripts.scripts:analyze_list",
        # Analysis Pipeline
        "monitor_aag = iqmon.devices.aag:monitor_aag",
        "monitor_weatherlink = iqmon.devices.weatherlink:monitor_davis_weather_link",
    ]
}


# modify the list of packages, to make sure that your package is defined correctly
setup(name=NAME,
      provides=NAME,
      version=VERSION,
      license='BSD2',
      description='Quick look imaging DRP.',
      long_description=open('readme.md').read(),
      author='Josh Walawender',
      author_email='jmwalawender@gmail.com',
      packages=['iqmon',],
      scripts=scripts,
      entry_points=entry_points,
      )
