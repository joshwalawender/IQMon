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
        "iqmon_ingest_monitor = iqmon.scripts.ingest_script:watch_directory",
        "iqmon_ingest_one = iqmon.scripts.ingest_script:analyze_one",
        "iqmon_ingest_cd = iqmon.scripts.ingest_script:change_directory",
        "iqmon_ingest_clear = iqmon.scripts.ingest_script:clear_queue",
        "iqmon_ingest_list = iqmon.scripts.ingest_script:list_queue",
        # Analysis Pipeline
        "iqmon_analyze_monitor = iqmon.scripts.analysis_script:watch_directory",
        "iqmon_analyze_one = iqmon.scripts.analysis_script:analyze_one",
        "iqmon_analyze_cd = iqmon.scripts.analysis_script:change_directory",
        "iqmon_analyze_clear = iqmon.scripts.analysis_script:clear_queue",
        "iqmon_analyze_list = iqmon.scripts.analysis_script:list_queue",
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
