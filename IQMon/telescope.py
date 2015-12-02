#!/usr/bin/env python
# encoding: utf-8
"""
IQMon.py

Created by Josh Walawender on 2013-07-22.
Copyright (c) 2013 . All rights reserved.
"""
from __future__ import division, print_function

## Import General Tools
import sys
import os
import re
import stat
import shutil
import datetime
if sys.version_info.major == 2:
    import subprocess32 as subprocess
elif sys.version_info.major == 3:
    import subprocess
import logging
import yaml
import math
import numpy as np
import pymongo
from pymongo import MongoClient


## Import Astronomy Specific Tools
import ephem
import astropy.units as u
import astropy.io.fits as fits
import astropy.coordinates as coords
import astropy.table as table
import astropy.wcs as wcs
import astropy.io.ascii as ascii


class TelescopeConfigError(Exception):
    pass

##-----------------------------------------------------------------------------
## Define Telescope object to hold telescope information
##-----------------------------------------------------------------------------
class Telescope(object):
    '''Object which contains information about the telescope which took the
    Image (see IQMon.Image object definition).
    
    Parameters
    ----------
    config_file : string
        Path to the configuration file to be read.  The telescope properties are
        definied when a configuration file is read.  The configuration file is a
        YAML formatted text file which can contain the following entries:

        name : string containing the name of the telescope

        logs_file_path : the path in to which log files will be written

        plot_file_path : the path in to which the plots (e.g. those from the
            image.make_PSF_plot() or image.make_zero_point_plot() methods) will
            be written.

        temp_file_path : the path used for temporary files created by the
            program

        mongo_address : string containing the address to connecto to the mongo
            server which, if used, will contain a database of image analysis
            results.

        mongo_port : integer containing the port number of the mongo server

        mongo_db : string containing the name of the mongo database to use

        mongo_collection : string containing the name of the mongo collection to
            use

        pixel_scale : float containing the pixel scale in arcseconds per pixel.
            If either the pixel_scale or both the focal_length and pixel_size
            are required to be in the configuration file.

        focal_length : interger containing the focal length of the telescope in
            mm.  This is used in estimating the pixel scale prior to plate
            solving.

        pixel_size : float containing the size of a pixel in microns.  This is
            used in estimating the pixel scale prior to plate solving.

        gain : float with the estimated gain (electroncs per ADU) of the
            detector.  This value is used by source extractor.  If it is not
            present a default value of 1.0 will be used.

        saturation : float with the saturation level of the detector in ADU.
            This is used in optionally marking saturated pixels in the jpegs
            made by the make_JPEG() method of the Image object.

        threshold_FWHM : float value (in units of pixels) with the threshold
            FWHM value.  If the image FWHM is above this value, the FWHM flag
            will be set.

        threshold_pointing_err : float value (in units of arcmin) with the
            threshold pointing error.  If the pointing error is above this
            value, the pointing error flag will be set.

        threshold_ellipticity : float value (unitless) with the threshold
            ellipticity.  If the ellipticity is above this value, the
            ellipticity flag will be set.

        threshold_zeropoint : float value (in magnitudes) with the threshold
            zero point.  If the zero point is above this value, the zero point
            flag will be set.

        units_for_FWHM : The units which will be used when displaying the FWHM
            value on the web page.  Also often used by customized scripts for
            displaying IQMon results.

        ROI : string representing the region of interest in pixels which the
            image should be cropped to.  Format is "[x1:x2,y1:y2]" where x1 is
            the minimum x pixel, x2 is the maximum x pixel, y1 is the minimum y
            pixel, and y2 is the maximum y pixel.  All pixel values will be
            forced to integers.

        PSF_measurement_radius : Radius in pixels of central region for which
            the image quality measurements of stars are used in determining the
            image FWHM and ellipticity.  This parameter allows the user to
            ignore the corners of the image if they wish to ignore the optical
            aberrations in that region.

        pointing_marker_size : Diameter (in arcminutes) of the pointing marker
            symbol on the image jpegs generated by the make_JPEG() method.

        SExtractor_params : Dictionary of source extractor parameters to be used
            when source extractor is called.  For example, when you want to set
            the CATALOG_TYPE to FITS_LDAC, use {'CATALOG_TYPE': 'FITS_LDAC'}.

        SCAMP_params : Dictionary of SCAMP parameters to be used when SCAMP is
            called.

        catalog : Dictionary containing information about the stellar catalog
            which is matched to the detected stars to determine the zero point.
            The get_catalog() method will query Vizier using these parameters.
        
            name : The name of the catalog which will be passed to a
                astroquery.vizier.Vizier() instance.
        
            columns : List of the column names to retrieve.  Defaults to
                ['_RAJ2000','_DEJ2000','UCAC4','Bmag','Vmag','gmag','rmag','imag']
                if the catalog is UCAC4.
        
            magmax : Maximum magnitude to retrieve.
        
            Remaining parameters are the dictionary which link the filter names
            in the fits header to the filter names in the Vizier catalog that is
            retrieved.  For example, if I want to compare the source extractor
            catalog with the UCAC4 catalog limited to stars brighter than
            magnitude 15 and I want to compare the 'rmag' filter magnitudes from
            the catalog to my images which have the FILTER header keyword set to
            'PSr', then I would have a configuration file which contains the
            following:

            catalog:
                name: 'UCAC4'
                magmax: 15.0
                PSr: 'rmag'
    '''
    def __init__(self, config_file):
        self.site = ephem.Observer()

        ## Read YAML Config File
        config_file = os.path.expanduser(config_file)
        if not os.path.exists(config_file):
            raise TelescopeConfigError('Configuration file {} not found'.format(config_file))
        with open(config_file, 'r') as yaml_string:
            config = yaml.load(yaml_string)
        if not isinstance(config, dict):
            raise TelescopeConfigError('Configuration file contents not parsed as dict')
        self.config = config

        ## Populate Configured Properties
        if 'name' in config.keys():
            self.name = str(config['name'])
        else:
            self.name = 'telescope'

        if 'temp_file_path' in config.keys():
            self.temp_file_path = os.path.expanduser(config['temp_file_path'])
        else:
            self.temp_file_path = os.path.join('/', 'tmp')

        if 'plot_file_path' in config.keys():
            self.plot_file_path = os.path.expanduser(config['plot_file_path'])
        else:
            self.plot_file_path = os.path.join('/', 'tmp')

        if 'logs_file_path' in config.keys():
            self.logs_file_path = os.path.expanduser(config['logs_file_path'])
        else:
            self.logs_file_path = os.path.join('/', 'tmp')

        if 'mongo_address' in config.keys():
            self.mongo_address = config['mongo_address']
        else:
            self.mongo_address = None

        if 'mongo_port' in config.keys():
            self.mongo_port = config['mongo_port']
        else:
            self.mongo_port = 27017 ## Use default mongo port

        if 'mongo_db' in config.keys():
            self.mongo_db = config['mongo_db']
        else:
            self.mongo_db = None

        if 'mongo_collection' in config.keys():
            self.mongo_collection = config['mongo_collection']
        else:
            self.mongo_collection = None

        ## Define astropy.units Equivalency for Arcseconds and Pixels
        self.pixel_scale_equivalency = [(u.pix, u.arcsec,
             lambda pix: (pix*u.radian.to(u.arcsec) * self.pixel_size\
                         / self.focal_length).decompose().value,
             lambda arcsec: (arcsec/u.radian.to(u.arcsec) * self.focal_length\
                            / self.pixel_size).decompose().value
             )]

        if not 'pixel_scale' in config.keys():
            if ('focal_length' in config.keys()) and ('pixel_size' in config.keys()):
                self.focal_length = config['focal_length'] * u.mm
                self.pixel_size = config['pixel_size'] * u.um
                self.pixel_scale = self.pixel_size.to(u.mm)\
                                   /self.focal_length.to(u.mm)\
                                   *u.radian.to(u.arcsec)*u.arcsec/u.pix
            else:
                raise TelescopeConfigError('Configuration file does not contain\
                                          information to determine pixel_scale')
        else:
            self.focal_length = None
            self.pixel_size = None
            self.pixel_scale = config['pixel_scale'] * u.arcsec/u.pix

        if 'latitude' in config.keys():
            self.latitude = float(config['latitude']) * u.deg
        else:
            self.latitude = None

        if 'longitude' in config.keys():
            self.longitude = float(config['longitude']) * u.deg
        else:
            self.longitude = None

        if 'altitude' in config.keys():
            self.altitude = float(config['altitude']) * u.meter
        else:
            self.altitude = None

        if 'gain' in config.keys():
            self.gain = config['gain'] / u.adu
        else:
            self.gain = 1.0 / u.adu

        if 'saturation' in config.keys():
            self.saturation = config['saturation'] * u.adu
        else:
            self.saturation = None

        if 'threshold_FWHM' in config.keys():
            self.threshold_FWHM = config['threshold_FWHM'] * u.pix
        else:
            self.threshold_FWHM = None

        if 'threshold_pointing_err' in config.keys():
            self.threshold_pointing_err = config['threshold_pointing_err'] * u.arcmin
        else:
            self.threshold_pointing_err = None

        if 'threshold_ellipticity' in config.keys():
            self.threshold_ellipticity = config['threshold_ellipticity']
        else:
            self.threshold_ellipticity = None

        if 'threshold_zeropoint' in config.keys():
            self.threshold_zeropoint = config['threshold_zeropoint'] * u.mag
        else:
            self.threshold_zeropoint = None

        if 'units_for_FWHM' in config.keys():
            self.units_for_FWHM = getattr(u, config['units_for_FWHM'])
        else:
            self.units_for_FWHM = u.pix

        if 'ROI' in config.keys():
            self.ROI = str(config['ROI'])
        else:
            self.ROI = None

        if 'PSF_measurement_radius' in config.keys():
            self.PSF_measurement_radius = config['PSF_measurement_radius'] * u.pix
        else:
            self.PSF_measurement_radius = None

        if 'pointing_marker_size' in config.keys():
            self.pointing_marker_size = config['pointing_marker_size'] * u.arcmin
        else:
            self.pointing_marker_size = 1 * u.arcmin

        if 'SExtractor_params' in config.keys():
            self.SExtractor_params = config['SExtractor_params']
        else:
            self.SExtractor_params = None

        if 'SCAMP_params' in config.keys():
            self.SCAMP_params = config['SCAMP_params']
        else:
            self.SCAMP_params = None

        if 'catalog' in config.keys():
            self.catalog_info = config['catalog']
        else:
            self.catalog_info = None

        ## create paths
        paths_to_check = []
        paths_to_create = []
        if self.temp_file_path: paths_to_check.append(self.temp_file_path)
        if self.plot_file_path: paths_to_check.append(self.plot_file_path)
        if self.logs_file_path: paths_to_check.append(self.logs_file_path)
        for path in paths_to_check:
#             print('Checking: {}'.format(path))
            while not os.path.exists(path):
#                 print('Need to create: {}'.format(path))
                paths_to_create.append(path)
                path = os.path.split(path)[0]
        while len(paths_to_create) > 0:
            new_path = paths_to_create.pop()
#             print('Creating: {}'.format(new_path))
            os.mkdir(new_path)


    def __del__(self):
        pass
#         print('Deleted telescope object')

    def __enter__(self):
        return self

    def __exit__(self ,type, value, traceback):
        self.__del__()
