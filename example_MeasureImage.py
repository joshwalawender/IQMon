#!/usr/bin/env python
'''
This is an example of how to create a script to analyze an image using the
IQMon module tools.  This script is not general and will need to be customized
to the particular telescope and camera system of interest.
'''

from __future__ import division, print_function

import sys
import os
from argparse import ArgumentParser
import re
import datetime
import math

import ephem
import astropy.units as u

import IQMon


##-------------------------------------------------------------------------
## Main Program
##-------------------------------------------------------------------------
def main():
    ##-------------------------------------------------------------------------
    ## Parse Command Line Arguments
    ##-------------------------------------------------------------------------
    ## create a parser object for understanding command-line arguments
    parser = ArgumentParser(description="Describe the script")
    ## add flags
    parser.add_argument("-v", "--verbose",
        action="store_true", dest="verbose",
        default=False, help="Be verbose! (default = False)")
    parser.add_argument("-c", "--clobber",
        action="store_true", dest="clobber",
        default=False, help="Delete previous logs and summary files for this image. (default = False)")
    ## add arguments
    parser.add_argument("filename",
        type=str,
        help="File Name of Input Image File")
    args = parser.parse_args()


    ##-------------------------------------------------------------------------
    ## Create Telescope Object
    ##-------------------------------------------------------------------------
    config_file = os.path.join(os.path.expanduser('~'), 'IQMon', 'config_VYSOS-5.yaml')
    tel = IQMon.Telescope(config_file)

    ##-------------------------------------------------------------------------
    ## Create IQMon.Image Object
    ##-------------------------------------------------------------------------
    image = IQMon.Image(args.filename, tel=tel)  ## Create image object

    ##-------------------------------------------------------------------------
    ## Create Filenames
    ##-------------------------------------------------------------------------
    logs_file = os.path.join(tel.logs_file_path, DataNightString+"_"+telescope+"_IQMonLog.txt")
    html_file = os.path.join(tel.logs_file_path, DataNightString+"_"+telescope+".html")
    yaml_file = os.path.join(tel.logs_file_path, DataNightString+"_"+telescope+"_Summary.txt")
    if args.clobber:
        if os.path.exists(logs_file): os.remove(logs_file)
        if os.path.exists(html_file): os.remove(html_file)
        if os.path.exists(yaml_file): os.remove(yaml_file)

    ##-------------------------------------------------------------------------
    ## Perform Actual Image Analysis
    ##-------------------------------------------------------------------------
    image.make_logger(logfile=logs_file, verbose=verbose)
    image.logger.info("###### Processing Image: {} ######".format(FitsFilename))
    image.read_image()           ## Create working copy of image (don't edit raw file!)
    image.read_header()           ## Extract values from header

    if not image.image_WCS:      ## If no WCS found in header ...
        image.solve_astrometry() ## Solve Astrometry
        image.read_header()       ## Refresh Header
    image.determine_pointing_error()            ## Calculate Pointing Error

    image.run_SExtractor()       ## Run SExtractor
    image.determine_FWHM()       ## Determine FWHM from SExtractor results

    image.run_SCAMP(catalog='UCAC-3')
    image.run_SWarp()
    image.read_header()           ## Extract values from header
    image.get_local_UCAC4(local_UCAC_command="/Users/joshw/Data/UCAC4/access/u4test", local_UCAC_data="/Users/joshw/Data/UCAC4/u4b")
    image.run_SExtractor(assoc=True)
    image.determine_FWHM()       ## Determine FWHM from SExtractor results
    image.measure_zero_point(plot=True)
    image.make_PSF_plot()

    small_JPEG = image.raw_file_basename+"_fullframe.jpg"
    image.make_JPEG(small_JPEG, binning=3,\
                    p1=0.15, p2=0.50,\
                    make_hist=False,\
                    mark_pointing=True,\
                    mark_detected_stars=False,\
                    mark_catalog_stars=False,\
                    mark_saturated=True,\
                    transform='rotate90')
                    )
    cropped_JPEG = image.raw_file_basename+"_crop.jpg"
    image.make_JPEG(cropped_JPEG,\
                    p1=0.15, p2=0.50,\
                    make_hist=False,\
                    mark_pointing=True,\
                    mark_detected_stars=True,\
                    mark_catalog_stars=False,\
                    mark_saturated=True,\
                    crop=(1024, 1024, 3072, 3072),\
                    transform='rotate90')
                    )
    
    image.clean_up()               ## Cleanup (delete) temporary files.
    image.calculate_process_time() ## Calculate how long it took to process this image
    fields=["Date and Time", "Filename", "Alt", "Az", "Airmass", "MoonSep", "MoonIllum", "FWHM", "ellipticity", "ZeroPoint", "PErr", "PosAng", "nStars", "ProcessTime"]
    image.add_web_log_entry(html_file, fields=fields)
    image.add_yaml_entry(yaml_file)
    

if __name__ == '__main__':
    main()





