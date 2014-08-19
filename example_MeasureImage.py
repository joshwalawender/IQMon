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
    path_temp = '/Users/vysosuser/IQMon/tmp'
    path_plots = '/Users/vysosuser/IQMon/Plots'
    tel = IQMon.Telescope(path_temp, path_plots)
    tel.name = 'SVQ100'
    tel.long_name = 'SVQ100'
    tel.focal_length = 580.*u.mm
    tel.pixel_size = 6.5*u.micron
    tel.aperture = 100.*u.mm
    tel.gain = 1.0 / u.adu
    tel.units_for_FWHM = 1.*u.pix
    tel.ROI = None
    tel.threshold_FWHM = 3.0*u.pix
    tel.threshold_pointing_err = 10.0*u.arcmin
    tel.threshold_ellipticity = 0.25*u.dimensionless_unscaled
    tel.pixel_scale = tel.pixel_size.to(u.mm)/tel.focal_length.to(u.mm)*u.radian.to(u.arcsec)*u.arcsec/u.pix
    tel.fRatio = tel.focal_length.to(u.mm)/tel.aperture.to(u.mm)
    tel.SExtractor_params = {'PHOT_APERTURES': '6.0',
                            'BACK_SIZE': '32',
                            'SEEING_FWHM': '2.5',
                            'SATUR_LEVEL': '50000',
                            'DETECT_MINAREA': '5',
                            'DETECT_THRESH': '5.0',
                            'ANALYSIS_THRESH': '5.0',
                            'FILTER': 'N',
                            }
    tel.pointing_marker_size = 4*u.arcmin
    ## Define Site (ephem site object)
    tel.site = ephem.Observer()
    tel.check_units()
    tel.define_pixel_scale()

    ##-------------------------------------------------------------------------
    ## Create IQMon.Image Object
    ##-------------------------------------------------------------------------
    image = IQMon.Image(args.filename, tel=tel)  ## Create image object

    ##-------------------------------------------------------------------------
    ## Create Filenames
    ##-------------------------------------------------------------------------
    path_log = os.path.join(os.path.expanduser('~'), 'IQMon', 'Logs')
    IQMonLogFileName = os.path.join(path_log, "IQMonLog.txt")
    htmlImageList = os.path.join(path_log, "IQMon.html")
    summaryFile = os.path.join(path_log, "IQMon.txt")
    if args.clobber:
        if os.path.exists(IQMonLogFileName): os.remove(IQMonLogFileName)
        if os.path.exists(htmlImageList): os.remove(htmlImageList)
        if os.path.exists(summaryFile): os.remove(summaryFile)

    ##-------------------------------------------------------------------------
    ## Perform Actual Image Analysis
    ##-------------------------------------------------------------------------
    image.make_logger(IQMonLogFileName, args.verbose)
    image.logger.info("###### Processing Image:  %s ######", FitsFilename)
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
    image.make_PSF_plot()
    image.measure_zero_point(plot=True)
    CatalogJPEG = image.rawFileBasename+"_catstars.jpg"
    image.make_JPEG(CatalogJPEG, mark_catalog_stars=True, mark_pointing=True, binning=2)
    
    image.clean_up()               ## Cleanup (delete) temporary files.
    image.calculate_process_time() ## Calculate how long it took to process this image
    fields=["Date and Time", "Filename", "Alt", "Az", "Airmass", "MoonSep", "MoonIllum", "FWHM", "ellipticity", "ZeroPoint", "PErr", "PosAng", "nStars", "ProcessTime"]
    image.add_web_log_entry(htmlImageList, fields=fields) ## Add line for this image to HTML table
    image.add_summary_entry(summaryFile)  ## Add line for this image to text table
    

if __name__ == '__main__':
    main()





