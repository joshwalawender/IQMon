# Readme: IQMon

Copyright © 2012-2013, Dr. Josh Walawender (email: jmwalawender@gmail.com). All rights reserved.


## Overview

IQMon is a python module which can be used to quickly analyze an image for on the fly reports of image quality (__I__mage __Q__uality __Mon__itor = __IQMon__).  Originally written to analyze data from the VYSOS telescopes, but also written with the intent of building a set of general tools to analyze any fits image.

The base finctionality is that it uses SExtractor to find stars in the image and report the typical Full Width at Half Max (FWHM) and ellipticity.  This allows quick and dirty evaluation of the image quality in near real time (a few to tens of seconds on modest hardware circa 2010).

If the image contains a WCS, the module can also compare the WCS coordinates of the central pixel to the pointing coordinates in the image header to determine the pointing error.  If no WCS is present, the module can also attempt to solve the astrometry in the image using the astrometry.net solver.  

## Requirements

* python2.7.X
* astropy (<http://www.astropy.org>)
* pyephem (<http://rhodesmill.org/pyephem/>)
* SExtractor (<http://www.astromatic.net/software/sextractor>)
* astrometry.net solver (<http://astrometry.net>)

* matplotlib (Should be bundled with most python installations)
* subprocess32

## Version History

* __v1.0__ (released on github.com 2013/08/??)
    * Rewritten as object oriented code.  Implements most capabilities of v0.1.
    * Does not color code HTML table with alerts to poor values of FWHM, pointing error, etc.
    * Runs roughly 2x faster than v0.X.
* __v0.X__ (frozen 2013/07/15)
    * Initial version, not under version control.
    * Works on VYSOS Data.
    * Implemented using functions instead of objects.
    * Reports image pointing info (alt, az, moon angle)
    * Analyzes image with SExtractor
        * Reports FWHM, ellipticity, background, number of stars detected
    * Can solve image with astrometry.net (poor error handling)
    * Reports pointing error and position angle
    * Makes jpegs of image and cropped version with circles overlayed on stars found by SExtractor
    * Makes HTML and text file versions of results with one line per image, usually one night of images per file.

## Code Structure

IQMon functionality centers around the use of three objects:  IQMon.Config, IQMon.Telescope, and IQMon.Image.  When used, you must create one of each of these objects.

* IQMon.Config is a singleton and holds general configuration information for the module.  It has several properties, all of which are paths to things like data directories or directories for temporary files.

* IQMon.Telescope object is also a singleton and holds detailed information on the telescope used to take the data.

* IQMon.Image object holds information about the image and contains the methods which do all of the image analysis which then fills in the object properties.

### Example Use

A typical code to use IQMon on a single image might be structured something like this:

    def listDarks():
        ## define a function to return a list of the paths
        ## to dark file(s) on your system
        return Darks

    def main():
        ## main program, possibly set up command line argument
        ## to accept the filename of the fits image you want to analyze
        FitsFile = "/path/to/fits/file"
        
        ## Establish IQMon Configuration
        config = IQMon.Config()  ## This reads configuration info in
                                 ## your $HOME/.IQMonConfig file
        ## Create Telescope Object
        tel = IQMon.Telescope()
        ## Define your telescope's properties
        tel.name = "MyTelescope"
        tel.longName = "MyReallyExcellentTelescope"
        tel.focalLength = 735.*u.mm
        tel.pixelSize = 9.0*u.micron
        tel.aperture = 135.*u.mm
        tel.gain = 1.6 / u.adu
        tel.unitsForFWHM = 1.*u.pix
        tel.ROI = "[1024:3072,1024:3072]"
        tel.thresholdFWHM = 2.5*u.pix
        tel.thresholdPointingErr = 5.0*u.arcmin
        tel.thresholdEllipticity = 0.30
        tel.pixelScale = tel.pixelSize.to(u.mm) / tel.focalLength.to(u.mm) * u.radian.to(u.arcsec) * u.arcsec / u.pix
        tel.fRatio = tel.focalLength.to(u.mm) / tel.aperture.to(u.mm)
        tel.SExtractorPhotAperture = 6.0*u.pix
        tel.SExtractorSeeing = 2.0*u.arcsec
        tel.site = ephem.Observer()
        
        ## Create IQMon.Image Object
        image = IQMon.Image(FitsFile)
        
        ## Create Logger Object
        IQMonLogFileName = "/path/to/my/log"
        image.htmlImageList = "/path/to/my/HTML/output"
        image.summaryFile = "/path/to/my/text/output"
        verbose = True
        logger = config.MakeLogger(IQMonLogFileName, verbose)
        
        ## Perform Actual Image Analysis
        tel.CheckUnits(logger)
        image.ReadImage(config)        ## Create working copy of image (don't edit raw file!)
        image.GetHeader(tel, logger)   ## Extract values from header
        image.MakeJPEG(image.rawFileBasename+"_full.jpg", tel, config, logger, rotate=True, binning=2)
        if not image.imageWCS:                          ## if no WCS found in header ...
            image.SolveAstrometry(tel, config, logger)  ## Solve Astrometry
            image.GetHeader(tel, logger)                ## Refresh Header
        image.DeterminePointingError(logger)            ## Calculate Pointing Error
        darks = ListDarks(image, tel, config, logger)   ## List dark files
        image.DarkSubtract(darks, tel, config, logger)  ## Dark Subtract Image
        image.Crop(tel, logger)                         ## Crop Image
        image.GetHeader(tel, logger)                    ## Refresh Header
        image.RunSExtractor(tel, config, logger)        ## Run SExtractor
        image.DetermineFWHM(logger)                     ## Determine FWHM from SExtractor Results
        image.CleanUp(logger)                           ## Delete temporary files
        image.CalculateProcessTime(logger)              ## Calculate how long it took to process this image
        image.AddWebLogEntry(tel, config, logger)       ## Add line for this image to your HTML table
        image.AddSummaryEntry(logger)                   ## Add line for this image to your text table





## Planned Features

… in no particular order:

* Add color coding of HTML output table to allow marking of poor FWHM, ellipticity, pointing error, moon separation, and airmass values.
* Robust error handling and timeout on astrometry.net call.
* Add method to write image analysis results to MySQL database
    * Add method to make HTML table from database of images rather than image by image.  This would not be a method of image, but either of telescope or config or perhaps a new object.
* Implement reading of raw DSLR images via dcraw.
* Refine WCS by adding distortion terms.
* Determine zero point of image by comparing SExtractor photometry with catalog magnitudes (using UCAC4).

## License Terms

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.