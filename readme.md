# Readme: IQMon

Copyright © 2012-2013, Dr. Josh Walawender (email: jmwalawender@gmail.com). All rights reserved.


## Overview

IQMon is a python module which can be used to quickly analyze an image for on the fly reports of image quality (**I**mage **Q**uality **Mon**itor = **IQMon**).  Originally written to analyze data from the VYSOS telescopes, but also written with the intent of building a set of general tools to analyze any fits image.

The base finctionality is that it uses SExtractor to find stars in the image and report the typical Full Width at Half Max (FWHM) and ellipticity.  This allows quick and dirty evaluation of the image quality in near real time (a few to tens of seconds on modest hardware circa 2010).

If the image contains a WCS, the module can also compare the WCS coordinates of the central pixel to the pointing coordinates in the image header to determine the pointing error.  If no WCS is present, the module can also attempt to solve the astrometry in the image using the astrometry.net solver.  

## Requirements

* python2.7.X
* astropy (<http://www.astropy.org>)
* pyephem (<http://rhodesmill.org/pyephem/>)
* SExtractor (<http://www.astromatic.net/software/sextractor>)
* astrometry.net solver (<http://astrometry.net>)

* matplotlib (Should be bundled with most python installations)
* subprocess

## Version History

* **v1.0.4**
	* MakeJPEG now marks the brightest 5000 stars rather than the first 5000 in the table.  Also annotates image to let viewer know more stars were detected.
	* added option to HTML output to choose which columns are displayed
	* other minor bug fixes and tweaks
	* fixed bug where axes in crop were reversed
* **v1.0.3**
    * Fixed handling of spaces in filenames by repalcing spaces in input file name with underscores when making working file.
    * Fixed bug in MakeJPEG which hard coded a binning of 2x2.
    * Added flag to ignore missing end card in fits header after running astrometry.net.
* **v1.0.2**
    * Fixed bug in jpeg creation which would crash when many (>~5000) stars were to be marked.  Now marks first 5000 stars.
    * Fixed bug in color coding of HTML.
* **v1.0.1**
    * Added astropy units converter between pixels and arcseconds to handle internal conversion of pixel scale.
    * Fixed bug in jpeg creation.  Will now handle rotation and marked stars.
    * Fixed bug when reading in summary text file.
* **v1.0** (released on github.com 2013/08/14)
    * Rewritten as object oriented code.  Implements most capabilities of v0.X.
    * Runs roughly 2x faster than v0.X.
* **v0.X** (frozen 2013/07/15)
    * Initial version, not under version control.
    * Deployed on live VYSOS data.
    * Analyzes image with SExtractor
        * Reports FWHM, ellipticity, background, number of stars detected
    * Reports image pointing info (alt, az, moon angle)
    * Can solve image with astrometry.net (with poor error handling)
    * Reports pointing error and position angle
    * Makes jpegs of image and cropped version with circles overlayed on stars found by SExtractor
    * Makes HTML and text file versions of results with one line per image, usually one night of images per file.


## Planned Features

… in no particular order:

* Implement reading of raw DSLR images via dcraw.
* Refine WCS by adding distortion terms.
* Determine zero point of image by comparing SExtractor photometry with catalog magnitudes (using UCAC4).


## Code Structure

IQMon functionality centers around the use of three objects:  IQMon.Config, IQMon.Telescope, and IQMon.Image.  When used, you must create one of each of these objects.

* IQMon.Config is a singleton and holds general configuration information for the module.  It has several properties, all of which are paths to things like data directories or directories for temporary files.

* IQMon.Telescope object is also a singleton and holds detailed information on the telescope used to take the data.

* IQMon.Image object holds information about the image and contains the methods which do all of the image analysis which then fills in the object properties.

### Example Use

A typical code to use IQMon on a single image might be structured something like this:

```
def listDarks():
    ## define a function to return a list of the paths
    ## to dark file(s) on your system
    return Darks

def main():
    ## main program
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
    image = IQMon.Image(FitsFile, tel, config)  ## Create image object
            
    ## Create Filenames and set verbosity
    IQMonLogFileName = "/path/to/my/log"
    htmlImageList = "/path/to/my/HTML/output"
    summaryFile = "/path/to/my/text/output"
    FullFrameJPEG = "/path/to/full/frame/jpeg"
    CropFrameJPEG = "/path/to/crop/frame/jpeg"
    verbose = True
    
    ## Perform Actual Image Analysis
    image.MakeLogger(IQMonLogFileName, verbose)
    image.logger.info("###### Processing Image:  %s ######", FitsFilename)
    image.logger.info("Setting telescope variable to %s", telescope)
    image.tel.CheckUnits()
    image.ReadImage()           ## Create working copy of image (don't edit raw file!)
    image.GetHeader()           ## Extract values from header
    image.MakeJPEG(FullFrameJPEG, rotate=True, binning=2)
    if not image.imageWCS:      ## If no WCS found in header ...
        image.SolveAstrometry() ## Solve Astrometry
        image.GetHeader()       ## Refresh Header
    image.DeterminePointingError() ## Calculate Pointing Error
    darks = ListDarks(image)    ## List dark files
    image.DarkSubtract(darks)   ## Dark Subtract Image
    image.Crop()                ## Crop Image
    image.GetHeader()           ## Refresh Header
    image.RunSExtractor()       ## Run SExtractor
    image.DetermineFWHM()       ## Determine FWHM from SExtractor results
    image.MakeJPEG(CropFrameJPEG, marked=True, binning=1)
    image.CleanUp()             ## Cleanup (delete) temporary files.
    image.CalculateProcessTime()## Calculate how long it took to process this image
    image.AddWebLogEntry(htmlImageList) ## Add line for this image to HTML table
    image.AddSummaryEntry(summaryFile)  ## Add line for this image to text table
```

## License Terms

Copyright (c) 2013, Josh Walawender
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.