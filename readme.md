# Readme: IQMon

Copyright © 2012-2013, Dr. Josh Walawender (email: jmwalawender@gmail.com). All rights reserved.


## Overview

IQMon is a python module which can be used to quickly analyze an image for on the fly reports of image quality (__I__mage __Quality__ __Mon__itor = __IQMon__).  Originally written to analyze data from the VYSOS telescopes, but also written with the intent of building a set of general tools to analyze any fits image.

The base finctionality is that it uses SExtractor to find stars in the image and report the typical Full Width at Half Max (FWHM) and ellipticity.  This allows quick and dirty evaluation of the image quality in near real time (a few to tens of seconds on modest hardware circa 2010).

If the image contains a WCS, the module can also compare the WCS coordinates of the central pixel to the pointing coordinates in the image header to determine the pointing error.  If no WCS is present, the module can also attempt to solve the astrometry in the image using the astrometry.net solver.  

## Requirements

* python2.7.X
* astropy (http://www.astropy.org)
* pyephem (http://rhodesmill.org/pyephem/)
* SExtractor (http://www.astromatic.net/software/sextractor)
* astrometry.net solver (http://astrometry.net)

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

## Use

IQMon functionality centers around the use of three objects:  IQMon.Config, IQMon.Telescope, and IQMon.Image.  When used, you must create one of each of these objects.

* IQMon.Config is a singleton and holds general configuration information for the module.  It has several properties, all of which are paths to things like data directories or directories for temporary files.

* IQMon.Telescope object is also a singleton and holds detailed information on the telescope used to take the data.

* IQMon.Image object holds information about the image and contains the methods which do all of the image analysis which then fills in the object properties.

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