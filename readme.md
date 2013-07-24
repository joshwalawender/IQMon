# Readme: IQMon

## Overview

IQMon is a python module which can be used to quickly analyze an image for on the fly reports on image quality (__I__mage __Quality__ __Mon__itor = __IQMon__).

Originally written to analyze data from the VYSO telescopes, but written with the intent of building a set of general tools to analyze any fits image.

## Version Map

* __v0.1__ Initial version, works on VYSOS Data.
	* Implemented using functions instead of objects.
	* Reports image pointing info (alt, az, moon angle)
	* Analyzes image with SExtractor
		* Reports FWHM, ellipticity, background, number of sars detected
	* Can solve image with astrometry.net (poor error handling)
	* Reports pointing error and position angle
	* Makes jpegs of image and cropped version with circles overlayed on stars found by SExtractor
	* Makes HTML and text file versions of results with one line per image, usually one night of images per file.
* __v1.0__ Object oriented code.  Implements capabilities of v0.1.
* __v1.1__ Robust error handling and timeout on astrometry.net call.
* __v1.2__ Implements reading of raw DSLR images via dcraw.
* __v1.3__ Refines WCS by adding distortion keywords.
* __v1.4__ Determines zero point of image by comparing SExtractor photometry with catalog magnitudes (using UCAC4).
