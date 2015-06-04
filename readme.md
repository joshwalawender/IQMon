# Readme: IQMon

Copyright © Dr. Josh Walawender (email: jmwalawender@gmail.com). All rights reserved.


## Overview

IQMon is a python module which can be used to quickly analyze an image for on the fly reports of image quality (Image Quality Monitor = IQMon).  It was originally written to provide a quick look analysis of data from robotic telescopes.

The base functionality is that it uses SExtractor to find stars in the image and report the typical Full Width at Half Max (FWHM) and ellipticity.  This allows quick and dirty evaluation of the image quality in near real time (a few to tens of seconds on modest hardware circa 2010).

If the image contains a WCS, the module can also compare the WCS coordinates of the central pixel to the pointing coordinates in the image header to determine the pointing error.  If no WCS is present, the module can also attempt to solve the astrometry in the image using the astrometry.net solver.  

In addition, the software can use SCAMP and SWarp to determine a plate solution which includes distortion parameters and then compare the instrumental magnitudes determined by SExtractor to catalog magnitudes to determine the photometric zero point of the image.


## Requirements

Python Modules:

* python2.7.X or python3.X
* astropy (<http://www.astropy.org>)
* pyephem (<http://rhodesmill.org/pyephem/>)
* astroquery
* numpy
* matplotlib
* subprocess32
* PyYAML
* PIL
* skimage

External Programs:

* SExtractor (<http://www.astromatic.net/software/sextractor>)
* SCAMP (<http://www.astromatic.net/software/scamp>)
* SWarp (<http://www.astromatic.net/software/swarp>)
* MissFITS (<http://www.astromatic.net/software/missfits>)
* astrometry.net solver (<http://astrometry.net>)


## Code Structure

IQMon functionality centers around the use of two objects:  IQMon.Telescope and IQMon.Image.  When used, you must create one of each of these objects.

* IQMon.Telescope object and holds detailed information on the telescope used to take the data and some custom configuration parameters.  The telescope object is passed to the image object on creation of an image object.

* IQMon.Image object holds information about the image and contains the methods which do all of the image analysis which then fills in the object properties.


### Example Use

See example_MeasureImage.py file in this repository.  This file may be slightly out of date relative to the main program as it is not updated as often because it is not in use with any real telescopes.

A better example is to look at the MeasureImage.py file in the VYSOStools repository (<https://github.com/joshwalawender/VYSOStools>).  This is a repository of programs for the VYSOS robotic telescopes and is always updated to use the very latest commit of IQMon.  It is a bit more complex than the example included in this repository because it handles two telescopes.


## License Terms

Please see LICENSE file.
