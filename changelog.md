## Version History

* __v1.5__
    * Add tornado web application for serving results stored in mongo
    * Additional doc strings
    * Add plot keyword to ``measure_FWHM`` method
    * Add keyword to pass parameters to ``run_SCAMP`` and ``run_SExtractor``
    * Get mogno address, port, database, and collection from config
    * Fix orientation of jpeg
    * Fix bug in ``crop`` method which would lose WCS info
    * Add ``read_header`` to ``solve_astrometry``, ``crop``, and ``run_SWarp``
* __v1.4__
    * Bug fixes
    * Add `add_mongo_entry()` method
    * Add `uncompress()` method to handle fpack compressed images
    * Use subprocess32 and timeout for external calls to SCAMP, SExtractor, astrometry.net
* __v1.3.4__
    * Various minor bug fixes.
    * Fixed bug in calculation of weights for FWHM and zero point calculations.  Very little actual impact on results.
* __v1.3.3__
    * Bug fixes.
    * Added version info for IQMon, SExtractor, SCAMP, and SWarp to log.
* __v1.3.2__
    * Reduced threshold number of stars to calculate zero point to 10.
    * Made downsampling factor for solve-field (from astrometry.net) a user settable argument in `solve_astrometry()`.  Defaults to 4.
    * Stored `__version__` property in IQMon module instead of in setup.py file.
    * Records IQMon version in the YAML entry for that image.
* __v1.3.1__
    * Bug fix to related to coordinate wrapping in `get_catalog` and `get_local_UCAC` methods.
    * Filtered out stars with magnitude of exactly 0.0 from zero point calculation as some catlaogs seem to return 0.0 in place of `nan` or `None`.
* __v1.3__
    * Can now query Vizier for USNO-B1.0 and UCAC4 catalogs using `get_catalog` method.
    * Calculation of FWHM, ellipticity, and zero point now all use a weighted average where the weight is equal to the signal to noise calculated from the source extractor photometry.
    * Improved configuration files to handle SCAMP configuration and photometric catalog information
    * Default behavior is now to write new IQMon log file for each image.
    * Minor improvements to PSF plot and zero point plot.
    * Added experimental `is_blank` method to try to guess whether image is blank.  Can be used to determine if time should be spent analyzing image.
    * Bug fixes
* __v1.2__
    * New configuration scheme.  The telescope configuration scheme is read from a YAML config file.  The config file also controls output directories for the logs, plots, and tmp files.
    * Flags if the image FWHM, ellipticity, pointing error, or zero point are above their respective threshold values are now recorded in the YAML file and are used to color code the HTML output.
    * New `add_yaml_entry` method (intended to replace `add_summary_entry`) stores image analysis results in YAML file.
    * Improvements to the PSF plots to show histograms and spatial distribution of high ellipticity and high FWHM stars.
    * The `make_JPEG` method can now flag saturated pixels
    * Image histogram is now optional output of `make_JPEG` method.
    * Numerous bug fixes.
* __v1.1__
    * Added ability to solve astrometry using SCAMP and rectify the image using SWarp
    * Added ability to solve for Zero Point using the SCAMP-solved data
    * Renamed many methods to better follow PEP8 conventions
    * Refactored `make_JPEG` to use matplotlib and PIL instead of ImageMagick.  The previous method is available as `make_JPEG_ImageMagick`.  This removes the 5000 star limit from marking, but the image annotations are no longer available.
* __v1.0.6__
    * Tweaks related to marking up of jpeg files.  Better markers, user settable maker size, and marks are labeled.
* __v1.0.5__
    * Bug fixes related to marking up of jpeg files.
* __v1.0.4__
    * MakeJPEG now marks the brightest 5000 stars rather than the first 5000 in the table.  Also annotates image to let viewer know more stars were detected.
    * added option to HTML output to choose which columns are displayed
    * other minor bug fixes and tweaks
    * fixed bug where axes in crop were reversed
* __v1.0.3__
    * Fixed handling of spaces in filenames by repalcing spaces in input file name with underscores when making working file.
    * Fixed bug in MakeJPEG which hard coded a binning of 2x2.
    * Added flag to ignore missing end card in fits header after running astrometry.net.
* __v1.0.2__
    * Fixed bug in jpeg creation which would crash when many (>~5000) stars were to be marked.  Now marks first 5000 stars.
    * Fixed bug in color coding of HTML.
* __v1.0.1__
    * Added astropy units converter between pixels and arcseconds to handle internal conversion of pixel scale.
    * Fixed bug in jpeg creation.  Will now handle rotation and marked stars.
    * Fixed bug when reading in summary text file.
* __v1.0__ (released on github.com 2013/08/14)
    * Rewritten as object oriented code.  Implements most capabilities of v0.X.
    * Runs roughly 2x faster than v0.X.
* __v0.X__ (frozen 2013/07/15)
    * Initial version, not under version control.
    * Deployed on live VYSOS data.
    * Analyzes image with SExtractor
        * Reports FWHM, ellipticity, background, number of stars detected
    * Reports image pointing info (alt, az, moon angle)
    * Can solve image with astrometry.net (with poor error handling)
    * Reports pointing error and position angle
    * Makes jpegs of image and cropped version with circles overlayed on stars found by SExtractor
    * Makes HTML and text file versions of results with one line per image, usually one night of images per file.
