This version of IQMon uses the Keck Data Reduction Framework (KDRPF; https://keckdrpframework.readthedocs.io/en/latest/)

```
> db.images.findOne( {'telescope': 'V5'} )
{
    "_id" : ObjectId("6161346bb780ebdf55874858"),
    "filename" : "V5_AS353A-PSr3-20211009at050106.fts",
    "telescope" : "V5",
    "compressed" : false,
    "target name" : "V5_AS353A",
    "exptime" : "100.0",
    "date" : ISODate("2021-10-09T05:01:11Z"),
    "az" : 216.0082,
    "alt" : 78.9609,
    "airmass" : 1.01918938497,
    "header_RA" : 291.18749999999994,
    "header_DEC" : 10.5,
    "moon_alt" : 14.023852963855044,
    "moon_separation" : 66.63240318988049,
    "moon_illumination" : 9.96456815169346,
    "FWHM_pix" : 5.006289586221885,
    "ellipticity" : 1.2798003112449319,
    "n_stars" : 5474,
    "zero point" : 22.74422572531131,
    "throughput" : 7.611131294360541,
    "sky background" : 2372.3722751214323,
    "perr_arcmin" : 52.13182009729838,
    "wcs" : "WCSAXES =                    2 / Number of coordinate axes                      CRPIX1  =        3003.41778564 / Pixel coordinate of reference point            CRPIX2  =        1310.08117676 / Pixel coordinate of reference point            PC1_1   =   -4.69284792725E-05 / Coordinate transformation matrix element       PC1_2   =    0.000695371629704 / Coordinate transformation matrix element       PC2_1   =    -0.00069458536881 / Coordinate transformation matrix element       PC2_2   =   -4.15398212941E-05 / Coordinate transformation matrix element       CDELT1  =                  1.0 / [deg] Coordinate increment at reference point  CDELT2  =                  1.0 / [deg] Coordinate increment at reference point  CUNIT1  = 'deg'                / Units of coordinate increment and value        CUNIT2  = 'deg'                / Units of coordinate increment and value        CTYPE1  = 'RA---TAN'           / TAN (gnomonic) projection + SIP distortions    CTYPE2  = 'DEC--TAN'           / TAN (gnomonic) projection + SIP distortions    CRVAL1  =         291.50558067 / [deg] Coordinate value at reference point      CRVAL2  =        9.81502190238 / [deg] Coordinate value at reference point      LONPOLE =                180.0 / [deg] Native longitude of celestial pole       LATPOLE =        9.81502190238 / [deg] Native latitude of celestial pole        RADESYS = 'FK5'                / Equatorial coordinate system                   EQUINOX =               2000.0 / [yr] Equinox of equatorial coordinates         END",
    "jpegs" : [
        "V5_AS353A-PSr3-20211009at050106.jpg"
    ]
}
```
