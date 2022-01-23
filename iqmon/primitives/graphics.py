from pathlib import Path
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from matplotlib.dates import HourLocator, MinuteLocator, DateFormatter
plt.style.use('classic')
import pymongo

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import stats

from keckdrpframework.primitives.base_primitive import BasePrimitive

from .utils import pre_condition, post_condition, get_memory_size, get_jpeg_dir



##-----------------------------------------------------------------------------
## Function: overlay_extracted_stars
##-----------------------------------------------------------------------------
def overlay_extracted_stars(ax, cfg, meta, objects, x0, nx, y0, ny):
    FWHM_pix = meta.get('fwhm', 8)
    binning = cfg['jpeg'].getint('binning', 1)
    ap_radius = cfg['Photometry'].getfloat('aperture_radius', 2)*FWHM_pix
    marker_size = cfg['jpeg'].getfloat('marker_size', 2)*FWHM_pix
    radius = marker_size/binning
    for star in objects:
        if star['x'] > x0 and star['x'] < x0+nx and star['y'] > y0 and star['y'] < y0+ny:
            x = star['x']-x0
            y = star['y']-y0
            ax.plot([x, x], [y+radius, y+2.5*radius], 'b', alpha=0.3)
            ax.plot([x, x], [y-radius, y-2.5*radius], 'b', alpha=0.3)
            ax.plot([x-radius, x-2.5*radius], [y, y], 'b', alpha=0.3)
            ax.plot([x+radius, x+2.5*radius], [y, y], 'b', alpha=0.3)


##-----------------------------------------------------------------------------
## Function: overlay_calibrator_stars
##-----------------------------------------------------------------------------
def overlay_calibrator_stars(ax, cfg, meta, calibrators, x0, nx, y0, ny):
    FWHM_pix = meta.get('fwhm', 8)
    binning = cfg['jpeg'].getint('binning', 1)
    ap_radius = cfg['Photometry'].getfloat('aperture_radius', 2)*FWHM_pix
    radius = ap_radius/binning
    for entry in calibrators:
        if entry['apflux'] > 0:
            xy = (entry['x'], entry['y'])
            ax.add_artist(plt.Circle(xy, radius=radius,
                                     edgecolor='r', facecolor='none', alpha=0.5))
            ax.add_artist(plt.Circle(xy, radius=int(np.ceil(1.5*radius)),
                                     edgecolor='r', facecolor='none', alpha=0.5))
            ax.add_artist(plt.Circle(xy, radius=int(np.ceil(2.0*radius)),
                                     edgecolor='r', facecolor='none', alpha=0.5))


##-----------------------------------------------------------------------------
## Function: overlay_catalog_stars
##-----------------------------------------------------------------------------
def overlay_catalog_stars(ax, cfg, meta, imagewcs, catalog, x0, nx, y0, ny):
    FWHM_pix = meta.get('fwhm', 8)
    binning = cfg['jpeg'].getint('binning', 1)
    marker_size = cfg['jpeg'].getfloat('marker_size', 2)*FWHM_pix
    radius = marker_size/binning
    for entry in catalog:
        x, y = imagewcs.all_world2pix(entry['raMean'], entry['decMean'], 1)
        ax.plot([x, x], [y+radius, y+2.5*radius], 'g', alpha=0.7)
        ax.plot([x, x], [y-radius, y-2.5*radius], 'g', alpha=0.7)
        ax.plot([x-radius, x-2.5*radius], [y, y], 'g', alpha=0.7)
        ax.plot([x+radius, x+2.5*radius], [y, y], 'g', alpha=0.7)


##-----------------------------------------------------------------------------
## Function: overlay_pointing
##-----------------------------------------------------------------------------
def overlay_pointing(ax, cfg, meta, imagewcs, header_pointing, x0, nx, y0, ny):
    radius = cfg['jpeg'].getfloat('pointing_radius', 40)
    x, y = imagewcs.all_world2pix(header_pointing.ra.degree,
                                  header_pointing.dec.degree, 1)
    ax.plot([nx/2-3*radius,nx/2+3*radius], [ny/2,ny/2], 'k-', alpha=0.5)
    ax.plot([nx/2, nx/2], [ny/2-3*radius,ny/2+3*radius], 'k-', alpha=0.5)
    # Draw crosshair on target
    ax.add_artist(plt.Circle((x, y), radius=radius, edgecolor='g', alpha=0.7,
                             facecolor='none'))
    ax.plot([x, x], [y+0.6*radius, y+1.4*radius], 'g', alpha=0.7)
    ax.plot([x, x], [y-0.6*radius, y-1.4*radius], 'g', alpha=0.7)
    ax.plot([x-0.6*radius, x-1.4*radius], [y, y], 'g', alpha=0.7)
    ax.plot([x+0.6*radius, x+1.4*radius], [y, y], 'g', alpha=0.7)


##-----------------------------------------------------------------------------
## Primitive: RenderJPEG
##-----------------------------------------------------------------------------
class RenderJPEG(BasePrimitive):
    """
    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.log = context.pipeline_logger
        self.cfg = self.context.config.instrument

    def _pre_condition(self):
        """Check for conditions necessary to run this process"""
        checks = [pre_condition(self, 'Skip image is not set',
                                not self.action.args.skip),
                 ]
        return np.all(checks)

    def _post_condition(self):
        """
        Check for conditions necessary to verify that the process ran
        correctly.
        """
        checks = []
        return np.all(checks)

    def _perform(self):
        """
        Returns an Argument() with the parameters that depend on this
        operation.
        """
        self.log.info(f"Running {self.__class__.__name__} action")

        im = self.action.args.ccddata.data

        plt.rcParams.update({'font.size': self.cfg['jpeg'].getint('font_size', 12)})
        binning = self.cfg['jpeg'].getint('binning', 1)
        vmin = np.percentile(im, self.cfg['jpeg'].getfloat('vmin_percent', 0.5))
        vmax = np.percentile(im, self.cfg['jpeg'].getfloat('vmax_percent', 99))
        dpi = self.cfg['jpeg'].getint('dpi', 72)
        nx, ny = im.shape
        sx = nx/dpi/binning
        sy = ny/dpi/binning
        FWHM_pix = self.action.args.meta.get('fwhm', 8)
        ap_radius = self.cfg['Photometry'].getfloat('aperture_radius', 2)*FWHM_pix
        marker_size = self.cfg['jpeg'].getfloat('marker_size', 2)*FWHM_pix

        if self.action.args.wcs is not None:
            pixel_scale = np.mean(proj_plane_pixel_scales(self.action.args.wcs))*60*60
        else:
            pixel_scale = self.cfg['Telescope'].getfloat('pixel_scale', 1)

        fig = plt.figure(figsize=(2*sx, 1*sy), dpi=dpi)

        plotpos = [ [ [0.010, 0.010, 0.550, 0.965], [0.560, 0.475, 0.380, 0.500] ],
                    [ None                        , [0.560, 0.260, 0.380, 0.200] ],
                    [ None                        , [0.560, 0.020, 0.380, 0.200] ],
                  ]
#         plotpos = [ [ [0.010, 0.010, 0.550, 0.965], [0.565, 0.775, 0.375, 0.200] ],
#                     [ None                        , [0.565, 0.550, 0.375, 0.200] ],
#                     [ None                        , [0.565, 0.010, 0.375, 0.500] ],
#                   ]

        ##-------------------------------------------------------------------------
        # Show JPEG of Image
        self.log.debug(f'  Rendering JPEG of image')
        jpeg_axes = plt.axes(plotpos[0][0])
        jpeg_axes.imshow(im, cmap=plt.cm.gray_r, vmin=vmin, vmax=vmax, origin='lower')
        jpeg_axes.set_xticks([])
        jpeg_axes.set_yticks([])
        titlestr = f'{self.action.args.meta.get("fitsfile")}: '

        # Overlay Extracted (blue)
        if self.cfg['jpeg'].getboolean('overplot_extracted', False) is True\
                and self.action.args.objects is not None:
            self.log.debug(f'  Overlay extracted stars')
            titlestr += f'blue=extracted({len(self.action.args.objects)}) '
            overlay_extracted_stars(jpeg_axes,
                                    self.cfg,
                                    self.action.args.meta,
                                    self.action.args.objects,
                                    0, nx, 0, ny)
            self.log.debug(f'  Done')

        # Overlay Calibrators (red)
        if self.cfg['jpeg'].getboolean('overplot_calibrators', False) is True\
                and self.action.args.associated_calibrators is not None\
                and self.action.args.wcs is not None:
            self.log.debug(f'  Overlay measured calibration stars')
            titlestr += f'red=calibrators({len(calibrators)}) '
            overlay_calibrator_stars(jpeg_axes,
                                     self.cfg,
                                     self.action.args.meta,
                                     self.action.args.associated_calibrators,
                                     0, nx, 0, ny)
            self.log.debug(f'  Done')

        # Overlay Catalog (green)
        if self.cfg['jpeg'].getboolean('overplot_catalog', False) is True\
                and self.action.args.calibration_catalog is not None\
                and self.action.args.wcs is not None:
            self.log.debug(f'  Overlay stars from catalog')
            titlestr += f'green=catalog({len(catalog)}) '
            overlay_catalog_stars(jpeg_axes,
                                  self.cfg,
                                  self.action.args.meta,
                                  self.action.args.wcs,
                                  self.action.args.calibration_catalog,
                                  0, nx, 0, ny)
            self.log.debug(f'  Done')

        # Overlay Pointing
        if self.cfg['jpeg'].getboolean('overplot_pointing', False) is True\
                and self.action.args.header_pointing is not None\
                and self.action.args.wcs is not None\
                and self.action.args.wcs_pointing is not None:
            self.log.debug(f'  Overlay pointing')
            overlay_pointing(jpeg_axes,
                             self.cfg,
                             self.action.args.meta,
                             self.action.args.wcs,
                             self.action.args.header_pointing,
                             0, nx, 0, ny)
            self.log.debug(f'  Done')

        # Wrap Up Full Image Display
        jpeg_axes.set_xlim(0,nx)
        jpeg_axes.set_ylim(0,ny)

        ##-------------------------------------------------------------------------
        # Show High Resolution JPEG of Center of Image
        self.log.debug(f'  Rendering central image')
        jpeg_axes2 = plt.axes(plotpos[0][1])
        jpeg_axes2.set_xticks([])
        jpeg_axes2.set_yticks([])
        dx = 1000
        dy = 650
        x0 = int(im.shape[1] / 2) - int(dx/2)
        y0 = int(im.shape[0] / 2) - int(dy/2)
        central_im = im[y0:y0+dy,x0:x0+dx]
        jpeg_axes2.set_title(f"Central {dx:d}x{dy:d} pixels")
        jpeg_axes2.imshow(central_im, cmap=plt.cm.gray_r, vmin=vmin, vmax=vmax, origin='lower')
        self.log.debug(f'  Done')

        # Overlay Extracted (blue)
        if self.cfg['jpeg'].getboolean('overplot_extracted', False) is True\
                and self.action.args.objects is not None:
            self.log.debug(f'  Overlay extracted stars')
            titlestr += f'blue=extracted({len(self.action.args.objects)}) '
            overlay_extracted_stars(jpeg_axes2,
                                    self.cfg,
                                    self.action.args.meta,
                                    self.action.args.objects,
                                    x0, dx, y0, dy)
            self.log.debug(f'  Done')

        # Overlay Calibrators (red)
        if self.cfg['jpeg'].getboolean('overplot_calibrators', False) is True\
                and self.action.args.associated_calibrators is not None\
                and self.action.args.wcs is not None:
            self.log.debug(f'  Overlay measured calibration stars')
            titlestr += f'red=calibrators({len(calibrators)}) '
            overlay_calibrator_stars(jpeg_axes2,
                                     self.cfg,
                                     self.action.args.meta,
                                     self.action.args.associated_calibrators,
                                     x0, dx, y0, dy)
            self.log.debug(f'  Done')

        # Wrap Up Central Image Display
        jpeg_axes2.set_xlim(0,dx)
        jpeg_axes2.set_ylim(0,dy)


        ##-------------------------------------------------------------------------
        ## If no stellar data is available
        ## plot histogram of pixel values
        if self.action.args.calibration_catalog is None and self.action.args.objects is None:
            self.log.debug(f'  Generating histogram of pixel values')
            pixel_axes = plt.axes(plotpos[1][1])
            mean, med, std = stats.sigma_clipped_stats(im)
            p1, p99 = np.percentile(im, 1), np.percentile(im, 99)
            pixel_axes.set_title(f"Histogram of Pixel Values (median = {med:.0f})")
            npix, bins, p = pixel_axes.hist(im.ravel(), color='b', alpha=0.5,
                                            bins=np.linspace(p1,p99,50))
            pixel_axes.plot([med, med], [0,max(npix)*1.2], 'r-', alpha=0.5)
            pixel_axes.set_xlabel('e-/s')
            pixel_axes.set_ylabel('N Pix')
            pixel_axes.set_ylim(0,max(npix)*1.2)

        ##-------------------------------------------------------------------------
        ## If extracted stars, but no photometry available
        if self.action.args.objects is not None and self.action.args.calibration_catalog is None:
            # Plot histogram of FWHM
            self.log.debug(f'  Generating histogram of FWHM values')
            fwhm_axes = plt.axes(plotpos[1][1])
            avg_fwhm = self.action.args.meta.get('fwhm')*pixel_scale
            minfwhm = np.percentile(self.action.args.objects['FWHM'], 0.1)*pixel_scale
            maxfwhm = np.percentile(self.action.args.objects['FWHM'], 90)*pixel_scale
            fwhm_axes.set_title(f"FWHM = {avg_fwhm:.1f} arcsec")
            nstars, bins, p = fwhm_axes.hist(self.action.args.objects['FWHM']*pixel_scale,
                                             bins=np.arange(minfwhm,maxfwhm,0.25*pixel_scale),
                                             color='b', alpha=0.5)
            fwhm_axes.plot([avg_fwhm, avg_fwhm], [0,max(nstars)*1.2], 'r-', alpha=0.5)
            fwhm_axes.set_ylabel('N stars')
            fwhm_axes.set_ylim(0,max(nstars)*1.2)
            fwhm_axes.set_xlabel(f'FWHM (arcsec) [1 pix = {pixel_scale:.2f} arcsec]')
            if self.action.args.associated_calibrators is not None:
                nstars, bins, p = fwhm_axes.hist(self.action.args.associated_calibrators['FWHM']*pixel_scale,
                                                 bins=np.arange(minfwhm,maxfwhm,0.25*pixel_scale),
                                                 color='r', alpha=0.5)
                fwhm_axes.plot([avg_fwhm, avg_fwhm], [0,max(nstars)*1.2], 'r-', alpha=0.5)
            self.log.debug(f'  Done')

            # Plot histogram of ellipticity
            self.log.debug(f'  Generating histogram of ellipticity values')
            elip_axes = plt.axes(plotpos[2][1])
            avg_elip = self.action.args.meta.get('ellipticity')
            elip_axes.set_title(f"ellipticity = {avg_elip:.2f} ")
            nstars, bins, p = elip_axes.hist(self.action.args.objects['ellipticity'],
                                             bins=np.arange(0.9,2.5,0.05),
                                             color='b', alpha=0.5)
            elip_axes.plot([avg_elip, avg_elip], [0,max(nstars)*1.2], 'r-', alpha=0.5)
            elip_axes.set_ylabel('N stars')
            elip_axes.set_ylim(0,max(nstars)*1.2)
            elip_axes.set_xlabel('ellipticity')
            if self.action.args.associated_calibrators is not None:
                nstars, bins, p = elip_axes.hist(self.action.args.associated_calibrators['ellipticity'],
                                                 bins=np.arange(0.9,2.5,0.05),
                                                 color='r', alpha=0.5)
            self.log.debug(f'  Done')

        ##-------------------------------------------------------------------------
        ## If photometry available
        if self.action.args.associated_calibrators is not None\
            and self.action.args.zero_point_f0 is not None:
            # Plot instrumental mags
            flux_axes = plt.axes(plotpos[1][1])
            self.log.debug(f'  Generating plot of flux vs catalog magnitude')
            flux_axes.plot(self.action.args.associated_calibrators['mag'],
                           self.action.args.associated_calibrators['apflux'], 'bo',
                           label='photutils', mec=None, alpha=0.6)
            flux_axes.set_xlim(min(self.action.args.associated_calibrators['mag']),
                               max(self.action.args.associated_calibrators['mag']))
            flux_axes.set_yscale('log')
            flux_axes.set_ylabel('Measured Flux (e-/s)')
            plt.grid()
            plt.legend(loc='best')
            self.log.debug(f'  Done')

            # Plot instrumental mag diffs
            diff_axes = plt.axes(plotpos[1][1])
            self.log.debug(f'  Generating plot of flux residual')
            deltamag = self.action.args.associated_calibrators['mag']\
                       - self.action.args.associated_calibrators['instmag']
            mean, med, std = stats.sigma_clipped_stats(deltamag)
            xmin = min(self.action.args.associated_calibrators['mag'])
            xmax = max(self.action.args.associated_calibrators['mag'])
            diff_axes.plot(self.action.args.associated_calibrators['mag'],
                           deltamag, 'bo',
                           label=f'Individual Zero Points',
                           mec=None, alpha=0.6)
            zp = self.action.args.zero_point
            label = f'Zero Point = {zp:.2f} (throughput={self.action.args.throughput:.3f})'
            diff_axes.plot([xmin, xmax], [zp,zp], 'k-', mec=None, alpha=0.6,
                           label=label)
            diff_axes.plot([xmin, xmax], [zp+0.1,zp+0.1], 'r-', mec=None, alpha=0.6)
            diff_axes.plot([xmin, xmax], [zp-0.1,zp-0.1], 'r-', mec=None, alpha=0.6,
                           label=f"+/-0.1 magnitude error (StdDev={std:.2f})")
            diff_axes.set_xlabel('Catalog Magnitude')
            diff_axes.set_xlim(xmin, xmax)
            diff_axes.set_ylabel('Catalog - Instrumental Mag')
            diff_axes.set_ylim(np.percentile(deltamag, 1)-0.25, np.percentile(deltamag, 99)+0.25)
            plt.legend(loc='best')
            plt.grid()
            self.log.debug(f'  Done')


        self.log.debug(f'  Saving JPEG file')
        jpeg_axes.set_title(titlestr)
        instrument = self.cfg['Telescope'].get('name')
        jpeg_filename = f'{self.action.args.meta.get("fitsfile").split(".")[0]}.jpg'
#         jpeg_file_path = Path('/var/www/plots/') / instrument / jpeg_filename
        jpeg_file_dir = Path(get_jpeg_dir(self.cfg, self.action.args.meta['date']))
        jpeg_file_dir.mkdir(parents=True, exist_ok=True)
        jpeg_file_path = jpeg_file_dir / jpeg_filename
        self.action.args.meta['jpegfile'] = '/'.join(jpeg_file_path.parts[-2:])
        plt.savefig(jpeg_file_path, dpi=dpi)
        self.log.debug(f'  Done')

#         self.log.info(f"Memory Size after {self.__class__.__name__} action: {get_memory_size(self):.1f} MB")
        return self.action.args

