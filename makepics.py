#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function

from astropy.io import fits
from astropy.wcs import WCS
import astropy
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.colors import LogNorm
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.convolution import convolve, convolve_fft, Box2DKernel
import numexpr as ne


with warnings.catch_warnings():
    warnings.simplefilter('ignore', astropy.wcs.FITSFixedWarning)
    warnings.simplefilter('ignore', astropy.io.fits.hdu.image.VerifyWarning)
    warnings.simplefilter('ignore', category=AstropyWarning)


def fitsopen(filename, chan, verbose=True):
    """
    Open the FITS file. Central 10x10 pixels are inspected for NaNs or all 0s (in case of flagging). If the selected channel is bad, the channel 10 pixels along is selected.

    Args:
        filename: Name of the file (str).
        chan: Channel index (int).
        verbose: Verbosity (bool).

    Returns:
        data: Astropy array of containing selected channel image.
        head: Astropy FITS header.

    """
    if verbose:
        print('Opening', filename)
    hdulist = fits.open(filename, memmap=True, mode='denywrite')
    hdu = hdulist[0]
    head = hdu.header
    size1 = head['NAXIS1']
    size2 = head['NAXIS2']
    size3 = head['NAXIS3']

    if chan > size3:
        raise Exception('Channel outside of range!')

    # Sniff around middle of image
    sniff = np.std(np.squeeze(hdu.data)[chan, size2 // 2 - 10:size2 // 2 + 10, size1 // 2 - 10:size1 // 2 + 10])
    if ~(sniff > 0) or ~np.isfinite(sniff):
        if chan > 0:
            chan += 10
        if chan < 0:
            chan -= 10
        if verbose:
            print('Bad channel! Trying N+10')

    if chan > size3:
        raise Exception('Channel outside of range!')

    else:
        if verbose:
            print('Good channel!')

    # Check for Stokes axis
    if head['NAXIS'] > 3:
        if verbose:
            print('Stokes axis present - [grumbles] - squishing it out!')
        data = np.squeeze(hdu.data)

    else:
        data = hdu.data
    data = ne.evaluate('sum(data**2, axis=0)')
    return data, head


def makeplot(filename, data, head, verbose=True, thumbnail=False, lims=None, dpi=1000):
    """
    Produce the plots. Will open an interactive matplotlib plot for inspection. After closing the image will be saved to a PNG.

    Args:
        filename: Name of the file (str).
        data: Astropy array of containing selected channel image.
        head: Astropy fits header.
        verbose: Verbosity (bool).
        thumbnail: EXPERIMENTAL -- Convolve the image with a 2D box and produce image (bool).
        lims: The vmin and vmax of the image.

    Returns:
        None

    """

    # Make thumbnail image - work in progress
    if thumbnail:
        if verbose:
            print('Making thumbnail')

        kernel = Box2DKernel(head['NAXIS1'] / 1000)
        astropy_conv = convolve_fft(data, kernel, allow_huge=True)
        image = np.power(astropy_conv, 2)
        if lims is None:
            rms = np.nanstd(image)
        proj = WCS(head).dropaxis(2).dropaxis(2)

        fig = plt.figure()
        fig.set_size_inches(8, 8)
        ax = fig.add_subplot(111, projection=proj)
        if lims is None:
            im = ax.imshow(image, origin='lower', cmap='cubehelix_r', vmax=10 * rms)
        if lims is not None:
            im = ax.imshow(image, origin='lower', cmap='cubehelix_r', vmin=lims[0], vmax=lims[1])
        lon = ax.coords[0]
        lon.set_ticklabel(size=8)
        lon.set_axislabel(r'$RA$')
        lon.display_minor_ticks(True)
        lat = ax.coords[1]
        lat.set_ticklabel(size=8)
        lat.display_minor_ticks(True)
        lat.set_axislabel(r'$DEC$')
        c = plt.colorbar(im)
        c.set_label('$S^2$')
        c.ax.tick_params(length=3)
        plt.title('Thumbnail')
        plt.show()

        outfile = re.sub('.fits', '.thumbnail.png', filename)
        fig.savefig(outfile, dpi=dpi)
        if verbose:
            print('Saved to', outfile)

    # Make main image
    if verbose:
        print('Making plot')
    image = data
    if lims is None:
            rms = np.nanstd(image)
    proj = WCS(head).dropaxis(2).dropaxis(2)
    print()
    fig = plt.figure()
    fig.set_size_inches(8, 8)
    ax = fig.add_subplot(111, projection=proj)
    if lims is None:
        im = ax.imshow(image, origin='lower', cmap='cubehelix_r', vmax=10 * rms)
    if lims is not None:
        im = ax.imshow(image, origin='lower', cmap='cubehelix_r', vmin=lims[0], vmax=lims[1])
    lon = ax.coords[0]
    lon.set_ticklabel(size=8)
    lon.set_axislabel(head['CTYPE1'])
    lon.display_minor_ticks(True)
    lat = ax.coords[1]
    lat.set_ticklabel(size=8)
    lat.display_minor_ticks(True)
    lat.set_axislabel(head['CTYPE2'])
    c = plt.colorbar(im)
    c.set_label('$S^2$')
    c.ax.tick_params(length=3)
    plt.title('Full image')
    plt.show()

    outfile = re.sub('.fits', '.medimage.png', filename)
    fig.savefig(outfile, dpi=dpi)
    if verbose:
        print('Saved to', outfile)


def main():
    import argparse

    # Parse the command line options
    descStr = """
    Open big file and make an average image of a selected channel.
    Image will be squared so that Stokes Q and U look reasonable.

    Saves a png with same name as FITS image + .medimage.png

    """

    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('fitsfile', metavar='fitsfile',
                        help='FITS file to open')
    parser.add_argument("chan", metavar="channel", type=int,
                        default=0, help="Channel to inspect")
    parser.add_argument("-v", dest="verbose", default=False,
                        action="store_true", help="Verbosity.")
    parser.add_argument("-t", dest="thumbnail", default=False,
                        action="store_true", help="Make thumbnail (takes extra time!)")
    parser.add_argument("-c", metavar=('vmin', 'vmax'), dest="lims", type=float, nargs=2, default=None, help="Limits for image (defaults to vmax=std(image))")
    parser.add_argument("-d", dest="dpi", type=float, default=1000, help="Resolution of saved images in DPI.")


    args = parser.parse_args()
    data, head = fitsopen(args.fitsfile, args.chan, args.verbose)
    makeplot(args.fitsfile, data, head, verbose=args.verbose, thumbnail=args.thumbnail, lims=args.lims, dpi=args.dpi)


if __name__ == "__main__":
    main()
