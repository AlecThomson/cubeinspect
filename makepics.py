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



with warnings.catch_warnings():
    warnings.simplefilter('ignore', astropy.wcs.FITSFixedWarning)
    warnings.simplefilter('ignore', category=AstropyWarning)


def fitsopen(filename, chan, verbose=True):
    if verbose:
        print('Opening', filename)
    hdulist = fits.open(filename, memmap=True, mode='denywrite')
    hdu = hdulist[0]
    head = hdu.header
    size1 = head['NAXIS1']
    size2 = head['NAXIS2']
    sniff = np.std(np.squeeze(hdu.data)[chan, size2 // 2 - 10:size2 // 2 + 10, size1 // 2 - 10:size1 // 2 + 10])
    if ~(sniff > 0) or ~np.isfinite(sniff):
        chan += 10
        if verbose:
            print('Bad channel! Trying N+10')
    else:
        if verbose:
            print('Good channel!')
    data = np.squeeze(hdu.data)
    data = data[chan]
    return data, head


def makeplot(filename, data, head, verbose=True, thumbnail=False, lims=None):
    if thumbnail:
    # Still in-progress
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
        fig.savefig(outfile, dpi=500)
        print('Saved to', outfile)
    if verbose:
        print('Making plot')
    image = np.power(data, 2)
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
    fig.savefig(outfile, dpi=500)
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
    parser.add_argument("-v", dest="verbose", default=True,
                        action="store_true", help="Verbosity.")
    parser.add_argument("-t", dest="thumbnail", default=False,
                        action="store_true", help="Make thumbnail (takes extra time!)")
    parser.add_argument("-c", metavar=('vmin', 'vmax'), dest="lims", type=float, nargs=2, default=None, help="Limits for image (defaults to vmax=std(image))")


    args = parser.parse_args()
    data, head = fitsopen(args.fitsfile, args.chan, args.verbose)
    makeplot(args.fitsfile, data, head, verbose=args.verbose, thumbnail=args.thumbnail, lims=args.lims)


if __name__ == "__main__":
    main()
