#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE

import numpy as np

from astropy.io import fits

c = 299792458   # [m/s] speed of light


def FTSfits(ftsname):

    if ftsname.endswith(".dat"):
        data = np.loadtxt(ftsname)
        w = data[:, 0]
        f = data[:, 1]
        f = f[::-1]
        w = 1e8 / w[::-1]
    elif ftsname.endswith(".fits"):
        hdu = fits.open(ftsname, ignore_blank=True, output_verify='silentfix')
      
        hdr = hdu[0].header        
        cdelt1 = hdr.get('CDELT1', 'none')
        
        if cdelt1 == 'none':
            wavetype = hdr.get('wavetype', 'none')
            unit = hdr.get('unit', 'none')
            w = hdu[1].data['wave']
            f = hdu[1].data['flux']#[::-1]
            
            if wavetype == 'wavenumber':  w = 1e8 / w[::-1]
            if unit == 'nm': w *= 10
        
        else:
            f = hdu[0].data[::-1]
            try:
                w = hdr['CRVAL1'] + hdr['CDELT1'] * (np.arange(f.size) + 1. - hdr['CRPIX1'])
            except:
                # for OES
                w = hdr['CRVAL1'] + hdr['CDELT1'] * (np.arange(f.size) + 1.)
            w = 1e8 / w[::-1]   # convert wavenumbers to wavelength [angstrom]

    return w, f


def resample(w, f, dv=100):
    '''
    dv: Sampling step for uniform log(lambda) [m/s]
    '''
    # define a supersampled log(wavelength) space with knot index j
    u = np.log(w)
    uj = np.arange(u[0], u[-1], dv/c)
    iod_j = np.interp(uj, u, f)

    return w, f, uj, iod_j


