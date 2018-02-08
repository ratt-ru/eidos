"""
Some utility functions to be used by other scripts
"""
import numpy as np
from numpy import pi, exp, sin, cos
from astropy.io import fits
import time

def write_fits(beam, freqs, diameter, filename):
    data = np.zeros((2,)+beam.shape)
    data[0,...] = beam.real
    data[1,...] = beam.imag
    
    # Create header
    hdr = fits.Header()
    fMHz = np.array(freqs)*1e6
    diam = float(diameter)
    ctypes = ['X', 'Y', 'FREQ', 'FEED1', 'FEED2', 'PART']
    crvals = [0.0, 0.0, fMHz[0], 0, 0, 0]
    cdelts = [diam/beam.shape[-2], diam/beam.shape[-1], fMHz[1]-fMHz[0], 1, 1, 1]
    cunits = ['deg', 'deg', 'Hz', '', '', '']
    nx, ny = beam.shape[-2], beam.shape[-1]
    if nx%2 == 0: crpixx, crpixy = nx/2+0.5, ny/2+0.5
    elif nx%2 == 1: crpixx, crpixy = nx/2+1, ny/2+1
    crpixs = [crpixx, crpixy, 1, 1, 1, 1]
    for i in range(len(data.shape)):
        ii = str(i+1)
        hdr['CTYPE'+ii] = ctypes[i]
        hdr['CRPIX'+ii] = crpixs[i]
        hdr['CRVAL'+ii] = crvals[i]
        hdr['CDELT'+ii] = cdelts[i]
        hdr['CUNIT'+ii] = cunits[i]
    hdr['TELESCOP'] = 'MeerKAT'
    hdr['DATE'] = time.ctime()
    
    # Write real and imag parts of data
    hdu = fits.PrimaryHDU(data, header=hdr)
    hdu.writeto(filename, overwrite=True)

def write_fits_single(beam, freqs, diameter, filename):
    data = np.zeros((2,)+beam.shape)
    data[0,...] = beam.real
    data[1,...] = beam.imag
    
    # Create header
    hdr = fits.Header()
    fMHz = np.array(freqs)*1e6
    diam = float(diameter)
    ctypes = ['X', 'Y', 'FEED1', 'FEED2', 'PART']
    crvals = [0.0, 0.0, 0, 0, 0]
    cdelts = [diam/beam.shape[-2], diam/beam.shape[-1], 1, 1, 1]
    cunits = ['deg', 'deg', '', '', '']
    nx, ny = beam.shape[-2], beam.shape[-1]
    if nx%2 == 0: crpixx, crpixy = nx/2+0.5, ny/2+0.5
    elif nx%2 == 1: crpixx, crpixy = nx/2+1, ny/2+1
    crpixs = [crpixx, crpixy, 1, 1, 1]
    for i in range(len(data.shape)):
        ii = str(i+1)
        hdr['CTYPE'+ii] = ctypes[i]
        hdr['CRPIX'+ii] = crpixs[i]
        hdr['CRVAL'+ii] = crvals[i]
        hdr['CDELT'+ii] = cdelts[i]
        hdr['CUNIT'+ii] = cunits[i]
    hdr['TELESCOP'] = 'MeerKAT'
    hdr['DATE'] = time.ctime()
    
    # Write real and imag parts of data
    hdu = fits.PrimaryHDU(data, header=hdr)
    hdu.writeto(filename, overwrite=True)

def split_into_eight(filename):
    d = fits.getdata(filename)
    d = np.nan_to_num(d)
    h = fits.getheader(filename)
    del h['C*4']
    del h['C*5']
    del h['C*6']
    P, C = ['R', 'I'], ['X', 'Y']
    for p in range(2):
        for i in range(2):
            for j in range(2):
                fits.writeto(filename[:-5]+'_%s%s_%s.fits'%(C[i],C[j],P[p]), d[p,i,j,...], h, overwrite=True)

def freq_to_idx(freqs=np.arange(857)+856, sb=[0,1300]):
    # convert a frequency to an index, default is MeerKAT L-band specifications
    start = (np.abs(freqs-sb[0])).argmin()
    end = (np.abs(freqs-sb[1])).argmin()
    return start, end
