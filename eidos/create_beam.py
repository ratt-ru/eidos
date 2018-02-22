#!/usr/bin/env python
from zernike import *
from spectrum import *
from util import *
import argparse
import sys
import os

def main(argv):
    parser=argparse.ArgumentParser(description='Create primary beam model of MeerKAT')
    parser.add_argument('-p', '--pixels', help='Number of pixels on one side', type=int, required=True)
    parser.add_argument('-d', '--diameter', help='Diameter of the required beam', default=6., type=float, required=False)
    parser.add_argument('-f', '--freq', help='A single freq, or the start, end freqs, and channel width in MHz', nargs='+', type=float, required=True)
    parser.add_argument('-c', '--coefficients-file', help='Coefficients file name', type=str, required=True)
    parser.add_argument('-P', '--prefix', help='Prefix of output beam beam file(s)', type=str, required=False)
    parser.add_argument('-o8', '--output-eight', help='Output complex volatge beams (8 files)', action='store_true')

    args = parser.parse_args(argv)

    if len(args.freq)==1: 
        nu = float(args.freq[0])
    elif len(args.freq)==2: 
        nu = np.arange(args.freq[0], args.freq[1], 1)
    elif len(args.freq)==3: 
        nu = np.arange(args.freq[0], args.freq[1], args.freq[2])

    coeffs = np.load(args.coefficients_file).item()
    
    mod = Zernike(coeffs, mode='recon', thresh=[15,8], Npix=args.pixels, freq=nu)

    if args.prefix:
        filename = args.prefix
    else:
        filename = 'meerkat_pb_jones_cube_%ichannels'%len(nu)

    if isinstance(nu, (int, float)):
        data = mod.recons
        write_fits_single(data, nu, args.diameter, filename)
        print 'Saved as %s'%filename
    else:
        if args.output_eight:
            data = mod.recons_all
            write_fits_eight(data, nu, args.diameter, filename)
            print "Saved as 8 files with prefix %s"%filename
        else:
            data = mod.recons_all
            write_fits(data, nu, args.diameter, filename)
            print "Saved as %s"%filename
