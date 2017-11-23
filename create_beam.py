#!/usr/bin/env python
from zernike import *
from spectrum import *
from util import *
import argparse, sys

def main():
    parser=argparse.ArgumentParser(description='Create primary beam model of MeerKAT')
    parser.add_argument('-p', help='Number of pixels on one side', type=int, required=True)
    parser.add_argument('-f', help='A single freq, or the start, end freqs, and channel width in MHz', nargs='+', type=float, required=True)
    parser.add_argument('-d', help='Diameter of the required beam', default=6., type=float, required=False)
    parser.add_argument('-c', help='Coefficients file name', type=str, default='meerkat_coeff_dict.npy', required=False)
    args=parser.parse_args()

    if len(args.f)==1: nu = float(args.f[0])
    elif len(args.f)==2: nu = np.arange(args.f[0], args.f[1], 1)
    elif len(args.f)==3: nu = np.arange(args.f[0], args.f[1], args.f[2])

    coeffs = np.load(args.c).item()
    
    mod = Zernike(coeffs, mode='recon', thresh=[15,8], Npix=args.p, freq=nu)

    if isinstance(nu, (int, float)):
        data = mod.recons
        filename = 'meerkat_pb_jones_%iMHz.fits'%int(nu)
        write_fits_single(data, nu, args.d, filename)
        print 'Saved as %s'%filename
    else:
        data = mod.recons_all
        filename = 'meerkat_pb_jones_cube_%ichannels.fits'%len(nu)
        write_fits(data, nu, args.d, filename)
        print "Saved as %s"%filename

if __name__=='__main__':
    main()
