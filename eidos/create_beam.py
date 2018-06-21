#!/usr/bin/env python
from zernike import *
from spectrum import *
from util import *
import argparse
import sys
import os

package_directory = os.path.dirname(__file__)

def main(argv):
    parser=argparse.ArgumentParser(description='Create primary beam model of MeerKAT')
    parser.add_argument('-p', '--pixels', help='Number of pixels on one side', type=int, required=True)
    parser.add_argument('-d', '--diameter', help='Diameter of the required beam', default=6., type=float, required=False)
    parser.add_argument('-f', '--freq', help='A single freq, or the start, end freqs, and channel width in MHz', nargs='+', type=float, required=True)
    parser.add_argument('-c', '--coefficients-file', help='Coefficients file name', type=str, \
        default=os.path.join(package_directory, "data", "meerkat_coeff_dict.npy") )
    parser.add_argument('-P', '--prefix', help='Prefix of output beam beam file(s)', type=str, required=False)
    parser.add_argument('-o8', '--output-eight', help='Output complex volatge beams (8 files)', action='store_true')
    parser.add_argument('-N', '--normalise', help='Normalise the E-Jones wrt central pixel to remove the bandpass', action='store_true')
    parser.add_argument('-S', '--Stokes', help="If given output will be in Stokes, i. e. Mueller, formalism, instead of the default Jones formalism. \
        Specify 'I', 'Q', 'U', 'V' for the Stokes beams, or 'M' to get the full Mueller matrix. \
        If you give 'IQ' instead, the leakage from Q to I will be provided, and so on for any other combination of Stokes parameters.", type=str, required=False)

    args = parser.parse_args(argv)

    if len(args.freq)==1:
        nu = float(args.freq[0])
    elif len(args.freq)==2: 
        nu = np.arange(args.freq[0], args.freq[1], 1)
    elif len(args.freq)==3: 
        nu = np.arange(args.freq[0], args.freq[1], args.freq[2])
    else: print "Do `eidos -h` to see how to input the frequency parameters"

    coeffs = np.load(args.coefficients_file).item()
    
    mod = Zernike(coeffs, mode='recon', thresh=[15,8], Npix=args.pixels, freq=nu)

    if args.prefix:
        filename = args.prefix
    else:
        try: chan = "%ichannels"%len(nu)
        except: chan = "%iMHz"%(int(nu))
        filename = 'meerkat_pb_cube_%s'%chan
    
    # get E-Jones
    if isinstance(nu, (int, float)): data = mod.recons
    else: data = mod.recons_all
    if len(data.shape)==4: data = np.expand_dims(data, axis=2)

    # Normalise
    if args.normalise: data = normalise_multifreq(data)

    # Convert to Stokes/Mueller formalism
    st = [['I', 'IQ', 'IU', 'IV'],
          ['QI', 'Q', 'QU', 'QV'],
          ['UI', 'UQ', 'U', 'UV'],
          ['VI', 'VQ', 'VU', 'V']]
    if args.Stokes!=None:
        m = args.Stokes
        data_M = jones_to_mueller_all(data)
        if m=='M': data = data_M
        else:
            ind = np.where(np.array(st)==m)
            data = np.zeros((1,1,)+data.shape[2:], dtype=np.complex)
            data[0,0,:,:,:] = data_M[ind]
        filename = filename+'_'+m

    save_fits(data, nu, args, filename)    

def save_fits(data, nu, args, filename):
    # Save as fits files
    if args.output_eight and args.Stokes==None:
        write_fits_eight(data, nu, args.diameter, filename)
        print "Saved as 8 files with prefix %s"%filename
    elif args.output_eight and args.Stokes!=None:
        print "8 output files can be created for Jones formalism only, not for Mueller"
        print "!WARNING: No output file created"
    else:
        write_fits(data.real, nu, args.diameter, filename+'_re.fits')
        print "Saved Real part as %s_re.fits"%filename
        if args.Stokes==None:
            write_fits(data.imag, nu, args.diameter, filename+'_im.fits')
            print "Saved Imaginary part as %s_im.fits"%filename
        
