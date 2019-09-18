from astropy.io import fits
from africanus.rime import zernike_dde
import numpy as np
import dask
import dask.array as da
import time
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# Import Eidos functions
from spectral import dct_recon, dct_recon_all
from spatial import recon_par, write_fits
from create_beam import zernike_parameters

def save_fits(data, nu, filename):
    write_fits(data.real, nu, 10, filename+'_re.fits')
    print("Saved Real part as %s_re.fits"%filename)
    write_fits(data.imag, nu, 10, filename+'_im.fits')
    print("Saved Imaginary part as %s_im.fits"%filename)

na = 1
nchan = 1
ntime = 1
npoly = 20

 # LB this is what we want to reconstruct
eidos_pb = fits.getdata("/home/landman/Students/Josh_v_Staden/Data/Beams/mk_1070_test_xx_re.fits") 

hdr = fits.getheader("/home/landman/Students/Josh_v_Staden/Data/Beams/mk_1070_test_xx_re.fits")

# Get center values, indices, deltas, and npix

# LB - these are always zero
# cr_l = hdr['CRVAL2']
# cr_m = hdr['CRVAL1']

# LB - crpix should alwys be in the centre
# crpix_l = hdr['CRPIX2']
# crpix_m = hdr['CRPIX1']

del_l = hdr['CDELT2']  # * np.pi/180
del_m = hdr['CDELT1']  # * np.pi/180

npix = hdr['NAXIS1']

# Generate an lm grid
l = np.arange(-(npix//2), npix//2) * del_l
m = np.arange(-(npix//2), npix//2) * del_m

# scale coordinates to unit disc
scaling = 5
l /= scaling
m /= scaling

ll, mm = np.meshgrid(l, m)

lm = np.vstack((ll.flatten(), mm.flatten()))

# Create arguments for Codex call. Separate arrays for real and imaginary coefficients and noll indices
coords = np.empty((3, npix * npix, ntime, na, nchan), dtype=np.float64)
coeffs_r = np.empty((na, nchan, 2, 2, npoly), dtype=np.float64)
coeffs_i = np.empty((na, nchan, 2, 2, npoly), dtype=np.float64)
noll_index_r = np.empty((na, nchan, 2, 2, npoly), dtype=np.int32)
noll_index_i = np.empty((na, nchan, 2, 2, npoly), dtype=np.int32)
parallactic_angles = np.zeros((ntime, nchan,), dtype=np.float64)  # LB - I think nchan should be na here
antenna_scaling = np.ones((na, nchan, 2))
frequency_scaling = np.ones((nchan,), dtype=np.float64)
pointing_errors = np.zeros((ntime, na, nchan, 2), dtype=np.float64)

coords[:2, :, 0, 0, 0], coords[2, :, 0, 0, 0] = lm[:, :], 0

# Get Zernike coefficients and noll indices (taken from create_beam.py in Eidos)
filename = "data/meerkat_beam_coeffs_ah_zp_dct.npy"
params, freqs = zernike_parameters(filename)

nu = 1070.

ch = abs(freqs-nu).argmin()
# B = recon_par(params[ch,:])

data = params[ch,:]

# Coeffs and nolls split into real and imaginary sections
coeffs, nolls, _, _ = data[:4]

# Asssign values to arrays
coeffs_r[0,0,:,:,:] = coeffs[0,:,:,:]
coeffs_i[0,0,:,:,:] = coeffs[1,:,:,:]
noll_index_r[0,0, :, :, :] = nolls[0, :, :, :]
noll_index_i[0,0, :, :, :] = nolls[1, :, :, :]

# Create a dictionary to write to zernike_real_imag_coeffs.npy for Codex Africanus
corr_letters = [b'x', b'y']
zernike_file_coeffs = {b'coeffs': {b'real':  {b'xx': None,
                                            b'xy': None,
                                            b'yx':None,
                                            b'yy':None},
                                    b'imag':  {b'xx': None,
                                            b'xy': None,
                                            b'yx':None,
                                            b'yy':None}
                                    },
                        b'noll_index': {b'real': {b'xx': None,
                                            b'xy': None,
                                            b'yx' : None,
                                            b'yy' : None},
                                        b'imag':  {b'xx': None,
                                            b'xy': None,
                                            b'yx':None,
                                            b'yy':None}
                        }
}
for i in range(2):
    for j in range(2):
        corr_index = corr_letters[i] + corr_letters[j]
        zernike_file_coeffs[b'coeffs'][b'real'][corr_index] = coeffs_r[0,0, i, j, :]
        zernike_file_coeffs[b'coeffs'][b'imag'][corr_index] = coeffs_i[0,0, i, j, :]
        zernike_file_coeffs[b'noll_index'][b'real'][corr_index] = noll_index_r[0,0, i, j, :]
        zernike_file_coeffs[b'noll_index'][b'imag'][corr_index] = noll_index_i[0,0, i, j, :]
# np.save("zernike_real_imag_coeffs.npy", zernike_file_coeffs)

# Calls to zernike_dde
zernike_real = zernike_dde(coords, coeffs_r, noll_index_r, parallactic_angles, frequency_scaling,  antenna_scaling, pointing_errors)
zernike_imag = 1j * zernike_dde(coords, coeffs_i, noll_index_i, parallactic_angles, frequency_scaling,  antenna_scaling, pointing_errors)
zernike_vals = (zernike_real + zernike_imag)[:, 0, 0, 0, 0, 0].reshape((1,1,npix,npix))

print(zernike_real.shape)
print(eidos_pb.shape)

plt.figure('africanus')
plt.imshow(zernike_real[:, 0, 0, 0, 0, 0].reshape(npix, npix))
plt.colorbar()

plt.figure('eidos')
plt.imshow(eidos_pb[0].real)
plt.colorbar()

plt.show()

# # Use save_fits function from Eidos
# save_fits(zernike_vals, nu, "save_fits_codex")