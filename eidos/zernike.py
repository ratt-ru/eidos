#!/usr/bin/env python
from util import *
from scipy.misc import factorial as fac

class Zernike(object):
    """
    Decompose and reconstruct a 2d image or 4d Jones matrix using Zernike polynomials
    """
    def __init__(self, data=[], Nmodes=50, threshold=None, Npix=None, m=0, n=0, mode='recon', idx=None, thresh=None, freq=None):
        self.Nmodes = Nmodes
        self.threshold = threshold
        self.thresh = thresh
        self.mode = mode
        self.npix=Npix
        self.idx = idx

        if isinstance(data, dict):
            if isinstance(freq, (float,int)):
                ch = freq_to_idx(sb=[0,freq])[1]
                self.dict_recon(data, ch)
            if isinstance(freq, (list,np.ndarray)):
                nchan = len(freq)
                self.recons_all = np.zeros((2,2,nchan,self.npix,self.npix), dtype=np.complex)
                for c in range(nchan):
                    print freq[c],
                    ch = freq_to_idx(sb=[0,freq[c]])[1]
                    self.dict_recon(data, ch)
                    self.recons_all[:,:,c,:,:] = self.recons

        else:
            if mode=='img':
                self.img = data
                self.reconstruct()

            if 'jones' in mode: self.jones_images(data)
            elif 'cube' in mode: self.all_freqs(data)

        if mode=='basis':
            if not Nmodes:
                self.unit_disk(npix=Npix)
                self.basis = self.zernike(m, n, self.grid_rho, self.grid_phi)*self.grid_mask
                self.basis = circular_mask(np.array(self.basis))
            if Nmodes:
                self.unit_disk(npix=Npix)
                self.basis = [self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for i in range(Nmodes)]
                self.basis = circular_mask(np.array(self.basis))

    def jones_images(self, data):
        self.coeffs_J = self.coeffs_trunc_J = np.zeros((data.shape[0], data.shape[1], self.Nmodes), dtype=np.complex)
        self.recon_full_J = self.recon_trunc_J = np.zeros(data.shape, dtype=np.complex)
        if 'recon' in self.mode: self.recon = np.zeros((2,2,self.npix,self.npix), dtype=np.complex)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                print "Fitting Zernike polynomials to Jones %i %i"%(i,j)
                if 'decom' in self.mode:
                    self.img = data[i,j,:,:]
                    self.decompose()
                    self.coeffs_J[i,j,:] = self.coeffs
                elif 'recon' in self.mode:
                    self.coeffs = data[i,j,:]
                    try: self.ind = self.idx[i,j,:]
                    except: self.ind = None
                    self.recon[i,j,:,:] = self.reconstruct()
                elif 'both' in self.mode:
                    self.img = data[i,j,:,:]
                    self.reconstruct()
                    self.coeffs_J[i,j,:] = self.coeffs
                    self.coeffs_trunc_J[i,j,:] = self.coeffs_trunc
                    self.recon_full_J[i,j,:,:] = self.recon_full
                    self.recon_trunc_J[i,j,:,:] = self.recon_trunc

    def dict_recon(self, data, ch):
        corrs = [['xx','xy'],['yx','yy']]
        self.recons = np.zeros((2,2,self.npix,self.npix), dtype=np.complex)
        for i in range(2):
            for j in range(2):
                corr = corrs[i][j]
                dat, ind = data[corr], data['basei'][corr]
                if i==j: thresh = self.thresh[0]
                elif i!=j: thresh = self.thresh[1]
                self.coeffs, self.ind = dat[ch,:thresh], ind[:thresh]
                self.reconstruct()
                self.recons[i,j,:,:] = self.recon

    def all_freqs(self, fits_file):
        start = time.time()
        d = fits.getdata(fits_file)
        h = fits.getheader(fits_file)
        dc = d[0,...] + 1j * d[1,...]
        dc = np.nan_to_num(dc)
        print dc.shape
        recons = np.zeros(dc.shape, dtype=np.complex)
        coeffs = np.zeros((2,2,dc.shape[2],self.Nmodes), dtype=np.complex)
        for f in range(dc.shape[2]):
            print '... Channel %i'%(f+1)
            self.jones_images(dc[:,:,f,:,:])
            if 'decom' in self.mode: coeffs[:,:,f,:] = self.coeffs_J
            elif 'recon' in self.mode: recons[:,:,f,:,:] = self.recon

        if 'decom' in self.mode: np.save(fits_file[:-5]+'_zp_%icoeffs.npy'%self.Nmodes, coeffs)

        elif 'recon' in self.mode:
            d[0,...], d[1,...] = recons.real, recons.imag
            fits.writeto(fits_file[:-5]+'_zp_recon_%icoeffs.fits'%self.Nmodes, circular_mask(d), h, overwrite=True)

        end = time.time()
        print "... Time taken: %.2f minutes"%((end-start)/60.)

    def zernike_rad(self, m, n, rho):
        """
        Calculate the radial component of Zernike polynomial (m, n) 
        given a grid of radial coordinates rho.
        """
        if (n < 0 or m < 0 or abs(m) > n):
            raise ValueError
        if ((n-m) % 2):
            return rho*0.0
        pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
        return sum(pre_fac(k) * rho**(n-2.0*k) for k in xrange((n-m)/2+1))

    def zernike(self, m, n, rho, phi):
        """
        Calculate Zernike polynomial (m, n) given a grid of radial
        coordinates rho and azimuthal coordinates phi.
        """
        if (m > 0): return self.zernike_rad(m, n, rho) * np.cos(m * phi)
        if (m < 0): return self.zernike_rad(-m, n, rho) * np.sin(-m * phi)
        return self.zernike_rad(0, n, rho)

    def noll_to_zern(self, j):
        """
        Convert linear Noll index to tuple of Zernike indices.
        j is the linear Noll coordinate, n is the radial Zernike index and m is the azimuthal Zernike index.
        @param [in] j Zernike mode Noll index
        @return (n, m) tuple of Zernike indices
        @see <https://oeis.org/A176988>.
        """
        # add 1 to start from 1
        j += 1

        n = 0
        j1 = j-1
        while (j1 > n):
            n += 1
            j1 -= n

        m = (-1)**j * ((n % 2) + 2 * int((j1+((n+1)%2)) / 2.0 ))
        return (n, m)

    def zernikel(self, j, rho, phi):
        """
        Calculate Zernike polynomial with Noll coordinate j given a grid of radial
        coordinates rho and azimuthal coordinates phi.
        """
        nm = self.noll_to_zern(j)
        m, n = nm[1], nm[0]
        return self.zernike(m, n, rho, phi)
    
    def unit_disk(self, npix=None):
        """Create an unit disk and convert to rho, phi"""
        if npix: nx, ny = npix, npix
        else: nx, ny = self.img.shape
        grid = (np.indices((nx, ny), dtype=np.float) - nx/2) / (nx*1./2) # create unit grid [-1,1]
        self.grid_rho = (grid**2.0).sum(0)**0.5 # rho = sqrt(x^2+y^2)
        self.grid_phi = np.arctan2(grid[0], grid[1]) # phi = itan(x/y)
        self.grid_mask = self.grid_rho <= 1 # boolean array specifying where rho<=1
    
    def decompose(self):
        """Decompose using SVD"""
        self.unit_disk()
        # Caculate Zernike bases given the maximum Noll index, N
        N = self.Nmodes
        basis = [self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for i in range(N)]

        # Calculate covariance between all Zernike polynomials
        self.cov_mat = np.array([[np.sum(zerni * zernj) for zerni in basis] for zernj in basis])

        # Invert covariance matrix using SVD (  A x = b  ==>  x =>  x= A^{pseudoinv} b)
        self.cov_mat_in = np.linalg.pinv(self.cov_mat)

        # Inner product between the img and the Zernike bases
        self.innerprod = np.array([np.sum(self.img * zerni) for zerni in basis])
        
        # Dot product between inverse covariance matrix and the innerprod to get the coeffs
        self.coeffs = np.dot(self.cov_mat_in, self.innerprod)

    def truncate(self):
        """Truncate the coefficients upto the given threshold"""
        sortedindex = np.argsort(np.abs(self.coeffs))[::-1]
        Ncoeff = self.coeffs.shape[-1]
        cutoff = np.int(np.round(Ncoeff*self.threshold/100.))
        
        print "Keeping %2.0f %% (N=%s) of the biggest coefficients"%(self.threshold,cutoff)

        self.coeffs_trunc = self.coeffs.copy() # copy of all coeff
        self.coeffs_trunc[sortedindex[cutoff:]] = 0 # put coeff below threshold to 0

    def reconstruct(self):
        """Reconstruct a model image from the coeffcicients"""
        if 'recon' in self.mode:
            self.unit_disk(self.npix)
            if self.threshold!=None:
                self.truncate()
                return np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs_trunc))
            else:
                self.recon = np.sum(self.coeffs[i] * self.zernikel(val, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.ind))

        if 'both' in self.mode:
            self.decompose()
            self.recon_full = np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs))
            self.res_full = (abs(self.img) - abs(self.recon_full))  * self.grid_mask
            self.truncate()
            self.recon_trunc = np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs_trunc))
            self.res_trunc = (abs(self.img) - abs(self.recon_trunc)) * self.grid_mask
            self.diff_full_trunc = (self.recon_full - self.recon_trunc) * self.grid_mask
