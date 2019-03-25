#!/usr/bin/env python
from util import *
from scipy.misc import factorial as fac

class Zernike(object):
    """
    Decompose and reconstruct a 2d image or 4d Jones matrix using Zernike polynomials
    """
    def __init__(self, data=[], Nmodes=50, threshold=None, Npix=None, m=0, n=0, mode='recon', idx=None, freq=None, fname=''):
        self.Nmodes = Nmodes
        self.threshold = threshold
        self.mode = mode
        self.npix=Npix
        self.idx = idx

        if mode=='img':
            self.img = data
            self.reconstruct()

        if 'jones' in mode: self.jones_images(data)

    def jones_images(self, data):
        self.coeffs_J = self.coeffs_trunc_J = np.zeros((data.shape[0], data.shape[1], self.Nmodes), dtype=np.complex)
        self.recon_full_J = self.recon_trunc_J = np.zeros(data.shape, dtype=np.complex)
        if 'recon' in self.mode: self.recon = np.zeros((2,2,self.npix,self.npix), dtype=np.complex)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                #print "Fitting Zernike polynomials to Jones %i %i"%(i,j)
                if self.threshold:
                    if i==j: self.thresh = self.threshold[0]
                    elif i!=j: self.thresh = self.threshold[1]
                if 'decom' in self.mode:
                    self.img = data[i,j,:,:]
                    self.decompose()
                    self.coeffs_J[i,j,:] = self.coeffs
                elif 'recon' in self.mode:
                    self.coeffs = data[i,j,:]
                    try: self.ind = list(self.idx[i,j,:])
                    except: self.ind = None
                    self.recon[i,j,:,:] = self.reconstruct()
                elif 'both' in self.mode:
                    self.img = data[i,j,:,:]
                    self.reconstruct()
                    self.coeffs_J[i,j,:] = self.coeffs
                    self.coeffs_trunc_J[i,j,:] = self.coeffs_trunc
                    self.recon_full_J[i,j,:,:] = self.recon_full
                    self.recon_trunc_J[i,j,:,:] = self.recon_trunc

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
        return sum(pre_fac(k) * rho**(n-2.0*k) for k in np.arange((n-m)/2+1))

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

    def truncate(self, thresh):
        """Truncate the coefficients upto the given threshold"""
        sortedindex = np.argsort(np.abs(self.coeffs))[::-1]
        Ncoeff = self.coeffs.shape[-1]
        cutoff = np.int(np.round(Ncoeff*thresh/100.))
        
        #print "Keeping %2.0f %% (N=%s) of the biggest coefficients"%(thresh,cutoff)

        self.coeffs_trunc = self.coeffs.copy() # copy of all coeff
        self.coeffs_trunc[sortedindex[cutoff:]] = 0 # put coeff below threshold to 0

    def best_coeffs(self, C, I):
        idx = np.argsort(np.abs(C))[::-1][:self.thresh]
        return C[idx], I[idx]

    def reconstruct(self):
        """Reconstruct a model image from the coeffcicients"""
        if 'recon' in self.mode:
            self.unit_disk(self.npix)
            if self.ind:
                if self.thresh: C, I = self.best_coeffs(self.coeffs, np.array(self.ind))
                else: C, I = self.coeffs, self.ind
                return np.sum(C[i] * self.zernikel(val, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(I))
            else:
                self.truncate(self.thresh)
                return np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs_trunc))

        if 'both' in self.mode:
            self.decompose()
            self.recon_full = np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs))
            self.res_full = (abs(self.img) - abs(self.recon_full))  * self.grid_mask
            self.truncate(self.thresh)
            self.recon_trunc = np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs_trunc))
            self.res_trunc = (abs(self.img) - abs(self.recon_trunc)) * self.grid_mask
            self.diff_full_trunc = (self.recon_full - self.recon_trunc) * self.grid_mask

def noll_to_zern(j):
        j += 1

        n = 0
        j1 = j-1
        while (j1 > n):
            n += 1
            j1 -= n

        m = (-1)**j * ((n % 2) + 2 * int((j1+((n+1)%2)) / 2.0 ))
        return (n, m)

def decom_par(data):
    Z = Zernike(data, mode='jones decom', Nmodes=300, threshold=[15,8])
    return Z.coeffs_J

def recon_par(data):
    coeffs, nolls, npix, thr = data[0], data[1], data[2], data[3]
    Zr = Zernike(coeffs[0,...], idx=nolls[0,...], threshold=[thr,thr], Npix=npix, mode='jones recon')
    Zi = Zernike(coeffs[1,...], idx=nolls[1,...], threshold=[thr,thr], Npix=npix, mode='jones recon')
    B = Zr.recon + 1j*Zi.recon
    return B
    