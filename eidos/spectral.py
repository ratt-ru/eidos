from scipy.fftpack import dct,idct
import numpy as np

def best_coeffs(d, thr):
    coeffs = {}
    modes = np.zeros((2,2,2,thr))
    coeffs = np.zeros((2,d.shape[0],2,2,thr))
    half = int(d.shape[0]/2)
    dat = [d.real, d.imag]
    for i in range(2):
        for j in range(2):
            for k,kk in enumerate(dat):
                avg = np.mean(abs(kk[:,i,j,:]), axis=0)
                idx = np.argsort(avg)[::-1][:thr]
                modes[k,i,j,:] = idx
                coeffs[k,:,i,j,:] = kk[:,i,j,idx]
                
    return coeffs, modes

def good_channels_sig(d, nu, c=[0,0], edge=1628):
    sig = np.array([np.std(abs(d[i,c[0],c[1],:,:])) for i in range(d.shape[0])])*1e3
    med = np.median(sig)
    e = abs(nu-edge).argmin()
    bads = np.where(np.logical_or(sig[:e]>med+1,sig[:e]<med-2))[0]
    goods = list(set(np.arange(len(nu)))-set(bads))
    return goods, sig

def good_channels_adjacent(C, I, ad=1):
    indices = []
    for i in range(I.shape[-1]):
        idx = int(I[i])
        d = C[:,idx]
        d1 = np.zeros(len(d))
        d1[ad:] = d[:-ad]
        res = abs(d-d1)
        ind = np.where(res>0.001)
        indices.extend(ind[0])
    bads = np.unique(indices)
    goods = list(set(np.arange(C.shape[0]))-set(bads))
    return goods, bads

def interpol(data, nu0, nu1):
    datao = np.zeros((2,len(nu1),2,2,data.shape[-1]))
    for i in range(2):
        for j in range(2):
            for k in range(data.shape[-1]):
                for p in range(2):
                    d = data[p,:,i,j,k]
                    datao[p,:,i,j,k] = np.interp(nu1, nu0[d!=0], d[d!=0])
    return datao


def dct_decom(x, y, thresh):
    dctcoeff = dct(y, norm='ortho')
    co = abs(dctcoeff)
    sorti = np.argsort(co)[::-1]
    if isinstance(thresh, float):
        noise = np.std(co[sorti][15:])
        ind = np.where(co>thresh*noise)
    elif isinstance(thresh, int): ind = sorti[:thresh]
    coeffs = dctcoeff[ind]
    return coeffs, ind

def dct_recon(coeffs, ind, nchan):
    co = np.zeros(nchan)
    co[ind] = coeffs
    return idct(co, norm='ortho')

def dct_recon_all(Co):
    F,C,I = Co['nu'], Co['dctc'], Co['dcti']
    recons = np.zeros((2,len(F),2,2,C.shape[-1]))
    for i in range(2):
        for j in range(2):
            for k in range(C.shape[-1]):
                for p in range(2):
                    recons[p,:,i,j,k] = dct_recon(C[p,:,i,j,k], list(map(int,I[p,:,i,j,k])), len(F))
    return recons

def dctise_all(data, nu, thresh):
    recons = np.zeros((2,len(nu),2,2,data.shape[-1]))
    coeffs = np.zeros((2,thresh,2,2,data.shape[-1]))
    inds = np.zeros((2,thresh,2,2,data.shape[-1]))
    for i in range(2):
        for j in range(2):
            for k in range(data.shape[-1]):
                for p in range(2):
                    d = data[p,:,i,j,k]
                    coeff, ind = dct_decom(nu, d, thresh)
                    recon = dct_recon(coeff, ind, len(nu))
                    coeffs[p,:,i,j,k] = coeff
                    inds[p,:,i,j,k] = ind
                    recons[p,:,i,j,k] = recon
    return recons, coeffs, inds
