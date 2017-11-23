#!/usr/bin/env python
import numpy as np, matplotlib.pyplot as plt, sys
from numpy import pi, exp, sin, cos, tan
from scipy import interpolate
from scipy.linalg import svd
from lmfit import Parameters, minimize, fit_report
from scipy.fftpack import dct,idct

def good_channels(nu=None, bads=None):
    if nu==None: nu = np.arange(856, 856+857)
    if bads==None: bads=([856,870], [920,960], [1125,1305], [1380,1384], [1463,1492], [1520,1630], [1670,1720])
    ind = []
    for bad in bads:
        ind.append(np.where((nu>=bad[0]) & (nu<=bad[1])))   
    flatten = lambda l: [item for sublist in l for item in sublist[0]]
    all = range(len(nu))
    bads = flatten(ind)
    goods = list(set(all)-set(bads))
    return goods, bads

def decompose(coeffs, nu, sb, bthresh=[10,5], dcthresh=[10,30], pars=[1.,5.,125.]):
    mod = {}
    corrs = [['xx','xy'],['yx','yy']]
    for i in range(2):
        for j in range (2):
            corr = corrs[i][j]
            mod[corr] = {}
            # take the unique coeffs in the 10 best from xx,yy and 5 best from xy,yx
            if i==j:
                ind = sortindex(coeffs[i,j,:,:], nu, sb, thresh=bthresh[0])
                thresh = dcthresh[0]
            if i!=j:
                ind = sortindex(coeffs[i,j,:,:], nu, sb, thresh=bthresh[1])
                thresh = dcthresh[1]
            lind, hind = list(ind[0]), list(ind[1])
            zind = np.unique(lind+hind)
            coeff = coeffs[i,j,:,zind].T

            # sort the best coeffs in decreasing energy order
            sind = sortindex(coeff, nu, sb=[[856,856+857]], thresh=None)[0]
            zind = zind[sind]
            coeff = coeffs[i,j,:,zind].T
            mod[corr]['basei'] = zind
            
            # calculate ampl and phase of the best coeffs
            A = np.abs(coeff)
            P = np.angle(coeff)
            
            # model the ampl and phase of each coeff
            nc = A.shape[1]
            dind = np.zeros((thresh, nc), dtype='int')
            dcoef = np.zeros((thresh, nc))
            phase = np.zeros((2,3, nc))
            for k in range(nc):
                # model amplitude using DCT
                y = A[:,k]
                A_intp = np.interp(nu, nu[y!=0.], y[y!=0.])
                dctcoeff = dct(A_intp, norm='ortho') # discrete Fourier transform, scipy.fftpack
                di = np.argsort(np.abs(dctcoeff))[::-1]
                dind[:,k] = di[:thresh] # find *dcthresh* strongest dct coeffs
                dcoef[:,k] = dctcoeff[dind[:,k]]
                
                # model the sin(phi) and cos(phi) using 3-parameter sine waves
                phase[0,:,k] = lmfit_sine(nu, np.sin(P[:,k]), a=pars[0],b=pars[1],c=pars[2]) # sin(phi)
                phase[1,:,k] = lmfit_sine(nu, np.cos(P[:,k]), a=pars[0],b=pars[1],c=pars[2]) # cos(phi)
                    
            # put the interpolated models in the dictionary: *mod*
            mod[corr]['dcti'] = dind
            mod[corr]['dctc'] = dcoef
            mod[corr]['sphi'] = phase[0,:,:]
            mod[corr]['cphi'] = phase[1,:,:]
    return mod

def reconstruct(mod, nu, mid=1300, data=None):
    recons = {}
    recons['data'] = {}
    recons['polar'] = {}
    recons['basei'] = {}
    corrs = [['xx','xy'],['yx','yy']]
    for i in range(2):
        for j in range(2):
            corr = corrs[i][j]
            midx = freq_to_idx(nu, [0,mid])[1]
            modij = mod[corr]
            basei = modij['basei']
            recons['data'][corr] = data[i,j,:,basei].T
            recons['basei'][corr] = basei
            nbase = len(basei)
            coeffs = np.zeros((3, len(nu), nbase))

            for k in range(nbase):
                # reconstruct the amplitude for all freqs
                ampl_dct = np.zeros(len(nu))
                ampl_dct[modij['dcti'][:,k]] = modij['dctc'][:,k]
                coeffs[0,:,k] = idct(ampl_dct, norm='ortho')

                # recon sin(phi) and cos(phi) using sine waves
                sphi = modij['sphi'][:,k]
                cphi = modij['cphi'][:,k]
                coeffs[1,:,k] = lmfit_res(sphi, nu)
                coeffs[2,:,k] = lmfit_res(cphi, nu)

            recons[corr] = coeffs[0,...] * (coeffs[2,...]+1j*coeffs[1,...])
            recons['polar'][corr] = coeffs
    return recons

def lmfit_res(pars, x, data=None):
    try:
        p = pars.valuesdict()
        a,b,c = p['a'], p['b'], p['c']
    except:
        a,b,c = pars[0], pars[1], pars[2]

    modul = a*np.cos(2*np.pi*(1./b)*x)
    model = sin(2.*pi* (1./c) *x + modul )
    
    if data is None: return model
    return (model - data)

def lmfit_sine(x, y, goods=good_channels()[0], a=1.,b=5.,c=125.):
    y[y==1.] = 0.
    #timestep = x[1]-x[0]
    #ks = abs(np.fft.fft(y))
    #fs = np.fft.fftfreq(x.size, timestep)
    
    p = Parameters()
    phases = np.random.random(3)
    
    p.add('a', value=a)
    p.add('b', value=x.size*b)
    p.add('c', value=c)
    
    out = minimize(lmfit_res, p, args=(x[goods],), kws={'data':y[goods]}, method='nelder')
    p = out.params.valuesdict()
    
    return [p['a'], p['b'], p['c']]

def inspect(coeffs, nu, idxs=None, original=[], show=True, Noll=True, ylim=[None,None], xlim=None, marker='-',\
    labels=None, colors=None, lw=1, msize=1.5, figsize=(20,10), share=False, part='imag', err=False):
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['font.size'] = 12
    if share==True: f, axes = plt.subplots(2,2, sharex=True, sharey=True)
    else: f, axes = plt.subplots(2,2, sharex=True, sharey=False)
    colrs = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'steelblue', 'gray', 'maroon']
    corrs = [['xx', 'xy'], ['yx', 'yy']]
    markers = ['.', '-']

    if isinstance(coeffs, dict):
        cdict = coeffs
        coeffs = dict_to_array(coeffs, part)
    if isinstance(original, dict):
        original = dict_to_array(original, part)

    if original==[]: original = coeffs

    errs = original-coeffs

    original[original==0.] = None
    original[original==1.] = None
    coeffs[coeffs==0.] = None
    coeffs[coeffs==1.] = None

    if lw<=2: alpha=None
    else: alpha = 0.4

    c = 0
    for i in range(2):
        for j in range(2):
            try:
                coeffs_avg = np.mean(coeffs[i,j,:,:], axis=0)
                sortidx = np.argsort(coeffs_avg)[::-1]
            except: None
            ax = axes[i,j]
            count = 0
            err = []
            if isinstance(labels, dict):
                labs = labels[corrs[i][j]]
            else: labs = None
            for idx in idxs:
                if colors==None: colr = colrs[count]
                else: colr = colors[count]
                if Noll==True:
                    ind = idx
                    try: label = labels[i][j][idx]
                    except: label=labels
                else:
                    ind = sortidx[idx]
                    label = sortidx[idx]+1

                if i==j: err = errs[i,j,ind]*1e2
                elif i!=j: err = errs[i,j,ind]*1e3
                lab = "%i"%label
                
                try: ax.plot(nu, coeffs[i,j,:,ind], marker, markersize=msize, color=colr, label=lab, lw=lw, alpha=alpha)
                except: ax.plot(nu, coeffs[i,j,:], marker, markersize=msize, color=colr, label=lab, lw=lw, alpha=alpha)

                try: ax.plot(nu, original[i,j,:,ind], '.', markersize=msize, color=colr)
                except: ax.plot(nu, original[i,j,:], '.', markersize=msize, color=colr)

                if i==j: ax.set_ylim(ylim[0])
                elif i!=j: ax.set_ylim(ylim[1])
                count += 1
            leg = ax.legend(ncol=3, loc='best')
            leg.get_frame().set_facecolor('none')
            leg.get_frame().set_linewidth('0.')
            if i==1: ax.set_xlabel('Frequency [MHz]')
            if j==0: ax.set_ylabel('Coefficient energy')
            #ax.set_title(corrs[i][j])
            if xlim==None: ax.set_xlim([min(nu), max(nu)+1])
            else: ax.set_xlim(xlim)
            if share==True: f.subplots_adjust(hspace=0.07, wspace=0.01)
            else: f.subplots_adjust(hspace=0.07, wspace=0.1)
    #plt.suptitle('Real (solid) and Imag (dashed) highest-energy coeffs: %i to %i'%(idxs[0], idxs[-1]), fontsize=14)
    
    if show: plt.show()

def dict_to_array(coeffs, part=None):
    corrs = [['xx','xy'],['yx','yy']]
    xx, xy, yx, yy = coeffs['xx'], coeffs['xy'], coeffs['yx'], coeffs['yy']
    nmax = np.max([xx.shape[1],xy.shape[1],yx.shape[1],yy.shape[1]])
    data = np.zeros((2,2,xx.shape[0],nmax), dtype='complex')
    for i in range(2):
        for j in range(2):
            corr = corrs[i][j]
            lim = coeffs[corr].shape[1]
            if part=='imag': data[i,j,:,:lim] = coeffs[corr].real
            if part=='real': data[i,j,:,:lim] = coeffs[corr].imag
            if part=='ampl':
                try: data[i,j,:,:lim] = coeffs['polar'][corr][0,...]
                except: data[i,j,:,:lim] = np.abs(coeffs[corr])
            if part=='sinp':
                try: data[i,j,:,:lim] = coeffs['polar'][corr][1,...]
                except: data[i,j,:,:lim] = np.sin(np.angle(coeffs[corr]))
            if part=='cosp':
                try: data[i,j,:,:lim] = coeffs['polar'][corr][2,...]
                except: data[i,j,:,:lim] = np.cos(np.angle(coeffs[corr]))
    return data

def reject_outliers(data, m=2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/(mdev if mdev else 1.)
    return data[s<m]

def flag(coeffs, freqs, fact=2):
    sb = map(int, np.linspace(0,coeffs.shape[2],10))
    for i in range(2):
        for j in range(2):
            for c in range(coeffs.shape[-1]):
                for s in range(len(sb)-1):
                    y = coeffs[i,j,sb[s]:sb[s+1],c]
                    med = np.median(y)
                    mad = np.median(abs(y-np.median(y)))
                    idx = np.where(abs(y) > (med+3*mad))
                    coeffs[i,j,idx,c] = 0.+1j*0.
    return coeffs

def flag_manual(coeffs, freqs, bads=([856,870], [920,960], [1125,1305], [1380,1384], [1463,1492], [1520,1630], [1670,1720])):
    for i in range(2):
        for j in range(2):
            for c in range(coeffs.shape[-1]):
                for bad in bads:
                    idx = np.where((freqs>=bad[0]) & (freqs<=bad[1]))
                    try: coeffs[i,j,idx,c] = 0.+1j*0.
                    except: coeffs[i,j,idx,c] = 0.
    return coeffs

def freq_to_idx(freqs=np.arange(857)+856, sb=[0,1300]):
    start = (np.abs(freqs-sb[0])).argmin()
    end = (np.abs(freqs-sb[1])).argmin()
    return start, end

def csort(coeffs, freqs, sb=None, thresh=20):
    if sb: s, e = freq_to_idx(freqs, sb)
    else: s, e = 0, None
    idxs = np.zeros((2,2,thresh), dtype=np.int)
    coeffs_sorted = np.zeros(coeffs.shape[:-1]+(thresh,), dtype=np.complex)
    for i in range(2):
        for j in range(2):
            coeffs_avg = np.mean(np.abs(coeffs[i,j,s:e,:]), axis=0)
            sortidx = np.argsort(coeffs_avg)[::-1]
            idxs[i,j,:] = sortidx[:thresh]
            coeffs_sorted[i,j,:,:] = coeffs[i,j,:,sortidx[:thresh]].T
    return coeffs_sorted, idxs

def sortindex(coeffs, freqs, sb=[[960,1125], [1305,1463]], thresh=10):
    ind = []
    for i in range(len(sb)):
        s, e = freq_to_idx(freqs, sb[i])
        coeffs_avg = np.mean(np.abs(coeffs[s:e,:]), axis=0)
        sortidx = np.argsort(coeffs_avg)[::-1]
        ind.append(sortidx[:thresh])

    return ind

def spline(coeffs, freqs, k=4, s=[8e-4,8e-4]):
    coeffs_interp = np.zeros(coeffs.shape)
    tck0, tck1 = [], []
    nc = coeffs.shape[-1]
    for i in range(2):
        for j in range(2):
            for c in range(nc):
                y = coeffs[i,j,:,c]
                if i==j:
                    tck = interpolate.splrep(freqs[y!=0.], y[y!=0.], k=k, s=s[0])
                    coeffs_interp[i,j,:,c] = interpolate.splev(freqs, tck, der=0)
                    if i==0: tck0.append(len(tck[0]))
                elif i!=j:
                    tck = interpolate.splrep(freqs[y!=0.], y[y!=0.], k=k, s=s[1])
                    coeffs_interp[i,j,:,c] = interpolate.splev(freqs, tck, der=0)
                    if i==0: tck1.append(len(tck[0]))
    print tck0
    print tck1
    return coeffs_interp

def linterpolate(coeffs, freqs):
        coeffs_interp = np.zeros(coeffs.shape)
        for i in range(2):
                for j in range(2):
                        for c in range(coeffs.shape[-1]):
                                y = coeffs[i,j,:,c]
                                coeffs_interp[i,j,:,c] = np.interp(freqs, freqs[y!=0.], y[y!=0.])
        return coeffs_interp

def fitter_phase(nu, phi, sb):
    fit = np.zeros((2,)+phi.shape)
    #sb = ([960,1125], [1305,1520]) # two good subbands of MeerKAT L-band
    s, e = freq_to_idx(nu, sb)
    for i in range(phi.shape[0]):
        for j in range(phi.shape[1]):
            for c in range(phi.shape[-1]):
                ph = phi[i,j,:,c]
                fit[0,i,j,:,c] = lmfit_sine(nu, np.sin(ph), s,e, wl=125.)
                fit[1,i,j,:,c] = lmfit_sine(nu, np.cos(ph), s,e, wl=125.)
    return fit

def dctransform(coeff, pcthresh=1e-3, mode='both', verbose=True):
    dctcoeff = dct(coeff, norm='ortho')
    sortedindex = np.argsort(np.abs(dctcoeff))[::-1]
    Ncoeff = dctcoeff.shape[-1]
    pcnumb = pcthresh # % of coefficients sorted in decreasing energy
    cutoff = np.int(np.round(Ncoeff*pcnumb))
    
    #if verbose: print "%2.2f %% (N=%s)"%(100.*pcnumb,cutoff),
    if verbose==True: print cutoff,
    approxcoeff=dctcoeff.copy() # copy of all coeff
    approxcoeff[sortedindex[cutoff:]]=0 # put coeff below threshold to 0

    recon = idct(approxcoeff, norm="ortho")
    return recon

def dctise(coeffs, thresh, mode='both'):
    recons = np.zeros(coeffs.shape, dtype=np.complex)
    for i in range(2):
        for j in range(2):
            for c in range(coeffs.shape[-1]):
                recons[i,j,:,c] = dctransform(coeffs[i,j,:,c], thresh, mode=mode)
    return recons
