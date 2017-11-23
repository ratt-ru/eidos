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
                fits.writeto(filename[:-5]+'_%s%s_%s.fits'%(C[i],C[j],P[p]), d[p,i,j,...], h, clobber=True)
