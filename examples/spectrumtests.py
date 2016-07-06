import scipy as sp
from scipy.interpolate import interp1d

import sys,os
from pyraf import iraf
from iraf import onedspec

import matplotlib.pyplot as plt

from astropy.io import ascii,fits
from astropy.table import Table

import mapspec as mps


def rebin2(ifile,xgrid,flag,mode='poly5'):
    onedspec.setParam('interp',mode)
    print onedspec.getParam('interp')
    xin = sp.genfromtxt(ifile,usecols=(0))
    y1 = xin[0]
    dy = xin[1] - xin[0]

    if os.path.isfile(ifile + '.fits'):
        os.remove(ifile + '.fits')
    onedspec.rspectext(input=ifile,output=ifile+'.fits',dtype='linear',
                       crval1 = y1,cdelt1 = dy)

    if os.path.isfile(ifile + '2.fits'):
        os.remove(ifile + '2.fits')
    onedspec.dispcor(input=ifile+'.fits',output=ifile+'2.fits',linearize='yes',
                     w1=xgrid[0] ,w2 = xgrid[-1],dw = xgrid[1]-xgrid[0])
   
    if os.path.isfile(ifile + '2.dat'):
        os.remove(ifile + '2.dat')
    onedspec.wspectext(input=ifile+'2.fits[1]',output=ifile+'2.dat')

    y = sp.genfromtxt(ifile+'2.dat',usecols=(1))
    xtest = sp.genfromtxt(ifile+'2.dat',usecols=(0))

    if flag==1:
        if os.path.isfile('rebin_temp'):
            os.remove('rebin_temp')
        os.system('awk \'{print $1\" \"$3}\' %s > rebin_temp'%ifile)

        if os.path.isfile('rebin_temp.fits'):
            os.remove('rebin_temp.fits')
        onedspec.rspectext(input='rebin_temp',output='rebin_temp.fits',dtype='linear',
                       crval1 = y1,cdelt1 = dy)

        if os.path.isfile('rebin_temp2.fits'):
            os.remove('rebin_temp2.fits')
        onedspec.dispcor(input='rebin_temp.fits',output='rebin_temp2.fits',linearize='yes',
                         w1=xgrid[0] ,w2 = xgrid[-1],dw = xgrid[1]-xgrid[0])

        if os.path.isfile('rebin_temp2.dat'):
            os.remove('rebin_temp2.dat')
        onedspec.wspectext(input='rebin_temp2.fits[1]',output='rebin_temp2.dat')
        
        z = sp.genfromtxt('rebin_temp2.dat',usecols=(1))
        return xgrid,y,z
    else:
        return xgrid,y




#Sp = Spectrum('ngc5548_140625_merge.dat')
#S1 = TextSpec('ngc5548.140625.merge.dat')
#S2 = MODSfitsSpec('ngc5548_140608_m1r.fits')
#S3 = MDMfitsSpec('n5548_001_final.fits')
#S4 = TextSpec('n5548_001_final.txt')

S1 = mps.spectrum.TextSpec('test.dat')

S1.set_interp()
xinterp = sp.r_[S1.wv.min():S1.wv.max():15000j]
yinterp = S1.interp(xinterp)

plt.plot(S1.wv,S1.f,'bo-',label='orig')
#xrebin = xinterp
xrebin = sp.r_[S1.wv.min():S1.wv.max():3]
#xrebin  = sp.r_[S1.wv.min() + 1000:S1.wv.max()-30:15]
#temp  = sp.r_[S1.wv.min() + 20:S1.wv.min()+200:15]
#xrebin = sp.r_[temp,xrebin]
print xrebin.min(),xrebin.max()
print S1.wv.min(),S1.wv.max()
S1.rebin(xrebin)


#print S1.wv.size

#
#insert = sp.searchsorted(S1.wv,xrebin)
#xinsert = sp.insert(S1.wv,insert,xrebin)
#xinsert = sp.unique(xinsert)
#
#xrebin2 = sp.zeros(xrebin.size)
#i = sp.digitize(xinsert,xrebin)
#for j in range(xrebin.size):
#    iuse = sp.where(i == j+1)[0]
#    xrebin2[j] = sp.mean(xinsert[iuse])
#
#print xrebin2

xr2,yr2,zr2 = rebin2('test.dat',xrebin,1,mode='spline3')
xr3,yr3,zr3 = rebin2('test.dat',xrebin,1,mode='sinc')
print xr2.shape,yr2.shape

#plt.plot(xinterp,yinterp,'ro-')
plt.plot(S1.wv,S1.f,'ro-',label='custom')
plt.plot(xr2,yr2,'co-',label='dispcor spline3')
plt.plot(xr3,yr3,'go-',label='dispcor sinc')
plt.legend(loc='upper right')
plt.show()

#print sp.c_[S2.wv,S2.f,S2.ef]
#print sp.c_[S4.wv,S4.f,S4.ef]
#print Sp
