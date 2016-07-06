import scipy as sp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import sys,os
sys.path.append('..')
from pyraf import iraf
from iraf import onedspec

from mapspec.spectrum import *

"""
This compares rebinned spectra using iraf (dispcor) and
mapspec.spectrum object
"""


#the old, hard way
def rebin_iraf(ifile,xgrid,flag,mode='poly5'):
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
                     w1=xgrid[0] ,w2 = xgrid[-1],dw = xgrid[1]-xgrid[0],flux='no')
   
    if os.path.isfile(ifile + '2.dat'):
        os.remove(ifile + '2.dat')
    onedspec.wspectext(input=ifile+'2.fits[1]',output=ifile+'2.dat',header='no')

    y = sp.genfromtxt(ifile+'2.dat',usecols=(1))
    xtest = sp.genfromtxt(ifile+'2.dat',usecols=(0))

    os.remove(ifile+'.fits')
    os.remove(ifile+'2.fits')
    os.remove(ifile+'2.dat')
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
                         w1=xgrid[0] ,w2 = xgrid[-1],dw = xgrid[1]-xgrid[0],flux='no')

        if os.path.isfile('rebin_temp2.dat'):
            os.remove('rebin_temp2.dat')
        onedspec.wspectext(input='rebin_temp2.fits[1]',output='rebin_temp2.dat',header='no')
        
        z = sp.genfromtxt('rebin_temp2.dat',usecols=(1))

        os.remove('rebin_temp.fits')
        os.remove('rebin_temp2.fits')
        os.remove('rebin_temp2.dat')
        return xgrid,y,z
    else:
        return xgrid,y

#sinc interpolation is slow, but is the most interesting comparison
S1 = TextSpec('test.dat',style='sinc')


plt.plot(S1.wv,S1.f,'bo-',label='orig')

xrebin = sp.r_[S1.wv.min():S1.wv.max():3]

S1.rebin(xrebin)


xr2,yr2,zr2 = rebin_iraf('test.dat',xrebin,1,mode='spline3')
xr3,yr3,zr3 = rebin_iraf('test.dat',xrebin,1,mode='sinc')


plt.plot(S1.wv,S1.f,'ro-',label='custom')
plt.plot(xr2,yr2,'co-',label='dispcor spline3')
plt.plot(xr3,yr3,'go-',label='dispcor sinc')
plt.legend(loc='upper right')
plt.show()

