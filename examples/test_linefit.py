import scipy as sp
import matplotlib.pyplot as plt
import sys

sys.path.append('..')

from spectrum import TextSpec,TextSpec_2c,EmissionLine,LineModel
"""
Examples for how to use the line fitting functionality, and to test
the install of spectrum.py.

These are for test.dat, which is a MODS 1 spectrum of NGC 5548 from
the LBT, taken in 2014.  Therefore, things are hard-coded.

You should get 3 plots---one shows results from fitting a single line,
one shows results for multiple components, and one shows results from
an empirical template.


"""

#Hbeta
window = [ [4821,5009],[4672,4711],[5167,5207] ]
#[OIII]lambda 5007
oiii   = [ [5066,5108],[5057,5066],[5108,5165] ]

s = TextSpec('test.dat',style='sinc')


l = EmissionLine(s,window[0],[window[1],window[2]])
plt.plot(l.wv,l.f,'k')

#The first experiment is to fit a single component with a different
#functional form.  Fitting is easy:
mg = LineModel(l,func = 'gaussian')
mh = LineModel(l,func = 'gauss-hermite')
ml = LineModel(l,func = 'lorentzian')
mv = LineModel(l,func = 'approx-voigt')
#if you wanted to see the chi^2 values.....
#print mg.chi2()
#print mh.chi2()
#print ml.chi2()
#print mv.chi2()

#new output grid
xgrid = sp.r_[l.wv.min():l.wv.max() + 1:1000j]

plt.plot(xgrid,mg(xgrid),'m',label='gauss')
plt.plot(xgrid,mh(xgrid),'orange',label='hermite')
plt.plot(xgrid,ml(xgrid),color='r',label='lorentz')
plt.plot(xgrid,mv(xgrid),color='purple',label='voigt')

plt.legend(loc='upper left')

#Now, we will remove narrow hbeta.  Will only do gauss-hermite for
#now, but in principle it is easy to do others 

#narrow line window---be careful that this is smaller than the
#templateline itself, so that there is space to shift around the
#center
#w = (l.wv >= 6651)*(l.wv <= 6680)
w = (l.wv >= 4925)*(l.wv <= 4960)

mgf = LineModel(l,func = 'gaussian',window = w,floating=1)
#mhf = LineModel(l,func = 'gauss-hermite',window = w,floating=1)
#mlf = LineModel(l,func = 'lorentzian',window = w,floating=1)
#mvf = LineModel(l,func = 'approx-voigt',window = w,floating=1)

#we will also try extracting narrow [OIII]5007 as a template for
#narrow hbeta
OIII = EmissionLine(s,oiii[0],[oiii[1],oiii[2]])
OIII.wv -= OIII.wv.mean()
mdat = LineModel(l,func = 'data',window = w, linedata = OIII)


plt.figure()
#original line in black
plt.plot(l.wv,l.f,'k')
low,u = plt.gca().get_ylim()


#for hbeta
xgrid2 = sp.r_[4920:4960:100j]

plt.plot(xgrid2, mgf(xgrid2),'orange',label='gauss-hermite')
#for empirical template, you have to add wavelength center back in)
plt.plot(OIII.wv + mdat.p[1], mdat(OIII.wv + mdat.p[1]),color='c',label='[OIII] rescaled')

#subtract out the narrow component, to see how good a model it is
mask = (l.wv >= OIII.wv.min()  + mdat.p[1])*(l.wv <= OIII.wv.max() + mdat.p[1] )

l.f[mask] -= mgf(l.wv[mask]) - mgf.p[-1]
plt.plot(l.wv,l.f,'r',label='minus Gauss-Hermite')
l.restore()

l.f[mask] -= mdat(l.wv[mask]) - mdat.p[-1]
plt.plot(l.wv,l.f,'b',label='minus [OIII] rescaled')
l.restore()

plt.legend(loc='upper right')
#Now fit multiple components

#Switch to Halpha, more impressive
window = [ [6336,6950],[6200,6336],[6950,7150] ]
l = EmissionLine(s, window[0],[window[1], window[2]])

#multi components
mg2 = LineModel(l,func = 'gaussian',nline = 5)

#at the moment, these others go negative, need to implement
#constrained fits
#mh2 = LineModel(l,func = 'gauss-hermite',nline = 3)
#ml2 = LineModel(l,func = 'lorentzian',nline = 5)
#mv2 = LineModel(l,func = 'approx-voigt',nline = 5)


#for halpha
xgrid3 = sp.r_[l.wv.min():l.wv.max() + 1:1000j]

plt.figure()
a = plt.gca()
a.plot(l.wv,l.f,'k')
a.plot(xgrid3,mg2(xgrid3),'r',label='gauss')
#a.plot(xgrid3,mh2(xgrid3),'orange',label='hermite')
#a.plot(xgrid,ml2(xgrid),color='r',label='lorentz')
#a.plot(xgrid,mv2(xgrid),color='purple',label='voigt')

mg2.plot_components(xgrid3,a,'r')
#mh2.plot_components(xgrid3,a,'orange')
#ml2.plot_components(xgrid,a,'r')
#mv2.plot_components(xgrid,a,'purple')
plt.legend(loc='upper right')
plt.show()
