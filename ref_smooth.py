import matplotlib
matplotlib.use('TkAgg') 

import scipy as sp
import matplotlib.pyplot as plt
from spectrum import *
from mapspec import *
from matplotlib.ticker import MultipleLocator
import sys


"""
This will calculate the resolution of a list of spectra and show the
distribution, in order to decide how much to smooth the reference.

Resolution comes from FWHM of Gaussian fits to line profile (assume
[OIII]lambda 5007).

Will do some basic sanity checks to make sure you believe the fits,
and let's you interactively choose how much to smooth the reference.
"""

if len(sys.argv) == 1:
    raise ValueError('need to enter redshift for FWHM calculation')
    

#list of spectra
speclist = sp.genfromtxt('speclist_use',dtype=str)
#line used to calculate resolution ([OIII]lambda 5007, again
window   = sp.genfromtxt('oiii.window')

#reference
sref = TextSpec('ref.txt',style='linear')

lref = EmissionLine(sref,window[0],[window[1],window[2]])

#model as a Gaussian, so that everything is measured in the same way
lmodelref = LineModel(lref,func='gaussian')


centdist = []
dispdist = []


plt.ion()
print 'name         FWHM, FWHM fit'
for spec in speclist:
    s = TextSpec(spec)
    shift0 = get_cc(s.f,sref.f,s.wv,sref.wv)

#    s.wv -= shift0
#    print shift0
    l = EmissionLine(s,window[0],[window[1],window[2]])
    lmodel = LineModel(l,func='gaussian')

    if lmodel.p[2] < 0: 
        lmodel.p[2] = -lmodel.p[2]

    dispdist.append(lmodel.p[2])
    centdist.append(lmodel.p[1])

    fwhm,dum1,dum2,dum3 = l.fwhm( (1. + float(sys.argv[1]))*5007.)
    print spec, fwhm,lmodel.p[2]*2.35


#worth checking if there are large wavelength shifts (where is the
#center of the fit?).  If so, consider uncommenting line 36
plt.hist(centdist)
plt.gca().set_title('Distribution of line centers')
plt.draw()
print 'Check that the center of the fits are close enough together that you trust the inferred line widths'
raw_input('Press enter to continue')

plt.clf()
#convert sigma to FWHM
res = lmodelref.p[2]*2.35
print 'parameters of referenc fit: (scale, center, width):',lmodelref.p

#check how good the gauasian fit is
plt.plot(lref.wv,lref.f,'ko-',label='reference data')
plt.plot(lref.wv,lmodelref(lref.wv),'bo-',label='model')
plt.legend(loc='upper left',frameon=0)
plt.gca().set_title('Reference line profile')
plt.draw()
print 'Check that you believe the fit the reference line profile--should only be good enough'
raw_input('Press enter to continue')

plt.clf()
dispdist = sp.array(dispdist)*2.35
plt.hist(dispdist)
l,u = plt.gca().get_ylim()
plt.plot([res,res],[l,u],'r-')
plt.gca().set_xlabel('FWHM (angstroms)')
plt.gca().set_title('Line Model FWHM')
plt.draw()

print 'red line shows native reference resolution:',res
cut = float(raw_input('Where to cut resolution distribution? (type "man" to enter smoothing width by hand)'))

plt.clf()

#smoothing width to get to max resolution below the cut
m = dispdist < cut
newres = sp.sqrt(dispdist[m].max()**2 - res**2 )/2.35
print 'max resolution above the cut:', dispdist[m].max()
print 'smoothing width (pixels):',newres*2.35, '('+str(2.35*newres/(s.wv[1] - s.wv[0]))+')'

#do the smoothing
sref.smooth(13,name = ('gaussian', newres/(s.wv[1] - s.wv[0]) ) )
lref = EmissionLine(sref,window[0],[window[1],window[2]])
lmodel = LineModel(lref,func='gaussian')
res2 = lmodel.p[2]*2.35

print 'new reference resolution:',res2
#plt.hist(dispdist,edgecolor='k',histtype='step')
plt.hist(dispdist)

l,u = plt.gca().get_ylim()
plt.plot([res,res],[l,u],'r',lw=2,label='Native Reference')
plt.plot([res2,res2],[l,u],'g',lw=2,label='Smoothed Reference')
plt.legend(loc='upper left',frameon=0)
plt.gca().set_xlabel('FWHM (angstroms)')
plt.gca().set_title('Line Model FWHM')
plt.draw()
raw_input('Press enter to to save reference, else Ctr+C to quit')




sp.savetxt('ref.smooth.txt',sp.c_[sp.real(sref.wv),sp.real(sref.f),sp.real(sref.ef)])

#print 'errors after here'
#print res,cut,dispdist[m].max(),newres,res2
fout = open('ref_resolution.dat','w')
fout.write('Ref native resolution:    % 2.4f\n'%res)
fout.write('Cut for worst resolution: % 2.4f\n'%cut)
fout.write('Worst object  below cut:  % 2.4f\n'%dispdist[m].max())
fout.write('Kernel width:             % 2.4f\n'%newres*2.35)
fout.write('New ref resolution:       % 2.4f\n'%res2)
fout.close()


