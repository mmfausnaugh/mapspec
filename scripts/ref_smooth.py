##import matplotlib
##matplotlib.use('TkAgg') 

import scipy as sp
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')

from mapspec.spectrum import *
from mapspec.mapspec import *


"""
This code will calculate the resolution of a list of spectra and show
the distribution, in order to decide how much to smooth a reference
spectrum to match some target value.

Resolution is measured by the FWHM of Gaussian fits to the line
profile.

'sref' is the spectrum you want to smooth, read in from a file
'infile.txt'.

'speclist' is the list of spectra to compare.  The code loops
through these, gets the FWHM of a user-specified emission line, and
shows the distribution of FWHMs.  

The user then enters a number, above which all FWHM values are
cut-----it is assumed these are outliers (no problem to set the cut
above the highest FWHM in the distribution).

The highest FWHM below the cut value is now the 'target'

Finally, the code will smooth 'sref' by a Gaussian kernel.  The width
of the kernel is:

kernel.FWHM**2 = target.FWHM**2 - sref.FWHM**2

The user must enter a file designating the window from which to
extract the emission line.  This follows the usual format, i.e., 3
lines specifying:

line_blue_edge     line_red_edge
bluecont_blue_edge bluecont_red_edge
redcont_blue_edge  redcont_red_edge

The code also does some basic sanity checks to make sure you believe
the fits, and let's you interactively specify how much to smooth
'sref', if you prefer.
"""

if len(sys.argv) == 1:
    print('Usage:')
    print('python  ref_smooth.py   infile.txt  speclist   window outfile.txt')
    print('infile.txt---  input spectrum to smooth (usually the reference for alignment)')
    print('speclist--  1 col ascii file with list of files to compare---determines target resolution')
    print('window---  window file designating wavelengths for the EmissionLine')
    print('outfile.txt---  output spectrum, after smoothing')
    sys.exit()

#reference
sref = TextSpec(sys.argv[1],style='linear')
#list of spectra
speclist = sp.genfromtxt(sys.argv[2],dtype=str)
#line used to calculate resolution
window   = sp.genfromtxt(sys.argv[3])
print(window)

lref = EmissionLine(sref,window[0],[window[1],window[2]])
lcenter = lref.wv_mean()

#model as a Gaussian, so that everything is measured in the same way
lmodelref = LineModel(lref,func='gaussian')


centdist = []
dispdist = []


plt.ion()
print('name         FWHM, FWHM fit')
for spec in speclist:
    s = TextSpec(spec)

#    shift0 = get_cc(s.f,sref.f,s.wv,sref.wv)
#    s.wv -= shift0
#    print shift0
    l = EmissionLine(s,window[0],[window[1],window[2]])
    lmodel = LineModel(l,func='gaussian')

    if lmodel.p[2] < 0: 
        lmodel.p[2] = -lmodel.p[2]

    dispdist.append(lmodel.p[2])
    centdist.append(lmodel.p[1])

    fwhm,dum1 = l.fwhm(lcenter)
    print(spec, fwhm,lmodel.p[2]*2.35)



#convert sigma to FWHM
res = lmodelref.p[2]*2.35
print('parameters of referenc fit: (scale, center, width):',lmodelref.p)

#check how good the gauasian fit is
plt.plot(lref.wv,lref.f,'ko-',label='reference data')
plt.plot(lref.wv,lmodelref(lref.wv),'bo-',label='model')
plt.legend(loc='upper left',frameon=0)
plt.gca().set_title('Reference line profile')
plt.draw()
print('Check that you believe the fit the reference line profile--need only be approximate')
input('Press enter to continue')

plt.clf()

#worth checking if there are large wavelength shifts (where is the
#center of the fit?).
plt.hist(centdist)
plt.gca().set_title('Distribution of line centers')
plt.draw()
print('Check that the center of the fits are close enough together that you trust the inferred line widths')
input('Press enter to continue')


plt.clf()

dispdist = sp.array(dispdist)*2.35
plt.hist(dispdist)
l,u = plt.gca().get_ylim()
plt.plot([res,res],[l,u],'r-')
plt.gca().set_xlabel('FWHM (angstroms)')
plt.gca().set_title('Line Model FWHM')
plt.draw()

print('red line shows native reference resolution:',res)

newres = None
cut = input('Where to cut resolution distribution? (type "m" for "manual" to enter smoothing width by hand)\n')
if 'm' in cut:
    #assumes the user will put in a FWHM
    newres = float(input('Enter smoothing width\n'))/2.35
else:
    cut = float(cut)

plt.clf()

#smoothing width to get to max resolution below the cut
if newres == None:
##    print 'dadada'
    m = dispdist < cut
    while len(dispdist[m]) == 0:
        cut = float(input( 'Try a different cut (cut is below the lowest FWHM)\n'))/2.35
        m = dispdist < cut

    newres = sp.sqrt(dispdist[m].max()**2 - res**2 )/2.35
    print('max resolution below the cut:', dispdist[m].max())
    print('smoothing width (pixels):',newres*2.35, '('+str(2.35*newres/(s.wv[1] - s.wv[0]))+')')

#do the smoothing
newres_pix =  newres/(s.wv[1] - s.wv[0])
sref.smooth(10.*newres_pix, name = ('gaussian', newres_pix) )
lref = EmissionLine(sref,window[0],[window[1],window[2]])
lmodel = LineModel(lref,func='gaussian')
res2 = lmodel.p[2]*2.35

print('new reference resolution:',res2)
#plt.hist(dispdist,edgecolor='k',histtype='step')
plt.hist(dispdist)

l,u = plt.gca().get_ylim()
plt.plot([res,res],[l,u],'r',lw=2,label='Native Reference')
plt.plot([res2,res2],[l,u],'g',lw=2,label='Smoothed Reference')
plt.legend(loc='upper left',frameon=0)
plt.gca().set_xlabel('FWHM (angstroms)')
plt.gca().set_title('Line Model FWHM')
plt.draw()
input('Press enter to to save reference, else Ctr+C to quit')




sp.savetxt(sys.argv[4],sp.c_[sp.real(sref.wv),sp.real(sref.f),sp.real(sref.ef)])

#print 'errors after here'
#print res,cut,dispdist[m].max(),newres,res2
fout = open('ref_resolution.dat','w')
fout.write('Ref native resolution:    % 2.4f\n'%res)
if 'm' in str(cut):
    fout.write('Cut for worst resolution: % manual\n')
    fout.write('Worst object  below cut:  % none\n')
else:
    fout.write('Cut for worst resolution: % 2.4f\n'%cut)
    fout.write('Worst object  below cut:  % 2.4f\n'%dispdist[m].max())

fout.write('Kernel width:             % 2.4f\n'%(newres*float(2.35)) )
fout.write('New ref resolution:       % 2.4f\n'%res2)
fout.close()


