import scipy as sp
from scipy.signal import correlate
import matplotlib.pyplot as plt
from spectrum import Spectrum,EmissionLine,TextSpec
"""
This is for making a reference image.  you provide a file with a list
of photometric nights, the program will align the spectra to a common
wavelength grid (based on wavelength shifts so that the [OIII]lambda
5007 line profiles match, and then average them.

"""

def get_chi2(s1,s2,shift):
    trim = int(shift/(s1.wv[1] - s1.wv[0])) + 1
    xnew = s2.wv[trim : -trim] - shift
    #resample the reference at the new wavelength grid
    y1,z1 = s1.interp(xnew)

    return sp.sum( 
        (y1 - s2.f[trim: - trim])**2/(z1**2 + s2.ef[trim: - trim]**2)
        )

def get_cc(y1,y2,x):
    cc = correlate(y1,y2,mode='same')
    i  = sp.where(cc == cc.max())[0]
    shift = (x[1] - x[0] )*(cc.size//2 - i)

    return shift

def tidy(xout,yout,zout):
    xmin = xout[0]
    for x in xout:
        if x.size < xmin.size:
            xmin = x

    w = sp.zeros(xmin.size)

    for i in range(sp.shape(xout)[0]):
        j = sp.in1d(xout[i],xmin)

        yout[i] = yout[i][j]
        zout[i] = zout[i][j]

    yout  = sp.array(yout)
    zout = sp.array(zout)
    ymean = sp.sum( yout/zout**2,axis = 0 )/sp.sum(1./zout**2, axis = 0)
    error = sp.sqrt( 
        1./sp.sum(1./zout**2,axis = 0) 
        )

    return xmin,ymean,error


def HM(ntrial,s1,s2,p):

    chi2 = get_chi2(s1,s2,p)
    chi2best = 1.e12

    pbest = p

    accept = 0

    for i in range(ntrial):
        if i%10 == 0 :
            print i,chi2best

        ptry = p + sp.randn()*0.1
 
        chi2try = get_chi2(s1,s2,ptry)


        if chi2try < chi2:
            
            p = ptry
            chi2 = chi2try

            accept += 1
            if chi2 < chi2best:
                pbest = ptry
                chi2best = chi2try


        else:
            prob = sp.exp(-chi2try/chi2)
            r = sp.rand()
            if r <= prob:
                p = ptry
                chi2 = chi2try
                accept += 1

    return chi2best,pbest,accept/float(ntrial)


#list of spectra for the reference
reflist = sp.genfromtxt('reflist',dtype='a')


S,L = [],[]
for ref in reflist:
    s = TextSpec(ref)
    s.set_interp(style='bspline')
    S.append(s)
    plt.plot(s.wv,s.f,'k')


trimmax = 0

xout = []
yout = []
zout = []

shiftout = []

#window for oxygen line, see do_map.py and run_map.sh for details
window = sp.genfromtxt('oiii.window')

lref = EmissionLine(S[0],window[0],[window[1],window[2]])
print lref.style

for s in S[1::]:
    
    shift0 = get_cc(S[0].f,s.f,S[0].wv)
    print shift0
    l = EmissionLine(s,window[0],[window[1],window[2]])
    chi,shiftuse,frac = HM(1000,lref,l,shift0)

    print chi,shiftuse,frac

    shiftout.append(shiftuse)
    s.wv -= shiftuse

    trim = int(shiftuse/(s.wv[1] - s.wv[0])) + 1
    if trim > trimmax: trimmax = trim
    print trimmax
    
    y1,z1 = s.interp(S[0].wv[trim:-trim])

    xout.append(S[0].wv[trim:-trim])
    yout.append(y1)
    zout.append(z1)


xout.append(S[0].wv)
yout.append(S[0].f)
zout.append(S[0].ef)

xref,yref,zref = tidy(xout,yout,zout)


sp.savetxt('ref.txt',sp.c_[xref,yref,zref])
sp.savetxt('shift.params',sp.c_[shiftout])
