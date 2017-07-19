import scipy as sp
from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.signal import get_window
from scipy import linalg,optimize
from numpy.polynomial.hermite import Hermite as H

from astropy.io import ascii,fits
from astropy.table import Table,Column

import re

from .sinc_interp import SincInterp
from .bsplines import Bspline

from copy import deepcopy

import matplotlib.pyplot as plt

__all__ = ['linear_interp_error','extinction','Spectrum','EmissionLine','LineModel','TextSpec','TextSpec_2c','FitsSpec']



class Spectrum(object):
    """
    A Spectrum class.

    In/out is handled by inheritance.  Ascii files (2 or 3 column) are
    quite general.  TextSpec uses scipy.genfromtxt.  fits files come
    in a wide variety, and so are defined individually.  fits files
    need astropy.io.fits.

    This object can do arbitrary rebinning, smoothing, and
    interpolation.  It also facilitate extinction, and line analysis.

    
    """
    def __init__(self,style = 'linear'):
        self._wv  = None
        self._f   = None
        self._ef  = None
        self.sky = None

        self.wv_orig = deepcopy(self.wv)
        self.f_orig  = deepcopy(self.f)
        self.ef_orig = deepcopy(self.ef)
        self.sky_orig= deepcopy(self.sky)
            
        #interpolator style
        self.style=style


    @property
    def wv(self):
        return self._wv

    #anytime an attribute is changed, reset the interpolation to match
    @wv.setter
    def wv(self,wvnew):
        self._wv  = wvnew
        if self.f is not None:
            if self.f.size == self._wv.size:
                self.set_interp(style = self.style)

    @property
    def f(self):
        return self._f

    @f.setter
    def f(self,fnew):
        self._f  = fnew
        self.set_interp(style = self.style)

    @property
    def ef(self):
        return self._ef

    @ef.setter
    def ef(self,efnew):
        self._ef  = efnew
        self.set_interp(style = self.style)



    def restore(self,verbose=True):
        if verbose:
            print('Restore data to default')
        #restore to original wv, f, and ef
        self.wv  = deepcopy(self.wv_orig)
        self.ef  = deepcopy(self.ef_orig)
        self.f   = deepcopy(self.f_orig)
#        self.sky = self.sky_orig

    def redistribute(self):
        """
        Perturb the spectrum by random gaussian deviates scaled to the
        error spectrum.
        """
        self.f += self.ef*sp.randn(self.ef.size)

    def interp(self,xnew):
        """
        Interpolate the spectrum at arbitrary points xnew
        """
        if self.ef is not None:
            zout = self._interpolator_error(xnew)
            if (zout <= 0).any():
                print('Warning:  some variances are negative, flagged with complex errors')
            return self._interpolator(xnew),sp.sqrt(zout)
        else:
            print('Warning:  No error spectrum, return array of ones')
            return self._interpolator(xnew),sp.ones(xnew.size)

    def set_interp(self,style='sinc',window1='lanczos',kw1=15,order1=3):
        """
        Choose which interpolation technique to use.  Relevent for
        self.interp() and self.rebin()


        Note that errors are only rigorous for linear interpolation.
        Otherwise, for now it will interpolate the squared error
        spectrum using the method of choice (i.e., splines, sinc,
        etc.).

        Uses scipy.interp1d---anything that can be pased to 'kind'
        keyword in this package is valid.
        """
        self.style = style
        if style == 'sinc':
            self._interpolator = SincInterp(self.wv,self.f, window=window1,kw=kw1)
            if self.ef is not None:  self._interpolator_error = SincInterp(self.wv,self.ef**2, window=window1,kw=kw1)

        elif style == 'bspline':
            self._interpolator = Bspline(self.wv,self.f, order=order1)
            if self.ef is not None:  self._interpolator_error = Bspline(self.wv,self.ef**2,order=order1)

        elif style == 'linear':
            self._interpolator = interp1d(self.wv,self.f,kind=style)
            if self.ef is not None:  self._interpolator_error = lambda xnew: linear_interp_error(self.wv,xnew, self.ef)

        else:
            self._interpolator = interp1d(self.wv,self.f,kind=style)
            if self.ef is not None:  self._interpolator_error = interp1d(self.wv,self.ef**2,kind=style)

    def rebin(self,xnew):
        """
        Rebin the spectrum on a new grid named xnew
        """

        #Does not need equal spaced bins, but why would you not?
        xnew.sort()

        fbin  = sp.zeros(xnew.size)
        efbin = sp.zeros(xnew.size)

        #up sampling is just interpolation
        m = (self.wv >= xnew[0])*(self.wv <= xnew[-1])
        if self.wv[m].size <= xnew.size - 1:
            fbin,efbin  = self.interp(xnew)
            
        else:
        #down sampling--
        #1) define bins so that xnew is at the center.
        #2) interpolate to account for fractional pixel weights
        #3) take the mean within each bin
            db  = 0.5*sp.diff(xnew)
            b2  = xnew[1::] - db
            b2  = sp.insert(b2,0,xnew[0])

            insert = sp.searchsorted(self.wv,b2)
            xinsert = sp.insert(self.wv,insert,xnew)
            xinsert = sp.unique(xinsert)
            yinsert,zinsert = self.interp(xinsert)

            i = sp.digitize(xinsert,b2)
            for j in range(b2.size):
                iuse = sp.where(i == j+1)[0]
                fbin[j]  = sp.mean(yinsert[iuse])
                efbin[j] = sp.mean(zinsert[iuse])

        self._wv = xnew
        if self.ef is not None:        
            self._ef = efbin            
        self.f = fbin
        assert self.wv.size == self.f.size

    def extinction_correct(self,E_BV,RV,u='microns'):
        """
        Uses CCM89 to extinction correct the spectrum.  You need to
        choose the reddening (which traces the hydrogen/dust column)
        and the extinction law (RV)
        """

        AV = RV*E_BV*extinction(self.wv,RV,unit=u)
        if self.ef is not None:
            self.ef /= 10**(-0.4*AV)
        self.f  /= 10**(-0.4*AV)

    def smooth(self,width,name='boxcar'):
        """
        Smooth a spectrum.

        width = total number of pixels in the kernel 

        name = a string used by sp.signal.get_window.  Things like
        ('gaussian' 1.5) are ok, see scipy documentation. 

        Because even widths are asymmetric, this only allows odd
        kernels.
        """
        if width%2 == 0:
            raise ValueError("Only allows odd widths (even kernels are asymetric)")
        W = get_window(name,width)
        W /= abs(sp.sum(W))

        fsmooth = sp.convolve(self.f,W,mode='same')
        #treat edge effects by replacing with original spectrum
        s1 = slice(0,width)
        s2 = slice(- (width ),None)

        fsmooth[s1] = self.f[s1]
        fsmooth[s2] = self.f[s2]

        if self.ef is not None:
            efsmooth = sp.sqrt(sp.convolve(self.ef**2,W**2,mode='same'))
            efsmooth[s1] = self.ef[s1]
            efsmooth[s2] = self.ef[s2]
            self.ef = efsmooth

        self.f = fsmooth



    def velocity_smooth(self,v_width):
        """ 
        Smooths with constant velcoity dispersion by resampling on
        even log intervals, smooth with gaussian of specific width,
        and then resampling to original wavelengths.

        Because log spacing is uneven, cannot use sinc interpolation
        (will use bsplines instead).

        Gaussian width is specified in km/s.
        """
        #get y for evenly spaced in log x
        xorig = deepcopy(self.wv)
        lwv = sp.log(self.wv)
        lognew = sp.r_[lwv.min():lwv.max():1j*lwv.size]
        xnew = sp.exp(lognew)
        xnew[0] = xorig[0]
        xnew[-1] = xorig[-1]

        #since bins will be uneven, cannot use sinc
        flag = 0
        if self.style=='sinc':
            flag = 1
            print('Warning:  Setting interp style to "bspline" since log spacing is uneven')
            self.set_interp(style='bspline')
            self.rebin(xnew)

        #assumes v_width in km/s
        dpix = v_width/2.998e5/(lognew[1] - lognew[0])
        #kernel width goes out to 5 sigma
        kw = round(dpix*10)
        if kw%2 == 0:
            kw += 1
        W = get_window(('gaussian', dpix),kw)
        W /= abs(sp.sum(W))

        fsmooth = sp.convolve(self.f,W,mode='same')

        s1 = slice(0,kw)
        s2 = slice(- (kw),None)
        fsmooth[s1] = self.f[s1]
        fsmooth[s2] = self.f[s2]


        if self.ef is not None:
            efsmooth = sp.sqrt(sp.convolve(self.ef**2,W**2,mode='same'))
            efsmooth[s1] = self.ef[s1]
            efsmooth[s2] = self.ef[s2]
            self.ef = efsmooth

        self.f = fsmooth
        
        self.rebin(xorig)

        if flag ==1:
            self.set_interp(style='sinc')


class EmissionLine(Spectrum):
    """
    Emission Line object.  Takes the flux of a Spectrum object within
    a given window to be an emission line.  Models the local continuum
    with a linear fit.  The local continuum is subtracted, but saved
    as an attribute.

    Standard error propogation is implemented, keeping track of the
    line flux error and continuum flux error.

    Parameters
    ----------
    Spec    = Spectrum Object

    window  = wavlength window of the emission line

    dwindow = extentions of window (in wavelengths), used to estimate
              the continuum

    """
    def  __init__(self,Spec,window,cwindow):
        super(EmissionLine,self).__init__()  
        self.style= Spec.style
        self.cwindow = cwindow
        ml = (Spec.wv >= window[0])*(Spec.wv <= window[1])
        mc  = (Spec.wv > cwindow[0][0])*(Spec.wv < cwindow[0][1])
        mc += (Spec.wv > cwindow[1][0])*(Spec.wv < cwindow[1][1])

        self.wv = Spec.wv[ml]

        wvfit = Spec.wv[mc] - Spec.wv[mc].mean()
        ffit  = Spec.f[mc]  - Spec.f[mc].mean()

        cfit,covar = linfit(wvfit,ffit,Spec.ef[mc])

        wvsub = self.wv - Spec.wv[mc].mean()

        self.cf  = cfit[0] + cfit[1]*wvsub
        self.cf += Spec.f[mc].mean()
        self.ecf = sp.sqrt(covar[0,0] + covar[1,1]*wvsub + 2*covar[0,1]*wvsub)

        self.f  = Spec.f[ml] - self.cf
        self.ef = sp.sqrt(Spec.ef[ml]**2 + self.ecf**2 )

        self.wv_orig = deepcopy(self.wv)
        self.f_orig  = deepcopy(self.f)
        self.ef_orig = deepcopy(self.ef)
#        self.sky_orig= deepcopy(self.sky)

        self._get_cumulative_flux()


    def _get_cumulative_flux(self):
        """
        Store a look-up table for the cumulative flux as a function of
        wavelength.  This seems to be an efficient way of getting
        things like velocity percentiles.
        """
        cflux = [0]
        for w in self.wv[1::]:
            m = self.wv <= w
            fint = simps(self.f[m], self.wv[m])
            cflux = sp.r_[cflux,fint]
        self._cumulative_flux = interp1d(cflux,self.wv)

    def integrate_line(self):
        """
        Integrate the line using Simpson's method.
        """
        return simps(self.f,self.wv)

    def equivalent_width(self):
        """
        Calculate the equivalent width, using the linear continuum fit
        under the line
        """
        ltot = simps(self.f,self.wv)
        cmean  = sp.sum(self.cf/self.ecf**2)/sp.sum(1./self.ecf**2)
        ecmean = sp.sqrt(1./sp.sum(1./self.ecf**2))
        return ltot/cmean
        
    def wv_mean(self):
        """
        First moment of the line, i.e., mean wavelength
        """
        moment = simps(self.wv*self.f,self.wv)
        norm = self.integrate_line()
        return moment/norm

    def wv_median(self):
        """
        50th percentile of the flux, i.e., median wavelength
        """
        ftarget = 0.5*self.integrate_line()
        return self._cumulative_flux(ftarget)


    def dispersion(self):
        """
        Second moment of the line, i.e., dispersion
        """
        l0 = self.wv_mean()
        moment = simps(self.wv**2*self.f,self.wv)
        norm   = self.integrate_line()

        return sp.sqrt(moment/norm - l0**2)

    def mad(self):
        """
        MAD = mean absolute deviation.  Defined same as dispersion,
        but uses absolute value instead of square.
        """
        wv_abs = sp.absolute(self.wv - self.wv_median())
        moment = simps(wv_abs*self.f,self.wv)
        norm = self.integrate_line()
        return moment/ norm

    def fwhm(self, center):
        """
        Full-width at half maximum of the line.  See Peterson et
        al. 2004 for a description of the algorith.

        center = center wavelength of the line.  This should match
        whatever data was read into the Spectrum/EmissionLine
        object---i.e., no redshifting is done, so you must feed it the
        observed/rest frame line wavelength as appropriate.

        Keeps track of red and blue side because of the possibility of
        double peaked lines.  A flag for double-peaked is also
        returned.
        """
        #need to give the center of the line, to decide if
        #double-peaked or not
        rmask = self.wv >=center
        bmask = self.wv < center

        #check double peaked
        bmax   = self.f[bmask].max()
        bmax_i = sp.where( self.f == bmax)[0]
        rmax   = self.f[rmask].max()
        rmax_i = sp.where( self.f == rmax)[0]

        wvmask = (self.wv > self.wv[bmax_i] )*(self.wv < self.wv[rmax_i] )
        if (self.f[wvmask].size > 0 ) and ( self.f[wvmask] < min(bmax, rmax)  ).any():
            doublepeaked = True
        else:
            doublepeaked = False

        if doublepeaked == False:            
            if rmax > bmax:
                bmax = rmax
                bmax_i = rmax_i
            else:
                rmax = bmax 
                rmax_i = bmax_i

        #blue side
        ftarget = 0.5*bmax
        for i in range(bmax_i + 1):
            if self.f[i] > ftarget:
                if i == 0:
                    print('Warning: edge of the blue wing > 0.5*blue_max')
                    b1 = self.wv[0]
                    break
                else:
                    slope = ( self.wv[i] - self.wv[i - 1] )/( self.f[i] - self.f[i - 1] )
                    b1 = self.wv[i - 1] + slope*(ftarget - self.f[i - 1])
                    break

        for i in range(bmax_i + 1):
            j = bmax_i - i
            if j == 0:
                print('Warning: hit edge of the blue wing')
                b2 = [self.wv[0]]
                break
            if self.f[j] < ftarget:
                slope = ( self.wv[j] - self.wv[j + 1] )/( self.f[j] - self.f[j + 1] )
                b2 = self.wv[j + 1] + slope*(ftarget - self.f[j + 1])
                break

        fwhm_blue = (b1 + b2)/2.

        #red side
        ftarget = 0.5*rmax
        for i in range(self.wv.size  - rmax_i):
            j = self.wv.size - 1 - i
            if self.f[j] > ftarget:
                if j == self.wv.size - 1:
                    print('Warning: edge of the red wing > 0.5*red_max')
                    r1 = self.wv[-1]
                    break
                else:
                    slope = ( self.wv[j] - self.wv[j + 1] )/( self.f[j] - self.f[j + 1] )
                    r1 = self.wv[j + 1] + slope*(ftarget - self.f[j + 1])
                    break

        for i in range(self.wv.size  - rmax_i):
            j = rmax_i + i
            if j == self.wv.size - 1:
                print('Warning: hit the edge of the red wing')
                r2 = [self.wv[-1]]
                break
            if self.f[j] < ftarget:
                slope = ( self.wv[j] - self.wv[j - 1] )/( self.f[j] - self.f[j - 1] )
                r2 = self.wv[j - 1] + slope*(ftarget - self.f[j - 1])
                break

        fwhm_red = (r1 + r2)/2.

        outdic = {'blue':(fwhm_blue - center)[0],
                  'red'  :(fwhm_red - center)[0],
                  'doublepeaked':doublepeaked
                      }
        return ( (fwhm_red - fwhm_blue)[0], outdic)

    def ip_wv(self,frac):
        """
        Returns interpercentile range (of flux) in terms of
        wavelengths.

        frac should be a decimal
        """
        f1 = 0.5*frac*self.integrate_line()
        f2 = (1. - 0.5*frac)*self.integrate_line()

        w1 = self._cumulative_flux(f1)
        w2 = self._cumulative_flux(f2)

        return (w2 - w1, w1, w2)

    def ipv(self,frac,center):
        """
        Returns interpercentile range (of flux) in terms of
        velocity.

        frac should be a decimal

        center is used to directly covner to velocity (returns km/s),
        and follows the same perscription as EmissionLine.fwhm()
        """
        dlam, lam1, lam2 = self.ip_wv(frac)
        #assumes km/s
        dvel = dlam/center*2.997925e5 
        dvel_blue = (lam1 - center)/center*2.997925e5 
        dvel_red  = (lam2 - center)/center*2.997925e5 
        return (dvel, dvel_blue, dvel_red)

class LineModel(EmissionLine):
    """
    Different analytic models, which are fit to a given emission line.

    An option is available to restrict the fit to a part of the line
    (window keyword).

    Can also do multiple components, although the model must be
    homogenous (lines keyword).

    Can 'float' the lines, by giving them a constant offset.  Say you
    only wanted to get rid of the narrow component, for example (set
    floating to True).

    One model is semi-empirical----you can give it an emission line
    object, which will be used as a template, and this will shift,
    rescale, and float the empircal line.

    Note that the functional models are defined after this class.
    """

    def __init__(self,EmLine,func = 'gaussian',window = None,nline=1,floating=False,linedata = None):
        super(EmissionLine,self).__init__()  
        self.wv   = EmLine.wv
        self.lf   = EmLine.f
        self.elf  = EmLine.ef

        self.floating = floating

        if window != None:
             x  = EmLine.wv[window]
             y  = EmLine.f[window] #- EmLine.lf[m].min()
             z  = EmLine.ef[window]
        else:
            x = EmLine.wv
            y = EmLine.f
            z = EmLine.ef

#        self.fuse,pinit = self._init_params(func,x,y,nline,linedata)
        #all the heavy lifting is in this step

        self.fuse,pinit,pbounds = self._init_params(func,x,y,nline,linedata)
#        print pbounds,pinit
        #all the heavy lifting is in this step
#        print pbounds
#        print pinit
        self.p,self.covar = fitfunc(self.fuse,pinit,x,y,z)
#        self.p,dummy = fitfunc_bound(self.fuse,sp.array(pinit),x,y,z,pbounds)
#        print self.p


    def _init_params(self,func,x,y,nline,*argv):
        self.func = func
        if func == 'gaussian':
            pinit = []
            for i in range(nline):
                pinit.append(y.mean())
                pinit.append(x.mean())
                pinit.append( sp.absolute(x.max() - x.mean())/3.)
                if self.floating == True:
                    pinit.append(0.0)

            if self.floating == True:
                fuse = lambda x,p: multiwrapper(floating(gauss),x,p,nparams[func] + 1)
                self.fuse2=floating(gauss)
            else:
                fuse = lambda x,p: multiwrapper(gauss,x,p,nparams[func])
                self.fuse2=gauss

        if func =='gauss-hermite':
            pinit = []
            for i in range(nline):
                pinit.append(y.mean())
                pinit.append(x.mean())
                pinit.append( sp.absolute(x.max() - x.mean())/3.)
                pinit.append(0.0)
                pinit.append(0.0)
                if self.floating == True:
                    pinit.append(0.0)

            if self.floating == True:
                fuse = lambda x,p: multiwrapper(floating(gauss_hermite),x,p,nparams[func] + 1)
                self.fuse2=floating(gauss_hermite)
            else:
                fuse = lambda x,p: multiwrapper(gauss_hermite,x,p,nparams[func])
                self.fuse2=gauss_hermite


        if func =='lorentzian':
            pinit = []
            for i in range(nline):
                pinit.append(y.mean())
                pinit.append(x.mean())
                pinit.append( sp.absolute(x.max() - x.mean())/3.)
                if self.floating == True:
                    pinit.append(0.0)

            if self.floating == True:
                fuse = lambda x,p: multiwrapper(floating(lorentzian),x,p,nparams[func] + 1)
                self.fuse2=floating(lorentzian)
            else:
                fuse = lambda x,p: multiwrapper(lorentzian,x,p,nparams[func])
                self.fuse2=lorentzian


        if func =='approx-voigt':
            pinit = []
            for i in range(nline):
                pinit.append(y.mean())
                pinit.append(x.mean())
                pinit.append( sp.absolute(x.max() - x.mean())/3.)
                pinit.append( sp.absolute(x.max() - x.mean())/5.)
                pinit.append(1.0)
                if self.floating == True:
                    pinit.append(0.0)

            if self.floating == True:
                fuse = lambda x,p: multiwrapper(floating(approx_voigt),x,p,nparams[func] + 1)
                self.fuse2=floating(approx_voigt)
            else:
                fuse = lambda x,p: multiwrapper(approx_voigt,x,p,nparams[func])
                self.fuse2=approx_voigt

        if func == 'data':
            #is there every a case not to give an offset for an
            #empirical line?
            pinit = [ 1.0,x.mean(),0.0 ]
            l = argv[0]
            fuse = lambda x,p : empiriline(x,p,l)

#        return fuse,pinit

#        if func =='voigt':
#            pinit = []
#            for i in range(nline):
#                pinit.append(y.mean())
#                pinit.append(x.mean())
#                pinit.append( sp.absolute(x.max() - x.mean())/3.)
#                pinit.append( sp.absolute(x.max() - x.mean())/5.)
#
#            fuse = MultiVoigt



#make sure there are no negative fluxes
#but can't get scipy constraint/bound system to work....        
        pbound = []
        for i in range(len(pinit)):
            #all the scales are above 0
            if i%nparams[func] == 0:
                pbound.append( (0,None) )
            else:
                pbound.append( (None,None))

        
        return fuse,pinit,tuple(pbound)
#        return fuse,pinit

 
    def __call__(self,x):
        return self.fuse(x,self.p)

    def chi2(self):
        return sp.sum( (self.lf - self.fuse(self.wv,self.p) )**2/self.elf**2 )

    def plot_components(self,x,ax,cuse):
        np = nparams[self.func]
        puse = []
        for i in range(len(self.p)):
            puse.append(self.p[i])
            if len(puse)%np == 0:
                y = self.fuse2(x, puse )
                ax.plot(x,y,color=cuse)
                puse = []

#stores number of paramters for each model
nparams = {'gaussian':3,
           'gauss-hermite':5,
           'lorentzian':3,
           'approx-voigt':5,
           'data':3}

def multiwrapper(func,x,p,np):
    """
    If you want several components, this will iteratively add them.
    """
    puse = []
    y = 0
    for i in range(len(p)):
        puse.append(p[i])
        if len(puse)%np == 0:
            y += func(x, puse )
            puse = []
    return y

def floating(func):
    """
    If you want a constant offset, this will modify any input function.
    """
    
    def floatfunc(x,p):
        return func(x,p[0:-1]) + p[-1]
    return floatfunc

#for feeding a line template
def empiriline(x,p,L):
    """
    Use the line L (which is an EmissionLine object) as a template.
    The line is shifted, then interpolated, then rescaled, and allowed
    to float.
    """
    xnew = x - p[1]
    m = (xnew >= L.wv.min())*(xnew <= L.wv.max() )
    ynew,znew = L.interp(xnew[m])
    return p[0]*ynew + p[2]

#functions for use with Line Model
def gauss(x,p):
    return p[0]*sp.exp(-0.5*(x - p[1])**2/p[2]**2)


def gauss_hermite(x,p):
    #The constants are only to make p[3] and p[4] correspond to a
    #rigorous definition of h3 and h4
    A = p[0]*(p[2]/sp.sqrt(2*sp.pi))
    w = (x - p[1])/p[2]
    h = H([ 1., 0., 0., p[3], p[4] ])
    
    return A*sp.exp(-0.5*w**2)*h(w)

def lorentzian(x,p):
    w = (x - p[1])/p[2]
    return p[0]/(1 + w**2)


def approx_voigt(x,p):
    #just sum the gaussian and lorentz
    w1 = (x - p[1])/p[2]
    w2 = (x - p[1])/p[3]
    g = sp.exp(-0.5*w1**2)
    l = 1./(1 + w2**2)

    return p[0]*(g + p[4]*l)



####testing
####might be too complicated for levenberg-marquadt

#def voigt(x,p):
#
#    print p
#    pwidth = x[1] - x[0]
#
#    prange1 = sp.r_[- (5*int(p[2]/pwidth)):5* int(p[2]/pwidth) + 1]
#    g = sp.exp(-0.5*prange1**2/pwidth**2)
#    g /= sp.sum(g)
#
#    prange2 = sp.r_[- abs(5*int(p[3]/pwidth)):5* abs(int(p[3]/pwidth)) + 1]
#    l = 1./(1 + prange2**2/(p[3]/pwidth)**2)
#    l /= sp.sum(l)
#    print g.size,l.size,prange2.size,prange1.size
#    v = sp.convolve(g,l,mode='same')
#    print v.size,prange2.size,prange1.size
#
#    if prange1.size > prange2.size:
#        puse = prange1
#    else:
#        puse = prange2
#    V = interp1d(puse,v)
#
#    w = (x-p[1])/pwidth
#
#    print puse,w
#    plt.plot(prange1*pwidth,g,'b')
#    plt.plot(prange2*pwidth,l,'r')
#    plt.plot(puse*pwidth,v,'k')
#    plt.plot( w*pwidth,V(w),'ko')
#    plt.show()
#
#    return p[0]*V(w)



class TextSpec(Spectrum):
    """
    A user level class, which will in intialize a Spectrum object from
    a 3 column ascii file----this is actually acomplished through
    inheritence.

    The user only needs to specify the file name, although they can
    change the interpolation method (linear, sinc, bsplines, etc.)
    here.
    """
    def __init__(self,ifile,style='linear',delimiter=None):
        super(TextSpec,self).__init__()  
        self.style=style
        x,y,z = sp.genfromtxt(ifile,unpack=1,usecols = (0,1,2), delimiter=None)
        self.wv = x
        self.f  = y
        self.ef = z

        self.wv_orig = deepcopy(self.wv)
        self.f_orig  = deepcopy(self.f)
        self.ef_orig = deepcopy(self.ef)
#            self.sky_orig= deepcopy(self.sky)
    def save(self,ofile):
        sp.savetxt(ofile, sp.c_[self.wv,self.f,self.ef])


class TextSpec_2c(Spectrum):
    """
    As TextSpec, but for 2 column ascii files.  The error spectrum is
    is set to a dummy array of ones.
    """

    def __init__(self,ifile,style='linear',delimiter=None):
        super(TextSpec_2c,self).__init__()  
        self.style=style
        x,y = sp.genfromtxt(ifile,unpack=1,usecols = (0,1),delimiter=None)
        self.wv = x
        self.f  = y
#            self.ef = sp.ones(y.size)

        self.wv_orig = deepcopy(self.wv)
        self.f_orig  = deepcopy(self.f)
        self.ef_orig = deepcopy(self.ef)

    def save(self,ofile):
        sp.savetxt(ofile, sp.c_[self.wv,self.f])

class FitsSpec(Spectrum):
    """    
    A user level class, which will in intialize a Spectrum object from
    a fits file.  Fits file is assumed to be a mutli-extension 1D
    array.  Uses astropy to read thing in.

    One should look in the fits header for how to define the
    wavelengths.  The python kewords give:

    x1    = the initial wavelength (this is assigned to pixel from p1key, see below).
    dxkey = Wavelength spacing, i.e., delta x
    p1key = The pixel where the wavelength solution starts

    The default keywords were chosen to match MDM spectra, which
    follow usual IRAF conventions.  However, these conventions are
    fairly loose, it's likely you'll have to change them for data from
    different spectrographs.
    """
    def __init__(self,ifile, extension = 0,data_axis = 1, error_axis =3, x1key ='CRVAL1', dxkey='CD1_1', p1key='CRPIX1',style='linear'):
        super(FitsSpec,self).__init__()  

        data  = fits.getdata(ifile,extension)
        self._f  = data[data_axis,0,:]
        self._ef = data[error_axis,0,:]

        x1 = fits.getval(ifile,x1key)
        p1 = fits.getval(ifile,p1key)
        dx = fits.getval(ifile,dxkey)

        self.wv = self._get_x(x1,p1,dx)
        self.f = self._f
        self.ef = self._ef

    def _get_x(self,x1,p1,dx):
        #helper for getting wavelengths from fits headers
        n = self._f.size
        xhigh = (n - p1 + 1)*dx + x1
        xlow  =  x1 - dx*(p1 - 1)
        return sp.r_[xlow:xhigh:dx]

    def save(self,ofile,head=None, x1key ='CRVAL1', dxkey='CD1_1', p1key='CRPIX1'):

        data = sp.array([
                [self.f],
                [self.ef]
                ])
    
        if head == None:
            fits.writeto(ofile,data,clobber=True)

        else:
            head[x1key] = self.wv[0]
            #assume even spacing......
            head[dxkey] = (self.wv[1] - self.wv[0])
            head[p1key] = 1

            head['COMMENT'] = 'Modified by mapspec on %s'%(time.strftime("%c"))

            fits.writeto(ofile,data,header=head,clobber=True)



def linfit(x,y,ey):
    """
    This function minimizes chi^2 for a least squares fit to a simple
    linear model.  The errors in x are assumed to be much smaller than
    the errors in y. 

    The function returns p, and covar.  p[0] is the zeropoint p[1] is
    the slope.  Covar is the covariance matrix---parameter errors are
    the sqrt of its diagonal.
    """
    C = sp.array([sp.sum(y/ey/ey),sp.sum(y*x/ey/ey)])
    A = sp.array([
            [sp.sum(1./ey/ey), sp.sum(x/ey/ey)],
            [sp.sum(x/ey/ey),  sp.sum(x*x/ey/ey)]
            ] )

    p     = linalg.solve(A,C)
    covar = linalg.inv(A)
    return p,covar

def fitfunc(func,pin,x,y,ey):
    """
    Non-linear least square fitter.  Assumes no errors in x.  Utilizes
    scipy.optimize.leastsq.  Uses a Levenberg-Marquardt algorithm
    (gradient chasing and Newton's method).  The user must provide an
    object that is the functional form to be fit, as well as initial
    guesses for the parameters.

    This function returns the parameters and the covariance matrix.
    Note that leastsq returns the jacobian, a message, and a flag as
    well.
    
    """
    merit = lambda params,func,x,y,ey:  (y - func(x,params))/ey
    p,covar,dum1,dum2,dum3 = optimize.leastsq(merit,pin,args=(func,x,y,ey),xtol = 1e-8,full_output=1)
    return p,covar

def fitfunc_bound(func,pin,x,y,ey,pbound):
    """
    Non-linear least square fitter, with an option of contraints

    ***********
    
    """
    merit = lambda params,func,x,y,ey:  sp.sum( 
        ( y - func(x,params) )**2/ey**2 
        )
    #return with weird object....
#########    out,dum1,dum2 = optimize.fmin_l_bfgs_b(merit,
#########                                           pin,
#########                                           args=(func,x,y,ey),
#########                                           #                                 tol = 1.e-15,
#########                                           #                            method='L-BFGS-B',
#########                                           bounds = pbound,
#########                                           approx_grad=True,
#########                                           epsilon=1.e-15,
##########                                           pgtol
#########                                           #                                 options={'disp':True}
#########                                           )
    print(pbound)
    out = optimize.minimize(merit,pin,args=(func,x,y,ey),method='TNC',bounds = pbound,options={'disp':True,'maxiter':10000})
    print('success:',out.success)
#    print out.message
#    errors = sp.sqrt(ey**2 * out.jac**2)
#    print out.keys()
#    print 'iterations:',out.nit,out.nfev#, out.njev, out.nhev
#    print 'WARNING!!!',dum2['warnflag'],dum2['task']
#    print out.x
    return out.x,out.jac

def linear_interp_error(x,xinterp,z):
    """
    Does the actual calculation for error propagation on linear interpolation.

    x = xold
    xinterp = grid for interpolatin x
    z = old error spectrum

    Note, returns variance, i.e. (new error spectrum)**2
    """

    i = sp.searchsorted(x,xinterp)
    f = (xinterp - x[i -1 ])/(x[i] - x[i -1])
    return f**2*z[i]**2 + (1 - f)**2*z[i-1]**2

def extinction(lambda1in,R,unit = 'microns'):
    """
    Calculates A(lambda)/A_V.  So, if we know E(B - V), we do
    A(lambda) = A(lambda)/A_V * E(B - V) * R.  

    R is alternatively R_V, usually 3.1---this parameterizes the extinction law, which you should know if you are using this function.

    This is the CCM89 extinction law, which assumes microns.
    """
    if 'ang' in unit:
        lambda1 = lambda1in/1.e4
    else:
        lambda1 = lambda1in

    if (lambda1 > 100).all():
        print("Check units!  This program assumes microns")

    if (lambda1 > 3.0).any():
        print ("Warning: extrapolating into the far IR (lambda > 3 microns)")
    if (lambda1 < 0.125).any():
        print('Warning:  extreme UV is an extrapolation')
    if (lambda1 < 0.1).any():
        print('warning: extrapolating into the extreme UV (lambda < 1000 A)')


    a = sp.zeros(lambda1.size)
    b = sp.zeros(lambda1.size)

    m = (lambda1 > 0.909)
    a[m] =  0.574*(1/lambda1[m])**(1.61)
    b[m] = -0.527*(1/lambda1[m])**(1.61)
        
    m = (lambda1 > 0.30303)*(lambda1 <= 0.909)
    x = 1/lambda1[m] - 1.82
    a[m] = 1 + 0.17699*x - 0.50447*x**2 - 0.02427*x**3 + 0.72085*x**4 + 0.01979*x**5 - 0.7753*x**6 + 0.32999*x**7
    b[m] =     1.41338*x + 2.28305*x**2 + 1.07233*x**3 - 5.38434*x**4 - 0.62251*x**5 + 5.3026*x**6 - 2.09002*x**7

    m = (lambda1 > 0.125)*(lambda1 <= 0.30303)
    x = 1/lambda1[m]
    a[m] =  1.752 - 0.316*x - 0.104/( (x - 4.67)**2 + 0.341) 
    b[m] = -3.090 + 1.825*x + 1.206/( (x - 4.62)**2 + 0.263) 

    m = (lambda1 > 0.125)*(lambda1 <= 0.1695)
    x = 1/lambda1[m]
    a[m] += -0.04473*(x - 5.9)**2 - 0.009779*(x-5.9)**3
    b[m] +=  0.21300*(x - 5.9)**2 + 0.120700*(x - 5.9)**3

    m = (lambda1 < 0.125)
    x = 1/lambda1[m]
    a[m] = -1.073 - 0.628*(x - 8.) + 0.137*(x - 8.)**2 - 0.070*(x - 8.)**3
    b[m] = 13.670 + 4.257*(x - 8.) - 0.420*(x - 8.)**2 + 0.374*(x - 8.)**3

    return a + b/R
