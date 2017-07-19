import scipy as sp
from scipy.signal import get_window 
from scipy.interpolate import interp1d 

__all__ = ['SincInterp']


def trim_sinc(kx):
    #equivalent to setting window to boxcar
    return sp.sinc(kx)

def lanczos(kx):
    w = 5.
    out = sp.sinc(kx/w)
    m = abs(kx) > w
    out[m] = 0
    return out*sp.sinc(kx)

def cubic_sinc_approx(kx):
    #see Park and Schowengerdt 1982
    alpha = -1.0
    out = sp.zeros(kx.size)

    m = abs(kx) < 1
    out[m] = (alpha + 2)*abs(kx[m])**3 - (alpha + 3)*abs(kx[m]**2) + 1
    m = (abs(kx) >= 1)*(abs(kx) < 2)
    out[m] = alpha*abs(kx[m])**3 - 5*alpha*abs(kx[m])**2 + 8*alpha*abs(kx[m]) - 4*alpha
    return out

func_dic = {"sinc":trim_sinc,
            "lanczos":lanczos,
            "cubic_conv":cubic_sinc_approx}



class SincInterp(object):
    """
    This class will construct an sinc kernel, and use convolvution
    methods to interpolate x,y pairs at arbitrary x.  The input x,y
    must be equally spaced, but the class will sort the data
    internally.  The interpolator does not handle extrapolation, but
    the the evaluted points xnew may be arbitrarily spaced, and in any
    order.

    sinc interpolation is the optimal choice, since this represents an
    ideal band-pass limited filter, which will preserve the power
    spectrum of the original data.  In practice, there is ringing
    because the sinc function must be evaluated on a finite interval
    (the exception is for periodic signals with no aliasing--- these
    are truly band-limited).  This can be alleviated with the choice
    of window.

    A good choice is the Lanczos window, which is not supported in
    scipy.signal.  This option is implemented internally, and set as a
    default with a 3 pixel half-width.  One can also choose to use a
    truncated sinc (seems not so bad if one chooses a wide window), a
    piecewise approximation (see Park & Schowengerdt 1982), or invoke
    sp.signal.get_window.

    scipy.sigal.resample also implements sinc interpolation, by
    increasing the sampling rate through zero-padding of the data's
    Fourier transform. This will result in ringing, but
    signal.resample can also include a window function to mitigate.
    The window is multiplied in Fourier space before zero-padding,
    which is different than the use of windows here.  signal.resample
    may be better if you want more samples, but not if you want to
    easily specify where to interpolate.  Some experiments also
    indicate that windowing the Fourier transform in this way results
    in larger errors than windowing the sinc function in real space.

    Only 1D for now.  Edge effects are handled by just replacing with
    linear interpolation.

    Parameters
    ----------
    x = independent variable
    y = dependent variable
    window = window for the sinc kernel
             (implemented in real space)

    Examples
    --------
    >>> si = SincInterp(x,y)
    >>> si.window ='lanczos'
    >>> yint = si(xnew)
    >>> si.window = 'boxcar' #from sp.signal.get_window
    >>> yin2 = si(xnew)
    """
    def __init__(self,x,y,window='lanczos',kw = 15):
        i = sp.argsort(x)
        self.x = x[i]
        self.y = y[i]
        #enforce equal spacing
        d = sp.diff(x)
        if (sp.absolute(d - d.max()) > 1.e-10 ).any():
            raise ValueError("abcissas are not evenly spaced")

        self.window=window
        #kernel half-width
        self.kw = kw

    def __call__(self,xnew):
        xnew.sort()
        #no extraplolation
        if (xnew < self.x[0]).any() or (xnew > self.x[-1]).any():
            raise ValueError("new abcissas are outside of original domain")

        self.out = sp.zeros(xnew.size)
        for ex in xnew:
            k = self._get_kernel(ex,self.window)
            yi = sp.convolve(self.y,k,mode='same')
            self._tidy(ex,xnew,yi)

        return self.out

    def _get_kernel(self,x1,func_name):
        #evaluate pixel shift
        if func_name in func_dic.keys():
            func = func_dic[func_name]
        else:
            func = lambda kx:  sp.sinc(kx)*get_window(func_name,kx.size)
        
        i = sp.searchsorted(self.x,x1)
        dpix = (x1 - self.x[i - 1])/(self.x[i] - self.x[i - 1])
        if dpix == 1.0: dpix = 0
        assert sp.absolute(dpix) < 1.0
        kx = sp.r_[-self.kw:self.kw + 1] + dpix
        k = func(kx)
        k = k/sp.sum(k)
        return k

    def _tidy(self,x1,xnew,yi):
        i = sp.searchsorted(self.x,x1) 
        #if x1 is in self.x, search.sorted returns the correct index.
        #Other wise, it matches the index for the next greater self.x,
        #but we want the index for the next lesser.
        if self.x[i] != x1:
            i -= 1
        assert i.size == 1

        #just do linear interpolation on the section spoiled by the
        #convolution
        if self.window =='lanczos':
            edge = 11.
        elif self.window =='cubic_conv':
            edge = 5.
        else:
            edge = 2*self.kw + 1
        if i < edge :
            z = interp1d(self.x,self.y)
            yi[i] = z(x1)

        elif  i > yi.size - edge:
            z = interp1d(self.x,self.y)
            yi[i] = z(x1)

        i2 = sp.where(xnew == x1)[0]
        assert i2.size == 1
        self.out[i2] = yi[i]

