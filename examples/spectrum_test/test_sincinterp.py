import scipy as sp
from scipy.signal import resample
import matplotlib.pyplot as plt
import sys
#sys.path.insert(0,'..')
#from sinc_interp import *
#from mapspec.sinc_interp2 import SincInterp
from mapspec.mapspec.sinc_interp import SincInterp

x = sp.r_[-1:0.85:51j]
dt = 0.5*(x[1] - x[0])
#try different functions
#y = sp.sin(2*sp.pi*x)
y = x**2
#scipy sinc resampling
y2,x2 = resample(y,500,t=x)#,window='han')

xa = x[0:-2] + dt
xb = sp.r_[-1:0.84:501j]

L = SincInterp(x,y)
L.window = 'lanczos'
ya = L(xa)
yb = L(xb)

plt.plot(x,y,'ko',ms=10,label='input data') #original
plt.plot(x2,y2,'r.',label='scipy.signal.resample') #original
plt.plot(xb,yb,'c.',label='mapspec.SincInterp at scipy spacing'.format(
    xb[1] - xb[0])) #sinc interp wavelengths 2
plt.plot(xa,ya,'m.',label='mapspec.SincInterp at 0.5 input') #sinc interp wavelengths 1

plt.gca().set_title('Comparison of mapspec.SincInterp and \n'
                    'scipy sinc interpolation on a parabola.  \nscipy '
                    'interpolates by filling in Fourier Domain, \n'
                    'and so has edge effects and subtle wiggles.',
                    fontsize=12)

plt.legend()
plt.figure()

#now try on some data, this is the 
x,y,z = sp.genfromtxt('test_data/test.LBT_MODS1.dat',unpack=1)
m = (x>4800)*(x<5200)
x = x[m]
y = y[m]
z = z[m]

L2 = SincInterp(x,y)
L2.window = 'lanczos'
x2 = sp.r_[x.min():x.max():1000j]
y2 = L2(x2)

x3 = sp.r_[x.min():x.max()-1:100j]
y3 = L2(x3)

a = sp.r_[4900:4950:200j]
b = sp.r_[5000:5100:3]
x4 = sp.r_[a,b]
y4 = L2(x4)

yi = L2(x)
#print(x2.shape,x3.shape,x.shape)


plt.plot(x,y,'b.-',label='orignal')  #original data
plt.plot(x2,y2,'r.',label='mapspec.SincInterp') #sinc interp
plt.plot(x3,y3,'c.',label='mapspec at 1/10 spacing') #sinc interp sparse
#plt.plot(x,yi,'ko')  #identity, i.e. original
plt.plot(x4,y4,'g.',label='mapspec with gap in array')  # gap

plt.legend(loc='upper left')

plt.gca().set_title('Comparison of some mapspec sinc interpolations on '
                    'NGC5548 Hbeta + [OIII] complex from LBT MODS1 \n'
                    'taken in 2014.  The sinc interpolation does not  '
                    'depend on spacing, and mapspec can handle gaps in \n'
                    'the requested wavlengths to interpolate.', 
                    fontsize=12)


plt.show()
