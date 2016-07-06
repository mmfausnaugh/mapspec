import scipy as sp
from scipy.signal import resample
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
#from sinc_interp import *
#from mapspec.sinc_interp2 import SincInterp
from mapspec.sinc_interp import SincInterp

x = sp.r_[-1:0.85:51j]
dt = 0.5*(x[1] - x[0])
#try different functions
#y = sp.sin(2*sp.pi*x)
y = x**2
#scipy sinc resampling
y2,x2 = resample(y,100,t=x)#,window='han')

xa = x[0:-2] + dt
xb = sp.r_[-1:0.84:501j]

L = SincInterp(x,y)
L.window = 'lanczos'
ya = L(xa)
yb = L(xb)

plt.plot(x2,y2,'ko') #original
plt.plot(xa,ya,'ro') #sinc interp wavelengths 1
plt.plot(xb,yb,'co') #sinc interp wavelengths 2


plt.figure()

#now try on some data, this is the 
x,y,z = sp.genfromtxt('test.dat',unpack=1)
m = (x>4800)*(x<5200)
x = x[m]
y = y[m]
z = z[m]

L2 = SincInterp(x,y)
L2.window = 'lanczos'
x2 = sp.r_[x.min():x.max():1000j]
y2 = L2(x2)

x3 = sp.r_[x.min():x.max()-1:100j]
print 'new grid'
y3 = L2(x3)

a = sp.r_[4900:4950:200j]
b = sp.r_[5000:5100:3]
x4 = sp.r_[a,b]
y4 = L2(x4)

yi = L2(x)
print x2.shape,x3.shape,x.shape


plt.plot(x,y,'b',label='orignal')  #original data
plt.plot(x2,y2,'ro',label='sinc interp') #sinc interp
plt.plot(x3,y3,'co',label='sparse interp') #sinc interp sparse
#plt.plot(x,yi,'ko')  #identity, i.e. original
plt.plot(x4,y4,'go',label='has a gap')  # gap



plt.show()
