import scipy as sp
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from mapspec.spectrum import *
import sys


"""This is for looking at the distribution of line widths, calculated
from the line dispersion.  The user chooses a cutoff, and the program will
calculate a smoothing width to try get the reference to match the worst spectrum.
"""

plt.ion()

speclist = sp.genfromtxt(sys.argv[1],dtype=str)
window   = sp.genfromtxt(sys.argv[2])


dispdist = []
for spec in speclist:
    s = TextSpec(spec)
    l = EmissionLine(s,window[0],[window[1],window[2]])

    dispdist.append(l.dispersion())


dispdist = sp.array(dispdist)


s = TextSpec(sys.argv[3])
l = EmissionLine(s,window[0],[window[1],window[2]])
res = l.dispersion()

plt.hist(dispdist)

l,u = plt.gca().get_ylim()
plt.plot([res,res],[l,u],'r--')

plt.draw()
cut = float(raw_input('Where to cut disp?'))

m = dispdist < cut
newres = sp.sqrt(dispdist[m].max()**2 - res**2 )

s.smooth(13,name = ('gaussian', newres/(s.wv[1] - s.wv[0]) ) )

l = EmissionLine(s,window[0],[window[1],window[2]])
res2 = l.dispersion()

plt.cla()
plt.hist(dispdist)
l,u = plt.gca().get_ylim()
plt.plot([res,res],[l,u],'r--')
plt.plot([res2,res2],[l,u],'g--')
plt.show()

sp.savetxt(sys.argv[3] + '.smooth.txt',sp.c_[sp.real(s.wv),sp.real(s.f),sp.real(s.ef)])

fout = open('ref_resolution.dat','w')
fout.write('Ref native resolution:    % 2.4f\n'%res)
fout.write('Cut for worst resolution: % 2.4f\n'%cut)
fout.write('Worst object  below cut:  % 2.4f\n'%dispdist[m].max())
fout.write('Kernel width:             % 2.4f\n'%newres)
fout.write('New ref resolution:       % 2.4f\n'%res2)
fout.close()


