import matplotlib
matplotlib.use('TkAgg') 

import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from spectrum import *
from mapspec import *
from copy import deepcopy

import sys


sref   = TextSpec(sys.argv[1])
window = sp.genfromtxt(sys.argv[2])

lref   = EmissionLine(sref,window[0],[ window[1],window[2] ] )

istyle = sys.argv[3]
lref.set_interp(style=istyle)

speclist = sp.genfromtxt(sys.argv[4],dtype=str)
fout = open('mapspec.params','a')

for spec in speclist:
    print spec
    s = TextSpec(spec)
    s.set_interp(style=istyle)

    s0 = get_cc(sref.f,s.f,sref.wv,s.wv)
    s.wv -= s0[0]

    l = EmissionLine(s,window[0],[ window[1],window[2] ])
    l.set_interp(style=istyle)


    f   = RescaleModel(lref,kernel="Delta")
    try: 
        chi2_delta,p_delta,frac_delta = metro_hast(1000,l,f,keep=0)
    except:
        chi2_delta,p_delta,frac_delta = 999,{'shift':-99, 'scale':-99}, 0 
        

    f    = RescaleModel(lref,kernel="Gauss")
 
    try: 
        chi2_gauss,p_gauss,frac_gauss,chain_gauss = metro_hast(5000,l,f,keep=1)
    except:
        chi2_gauss,p_gauss,frac_gauss = 999,{'shift':-99, 'scale':-99, 'width':-99}, 0 
        
    if chi2_delta < chi2_gauss:
        f.p = {'shift':p_delta['shift'], 'scale':p_delta['scale'], 'width': 0.001 }

    sout,dummy = f.output(s)
    sp.savetxt('scale_'+spec,sp.c_[sout.wv,sout.f,sout.ef],fmt='% 6.2f % 4.4e % 4.4e')

    f    = RescaleModel(lref,kernel="Hermite")
    f.make_dist_prior(chain_gauss,'width')


    try:
        chi2_herm,p_herm,frac_herm,chain_herm = metro_hast(50000,l,f,keep=1)
    except:
        chi2_herm,p_herm,frac_herm = 999, {'shift':99,'scale':-99,'width':-99,'h3':-99,'h4':-99}, 0 

    if chi2_delta < chi2_herm:
        print 'Alt!'
        f.p = {'shift':p_delta['shift'], 'scale':p_delta['scale'], 'width': 0.001 , 'h3':0.0, 'h4':0.0}

    sout,dummy = f.output(s)
    sp.savetxt('scale.h._'+spec,sp.c_[sout.wv,sout.f,sout.ef],fmt='% 6.2f % 4.4e % 4.4e')

    fout.write(
        "%15s %10.2f % 8.4f % 8.4f % 5.2f %10.2f % 8.4f % 8.4f % 8.4f % 5.2f % 10.2f % 8.4f % 8.4f % 8.4f % 5.4e % 5.4e %8.4f\n"%
        (spec,
         chi2_delta,p_delta['shift'],p_delta['scale'], frac_delta,
         chi2_gauss,p_gauss['shift'],p_gauss['scale'],p_gauss['width'],frac_gauss,
         chi2_herm,p_herm['shift'],p_herm['scale'],p_herm['width'],p_herm['h3'],p_herm['h4'],frac_herm)
        )
    fout.flush()
#    chain_gauss.save(spec+'.chain.gauss')
#    chain_herm.save(spec+'.chain.herm')
        
    plt.close('all')
fout.close()
