import matplotlib
matplotlib.use('TkAgg') 

import scipy as sp
import matplotlib.pyplot as plt
from mapspec.spectrum import *
from mapspec.mapspec import *
from copy import deepcopy

import sys

sref   = TextSpec(sys.argv[1])
window = sp.genfromtxt(sys.argv[2])

lref   = EmissionLine(sref,window[0],[ window[1],window[2] ] )

istyle = sys.argv[3]
lref.set_interp(style=istyle)

speclist = sp.genfromtxt(sys.argv[4],dtype=str)
fout = open('mapspec.params2','a')

for spec in speclist:
    print spec
    s = TextSpec(spec)
    s.set_interp(style=istyle)
    l = EmissionLine(s,window[0],[ window[1],window[2] ])
    l.set_interp(style=istyle)

    s0 = get_cc(sref.f,s.f,sref.wv,s.wv)
    f    = Model(lref,kernel="Gauss")
    f.p[0] = s0[0]

    f2    = Model(l,kernel="Gauss")
    f2.p[0] = -s0[0]

    chi2_gauss,p_gauss,frac_gauss,chain_gauss = metro_hast(5000,l,f,keep=1)
    chi2_ref,p_ref,frac_ref = metro_hast(5000,lref,f2,keep = 0)
    
    if chi2_ref < chi2_gauss:
        f.p = [-p_ref[0],1./p_ref[1],0.001]
        xout,yout,zout,dummy  = f.output(s)
        p_gauss = [ -p_ref[0],1./p_ref[1],-p_ref[2] ]
        chi2_gauss = chi2_ref
    else:
        f.p = p_gauss
        xout,yout,zout,dummy = f.output(s)

    sp.savetxt('scale_'+spec,sp.c_[xout,yout,zout],fmt='% 6.2f % 4.4e % 4.4e')

    f    = Model(lref,kernel="Hermite")
    f.p[0] = s0[0]

    f2      = Model(l,kernel="Hermite")
    f2.p[0] = -s0[0]

    chi2_herm,p_herm,frac_herm,chain_herm = metro_hast(10000,l,f,keep=1)
    chi2_ref,p_ref,frac_ref = metro_hast(5000,lref,f2,keep = 0)

    if chi2_ref < chi2_herm:
        f.p = [-p_ref[0],1./p_ref[1],0.001,0.0,0.0]
        xout,yout,zout,dummy  = f.output(s)
        p_herm = [ -p_ref[0],1./p_ref[1],-p_ref[2],p_ref[3],p_ref[4] ]
        chi2_herm = chi2_ref
    else:
        f.p = p_herm
        xout,yout,zout,dummy = f.output(s)

    try:
        sp.savetxt('scale.h._'+spec,sp.c_[xout,yout,zout],fmt='% 6.2f % 4.4e % 4.4e')
    except:
        continue

#    sp.savetxt(spec+'.chain.gauss',chain_gauss.pchain)
#    sp.savetxt(spec+'.chain.herm',chain_herm.pchain)

    
    fout.write(
        "%15s %10.2f % 8.4f % 8.4f % 8.4f % 5.2f % 10.2f % 8.4f % 8.4f % 8.4f % 5.4e % 5.4e %8.4f\n"%
        (spec,chi2_gauss,p_gauss[0],p_gauss[1],p_gauss[2],frac_gauss,
         chi2_herm,p_herm[0],p_herm[1],p_herm[2],p_herm[3],p_herm[4],frac_herm)
        )
    fout.flush()
    
    plt.close('all')
fout.close()
