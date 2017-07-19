import scipy as sp
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')


from mapspec.spectrum import *
from mapspec.mapspec import *

from copy import deepcopy



sref   = TextSpec(sys.argv[1])
window = sp.genfromtxt(sys.argv[2])

lref   = EmissionLine(sref,window[0],[ window[1],window[2] ] )

istyle = sys.argv[3]
lref.set_interp(style=istyle)

speclist = sp.genfromtxt(sys.argv[4],dtype=str)
fout = open('mapspec.params','a')

for spec in speclist:
    print(spec)
    s = TextSpec(spec)
    s.set_interp(style=istyle)
    l = EmissionLine(s,window[0],[ window[1],window[2] ])
    l.set_interp(style=istyle)

    s0 = get_cc(sref.f,s.f,sref.wv,s.wv)        

    f    = RescaleModel(lref,kernel="Gauss")
    f.p['shift'] = s0[0]
    f2    = RescaleModel(l,kernel="Gauss")
    f2.p['shift'] = s0[0]
 
    try: 
        chi2_gauss,p_gauss,frac_gauss,chain_gauss = metro_hast(5000,l,f,keep=1)
    except:
        chi2_gauss,p_gauss,frac_gauss = 999,{'shift':-99, 'scale':-99, 'width':-99}, 0 
    try: 
        chi2_ref,p_ref,frac_ref = metro_hast(5000,lref,f2,keep=0)
    except:
        chi2_ref,p_ref,frac_ref = 999,{'shift':-99, 'scale':-99, 'width':-99}, 0 
        
    if chi2_ref < chi2_gauss:
        f.p = {'shift':-p_ref['shift'], 'scale':1./p_ref['scale'], 'width': 0.001 }
        p_gauss = {'shift':-p_ref['shift'], 'scale':1./p_ref['scale'], 'width':-p_ref['width']}
        chi2_gauss = chi2_ref
    sout,dummy = f.output(s)
    sp.savetxt('scale_'+spec,sp.c_[sout.wv,sout.f,sout.ef],fmt='% 6.2f % 4.4e % 4.4e')

    f    = RescaleModel(lref,kernel="Hermite")
    f.p['shift'] = s0[0]
    if chi2_gauss < chi2_ref:
        f.make_dist_prior(chain_gauss,'width')
    f2    = RescaleModel(l,kernel="Hermite")
    f2.p['shift'] = s0[0]


    try:
        chi2_herm,p_herm,frac_herm,chain_herm = metro_hast(10000,l,f,keep=1)
    except:
        chi2_herm,p_herm,frac_herm = 999, {'shift':99,'scale':-99,'width':-99,'h3':-99,'h4':-99}, 0 
    try:
        chi2_ref,p_ref,frac_ref,chain_ref = metro_hast(10000,lref,f2,keep=0)
    except:
        chi2_ref,p_ref,frac_ref = 999, {'shift':99,'scale':-99,'width':-99,'h3':-99,'h4':-99}, 0 

    if chi2_ref < chi2_herm:
        f.p = {'shift':-p_ref['shift'], 'scale':1./p_ref['scale'], 'width': 0.001 , 'h3':0.0, 'h4':0.0}
        p_herm = {'shift':-p_ref['shift'], 'scale':1./p_ref['scale'], 'width': -p_ref['width'] , 'h3':0.0, 'h4':0.0}
        chi2_herm = chi2_ref

    sout,dummy = f.output(s)
    sp.savetxt('scale.h._'+spec,sp.c_[sout.wv,sout.f,sout.ef],fmt='% 6.2f % 4.4e % 4.4e')

    fout.write(
        "%15s %10.2f % 8.4f % 8.4f % 8.4f % 5.2f % 10.2f % 8.4f % 8.4f % 8.4f % 5.4e % 5.4e %8.4f\n"%
        (spec,
#         chi2_delta,p_delta['shift'],p_delta['scale'], frac_delta,
         chi2_gauss,p_gauss['shift'],p_gauss['scale'],p_gauss['width'],frac_gauss,
         chi2_herm,p_herm['shift'],p_herm['scale'],p_herm['width'],p_herm['h3'],p_herm['h4'],frac_herm)
        )
    fout.flush()
    chain_gauss.save(spec+'.chain.gauss')
    chain_herm.save(spec+'.chain.herm')
        
    plt.close('all')
fout.close()
