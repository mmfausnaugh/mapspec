import scipy as sp
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')


from mapspec.spectrum import *
from mapspec.mapspec import *

from copy import deepcopy


outgrid = sp.genfromtxt('outputgrid.dat')
def sub(s1,s2,params):
#    s1.wv += params['shift']
    print(params)
    f1  = s1.f/params['scale']
    ef1 = s1.ef/params['scale']

    m1 = (s1.wv > outgrid[0])*(s1.wv < outgrid[1])
    m2 = (s2.wv > outgrid[0])*(s2.wv < outgrid[1])

    f = f1[m1] - s2.f[m2]
    ef = sp.sqrt(ef1[m1]**2 + s2.ef[m2]**2)
    return s1.wv[m1],f, ef



sref   = TextSpec(sys.argv[1])
window = sp.genfromtxt(sys.argv[2])

lref   = EmissionLine(sref,window[0],[ window[1],window[2] ] )

istyle = sys.argv[3]
lref.set_interp(style=istyle)

speclist = sp.genfromtxt(sys.argv[4],dtype=str,usecols=(0))
fout = open('mapspec.diff.params','a')

for spec in speclist:
    print(spec)
    s = TextSpec(spec)
    s.set_interp(style=istyle)

    s0 = get_cc(sref.f,s.f,sref.wv,s.wv)
    s.wv -= s0[0]

    l = EmissionLine(s,window[0],[ window[1],window[2] ])
    l.set_interp(style=istyle)


    f   = RescaleModel(l,kernel="Delta")
    try: 
        chi2_delta,p_delta,frac_delta = metro_hast(1000,lref,f,keep=0)
    except:
        chi2_delta,p_delta,frac_delta = 999,{'shift':-99, 'scale':-99}, 0 
        

    f    = RescaleModel(l,kernel="Gauss")
 
    try: 
        chi2_gauss,p_gauss,frac_gauss,chain_gauss = metro_hast(5000,lref,f,keep=1,plot=0)
    except:
        chi2_gauss,p_gauss,frac_gauss = 999,{'shift':-99, 'scale':-99, 'width':-99}, 0 
        
    if chi2_delta < chi2_gauss:
        print('prefers delta')
        f.p = {'shift':p_delta['shift'], 'scale':1.0, 'width': 0.001 }
        p_gauss['scale'] = p_delta['scale']
    else:
        f.p = {'shift':p_gauss['shift'], 'scale':1.0, 'width': p_gauss['width'] }

    sout,lost_pix = f.output(sref)
    xout, fdiff, efdiff = sub(s,sout,p_gauss)
    sp.savetxt('diff_'+spec,sp.c_[xout,fdiff,efdiff],fmt='% 6.2f % 4.4e % 4.4e')

    
    f    = RescaleModel(l,kernel="Hermite")
#    f.make_dist_prior(chain_gauss,'width')



    try:
        chi2_herm,p_herm,frac_herm,chain_herm = metro_hast(50000,lref,f,keep=1,plot=0)
    except:
        chi2_herm,p_herm,frac_herm = 999, {'shift':99,'scale':-99,'width':-99,'h3':-99,'h4':-99}, 0 

    if chi2_delta < chi2_herm:
        print('prefers delta')
        f.p = {'shift':p_delta['shift'], 'scale':1.0, 'width': 0.001 , 'h3':0.0, 'h4':0.0}
        p_herm['scale'] = p_delta['scale']
    else:
        f.p = {'shift':p_herm['shift'], 'scale':1.0, 'width': p_herm['width'] , 'h3':p_herm['h3'], 'h4':p_herm['h4']}

    sout,lost_pix = f.output(sref)
    xout, fdiff, efdiff = sub(s,sout,p_herm)
    sp.savetxt('diff.h._'+spec,sp.c_[xout,fdiff,efdiff],fmt='% 6.2f % 4.4e % 4.4e')

    fout.write(
        "%15s % 8.4f  %10.2f % 8.4f % 8.4f % 5.2f %10.2f % 8.4f % 8.4f % 8.4f % 5.2f % 10.2f % 8.4f % 8.4f % 8.4f % 5.4e % 5.4e %8.4f\n"%
        (spec,s0[0],
         chi2_delta,p_delta['shift'],p_delta['scale'], frac_delta,
         chi2_gauss,p_gauss['shift'],p_gauss['scale'],p_gauss['width'],frac_gauss,
         chi2_herm,p_herm['shift'],p_herm['scale'],p_herm['width'],p_herm['h3'],p_herm['h4'],frac_herm)
        )
    fout.flush()
    chain_gauss.save(spec+'.chain.dgauss')
    try:
        chain_herm.save(spec+'.chain.dherm')
    except:
        continue
    plt.close('all')
fout.close()
