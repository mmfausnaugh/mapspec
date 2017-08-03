import scipy as sp
from scipy import linalg
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')

from mapspec.spectrum import *
from mapspec.mapspec  import *


sref = TextSpec('ref.smooth.txt',style='linear')

speclist = sys.argv[2:-2]
date     = sys.argv[-2]
ofile    = sys.argv[-1]

print(date)

outputgrid = sp.genfromtxt('outputgrid.dat')

if sys.argv[1] == 'gauss':
    s0 = sp.genfromtxt('mapspec.params',usecols=(1))
    fnames = sp.genfromtxt('mapspec.params',usecols=(0),dtype=str)
    chi1,chi2 = sp.genfromtxt('mapspec.params',unpack=1,usecols=(2,6))
    mask = chi1 < chi2

    shift, scale, width = sp.genfromtxt('mapspec.params',unpack=1,usecols=(7,8,9))
    dshift, dscale = sp.genfromtxt('mapspec.params',unpack=1,usecols=(3,4))
    shift[mask] = dshift[mask]
    scale[mask] = dscale[mask]
    width[mask] = 0.0001

    rescale = RescaleModel(sref,kernel = 'Gauss')


if sys.argv[1] == 'gh':
    s0 = sp.genfromtxt('mapspec.params',usecols=(1))
    fnames = sp.genfromtxt('mapspec.params',usecols=(0),dtype=str)
    chi1,chi2 = sp.genfromtxt('mapspec.params',unpack=1,usecols=(2,11))
    mask = chi1 < chi2

    shift, scale, width,h3,h4 = sp.genfromtxt('mapspec.params',unpack=1,usecols=(12,13,14,15,16))
    dshift, dscale = sp.genfromtxt('mapspec.params',unpack=1,usecols=(3,4))
    shift[mask] = dshift[mask]
    scale[mask] = dscale[mask]
    width[mask] = 0.0001
    h3[mask] = 0.0
    h4[mask] = 0.0

    rescale = RescaleModel(sref,kernel = 'Hermite')

yaverage = 0
weights  = 0
Covar = 0
#read in everything that is needed
inSpecs,outSpecs,Masks,Chains = [],[],[],[]
for i in range(len(speclist)):

    j = [ k for k in range(fnames.size) if fnames[k] in speclist[i] ]
    print(j,  fnames[j],len(j))
    assert len(j) == 1
    fuse = fnames[j][0]

    if sys.argv[1] == 'gauss':
        p = {'shift':shift[j], 'scale':scale[j], 'width' :width[j]}
    if sys.argv[1] == 'gh':
        p = {'shift':shift[j], 'scale':scale[j], 'width' :width[j],'h3':h3[j],'h4':h4[j]}
    rescale.p = p

    s = TextSpec(fuse,style='linear')
    s.wv -= s0[j]
    inSpecs.append(s)

    sout,lost_pix= rescale.output(s,getcovar=False)
    m = (sout.wv >= outputgrid[0])*(sout.wv <= outputgrid[1])
    sout.wv = sout.wv[m]
    sout.f  = sout.f[m]
    sout.ef  = sout.ef[m]
    outSpecs.append(sout)
    
    c = Chain()
    if sys.argv[1] == 'gauss':
        c.read('chains/'+fuse+'.chain.gauss')
    if sys.argv[1] == 'gh':
        c.read('chains/'+fuse+'.chain.herm')
        
    c.burn(0.5)
    Chains.append(c)

    if sys.argv[1] == 'gauss':
        covar = sp.genfromtxt('covar_matrices/covar_'+fuse)
    if sys.argv[1] == 'gh':
        covar = sp.genfromtxt('covar_matrices/covar.h._'+ fuse)
     
    covar = covar[m,:]
    covar = covar[:,m]


    #get the average, measurement errors
    yaverage += sout.f/sout.ef**2
    weights  += 1./sout.ef**2
    Covar    += linalg.inv(covar)

yaverage /= weights
measure_error = sp.sqrt(1./weights)

#bootstrap the average to get model errors
if sys.argv[1] == 'gauss':
    pboot = {'shift':0.0, 'scale':1.0, 'width' : 1.0}
if sys.argv[1] == 'gh':
    pboot = {'shift':0.0, 'scale':1.0, 'width' :1.0,'h3':0.0,'h4':0.0}

averages = []
for i in range(1000):
    if i %100 == 0: print('iteration:  ',i)
    
    yboot      = 0
    weightboot = 0

    for j in range(len(inSpecs)):
        suse = inSpecs[j]
        cuse = Chains[j]


        params = sp.transpose(cuse.pchain)
        indexrandom = (sp.rand(params.shape[0])*params.shape[1]).astype(int)
        for k, key in enumerate(pboot.keys()):
            pboot[key] = params[ cuse.index[key],[indexrandom[k]] ]
#        pboot['shift'] = 
        rescale.p = pboot

        
        soutboot, lost_pix = rescale.output(suse,getcovar=False)
        m = (soutboot.wv >= outputgrid[0])*(soutboot.wv <= outputgrid[1])
        soutboot.wv = soutboot.wv[m]
        soutboot.f  = soutboot.f[m]
        soutboot.ef = soutboot.ef[m]

        yboot += soutboot.f/soutboot.ef**2
        weightboot += 1./soutboot.ef**2                
    averages.append(yboot/weightboot)

l,m,u = sp.percentile(averages,[16,50,84],axis=0)
model_error = (u - l)/2.
error = sp.sqrt(model_error**2 + measure_error**2)

Covar = linalg.inv(Covar)
#plt.plot(sp.sqrt(Covar[sp.diag_indices_from(Covar)]/model_error**2))
#plt.show()
Covar[sp.diag_indices_from(Covar)] += model_error**2

sp.savetxt(ofile, sp.c_[sout.wv, yaverage, error ],fmt='% 6.2f % 4.4e % 4.4e')
sp.savetxt('covar_matrices/'+ofile+'covar',Covar)

if (Covar.diagonal() <0).any():
    f = open('covarlog','a')
    f.write(date+'\n')
    f.close()

#plt.imshow(averages,interpolation='nearest')
#print sp.shape(averages)
#print sp.c_[model_error,measure_error,error]
#plt.figure()
#plt.plot(sout.wv,model_error/yaverage,'c')
#plt.plot(sout.wv,measure_error/yaverage,'r')
#plt.plot(sout.wv,error/yaverage,'k')
#
#F,(ax1,ax2) = plt.subplots(2,1,sharex='col')
#
#ax1.plot(sout.wv,yaverage,'k')
#for i in range(3):
#    ax1.plot(inSpecs[i].wv,inSpecs[i].f,'r')
#    ax1.plot(outSpecs[i].wv,outSpecs[i].f,'b')
##plt.figure()
#ax2.plot(sout.wv,m-yaverage,'b')
#ax2.fill_between(sout.wv, l-yaverage,u-yaverage,facecolor='b',alpha=0.25)
#ax2.set_ylim([-0.2e-15,0.2e-15])
#plt.show()


#to do, covariances......
