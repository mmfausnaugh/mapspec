import matplotlib
matplotlib.use('TkAgg') 

import scipy as sp
from scipy.integrate import simps
import matplotlib.pyplot as plt
from spectrum import TextSpec,EmissionLine
from copy import deepcopy
#these are 'probalists', turns out to matter in order to match van der
#Marel & Franx 1993
from numpy.polynomial.hermite_e import HermiteE as H

__all__ = ["Model","Chain","get_cc","metro_hast","gibbs"]

debug = True

class Model(object):
    """
    This model stores parameters (shift, scale, and convolution
    kernel), performs the operation, and evaluates the likelihood.
    """

    def __init__(self,Lref,kernel='Hermite'):

        #this is an emission line object, to which we scale to 
        self.Lref = Lref  
        self.count = 0

        #empty list to add prior dist, if wanted
        self.prior_indices = []
        self.bins    = []
        self.prior_prob    = []

        if kernel == 'Delta':
            self._get_kernel = lambda x: [1.0]
            self.p    = [0.0,1.0]
            self.step_scale = [0.03,0.01] 

        elif kernel == 'Gauss':
            self._get_kernel = self._Gauss

            #params are shift, scale, and width of gaussian
            #convolution
            self.p    = [0.0,2.5,0.01]
            self.step_scale = [0.05,0.02,0.750] 


        if kernel == 'Hermite':
            self._get_kernel = self._Hermite

            #params are shift, scale, width, h3, and h4 of
            #gauss-hermite polynomial

            self.p    = [0.0,1.0,0.01,0.0,0.0]
            self.step_scale = [0.05,0.02,0.30,0.03,0.03] 


    def __call__(self,L):
        y,var,mask  = self._get_yz(L)
        lnlikely  = self._get_chi2(y,var,mask)
        #add in priors.
        if len(self.prior_indices) > 0:
            lnlikely += self._add_priors()
#        lnlikely += sp.log10(1./self.p[1])
        return lnlikely

    def _add_priors(self):
        prior = 0
        for k,i in enumerate(self.prior_indices):
            prior += -2.*sp.log( self.prior_prob[k](self.p[i]) )
        return prior


    def make_prior_from_dist(self,C,i,burn = 0.5):
        params = sp.transpose(C.pchain)
        self.prior_indices.append(i)
        prior_dist = params[i]
        icut = prior_dist.size*burn
        prior_dist = prior_dist[icut::]
        c1,m,c2 = sp.percentile(prior_dist,[16,50,84])
        print c1,m,c2

        self.prior_prob.append(
            lambda x: sp.exp(-0.5*(x - m)**2/ ( (c2-c1)/2 )**2 )
            )
        self.p[i] = m

    def make_prior_from_param(self,i,p):
        #assume gaussian, with width of n 1% of fitted value
        self.prior_indices.append(i)
        sig = 0.05*p
        self.prior_prob.append(
            lambda x: sp.exp(-0.5*(x - p)**2/sig**2)
            )


    def output(self,S):
        s = deepcopy(S)
        #convolve
        k = self._get_kernel(s)
        s.ef = sp.convolve(s.ef**2,k,mode='same')
        s.f = sp.convolve(s.f,k,mode='same')
        #shift
        s.wv -= self.p[0]
        #scale
        s.ef *= self.p[1]
        s.f  *= self.p[1]

        m = (S.wv >= s.wv.min()  )*(S.wv <= s.wv.max()  )
        y,z = s.interp(S.wv[m])

        #edge effects?????********

        return S.wv[m],y,z,m


    def _get_yz(self,L):
        l = deepcopy(L)
         #convolve
        k = self._get_kernel(l)
        l.ef = sp.sqrt(sp.convolve(l.ef**2,k,mode='same'))
        l.f = sp.convolve(l.f,k,mode='same')
        #shift
        l.wv -= self.p[0]
        #scale
        l.ef *= self.p[1]
        l.f  *= self.p[1]


        #trim 10% of data to help with edge effects and shifting the data
        trim = round(0.05*self.Lref.wv.size)
        m = (self.Lref.wv >= l.wv[trim] )*(self.Lref.wv <= l.wv[-trim] )

        y,z = l.interp(self.Lref.wv[m])

        if (z**2 < 0).any(): 
            raise ValueError('Weights became negative')

        return y,z**2,m

    def _Gauss(self,line):
        dlambda = line.wv[1] - line.wv[0]
        pixwidth = self.p[2]/dlambda
        prange = sp.r_[ -(line.wv.size //2) + 1 : (line.wv.size)//2  ]
        assert prange.size %2 == 1

        k = sp.exp(-0.5* prange**2/pixwidth**2)
        k /= abs(sp.sum(k))
        return k

    def _Hermite(self,line):
        dlambda = line.wv[1] - line.wv[0]
        pixwidth = self.p[2]/dlambda
        prange = sp.r_[ -(line.wv.size //2) + 1 : (line.wv.size)//2  ]
        
        h = H([ 1.0, 0.0, 0.0,self.p[3], self.p[4] ])
        #although the constants will divide out when normalizing the
        #kernel, they are important for making sure that h3 and h4 are
        #defined correctly.  Technically, this only matters for
        #choosing good intervals of h3 and h4, otherwise it is just a
        #mismatch of units

        #Equations are defined in van der Marel & Franx 1993
        k = 1./pixwidth/sp.sqrt(2*sp.pi)*sp.exp(-0.5*prange**2/pixwidth**2)*h(prange/pixwidth)
        k /= abs(sp.sum(k))


        return k


    def _get_chi2(self,y,v,m):
        return sp.sum(    (self.Lref.f[m] - y)**2/(self.Lref.ef[m]**2 + v))


    def step(self):
        self.count = self.count + 1
        pout = self.step_scale*sp.randn(len(self.p)) + self.p


        if len(pout) > 2:
            dlambda = self.Lref.wv[1] - self.Lref.wv[0]
            #if width is too small, the kernel is
            #undersampled.  Weird things will happen, like flux that
            #goes negative
            if pout[2]/dlambda <0.5:
                pout[2] = 0.5*dlambda
            if pout[2] > 5.0:
                pout[2] = 5.0

        if len(pout) > 3:
            #experiments have found that h3 and h4 between -0.3 and
            #0.3 should be adequate (very diverse line shapes appear)
            if pout[3] < -0.3:
                pout[3] = -0.3
            elif pout[3] > 0.3:
                pout[3] = 0.3
            if pout[4] < -0.3:
                pout[4] = -0.3
            elif pout[4] > 0.3:
                pout[4] = 0.3

        return pout




class Chain(object):
    def __init__(self,M,chi2):
        self.pchain   = [deepcopy(M.p)]
        self.lnlikely = [deepcopy(chi2)]

        self.figure,(self.axes) = plt.subplots(sp.transpose(self.pchain).shape[0] + 1,1)

    def add(self,M,chi2):
        self.pchain.append(deepcopy(M.p))
        self.lnlikely.append(chi2)



    def plot(self,interact = 1):
        self.plotp = sp.transpose(self.pchain)
        for ax in self.axes: ax.cla()

        self.axes[0].plot(self.lnlikely)
        self.axes[0].set_ylabel('ln likelihood')

        for i in range(self.plotp.shape[0]):
            self.axes[i + 1].plot(self.plotp[i])

        self.axes[1].set_ylabel('shift')
        self.axes[2].set_ylabel('scale')

        for i in range(len(self.axes[3::])):
            self.axes[i+3].set_ylabel('kernel p%i'%i)
        for ax in self.axes[0:-1]:
            ax.set_xticks([])


        self.figure.subplots_adjust(hspace=0)
        if interact == 1:  
            plt.draw()



def get_cc(y1,y2,x1,x2):
    if x1.size > x2.size:
        m = (x1 > x2.min())*(x1 < x2.max())
        yuse = [y1[m],y2]
    else:
        m = (x2 > x1.min())*(x2 < x1.max())
        yuse = [y1,y2[m]]


    cc = sp.correlate(yuse[0],yuse[1],mode='same')
    i  = sp.where(cc == cc.max())[0]
    shift = (x1[1] - x1[0] )*(cc.size//2 - i)
    return shift



def metro_hast(ntrial,D,M,plot=0,keep=0):
    Mtry= deepcopy(M)
    chi2 = 1.e12
    chi2best = 1.e12

    pbest = deepcopy(M.p)
    accept = 0

    c = Chain(M,M(D))
    if plot ==1:
        plt.ion()
        c.plot()

    for i in range(ntrial):
        Mtry.p = M.step()
        chi2try = Mtry(D)

        if chi2try < chi2:
            
            M.p = deepcopy(Mtry.p)
            chi2 = deepcopy(chi2try)

            accept += 1
            c.add(M,chi2)
                
            if chi2 < chi2best:
                chi2best = deepcopy(chi2)
                pbest = deepcopy(M.p)

        else:
            prob = sp.exp(-chi2try/chi2)
            r = sp.rand()
            if r <= prob:
                M.p = deepcopy(Mtry.p)
                chi2 = deepcopy(chi2try)
                accept += 1
            c.add(M,chi2)
                
        if i%100 == 0 :
            print i,chi2best,chi2try
            if plot ==1:
                c.plot()

    if keep == 1:
        return chi2best,pbest,accept/float(ntrial),c
    else:
        return chi2best,pbest,accept/float(ntrial)


def gibbs(ntrial,D,M,keep=0,plot=0):
    #shift is help constant
    #scale ????
    #width 0 to some max
    #h3 -1 to 1
    #h4 -1 to 1


    c = Chain(M,M(D))

    chi2best = 1.e12
    pbest    = deepcopy(M.p)

    grid_dense1 = 500
    grid_dense2 = 500
    for j in range(ntrial):
 #       plt.figure()
        wgrid  = sp.r_[0.001:3:100j]
        h3grid = sp.r_[-1:1:1j*grid_dense1]
        h4grid = sp.r_[-1:1:1j*grid_dense2]

        pgrid = [wgrid,h3grid,h4grid]



        for i in range( sp.shape(pgrid)[0] ):
            z = []
            for p in pgrid[i]:
                M.p[i + 2] = p
                z.append(sp.exp(-0.5*M(D)) )
            k = simps(z,x = pgrid[i])
            M.p[i + 2] = sp.rand()*(pgrid[i][-1] - pgrid[i][0]) + pgrid[i][0]
            r = sp.rand()
        #I call this the "throwing darts" method, it's nice because it
        #doesn't depend on the grid (as the inverse transform does),
        #but it is not the most efficient and can be very slow if the
        #PDF is narrow
            count = 0
            while r > sp.exp(-0.5*M(D))/k:
                count += 1
                M.p[i + 2] = sp.rand()*(pgrid[i][-1] - pgrid[i][0]) + pgrid[i][0]
                r = sp.rand()
            if i == 1:  grid_dense1 = max(1500./len(sp.array(z)[z/k > 0.001]), 50)
            if i == 2:  grid_dense2 = max(1500./len(sp.array(z)[z/k > 0.001]), 50)
#            print i,count,len(sp.array(z)[z/k > 0.001])
#            if i == 1 or i ==2: print grid_dense1,grid_dense2
#            plt.plot(pgrid[i],z/k,label='p %i'%i)
#        plt.gca().legend(loc='upper right')
#        plt.show()
        c.add(M,M(D))
        if M(D) < chi2best:
            chi2best = M(D)
            pbest = M.p

#        if j%100 == 0:
        print j,chi2best,M(D),grid_dense1, grid_dense2

        if plot == 1:
            c.plot()

    return chi2best,pbest
