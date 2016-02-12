import scipy as sp
from scipy.integrate import simps
import matplotlib.pyplot as plt
from spectrum import *
from copy import deepcopy
#these are 'probalists', turns out to matter in order to match van der
#Marel & Franx 1993
from numpy.polynomial.hermite_e import HermiteE as H
import re


__all__ = ["RescaleModel","Chain","get_cc","metro_hast"]

debug = True

class RescaleModel(object):
    """
    This model stores parameters (shift, scale, and convolution
    kernel), performs the operation, and evaluates the likelihood.
    """

    def __init__(self,Lref,kernel='Hermite'):

        #this is an emission line object, to which we scale to 
        self.Lref = Lref  

        #empty dictionary that will hold prior distributions.  keys
        #will be same as parameters, value will be a function that
        #evaluates the prior probability at the input parameter.
        self.prior_prob = {}


        if kernel == 'Delta':
            self._get_kernel = lambda x: [1.0]
            self.p = {'shift':0.0, 'scale':1.0}
            self.set_scale = {'shift':0.05, 'scale':0.02}

        elif kernel == 'Gauss':
            self._get_kernel = self._Gauss

            #params are shift, scale, and width of gaussian
            #convolution
            width_min = 0.51*(self.Lref.wv[1] - self.Lref.wv[0])
            self.p         = {'shift':0.00, 'scale':1.00, 'width': width_min}
            self.set_scale = {'shift':0.05, 'scale':0.02, 'width':0.30}


        if kernel == 'Hermite':
            self._get_kernel = self._Hermite

            #params are shift, scale, width, h3, and h4 of
            #gauss-hermite polynomial
            width_min = 0.46*(self.Lref.wv[1] - self.Lref.wv[0])
            self.p         = {'shift':0.00, 'scale':1.00, 'width': width_min,
                              'h3':0.0, 'h4':0.0}
            self.set_scale = {'shift':0.05, 'scale':0.02, 'width':0.30,
                              'h3':0.03, 'h4':0.03}

    def __call__(self,L):
        y,var,mask  = self._get_yz(L)
        lnlikely  = self._get_chi2(y,var,mask)
        lnlikely += self._add_priors()

        #take care of limited parameter space by setting prior
        #probability to 0 [ -ln(prob) = inf]
        lnlikely += self._prior_limits()

        return lnlikely


    def _get_chi2(self,y,v,m):
        return sp.sum(    (self.Lref.f[m] - y)**2/(self.Lref.ef[m]**2 + v))


    def _add_priors(self):
        prior = 0
        for key in self.prior_prob.keys():
            prior += -2.*sp.log(
                self.prior_prob[key](self.p[key]) 
                )

        return prior


    def make_dist_prior(self,C,pname,burn = 0.5):
        params = sp.transpose(C.pchain)
        prior_dist = params[  C.index[pname]   ]
        icut = prior_dist.size*burn
        prior_dist = prior_dist[icut::]
        c1,m,c2 = sp.percentile(prior_dist,[16,50,84])
        print c1,m,c2

        #will model the distribution as a  gaussian, for now....
        self.prior_prob[pname] = lambda x: sp.exp(-0.5*(x - m)**2/ ( (c2-c1)/2 )**2 )
        self.p[pname] = m
        



    def make_func_prior(self,pname,func,params):
        self.prior_prob[pname] = lambda x: func(x,params) 


    def output(self,S):
        s = deepcopy(S)
        #convolve
        k = self._get_kernel(s)
        s.ef = sp.sqrt(sp.convolve(s.ef**2,k,mode='same'))
        s.f = sp.convolve(s.f,k,mode='same')
        #shift
        s.wv -= self.p['shift']
        #scale
        s.ef *= self.p['scale']
        s.f  *= self.p['scale']

        m = (S.wv >= s.wv.min()  )*(S.wv <= s.wv.max()  )
        y,z = s.interp(S.wv[m])

        #edge effects?????********
        sout = Spectrum()
        sout.wv = deepcopy(S.wv[m])
        sout.f  = y
        sout.ef = z

        return sout,m


    def _get_yz(self,L):
        l = deepcopy(L)
         #convolve
        k = self._get_kernel(l)
        l.ef = sp.sqrt(sp.convolve(l.ef**2,k,mode='same'))
        l.f = sp.convolve(l.f,k,mode='same')
        #shift
        l.wv -= self.p['shift']
        #scale
        l.ef *= self.p['scale']
        l.f  *= self.p['scale']


        #trim 10% of data to help with edge effects and shifting the
        #data.  This number is hard-coded so that the degrees of
        #freedom are kept fixed during the fit.
        trim = round(0.05*self.Lref.wv.size)
        m = (self.Lref.wv >= l.wv[trim] )*(self.Lref.wv <= l.wv[-trim] )

        y,z = l.interp(self.Lref.wv[m])

        if (z**2 < 0).any(): 
            raise ValueError('Weights (variances) became negative')

        return y,z**2,m

    def _Gauss(self,line):
        dlambda = line.wv[1] - line.wv[0]
        pixwidth = self.p['width']/dlambda
        prange = sp.r_[ -(line.wv.size //2) + 1 : (line.wv.size)//2  ]
        assert prange.size %2 == 1

        k = sp.exp(-0.5* prange**2/pixwidth**2)
        k /= abs(sp.sum(k))
        return k

    def _Hermite(self,line):
        dlambda = line.wv[1] - line.wv[0]
        pixwidth = self.p['width']/dlambda
        prange = sp.r_[ -(line.wv.size //2) + 1 : (line.wv.size)//2  ]
        
        h = H([ 1.0, 0.0, 0.0,self.p['h3'], self.p['h4'] ])
        #although the constants will divide out when normalizing the
        #kernel, they are important for making sure that h3 and h4 are
        #defined correctly.  Technically, this only matters for
        #choosing good intervals of h3 and h4, otherwise it is just a
        #mismatch of units

        #Equations are defined in van der Marel & Franx 1993
        k = 1./pixwidth/sp.sqrt(2*sp.pi)*sp.exp(-0.5*prange**2/pixwidth**2)*h(prange/pixwidth)
        k /= abs(sp.sum(k))


        return k


    def step(self):
        pout = {}
        for key in self.p.keys():
            pout[key] = self.p[key] + self.set_scale[key]*sp.randn()

        #there is a peculiarity that the variances go negative on the
        #GH model if the kernel is too narrow.  Although the
        #ln-probability will be -inf, this step prevents evaluation of
        #negative variances, which would otherwise force the program to quit
        if len(pout) > 2:
            if pout['width']/(self.Lref.wv[1] - self.Lref.wv[0]) < 0.45:
                pout['width'] = 0.45*(self.Lref.wv[1] - self.Lref.wv[0])

        return pout

    def _prior_limits(self):
        prior = 0
        if len(self.p) > 2:
            dlambda = self.Lref.wv[1] - self.Lref.wv[0]
            #if width is too small, the kernel is undersampled.  Weird
            #things will happen, so this represents a lower limit.
            if self.p['width']/dlambda <0.5:
                prior = sp.inf

        if len(self.p) > 3:
            #experiments have found that h3 and h4 between -0.3 and
            #0.3 should be adequate (very diverse line shapes appear)
            if self.p['h3'] < -0.3:
                prior = sp.inf
            elif self.p['h3'] > 0.3:
                prior = sp.inf
            if self.p['h4'] < -0.3:
                prior = sp.inf
            elif self.p['h4'] > 0.3:
                prior = sp.inf

        return prior




class Chain(object):
    def __init__(self):
        self.pchain   = []
        self.lnlikely = []

        self.index   = {}

#        self.figure,(self.axes) = plt.subplots(len(pnames) + 1,1)

    def add(self,M,chi2):
        if len(self.pchain) == 0:
            for i,k in enumerate(M.p.keys() ):            
                self.index[k] = i

            self.figure,(self.axes) = plt.subplots( len(M.p.keys()) + 1,1)

        self.pchain.append(deepcopy( M.p.values() ))
        self.lnlikely.append(chi2)

    def save(self,ofile):
        head = 'lnlikely   '
        outindex = []
        for key in self.index.keys():
            head += key+'   '
            outindex.append(self.index[key])
        sp.savetxt(ofile,sp.c_[self.lnlikely,sp.array(self.pchain)[:,outindex]],header=head)

    def read(self,ifile):
        fin = open(ifile,'r')
        line = fin.readline()
        pname = re.split('   ',line)
        if pname[0] != '# lnlikely':
            raise ValueError('Note a mapspec chain file! (must begin with lnlikely)')

        input_chain = sp.genfromtxt(ifile)
        self.lnlikely = input_chain[:,0]
        self.pchain   = input_chain[:,1::]
        for i,p in enumerate(pname[1:-1]):
            self.index[p] = i

        assert max(self.index.values()) == sp.transpose(self.pchain).shape[0] - 1

        self.figure,(self.axes) = plt.subplots( sp.transpose(self.pchain).shape[0] + 1,1)


    def plot(self,interact = 1):
        plotp = sp.transpose(self.pchain)
        for ax in self.axes: ax.cla()

        self.axes[0].plot(self.lnlikely)
        self.axes[0].set_ylabel('ln likelihood')

        for key in self.index.keys():
            self.axes[self.index[key] + 1].plot(
                plotp[ self.index[key] ] 
                )
            self.axes[self.index[key] + 1].set_ylabel(key)

        for ax in self.axes[0:-1]:
            ax.set_xticks([])
        self.figure.subplots_adjust(hspace=0)
        if interact == 1:  
            plt.draw()
        return

    def plot_hist(self):
        plotp = sp.transpose(self.pchain)
        for ax in self.axes: ax.cla()
        
        self.axes[0].hist(self.lnlikely,bins = 0.01*len(self.lnlikely))
        self.axes[0].set_xlabel('ln likelihood')
        for key in self.index.keys():
            self.axes[self.index[key] + 1].hist(
                plotp[ self.index[key] ] , bins = 0.01*len(plotp[self.index[key]])
                )
            self.axes[self.index[key] + 1].set_xlabel(key)

#        self.figure.set_size_inches(8,5)
        self.figure.tight_layout()




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

    c = Chain()
    c.add(M,M(D))
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
                
        if i%500 == 0 :
            print i,chi2best,chi2try
            if plot ==1:
                c.plot()

    if keep == 1:
        return chi2best,pbest,accept/float(ntrial),c
    else:
        return chi2best,pbest,accept/float(ntrial)


