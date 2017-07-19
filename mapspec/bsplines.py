import scipy as sp
from scipy.interpolate import splrep,splev

class Bspline(object):
    """A class to wrap around scipy bspline fits
    """
    def __init__(self,x,y,order=3):
        self.x = x
        self.y = y
        self._order = order
        
        self.tck = splrep(self.x,self.y, k = self._order)

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self,o_new):
        self._order = o_new
        self.tck = splrep(self.x,self.y, k = self._order)


    def __call__(self,xnew):
        return splev(xnew,self.tck,ext=2)
