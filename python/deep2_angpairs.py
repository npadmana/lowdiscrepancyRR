#!/usr/bin/env python
"""Generate random angular points in the deep2 mask."""

import numpy as np
import pyfits

class Deep2Pairs():
    def __init__(self):
        data = pyfits.open('../data/deep2/windowf.31.fits')
        header = data[0].header
        data = data[0].data
        self.data = data

        self.RAmin = header['RALIM0']
        self.RAmax = header['RALIM1']
        self.nRA = data.shape[1]
        self.dRA = (self.RAmax - self.RAmin)/self.nRA

        self.DECmin = header['DECLIM0']
        self.DECmax = header['DECLIM1']
        self.nDEC = data.shape[0]
        self.dDEC = (self.DECmax - self.DECmin)/self.nDEC
    #...

    def __call__(self,N):
        """Return N points in the deep2 mask, with their weights."""
        RA = np.random.uniform(self.RAmin,self.RAmax,N)
        zmin = np.sin(np.radians(self.DECmin))
        zmax = np.sin(np.radians(self.DECmax))
        z = np.random.uniform(zmin,zmax,N)
        Dec = np.degrees(np.arcsin(z))
        weights = self.get_weights(RA,Dec)
        return np.array(zip(RA,Dec,weights))
    #...

    def get_weights(self,RA,Dec):
        """Find the weights associated with the given RA/Dec pairs."""
        iRA = np.array((RA - self.RAmin)/self.dRA,'i4')
        iDEC = np.array((Dec - self.DECmin)/self.dDEC,'i4')
        return self.data[iDEC,iRA]
    #...
#...

np.random.seed(100)
format = '%14.10f'
dp = Deep2Pairs()
Npoints = 1000
filename = '../data/deep2/deep2_angular_%d.txt'
np.savetxt(filename%Npoints,dp(Npoints),fmt=format)
Npoints = 10000
np.savetxt(filename%Npoints,dp(Npoints),fmt=format)
Npoints = 100000
np.savetxt(filename%Npoints,dp(Npoints),fmt=format)
