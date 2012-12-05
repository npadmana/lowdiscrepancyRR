#!/usr/bin/env python
"""Generate some real (random) data in the DEEP2 mask."""
import numpy as np
import pyfits
from deep2_angpairs import Deep2Pairs

def make_data(d2p,N):
    """Save a data file with roughly N "real" objects in it, all in non-zero
    parts of the deep2 mask."""
    # increase to cover those points that land in weight=0 regions.
    points = d2p(N*1.08)
    points = points[points[:,2] != 0]
    outfile = '../data/deep2/deep2_data_%d.dat'%N
    out = open(outfile,'w')
    out.write('#RA Dec weight\n')
    np.savetxt(out,points,fmt='%12.10f')
    out.close()
#...

d2p = Deep2Pairs()
make_data(d2p,1000)
make_data(d2p,20000)
