#!/usr/bin/env python
"""
Compute pair counts on N RR pairs, in either a periodic box, or angular pairs 
in the DEEP2 mask.
"""
import numpy as np
from pairs import mr_wpairs,mr_angpairs
import os
import sys
import glob
import time

from mpi4py import MPI

import deep2_angpairs

def do_counting(comm,rank,bins,data1,data2,angular=False):
    """
    Count pairs in bins between data1 and data2.
    To count a dataset against itself, you must pass data.copy() as data2.
    If angular is set, assume data[:,2] is the weights.
    """
    t1 = time.time()
    print '%d of %d: Setting up...'%(comm.rank,comm.size)
    if angular:
        wpairs = mr_angpairs.angular_wpairs(comm,data1[:,:2],data2[:,:2],w1=data1[:,2],w2=data2[:,2])
    else:
        wpairs = mr_wpairs.radial_wpairs(comm,data1,data2,minpart=200,ntrunc=1000)
    t2 = time.time()
    print '%d of %d: Setup time = %8.5f'%(comm.rank,comm.size,t2-t1)
    sys.stdout.flush()

    t1 = time.time()
    print '%d of %d: Computing counts in bins...'%(comm.rank,comm.size)
    counts = wpairs(bins)
    t2 = time.time()
    print '%d of %d: Processing time = %8.5f'%(comm.rank,comm.size,t2-t1)
    sys.stdout.flush()

    return counts
#...

def write_1d(bins,counts,outfilename,weightsum=None):
    """Write the (1-dimensional) results to outfilename."""
    outfile = file(outfilename,'wb')
    if weightsum:
        outfile.write(str(weightsum)+'\n')
    header = ('Start','Stop','Count\n')
    outfile.write('\t'.join(header))
    for x in zip(bins[:-1],bins[1:],counts):
        outfile.write(('%10f\t%10f\t%10f\n')%x)
    print 'Wrote results to:',outfilename
    sys.stdout.flush()
#...

def do_one(comm,rank,bins,Nr,outname,d2p=None,data=None,lowdiscr=None):
    """
    Generate one random realization, and compute either its
    auto-correlation, if data == None
    cross-correlation if data is an array.
    """
    # careful! want to only have one of these floating around...
    if rank == 0:
        if d2p is not None:
            if lowdiscr is not None:
                points = d2p(Nr,vectorShift=True)
            else:
                points = d2p(Nr)
            points2 = points.copy()
            weightsum = points[:,2].sum()
        else:
            points = np.random.random((Nr,3))
            points2 = points.copy()
            weightsum = None
    else:
        points = None
        points2 = None
    points = comm.bcast(points,root=0)
    if data is not None:
        points2 = data.copy()
    else:
        points2 = comm.bcast(points2,root=0)
    if rank == 0:
        print '--------- Processing:',outname

    angular = (d2p is not None)

    counts = do_counting(comm,rank,bins,points,points2,angular=angular)
    if rank == 0:
        write_1d(bins,counts,outname,weightsum=weightsum)
#...

def main(argv=None):
    if argv is None: argv = sys.argv[1:]
    from optparse import OptionParser, OptionGroup
    
    usage = '%prog [OPTION]'
    parser = OptionParser(usage)
    parser.add_option('-n',dest='n',type=int,default=10,
                      help='Number of random catalogs of each size to create (%default)')
    parser.add_option('--Nr',dest='Nr',type=int,default=1000,
                      help='Number of random points (%default)')
    parser.add_option('-b','--bins',dest='binfile',default='distancebins_small_log2.dat',
                      help='file containing the bins to compute pairs in (%default)')
    parser.add_option('-s','--seed',dest='seed',default='1',type=int,
                      help='random number seed for numpy (%default).')
    parser.add_option('--deep2',dest='deep2',default=False,action='store_true',
                      help='Use the DEEP2 mask for angular pairs (%default).')
    parser.add_option('--data',dest='data',default='',
                      help='Use this file to calculate DR, instead of RR. Expected to have 2 or 3 position columns + 1 weight column, and the header marked with "#".')
    parser.add_option('--lowdiscr',dest='lowdiscr',default='',
                      help='Use this low discrepency sequence when calculating DR.')
    (opts,args) = parser.parse_args(argv)

    bins = np.loadtxt(opts.binfile)

    np.random.seed(opts.seed)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if opts.lowdiscr != '':
        lowdiscr = np.loadtxt(opts.lowdiscr)
        filename = 'NN_nn-lowdiscr.dat'
    else:
        lowdiscr = None
        filename = 'NN_nn.dat'

    if opts.deep2:
        d2p = deep2_angpairs.Deep2Pairs(randarr=lowdiscr)
        outfilename = '../data/deep2/counts/'+filename
    else:
        d2p = None
        outfilename = '../data/counts/'+filename

    if opts.data != '':
        data = np.loadtxt(opts.data)
        outfilename = outfilename.replace('nn','nn-DR')
        Ndata = opts.data.split('_')[-1].split('.')[0]
        outfilename = outfilename.replace('counts','counts_D'+Ndata)
    else:
        data = None

    if rank == 0:
        if not os.path.exists(os.path.dirname(outfilename)):
            os.mkdir(os.path.dirname(outfilename))

    for n in range(opts.n):
        outname = outfilename.replace('NN',str(opts.Nr)).replace('nn',str(n))
        do_one(comm,rank,bins,opts.Nr,outname,d2p=d2p,data=data,lowdiscr=lowdiscr)
#...
            
if __name__ == "__main__":
    sys.exit(main())
