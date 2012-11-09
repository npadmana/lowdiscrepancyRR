#!/usr/bin/env python
"""
Compute pair counts on N RR pairs.
"""
import numpy as np
from pairs import mr_wpairs
import sys
import time

from mpi4py import MPI

def do_counting(comm,rank,bins,data1,data2):
    """
    Count pairs in bins between data1 and data2.
    To count a dataset against itself, you must pass data.copy() as data2.
    """
    t1 = time.time()
    print '%d of %d: Setting up...'%(comm.rank,comm.size)
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

def write_1d(bins,counts,outfilename):
    """Write the (1-dimensional) results to outfilename."""
    outfile = file(outfilename,'wb')
    header = ('Start','Stop','Count\n')
    outfile.write('\t'.join(header))
    for x in zip(bins[:-1],bins[1:],counts):
        outfile.write(('%10f\t%10f\t%10f\n')%x)
    print 'Wrote results to:',outfilename
    sys.stdout.flush()
#...

def do_one(comm,rank,bins,Nr,outname):
    points = np.random.random((Nr,3))
    points2 = points.copy()
    if rank == 0:
        print '--------- Processing:',outname
    counts = do_counting(comm,rank,bins,points,points2)
    if rank == 0:
        write_1d(bins,counts,outname)
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
    (opts,args) = parser.parse_args(argv)

    bins = np.loadtxt(opts.binfile)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    outfilename = '../data/NN_nn.dat'

    for n in range(opts.n):
        outname = outfilename.replace('NN',str(opts.Nr)).replace('nn',str(n))
        do_one(comm,rank,bins,opts.Nr,outname)
#...
            
if __name__ == "__main__":
    sys.exit(main())
