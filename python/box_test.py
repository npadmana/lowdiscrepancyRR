#!/usr/bin/env python
"""
Generate random catalogs in a 0-1 box and computer their RR pair-counts.
"""
import numpy as np
from pairs import mr_wpairs
import sys
import subprocess

from mpi4py import MPI

script = """#!/bin/bash
#PBS -q long
#PBS -l nodes=NODES:ppn=8
#PBS -l walltime=48:0:00
#PBS -l mem=RAMgb
#PBS -o temp/NN_nn.out
#PBS -j oe

cd astronomy/fastRR/

mpirun -np NCPU python mr_2pc_bin.py -b distancebins_small_log2.dat -o OUT DATA
"""
bashfilename = 'temp/NN_nn.sh'
datafilename = 'temp/randoms_NN_nn.txt'

def do_one(Nr,n):
    """Generate a random file with Nr points, labeled n, and compute RR on it."""
    if Nr > 5e4:
        NODES = 1
    else:
        NODES = 2
    NCPU = 8*NODES
    RAM = 40*NODES
    points = np.random.random((3,Nr))
    dataname = datafilename.replace('NN',str(Nr)).replace('nn',str(n))
    np.savetxt(dataname)

    bashname = bashfilename.replace('NN',str(Nr)).replace('nn',str(n))
    bash = script.replace('NODES',str(NODES)).replace('RAM',str(RAM)).replace('NCPU',str(NCPU))
    bash = bash.replace('OUT',bashname.replace('sh','dat')).replace('DATA',dataname)
    bash = bash.replace('NN',str(Nr)).replace('nn',str(n))
    outfile = open(bashname,'w')
    outfile.write(bash)
    outfile.close()
    print 'Started:',bashname
    #subprocess.call('qsub '+bashname,shell=True)
#...

def main(argv=None):
    if argv is None: argv = sys.argv[1:]
    from optparse import OptionParser, OptionGroup
    
    import pdb
    pdb.set_trace()
    usage = '%prog [OPTION]'
    parser = OptionParser(usage)
    parser.add_option('-n',dest='n',type=int,default=3,
                      help='Number of random catalogs of each size to create (%default)')

    (opts,args) = parser.parse_args(argv)
    
    Nrands = 1000*(2**np.arange(0,8,dtype=np.int))

    for Nr in Nrands:
        print 'Doing:',Nr
        for n in range(opts.n):
            print 'Number:',n
            do_one(Nr,n)
#...

if __name__ == "__main__":
    sys.exit(main())
