#!/usr/bin/env python
"""
Generate random catalogs in a 0-1 box and computer their RR pair-counts.
"""
import numpy as np
import sys
import subprocess

from mpi4py import MPI

script = """#!/bin/bash
#PBS -q long
#PBS -l nodes=NODES:ppn=8
#PBS -l walltime=48:0:00
#PBS -l mem=RAMgb
#PBS -o ../temp/NN_nn.out
#PBS -j oe

cd astronomy/lowdiscrRR/python/

mpirun -np NCPU python box_test_N.py -s SEED -b distancebins_small_log2.dat --Nr NN -n nn
"""
bashfilename = '../temp/NN_nn.sh'

def do_one(Nr,n):
    """Generate a random file with Nr points, labeled n, and compute RR on it."""
    if Nr < 5e4:
        NODES = 1
    elif Nr < 2e5:
        NODES = 2
    else:
        NODES = 4
    NCPU = 8*NODES
    RAM = 40*NODES

    bashname = bashfilename.replace('NN',str(Nr)).replace('nn',str(n))
    bash = script.replace('NODES',str(NODES)).replace('RAM',str(RAM)).replace('NCPU',str(NCPU))
    bash = bash.replace('NN',str(Nr)).replace('nn',str(n))
    bash = bash.replace('SEED',str(Nr))
    outfile = open(bashname,'w')
    outfile.write(bash)
    outfile.close()
    print 'Started:',bashname
    subprocess.call('qsub '+bashname,shell=True)
#...

def main(argv=None):
    if argv is None: argv = sys.argv[1:]
    from optparse import OptionParser, OptionGroup
    
    usage = '%prog [OPTION]'
    usage += '\n\nGenerate (min,...,min*k^maxpow) random points and calculate RR.'
    parser = OptionParser(usage)
    parser.add_option('-n',dest='n',type=int,default=10,
                      help='Number of random catalogs of each size to create (%default)')
    parser.add_option('--min',dest='min',type=int,default=1000,
                      help='Minimum number of randoms (%default)')
    parser.add_option('--maxpow',dest='maxpow',type=int,default=5,
                      help='Maximum power of k to multiply min by (%default)')
    parser.add_option('--k',dest='k',type=int,default=2,
                      help='Power to raise to when multiplying by min (%default)')
    (opts,args) = parser.parse_args(argv)

    Nrands = opts.min*(opts.k**np.arange(0,opts.maxpow,dtype=np.int))

    for Nr in Nrands:
        print 'Doing:',Nr
        do_one(Nr,opts.n)
#...

if __name__ == "__main__":
    sys.exit(main())
