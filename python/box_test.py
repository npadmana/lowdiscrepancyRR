#!/usr/bin/env python
"""
Generate random catalogs in a 0-1 box, or the DEEP2 mask,
and computer their RR pair-counts.
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

mpirun -np NCPU python box_test_N.py MASK -s SEED -b BINS --Nr NN -n nn DATA
"""
bashfilename = '../temp/NN_nn.sh'

def do_one(Nr,n,deep2=False,data=''):
    """Generate a random file with Nr points, labeled n, and compute RR on it."""
    if Nr < 5e4:
        NODES = 1
    elif Nr < 2e5:
        NODES = 2
    else:
        NODES = 4

    if deep2:
        binfile = 'distancebins_angular_0.25.dat'
        mask = '--deep2'
        NODES *= 8
    else:
        binfile = 'distancebins_small_log2.dat'
        mask = ''

    NCPU = 8*NODES
    RAM = 40*NODES

    bashname = bashfilename.replace('NN',str(Nr)).replace('nn',str(n))
    bash = script.replace('NODES',str(NODES)).replace('RAM',str(RAM)).replace('NCPU',str(NCPU))
    bash = bash.replace('NN',str(Nr)).replace('nn',str(n))
    bash = bash.replace('SEED',str(Nr)).replace('BINS',binfile).replace('MASK',mask)
    bash = bash.replace('DATA',data)
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
    parser.add_option('--deep2',dest='deep2',default=False,action='store_true',
                      help='Use the DEEP2 mask for angular pairs (%default).')
    parser.add_option('--data',dest='data',default='',
                      help='Use this file to calculate DR, instead of RR.')
    (opts,args) = parser.parse_args(argv)

    Nrands = opts.min*(opts.k**np.arange(0,opts.maxpow,dtype=np.int))

    for Nr in Nrands:
        print 'Doing:',Nr
        do_one(Nr,opts.n,deep2=opts.deep2,data=opts.data)
#...

if __name__ == "__main__":
    sys.exit(main())
