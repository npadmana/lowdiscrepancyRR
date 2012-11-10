#!/bin/bash
#PBS -q long
#PBS -l nodes=16:ppn=8
#PBS -l walltime=120:00:00
#PBS -l mem=640gb
#PBS -o do_2point_ang.out
#PBS -j oe
N_CPU=128    #TS:nodes*ppn

cd ~/astronomy/fastRR/

# 1d redshift space correlation function
BINS_1D="-b distancebins_angular_10.dat"

# angular correlation function
RANDOMS=boss_geometry_2011_06_10-randoms-1000000.txt

echo ./mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D $RANDOMS
./mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D $RANDOMS
