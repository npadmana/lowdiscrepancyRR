#!/bin/bash
#PBS -q long
#PBS -l nodes=16:ppn=8
#PBS -l walltime=120:00:00
#PBS -l mem=640gb
#PBS -o do_2point_ang.out
#PBS -j oe
N_CPU=128    #TS:nodes*ppn

cd ~/astronomy/lowdiscrRR/python/

#BINS_1D="-b distancebins_angular_10.dat"

# angular correlation function
#RANDOMS=../data/boss_geometry_2011_06_10-randoms-100000.txt
#echo mpirun -np $N_CPU python mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D -o angpairs.dat $RANDOMS
#mpirun -np $N_CPU python mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D -o angpairs_100000.dat $RANDOMS

#RANDOMS=../data/boss_geometry_2011_06_10-randoms-200000.txt
#echo mpirun -np $N_CPU python mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D -o angpairs.dat $RANDOMS
#mpirun -np $N_CPU python mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D -o angpairs_200000.dat $RANDOMS

#RANDOMS=../data/boss_geometry_2011_06_10-randoms-1000000.txt
#echo mpirun -np $N_CPU python mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D -o angpairs.dat $RANDOMS
#mpirun -np $N_CPU python mr_2pc_bin.py --counter angular --skiprows 1 $BINS_1D -o angpairs_1000000.dat $RANDOMS

BINS_1D="-b distancebins_angular_0.25.dat"

RANDOMS=../data/deep2/deep2_angular_1000.txt
echo mpirun -np $N_CPU python mr_2pc_bin.py -w --counter angular -o angpairs_1000.dat $BINS_1D $RANDOMS
mpirun -np $N_CPU python mr_2pc_bin.py -w --counter angular -o angpairs_1000.dat $BINS_1D $RANDOMS

RANDOMS=../data/deep2/deep2_angular_10000.txt
echo mpirun -np $N_CPU python mr_2pc_bin.py -w --counter angular -o angpairs_10000.dat $BINS_1D $RANDOMS
mpirun -np $N_CPU python mr_2pc_bin.py -w --counter angular  -o angpairs_10000.dat $BINS_1D $RANDOMS

RANDOMS=../data/deep2/deep2_angular_100000.txt
echo mpirun -np $N_CPU python mr_2pc_bin.py -w --counter angular -o angpairs_1000000.dat $BINS_1D $RANDOMS
mpirun -np $N_CPU python mr_2pc_bin.py -w --counter angular -o angpairs_100000.dat $BINS_1D $RANDOMS
