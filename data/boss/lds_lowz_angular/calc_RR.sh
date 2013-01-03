#!/bin/bash
#PBS -q long
#PBS -l nodes=8:ppn=8
#PBS -l mem=256gb
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -o calc_RR.out

cd $PBS_O_WORKDIR
ulimit -l unlimited
mpirun -np 64 ../../../cpp/build/boss/bossRR lowz_input.cfg 

