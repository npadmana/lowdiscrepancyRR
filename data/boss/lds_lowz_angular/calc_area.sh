#!/bin/bash
#PBS -q long
#PBS -l nodes=8:ppn=8
#PBS -l mem=256gb
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -o calc_area.out

cd $PBS_O_WORKDIR
ulimit -l unlimited
mpirun -np 64 ../../../cpp/build/boss/bossmask_area lowz_area_input.cfg 

