#!/bin/bash
#PBS -q long
#PBS -l nodes=8:ppn=8
#PBS -l mem=256gb
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -o calc_RR_convergence.out

cd $PBS_O_WORKDIR
ulimit -l unlimited
mpirun -np 64 ../../../cpp/build/boss/bossRR lowz_convergence_0_1.cfg
mpirun -np 64 ../../../cpp/build/boss/bossRR lowz_convergence_1_0.cfg
mpirun -np 64 ../../../cpp/build/boss/bossRR lowz_convergence_10_0.cfg
mpirun -np 64 ../../../cpp/build/boss/bossRR lowz_convergence_25_0.cfg
mpirun -np 64 ../../../cpp/build/boss/bossRR lowz_convergence_50_0.cfg

