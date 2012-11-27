#!/bin/bash
#PBS -q long
#PBS -l nodes=1
#PBS -l walltime=20:00:00
#PBS -j oe

cd $PBS_O_WORKDIR

python runner.py 
