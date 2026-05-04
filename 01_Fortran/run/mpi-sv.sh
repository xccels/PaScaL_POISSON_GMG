#!/bin/bash
#SBATCH -J TDM_core64_run10
#SBATCH -p batch
##SBATCH -w cpu05
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH --time 02:30:00
##SBATCH --exclude=gpu02

# make clean
# make all

nproc=16
echo "mpirun -np $nproc ./a.out"
mpirun -np $nproc ./a.out