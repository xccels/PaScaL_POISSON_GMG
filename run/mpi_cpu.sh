#!/bin/bash
#SBATCH -J mpi_cpu_job
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH --time 00:30:00

module purge

#NVIDIA Compiler
module load nvhpc/23.7 

#Intel Compiler
#module load mpi

mpirun ./a.out