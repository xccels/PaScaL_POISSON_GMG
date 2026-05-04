#!/bin/bash
#SBATCH -J GMG_C_core2
#SBATCH -p batch
#SBATCH -w cpu05
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH -o run/%x.o%j
#SBATCH -e run/%x.e%j
#SBATCH --time 02:30:00
##SBATCH --exclude=gpu02

cd "$SLURM_SUBMIT_DIR"


export OMP_NUM_THREADS=1

nproc=2
EXE="run/poisson"
INPUT="run/PARA_INPUT.inp"

# echo "Working directory: $(pwd)"
# echo "Executable:        ${EXE}"
# echo "Input file:        ${INPUT}"
echo "mpirun -np ${nproc} ${EXE} ${INPUT}"

mpirun -np ${nproc} ${EXE} ${INPUT}
