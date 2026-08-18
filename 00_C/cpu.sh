#!/bin/bash
#SBATCH -J GMG_CPU_test
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --comment="field=mae;appl=in_house"
#SBATCH -o run/%x.o%j
#SBATCH -e run/%x.e%j
#SBATCH --time=02:00:00

cd "$SLURM_SUBMIT_DIR"


export OMP_NUM_THREADS=1

nproc=4
EXE="run/poisson"
INPUT="run/PARA_INPUT.inp"

# echo "Working directory: $(pwd)"
# echo "Executable:        ${EXE}"
# echo "Input file:        ${INPUT}"
echo "mpirun -np ${nproc} ${EXE} ${INPUT}"

# mpirun -np ${nproc} ${EXE} ${INPUT}

srun --ntasks="$SLURM_NTASKS" --cpus-per-task=1 \
     "$EXE" "$INPUT"