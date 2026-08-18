#!/bin/bash
#SBATCH -J GMG_GPU_test
#SBATCH --comment="field=mae;appl=in_house"
#SBATCH --partition=amd_a100nv_8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH -o run/%x.o%j
#SBATCH -e run/%x.e%j
#SBATCH --time=02:00:00

cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc/15.2.0
module load mpi/openmpi-4.1.8
module load nvhpc/25.11_

export OMP_NUM_THREADS=1

echo "Running on GPU..."
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_NODELIST = $SLURM_NODELIST"
echo "SLURM_NTASKS = $SLURM_NTASKS"


mkdir -p run/result

srun -n 2 ./run/test_poisson