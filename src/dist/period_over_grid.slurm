#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --partition=sapphire
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-4920
#SBATCH --job-name=juliamulti
#SBATCH -o ./slurm_out/period_grid_%a.out # STDOUT

cd /data/gpfs/projects/punim1541/ode_immunity/

module load Julia/1.9.3-linux-x86_64 
module load SUNDIALS/6.3.0

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia src/dist/period_over_grid.jl ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_MAX}