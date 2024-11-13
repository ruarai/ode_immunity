#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --partition=sapphire
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --array=1-8350
#;;;;SBATCH --array=981,1019,1168,1235,1258,1683,1712,1935,2046,2256,2257,2352,2415,2461,2818,3231,3344,3459,3611,3703,4150,4590,4592,4848,4856,5006,5033,5079,5243,5499,5525,5603,5791,6023,6047,6059,6600,6742,7060,7061,7173,7377,7410,7476,7521,7800,7805,7809,7820,7833,7847,8057
#SBATCH --job-name=juliamulti
#SBATCH -o ./slurm_out/period_grid_%a.out # STDOUT

cd /data/gpfs/projects/punim1541/ode_immunity/

module load Julia/1.9.3-linux-x86_64 
module load SUNDIALS/6.3.0

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia src/dist/period_over_grid.jl ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_MAX}
#julia src/dist/period_over_grid.jl ${SLURM_ARRAY_TASK_ID} 8350
