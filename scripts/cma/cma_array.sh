#!/bin/bash

# Request resources (per task):
#SBATCH -c 1           # 1 CPU core
#SBATCH --mem=8G       # 1 GB RAM
#SBATCH --time=3:0:0   # 6 hours (hours:minutes:seconds)
#SBATCH --output=./slurm_output/optimise_%A_%a.out

# Run on the shared queue
#SBATCH -p shared

#SBATCH --array=0-15

# Run program:
module load python/3.13.9
module load openmpi
source .venv/bin/activate
python scripts/cma/test_predictions_sweep.py $1 ${SLURM_ARRAY_TASK_ID}

