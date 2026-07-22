#!/bin/bash

# Request resources (per task):
#SBATCH -c 1           # 1 CPU core
#SBATCH --mem=8G       # 1 GB RAM
#SBATCH --time=70:0:0   # 6 hours (hours:minutes:seconds)

# Run on the shared queue
#SBATCH -p shared

# Each separate task can be identified based on the SLURM_ARRAY_TASK_ID
# environment variable:

# Run program:
module load python/3.13.9
module load openmpi
source .venv/bin/activate
python scripts/cma/run_cma.py $1
