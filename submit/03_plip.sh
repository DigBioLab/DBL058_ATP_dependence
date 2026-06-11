#!/bin/bash
#BSUB -J atp_plip
#BSUB -q hpc
#BSUB -n 16
#BSUB -W 12:00
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o logs/plip_%J.stdout
#BSUB -B
#BSUB -N

source /work3/s212001/miniforge3/etc/profile.d/conda.sh
conda activate atp_project_env

cd /work3/s212001/ATP_dependent
python scripts/03_run_plip.py --n-workers 16
