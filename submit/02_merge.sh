#!/bin/bash
#BSUB -J rmsd_merge
#BSUB -q hpc
#BSUB -n 1
#BSUB -W 1:00
#BSUB -R "rusage[mem=32GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o logs/merge_%J.stdout
#BSUB -B
#BSUB -N

source /work3/s212001/miniforge3/etc/profile.d/conda.sh
conda activate atp_project_env

cd /work3/s212001/ATP_dependent
python scripts/02_merge_rmsd_bands.py
