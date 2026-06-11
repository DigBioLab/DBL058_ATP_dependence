#!/bin/bash
#BSUB -J atp_download
#BSUB -q hpc
#BSUB -n 4
#BSUB -W 4:00
#BSUB -R "rusage[mem=16GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o logs/download_%J.stdout
#BSUB -B
#BSUB -N

source /work3/s212001/miniforge3/etc/profile.d/conda.sh
conda activate atp_project_env

cd /work3/s212001/ATP_dependent
python scripts/01_download_and_preprocess.py
