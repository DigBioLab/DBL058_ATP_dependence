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

CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate atp_project_env

cd "$(dirname "$(realpath "$0")")/.."
python scripts/01_download_and_preprocess.py
