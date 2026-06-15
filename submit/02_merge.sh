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

CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate atp_project_env

cd "$(dirname "$(realpath "$0")")/.."
python scripts/02_merge_rmsd_bands.py
