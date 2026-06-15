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

CONDA_BASE="${CONDA_BASE:-$HOME/miniforge3}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate atp_project_env

cd "$(dirname "$(realpath "$0")")/.."
python scripts/03_run_plip.py --n-workers 16
