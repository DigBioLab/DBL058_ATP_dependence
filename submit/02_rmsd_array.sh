#!/bin/bash
#BSUB -J rmsd[1-40]
#BSUB -q hpc
#BSUB -n 1
#BSUB -W 12:00
#BSUB -R "rusage[mem=16GB]"
#BSUB -R "span[hosts=1]"
#BSUB -o logs/rmsd_%J_%I.stdout
#BSUB -B
#BSUB -N

source /work3/s212001/miniforge3/etc/profile.d/conda.sh
conda activate atp_project_env

cd /work3/s212001/ATP_dependent
python scripts/02_compute_rmsd_band.py --band-idx $LSB_JOBINDEX --n-bands 40
