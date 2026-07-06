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
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
	source "$CONDA_BASE/etc/profile.d/conda.sh"
	conda activate atp_project_env || echo "Warning: conda activate failed"
else
	echo "Warning: conda.sh not found at $CONDA_BASE/etc/profile.d/conda.sh"
fi

# change to repo root (resolve even when sourced/run from different cwd)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir/.." || exit 1

# Allow running outside LSF for testing; prefer LSF index when available
BAND_IDX="${LSB_JOBINDEX:-${LSB_ARRAY_INDEX:-}}"
if [ -z "$BAND_IDX" ]; then
	echo "Warning: LSF array index not set; defaulting to band 1"
	BAND_IDX=1
fi
python scripts/02_merge_rmsd_bands.py
