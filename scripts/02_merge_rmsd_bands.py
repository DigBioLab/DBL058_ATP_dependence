"""
Merge all RMSD band files produced by 02_compute_rmsd_band.py into a single
symmetric N×N matrix.

Run after the array job completes (submit/02_merge.sh).
Saves: output/rmsd_matrix_filtered.npy
"""

from pathlib import Path
import numpy as np

BANDS_DIR  = Path('output/rmsd_bands')
OUT_PATH   = Path('output/rmsd_matrix_filtered.npy')
SDF_PATH   = Path('output/atp_conformers_raw.sdf')


def main():
    band_files = sorted(BANDS_DIR.glob('band_*.npy'))
    if not band_files:
        raise FileNotFoundError(f'No band files found in {BANDS_DIR}')

    print(f'Found {len(band_files)} band files.')

    bands = [np.load(f) for f in band_files]
    bands = [b for b in bands if b.size > 0]

    full_mat = np.vstack(bands).astype(np.float32)
    N = full_mat.shape[1]

    if full_mat.shape[0] != N:
        raise ValueError(
            f'Assembled matrix is {full_mat.shape} — expected ({N}, {N}). '
            f'Some bands may be missing or N changed since bands were computed.'
        )

    print(f'Assembled upper triangle: {full_mat.shape}')

    full_mat = full_mat + full_mat.T

    print(f'Symmetrised. max RMSD = {full_mat.max():.3f} Å')
    np.save(OUT_PATH, full_mat)
    print(f'Saved to {OUT_PATH}')


if __name__ == '__main__':
    main()
