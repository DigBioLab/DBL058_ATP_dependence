"""
Compute one horizontal band of the pairwise RMSD matrix.
Run as an LSF array job (submit/02_rmsd_array.sh); each task handles
rows [start, end) of the full N×N matrix and only computes the upper
triangle (j > i). The merge script symmetrises the result.

Saves: output/rmsd_bands/band_{band_idx:04d}.npy  shape: (band_size, N)
"""

import argparse, math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

SDF_PATH   = Path('output/atp_conformers_raw.sdf')
BANDS_DIR  = Path('output/rmsd_bands')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--band-idx', type=int, required=True,
                        help='1-based band index (set from $LSB_JOBINDEX)')
    parser.add_argument('--n-bands',  type=int, default=40,
                        help='Total number of bands (must match array job size)')
    args = parser.parse_args()

    BANDS_DIR.mkdir(parents=True, exist_ok=True)

    out_path = BANDS_DIR / f'band_{args.band_idx:04d}.npy'
    if out_path.exists():
        print(f'Band {args.band_idx} already computed — skipping.')
        return

    print(f'Loading conformers from {SDF_PATH}...')
    suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False, sanitize=False)
    mols  = [m for m in suppl if m is not None]
    N     = len(mols)
    print(f'  {N} molecules loaded.')

    band_size = math.ceil(N / args.n_bands)
    start     = (args.band_idx - 1) * band_size
    end       = min(args.band_idx * band_size, N)

    if start >= N:
        print(f'Band {args.band_idx} is out of range (N={N}) — saving empty.')
        np.save(out_path, np.zeros((0, N), dtype=np.float32))
        return

    actual_size = end - start
    print(f'Band {args.band_idx}/{args.n_bands}: rows [{start}, {end})  '
          f'({actual_size} rows, {N} cols)')

    band_mat = np.zeros((actual_size, N), dtype=np.float32)

    for local_i, global_i in enumerate(range(start, end)):
        if local_i % 20 == 0:
            print(f'  row {global_i}/{end-1}', flush=True)
        for j in range(global_i + 1, N):
            try:
                rms = rdMolAlign.GetBestRMS(mols[global_i], mols[j])
            except Exception:
                rms = 99.0
            band_mat[local_i, j] = rms

    np.save(out_path, band_mat)
    print(f'Saved {out_path}  shape={band_mat.shape}')


if __name__ == '__main__':
    main()
