#!/usr/bin/env python3
"""
Grid cavity scan — the measurement behind section 1.3 of the Fv plan, made
reproducible.

Section 1.3 reported a threshold sweep over the VH/VL interface of apo 1SY6
(0 / 49 / 208 / 578 A^3 at >=90 / 80 / 70 / 60 % buriedness) from an ad-hoc
script that was never committed. This reimplements the stated method so the
same number can be computed on RFD3 outputs and compared like for like:

    1 A cubic voxels, probe-accessible, ray-occlusion buriedness.

Method, per voxel:
  1. Keep voxels inside the region of interest and no closer to any protein
     heavy atom than PROBE + VDW (i.e. a probe sphere fits there).
  2. Discard voxels further than SHELL_MAX from the protein — bulk solvent.
  3. Buriedness = fraction of RAY_DIRS evenly spread rays from the voxel that
     hit a protein heavy atom within RAY_MAX.
  4. Group surviving voxels into connected cavities (6-connectivity) and report
     total and largest volume at each buriedness threshold.

Buriedness by ray occlusion is a crude proxy — it is sensitive to RAY_DIRS,
RAY_MAX and the vdW radii, so the sweep across thresholds carries the signal,
not any single number. Absolute volumes are not comparable to numbers from a
different cavity finder.

Usage
    # one structure, region = whole VH/VL interface
    python3 scripts/cavity_scan.py --struct inputs/okt3_fv.cif

    # around the ligand in RFD3 outputs, as a sweep summary
    python3 scripts/cavity_scan.py --glob 'diffusion_out/*/*.cif.gz' \\
        --around-ligand ATP --out cavity.csv
"""
import argparse
import glob as globmod
import gzip
import io
import json
import os
import sys

import numpy as np

# Heavy-atom vdW radii (A). Missing elements fall back to CARBON.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
VDW_DEFAULT = 1.70

PROBE = 1.4         # water probe
VOXEL = 1.0         # grid spacing, per section 1.3
SHELL_MAX = 8.0     # voxels further than this from protein are bulk solvent
RAY_MAX = 12.0      # ray hit search distance
RAY_DIRS = 42       # evenly spread directions per voxel
THRESHOLDS = (0.60, 0.70, 0.80, 0.90)


def load_structure(path):
    """(coords, elements, res_names) heavy atoms only. Handles .cif / .cif.gz."""
    from biotite.structure.io.pdbx import CIFFile, get_structure
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            cif = CIFFile.read(io.StringIO(fh.read()))
    else:
        cif = CIFFile.read(path)
    arr = get_structure(cif, model=1)
    arr = arr[arr.element != "H"]
    return arr.coord.astype(float), np.asarray(arr.element), np.asarray(arr.res_name)


def fibonacci_directions(n):
    """n roughly evenly spread unit vectors."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def connected_components(voxels):
    """6-connectivity grouping over integer voxel indices. Returns sizes."""
    index = {tuple(v): i for i, v in enumerate(voxels)}
    seen, sizes = set(), []
    nbrs = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    for start in index:
        if start in seen:
            continue
        stack, size = [start], 0
        seen.add(start)
        while stack:
            cur = stack.pop()
            size += 1
            for d in nbrs:
                nb = (cur[0] + d[0], cur[1] + d[1], cur[2] + d[2])
                if nb in index and nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def scan(coords, elements, centre, half_extent):
    """Cavity volumes by buriedness threshold within a box around `centre`."""
    radii = np.array([VDW.get(e, VDW_DEFAULT) for e in elements])

    lo = centre - half_extent
    hi = centre + half_extent
    axes = [np.arange(lo[i], hi[i] + VOXEL, VOXEL) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)

    # Chunked distance evaluation: the full grid x atoms matrix is too large.
    keep_pts, keep_idx = [], []
    ijk = np.stack(np.meshgrid(*[np.arange(len(a)) for a in axes],
                               indexing="ij"), axis=-1).reshape(-1, 3)
    for s in range(0, len(grid), 4000):
        block = grid[s:s + 4000]
        d = np.linalg.norm(block[:, None, :] - coords[None, :, :], axis=-1)
        clearance = d - radii[None, :]
        dmin = clearance.min(axis=1)
        ok = (dmin >= PROBE) & (dmin <= SHELL_MAX)
        keep_pts.append(block[ok])
        keep_idx.append(ijk[s:s + 4000][ok])
    if not any(len(k) for k in keep_pts):
        return {f"{t:.0%}": (0.0, 0.0) for t in THRESHOLDS}, 0
    pts = np.concatenate(keep_pts)
    idx = np.concatenate(keep_idx)

    # Ray-occlusion buriedness.
    dirs = fibonacci_directions(RAY_DIRS)
    steps = np.arange(1.0, RAY_MAX + 1.0, 1.0)
    buried = np.zeros(len(pts))
    for s in range(0, len(pts), 200):
        block = pts[s:s + 200]
        # sample points along every ray: (block, dirs, steps, 3)
        samples = (block[:, None, None, :]
                   + dirs[None, :, None, :] * steps[None, None, :, None])
        flat = samples.reshape(-1, 3)
        hit = np.zeros(len(flat), dtype=bool)
        for c in range(0, len(flat), 20000):
            sub = flat[c:c + 20000]
            d = np.linalg.norm(sub[:, None, :] - coords[None, :, :], axis=-1)
            hit[c:c + 20000] = (d - radii[None, :] < 1.0).any(axis=1)
        hit = hit.reshape(len(block), RAY_DIRS, len(steps)).any(axis=2)
        buried[s:s + 200] = hit.mean(axis=1)

    vox_vol = VOXEL ** 3
    result = {}
    for t in THRESHOLDS:
        sel = idx[buried >= t]
        if len(sel) == 0:
            result[f"{t:.0%}"] = (0.0, 0.0)
            continue
        sizes = connected_components(sel)
        result[f"{t:.0%}"] = (len(sel) * vox_vol, sizes[0] * vox_vol)
    return result, len(pts)


def region_for(coords, res_names, around_ligand, half_extent):
    """Box centre: ligand centroid if present, else protein centroid."""
    if around_ligand:
        mask = res_names == around_ligand
        if mask.any():
            return coords[mask].mean(axis=0), np.full(3, half_extent)
    return coords.mean(axis=0), np.full(3, half_extent)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--struct", help="single structure (.cif/.cif.gz)")
    ap.add_argument("--glob", help="glob of structures to scan")
    ap.add_argument("--around-ligand", default=None,
                    help="centre the box on this residue name (e.g. ATP)")
    ap.add_argument("--half-extent", type=float, default=12.0,
                    help="half box size in A (default 12 -> 24 A cube)")
    ap.add_argument("--out", default=None, help="write per-structure CSV here")
    ap.add_argument("--file-list", default=None,
                    help="text file with one structure path per line")
    ap.add_argument("--shard", default=None, metavar="I/N",
                    help="process only shard I of N (1-based), for LSF arrays")
    args = ap.parse_args()

    paths = []
    if args.struct:
        paths.append(args.struct)
    if args.glob:
        paths.extend(sorted(globmod.glob(args.glob)))
    if args.file_list:
        with open(args.file_list) as fh:
            paths.extend(l.strip() for l in fh if l.strip())
    if not paths:
        sys.exit("nothing to scan: pass --struct, --glob and/or --file-list")

    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        if not 1 <= i <= n:
            sys.exit(f"bad --shard {args.shard}")
        paths = paths[i - 1::n]
        print(f"shard {i}/{n}: {len(paths)} structures", flush=True)
        if not paths:
            return

    rows = []
    for i, p in enumerate(paths, 1):
        coords, elements, res_names = load_structure(p)
        # Cavities are measured in the protein only; the ligand must not fill
        # the very cavity being measured.
        prot = res_names != (args.around_ligand or "\0")
        centre, half = region_for(coords, res_names, args.around_ligand,
                                  args.half_extent)
        res, n_probe = scan(coords[prot], elements[prot], centre, half)
        row = {"structure": os.path.basename(p), "n_probe_voxels": n_probe}
        for k, (total, largest) in res.items():
            row[f"total_{k}"] = total
            row[f"largest_{k}"] = largest
        rows.append(row)
        print(f"[{i}/{len(paths)}] {os.path.basename(p)}  "
              + "  ".join(f"{k}:{v[0]:.0f}/{v[1]:.0f}" for k, v in res.items()),
              flush=True)

    if args.out:
        import pandas as pd
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\nwrote {args.out}  ({len(rows)} structures)")
        print("columns: total_X / largest_X = total and largest single cavity "
              "volume (A^3) at buriedness >= X")


if __name__ == "__main__":
    main()
