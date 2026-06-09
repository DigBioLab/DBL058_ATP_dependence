#!/usr/bin/env python3
"""
Generate per-cluster SDF files with up to N random representative conformers.
Re-extracts molecules from local PDB cache (no re-downloads).

Usage:
    python scripts/make_cluster_sdfs.py [--output_dir output] [--n_reps 10] [--seed 42]
"""

import os, json, argparse, random
import numpy as np
import requests
from pathlib import Path
from collections import defaultdict
from rdkit import Chem
from rdkit.ML.Cluster import Butina


def get_atp_pdb_ids(max_results=9999):
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb/ATP"
    print("  Fetching ATP PDB ID list from PDBe (for ordering)...")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    atp_data = data.get("ATP", [])
    ids = []
    for item in atp_data:
        if isinstance(item, str):
            ids.append(item.upper())
        elif isinstance(item, dict):
            for key in ("pdb_id", "entry_id", "id", "pdbId"):
                if key in item:
                    ids.append(item[key].upper())
                    break
    return ids[:max_results]


def extract_atp_from_pdb_text(pdb_text, pdb_id):
    atp_lines = [l for l in pdb_text.split("\n")
                 if (l.startswith("HETATM") and l[17:20].strip() == "ATP")
                 or l.startswith("CONECT") or l.startswith("END")]
    if not atp_lines:
        return []
    groups = defaultdict(list)
    for l in atp_lines:
        if l.startswith("HETATM"):
            chain = l[21]
            resnum = l[22:26].strip()
            groups[f"{chain}_{resnum}"].append(l)
    mols = []
    for key, lines in groups.items():
        tmp = "\n".join(lines + ["END"])
        mol = Chem.MolFromPDBBlock(tmp, sanitize=False, removeHs=False)
        if mol and mol.GetNumAtoms() > 5:
            try:
                Chem.SanitizeMol(mol)
                mol.SetProp("_Name", f"{pdb_id}_{key}")
                mols.append(mol)
            except Exception:
                pass
    return mols


def verify_order(mols, cluster_meta_list):
    """Check that centroid indices match centroid names from clusters.json."""
    mismatches = 0
    for cm in cluster_meta_list[:10]:
        idx = cm["centroid_idx"]
        expected = cm["centroid_name"]
        if idx >= len(mols):
            print(f"  [!] rank {cm['rank']}: index {idx} out of range (only {len(mols)} mols)")
            mismatches += 1
            continue
        actual = mols[idx].GetProp("_Name") if mols[idx].HasProp("_Name") else ""
        if actual != expected:
            print(f"  [!] rank {cm['rank']}: expected '{expected}', got '{actual}'")
            mismatches += 1
    return mismatches


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", default="output")
    parser.add_argument("--n_reps", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    d = args.output_dir
    random.seed(args.seed)

    # ── Load metadata ─────────────────────────────────────────────────────────
    with open(f"{d}/clusters.json") as f:
        meta = json.load(f)
    cutoff = meta["cutoff_A"]
    N_expected = meta["n_conformers"]

    # ── Re-cluster from saved matrix (fast) ───────────────────────────────────
    print("Loading RMSD matrix and re-clustering...")
    rmsd_mat = np.load(f"{d}/rmsd_matrix.npy")
    N = rmsd_mat.shape[0]
    dist_list = [float(rmsd_mat[i, j]) for i in range(1, N) for j in range(i)]
    clusters = Butina.ClusterData(dist_list, N, cutoff, isDistData=True, reordering=True)
    clusters = sorted(clusters, key=len, reverse=True)
    print(f"  {len(clusters)} clusters, {N} conformers total")

    # ── Re-extract molecules from local cache ─────────────────────────────────
    print("\nRe-extracting molecules from local PDB cache (no downloads)...")
    pdb_ids = get_atp_pdb_ids(9999)
    cache_dir = Path("pdb_cache")
    mols = []
    skipped = 0
    for i, pdb_id in enumerate(pdb_ids):
        if i % 500 == 0:
            print(f"  {i}/{len(pdb_ids)} PDB IDs processed, {len(mols)} molecules extracted...")
        local = cache_dir / f"{pdb_id.lower()}.pdb"
        if not local.exists():
            skipped += 1
            continue
        extracted = extract_atp_from_pdb_text(local.read_text(), pdb_id)
        mols.extend(extracted)

    print(f"  Extracted {len(mols)} molecules ({skipped} PDB IDs not in cache)")
    if len(mols) != N_expected:
        print(f"  WARNING: got {len(mols)}, expected {N_expected} — order may differ!")

    # ── Verify ordering against saved centroid names ──────────────────────────
    print("\nVerifying molecule ordering against clusters.json centroid names...")
    mismatches = verify_order(mols, meta["clusters"])
    if mismatches == 0:
        print("  OK — molecule order matches original run")
    else:
        print(f"  {mismatches} mismatches found. Output may be incorrect. Aborting.")
        return

    # ── Write per-cluster SDFs ────────────────────────────────────────────────
    out_dir = Path(d) / "cluster_representatives"
    out_dir.mkdir(exist_ok=True)
    print(f"\nWriting representative SDFs to {out_dir}/")

    for ci, cl in enumerate(clusters):
        rank = ci + 1
        size = len(cl)
        n_sample = min(args.n_reps, size)

        # Always include centroid first, then random sample from the rest
        centroid_idx = cl[0]
        rest = list(cl[1:])
        random.shuffle(rest)
        sampled = [centroid_idx] + rest[:n_sample - 1]

        sdf_path = out_dir / f"cluster_{rank:02d}_representatives.sdf"
        with Chem.SDWriter(str(sdf_path)) as w:
            w.SetKekulize(False)
            for pos, idx in enumerate(sampled):
                mol = mols[idx]
                mol.SetProp("Cluster_Rank", str(rank))
                mol.SetProp("Cluster_Size", str(size))
                mol.SetProp("Is_Centroid", "1" if pos == 0 else "0")
                mol.SetProp("Conformer_Index", str(idx))
                w.write(mol)

        print(f"  cluster_{rank:02d}: {n_sample} reps  (cluster size {size:,})")

    print(f"\nDone. {len(clusters)} SDF files written to {out_dir}/")


if __name__ == "__main__":
    main()
