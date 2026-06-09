#!/usr/bin/env python3
"""
ATP Conformer Analysis from PDB
================================
Downloads all ATP-containing structures from RCSB and clusters their conformations.

Usage:
    python scripts/atp_conformer_analysis.py [--test] [--cutoff 1.5] [--max_structs 3500]

Requirements:
    pip install rdkit biopython rcsbsearchapi numpy scipy matplotlib scikit-learn
"""

import os, sys, time, argparse, json, gzip
import numpy as np
import requests
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.ML.Cluster import Butina

# ─────────────────────────────────────────────────────────────────────
# 1. FETCH ATP-CONTAINING PDB IDs FROM RCSB
# ─────────────────────────────────────────────────────────────────────
def get_atp_pdb_ids(max_results=9999):
    """Get all PDB IDs containing ATP from the PDBe REST API."""
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb/ATP"
    print("Fetching ATP entry list from PDBe...")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    
    # Debug: print the actual structure so we can see what we're getting
    data = resp.json()
    sample = str(data)[:300]
    print(f"Response sample: {sample}")
    
    atp_data = data.get("ATP", [])
    
    # Handle both formats: list of strings OR list of dicts
    ids = []
    for item in atp_data:
        if isinstance(item, str):
            ids.append(item.upper())
        elif isinstance(item, dict):
            # Try common key names
            for key in ("pdb_id", "entry_id", "id", "pdbId"):
                if key in item:
                    ids.append(item[key].upper())
                    break
    
    print(f"Found {len(ids)} ATP-containing structures.")
    return ids[:max_results]

# ─────────────────────────────────────────────────────────────────────
# 2. DOWNLOAD PDB AND EXTRACT ATP HETATM RECORDS
# ─────────────────────────────────────────────────────────────────────
def download_pdb(pdb_id, cache_dir="pdb_cache"):
    """Download a PDB file, using local cache if available."""
    os.makedirs(cache_dir, exist_ok=True)
    local = Path(cache_dir) / f"{pdb_id.lower()}.pdb"
    if local.exists():
        return local.read_text()

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    resp = requests.get(url, timeout=30)
    if resp.status_code == 200:
        local.write_text(resp.text)
        return resp.text
    return None

def extract_atp_from_pdb_text(pdb_text, pdb_id):
    """Extract ATP molecule(s) from PDB text using RDKit."""
    # Write ATP HETATM lines to a temp file
    atp_lines = [l for l in pdb_text.split("\n") 
                 if (l.startswith("HETATM") and l[17:20].strip() == "ATP")
                 or (l.startswith("CONECT"))
                 or l.startswith("END")]
    if not atp_lines:
        return []
    
    # Group by chain/residue (multiple ATPs per structure possible)
    from collections import defaultdict
    groups = defaultdict(list)
    for l in atp_lines:
        if l.startswith("HETATM"):
            chain = l[21]
            resnum = l[22:26].strip()
            key = f"{chain}_{resnum}"
            groups[key].append(l)
    
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

# ─────────────────────────────────────────────────────────────────────
# 3. ALIGN ALL ATP MOLECULES TO A REFERENCE
# ─────────────────────────────────────────────────────────────────────
def build_reference_atp():
    """Generate a clean ATP reference conformer from SMILES."""
    atp_smiles = "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
    mol = Chem.MolFromSmiles(atp_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

def align_to_reference(mols, ref_mol):
    """Align all mols to ref_mol using MCS (maximum common substructure)."""
    from rdkit.Chem import rdFMCS
    aligned = []
    mcs = rdFMCS.FindMCS([mols[0], ref_mol], 
                          atomCompare=rdFMCS.AtomCompare.CompareElements,
                          bondCompare=rdFMCS.BondCompare.CompareAny,
                          timeout=5)
    ref_pattern = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref_mol.GetSubstructMatch(ref_pattern)
    
    for mol in mols:
        probe_match = mol.GetSubstructMatch(ref_pattern)
        if probe_match and ref_match:
            atom_map = list(zip(probe_match, ref_match))
            try:
                rdMolAlign.AlignMol(mol, ref_mol, atomMap=atom_map)
                aligned.append(mol)
            except Exception:
                aligned.append(mol)  # keep unaligned
        else:
            aligned.append(mol)
    return aligned

# ─────────────────────────────────────────────────────────────────────
# 4. COMPUTE PAIRWISE RMSD & CLUSTER
# ─────────────────────────────────────────────────────────────────────
def compute_rmsd_matrix(mols):
    """Compute pairwise heavy-atom RMSD matrix."""
    N = len(mols)
    mat = np.zeros((N, N))
    for i in range(N):
        for j in range(i+1, N):
            try:
                rms = rdMolAlign.GetBestRMS(mols[i], mols[j])
            except Exception:
                rms = 99.0
            mat[i,j] = mat[j,i] = rms
    return mat

def cluster_butina(rmsd_matrix, cutoff=1.5):
    """Butina clustering on RMSD matrix."""
    N = rmsd_matrix.shape[0]
    dist_list = [rmsd_matrix[i,j] for i in range(1,N) for j in range(i)]
    clusters = Butina.ClusterData(dist_list, N, cutoff, isDistData=True, reordering=True)
    return sorted(clusters, key=len, reverse=True)

# ─────────────────────────────────────────────────────────────────────
# 5. MAIN
# ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true", help="Use 50 synthetic conformers")
    parser.add_argument("--cutoff", type=float, default=1.5, help="RMSD cutoff in Angstrom")
    parser.add_argument("--max_structs", type=int, default=3500)
    parser.add_argument("--output_dir", default="output")
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.test:
        print("=== TEST MODE: Generating synthetic ATP conformers ===")
        atp_smiles = "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
        mol = Chem.MolFromSmiles(atp_smiles)
        mol_h = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=50, randomSeed=42)
        with Chem.SDWriter(f"{args.output_dir}/atp_conformers.sdf") as w:
            w.SetKekulize(False)
            for cid in cids: w.write(mol_h, confId=cid)
        mols = [m for m in Chem.SDMolSupplier(f"{args.output_dir}/atp_conformers.sdf", removeHs=False) if m]
    else:
        print("=== PRODUCTION MODE: Downloading from PDB ===")
        pdb_ids = get_atp_pdb_ids(args.max_structs)
        mols = []
        for i, pdb_id in enumerate(pdb_ids):
            if i % 100 == 0: print(f"  Processing {i}/{len(pdb_ids)}...")
            pdb_text = download_pdb(pdb_id)
            if pdb_text:
                extracted = extract_atp_from_pdb_text(pdb_text, pdb_id)
                mols.extend(extracted)
            time.sleep(0.05)  # be nice to RCSB
        print(f"Extracted {len(mols)} ATP instances from {len(pdb_ids)} structures")

    # Filter: keep only conformers with all heavy atoms present
    _ref_ha = Chem.MolFromSmiles(
        "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)"
        "COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
    ).GetNumHeavyAtoms()
    before = len(mols)
    mols = [m for m in mols if m.GetNumHeavyAtoms() == _ref_ha]
    print(f"Heavy atom filter ({_ref_ha} expected): "
          f"{len(mols)} kept, {before - len(mols)} excluded")

    print(f"\nComputing {len(mols)}x{len(mols)} RMSD matrix...")
    rmsd_mat = compute_rmsd_matrix(mols)
    np.save(f"{args.output_dir}/rmsd_matrix.npy", rmsd_mat)
    
    print(f"Clustering at {args.cutoff} Å cutoff...")
    clusters = cluster_butina(rmsd_mat, args.cutoff)
    
    # Save results
    results = {
        "n_conformers": len(mols),
        "n_clusters": len(clusters),
        "cutoff_A": args.cutoff,
        "clusters": [
            {
                "rank": i+1,
                "size": len(cl),
                "pct": round(len(cl)/len(mols)*100, 1),
                "centroid_idx": cl[0],
                "members": list(cl),
                "centroid_name": mols[cl[0]].GetProp("_Name") if mols[cl[0]].HasProp("_Name") else str(cl[0]),
            }
            for i, cl in enumerate(clusters[:20])
        ]
    }
    with open(f"{args.output_dir}/clusters.json", "w") as f:
        json.dump(results, f, indent=2)
    
    # Write cluster centroids as SDF
    with Chem.SDWriter(f"{args.output_dir}/cluster_centroids.sdf") as w:
        w.SetKekulize(False)
        for i, cl in enumerate(clusters[:10]):
            mol = mols[cl[0]]
            mol.SetProp("Cluster_Rank", str(i+1))
            mol.SetProp("Cluster_Size", str(len(cl)))
            w.write(mol)
    
    print(f"\nDone! Results in {args.output_dir}/")
    print(f"  - {len(clusters)} clusters found")
    print(f"  - Largest: {len(clusters[0])} conformers ({len(clusters[0])/len(mols)*100:.1f}%)")

if __name__ == "__main__":
    main()
