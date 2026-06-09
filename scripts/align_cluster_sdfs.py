#!/usr/bin/env python3
"""
Align all cluster representative SDF files to a reference ATP conformer.

Usage:
    python scripts/align_cluster_sdfs.py [--output_dir output]
"""

import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdFMCS


def build_reference_atp():
    atp_smiles = "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
    mol = Chem.MolFromSmiles(atp_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def align_mol(mol, ref_mol, ref_pattern, ref_match):
    """Align mol onto ref_mol using a precomputed MCS pattern."""
    probe_match = mol.GetSubstructMatch(ref_pattern)
    if not probe_match or not ref_match:
        return mol, False
    atom_map = list(zip(probe_match, ref_match))
    try:
        rdMolAlign.AlignMol(mol, ref_mol, atomMap=atom_map)
        return mol, True
    except Exception:
        return mol, False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", default="output")
    args = parser.parse_args()

    rep_dir = Path(args.output_dir) / "cluster_representatives"
    out_dir = Path(args.output_dir) / "cluster_representatives_aligned"
    out_dir.mkdir(exist_ok=True)

    sdf_files = sorted(rep_dir.glob("cluster_*_representatives.sdf"))
    if not sdf_files:
        print(f"No SDF files found in {rep_dir}")
        return

    # ── Build reference and MCS pattern ──────────────────────────────────────
    print("Building reference ATP conformer...")
    ref_mol = build_reference_atp()

    # Load one probe molecule to compute MCS once
    first_mol = next(Chem.SDMolSupplier(str(sdf_files[0]), removeHs=False))
    print("Computing MCS between reference and probe...")
    mcs = rdFMCS.FindMCS(
        [first_mol, ref_mol],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=10,
    )
    ref_pattern = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref_mol.GetSubstructMatch(ref_pattern)
    print(f"  MCS: {mcs.numAtoms} atoms, {mcs.numBonds} bonds")

    # ── Align each SDF file ───────────────────────────────────────────────────
    for sdf_path in sdf_files:
        out_path = out_dir / sdf_path.name
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mols = [m for m in supplier if m is not None]

        n_ok, n_fail = 0, 0
        with Chem.SDWriter(str(out_path)) as w:
            w.SetKekulize(False)
            for mol in mols:
                mol, ok = align_mol(mol, ref_mol, ref_pattern, ref_match)
                if ok:
                    n_ok += 1
                else:
                    n_fail += 1
                w.write(mol)

        print(f"  {sdf_path.name}: {n_ok} aligned, {n_fail} failed")

    print(f"\nDone. Aligned SDFs in {out_dir}/")


if __name__ == "__main__":
    main()
