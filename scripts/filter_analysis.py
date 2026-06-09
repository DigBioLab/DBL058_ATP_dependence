#!/usr/bin/env python3
"""
Filter ATP conformers by heavy atom completeness, align to a reference conformer,
plot per-cluster exclusions, and regenerate filtered+aligned cluster representative SDFs.

Usage:
    python scripts/filter_analysis.py [--output_dir output] [--plot_dir pictures] [--n_reps 10] [--seed 42]
"""

import json, argparse, random
from collections import defaultdict, Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import requests
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdFMCS
from rdkit.ML.Cluster import Butina

ATP_SMILES = "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
EXPECTED_HA = Chem.MolFromSmiles(ATP_SMILES).GetNumHeavyAtoms()
SMALL_THRESHOLD = 10


def get_atp_pdb_ids(max_results=9999):
    url = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb/ATP"
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
            chain, resnum = l[21], l[22:26].strip()
            groups[f"{chain}_{resnum}"].append(l)
    mols = []
    for key, lines in groups.items():
        mol = Chem.MolFromPDBBlock("\n".join(lines + ["END"]),
                                   sanitize=False, removeHs=False)
        if mol and mol.GetNumAtoms() > 5:
            try:
                Chem.SanitizeMol(mol)
                mol.SetProp("_Name", f"{pdb_id}_{key}")
                mols.append(mol)
            except Exception:
                pass
    return mols


def build_alignment_pattern():
    """Build reference ATP conformer (heavy atoms only) for alignment."""
    ref_mol = Chem.MolFromSmiles(ATP_SMILES)
    ref_mol = Chem.AddHs(ref_mol)
    AllChem.EmbedMolecule(ref_mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(ref_mol)
    return Chem.RemoveHs(ref_mol)  # align on heavy atoms only


def align_mol(mol, ref_mol, ref_pattern, ref_match):
    probe_match = mol.GetSubstructMatch(ref_pattern)
    if not probe_match or not ref_match:
        return False
    try:
        rdMolAlign.AlignMol(mol, ref_mol, atomMap=list(zip(probe_match, ref_match)))
        return True
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", default="output")
    parser.add_argument("--plot_dir",   default="pictures")
    parser.add_argument("--n_reps", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    d = args.output_dir
    Path(args.plot_dir).mkdir(parents=True, exist_ok=True)
    random.seed(args.seed)

    print(f"Expected ATP heavy atoms: {EXPECTED_HA}")

    # ── Re-cluster from saved RMSD matrix ─────────────────────────────────────
    with open(f"{d}/clusters.json") as f:
        meta = json.load(f)
    cutoff = meta["cutoff_A"]

    print("Loading RMSD matrix and re-clustering...")
    rmsd_mat = np.load(f"{d}/rmsd_matrix.npy")
    N = rmsd_mat.shape[0]
    dist_list = [float(rmsd_mat[i, j]) for i in range(1, N) for j in range(i)]
    clusters = Butina.ClusterData(dist_list, N, cutoff, isDistData=True, reordering=True)
    clusters = sorted(clusters, key=len, reverse=True)
    print(f"  {len(clusters)} clusters, {N} conformers")

    # ── Re-extract molecules from local cache ─────────────────────────────────
    print("\nRe-extracting molecules from PDB cache (no downloads)...")
    pdb_ids = get_atp_pdb_ids(9999)
    cache_dir = Path("pdb_cache")
    mols = []
    for i, pdb_id in enumerate(pdb_ids):
        if i % 500 == 0:
            print(f"  {i}/{len(pdb_ids)} processed, {len(mols)} molecules...")
        local = cache_dir / f"{pdb_id.lower()}.pdb"
        if local.exists():
            mols.extend(extract_atp_from_pdb_text(local.read_text(), pdb_id))
    print(f"  {len(mols)} molecules extracted")

    # Verify ordering
    mismatches = 0
    for cm in meta["clusters"][:10]:
        idx, expected = cm["centroid_idx"], cm["centroid_name"]
        actual = mols[idx].GetProp("_Name") if mols[idx].HasProp("_Name") else ""
        if actual != expected:
            mismatches += 1
            print(f"  [!] rank {cm['rank']}: expected '{expected}', got '{actual}'")
    if mismatches:
        print(f"  {mismatches} ordering mismatches — aborting.")
        return
    print("  Molecule order OK")

    # ── Build alignment reference ─────────────────────────────────────────────
    print("\nBuilding alignment reference and MCS pattern...")
    ref_mol = build_alignment_pattern()
    # Use first complete molecule to compute MCS (more representative than SMILES ref alone)
    first_complete = next(m for m in mols if m.GetNumHeavyAtoms() == EXPECTED_HA)
    mcs = rdFMCS.FindMCS(
        [first_complete, ref_mol],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=10,
    )
    ref_pattern = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref_mol.GetSubstructMatch(ref_pattern)
    print(f"  MCS: {mcs.numAtoms} atoms, {mcs.numBonds} bonds")

    # ── Heavy atom filter ─────────────────────────────────────────────────────
    ha_counts = np.array([m.GetNumHeavyAtoms() for m in mols])
    passes = ha_counts == EXPECTED_HA
    n_excl = int((~passes).sum())

    print(f"\nHeavy atom filter (expected = {EXPECTED_HA}):")
    print(f"  Kept:     {passes.sum():,}")
    print(f"  Excluded: {n_excl:,} ({n_excl / N * 100:.1f}%)")
    print(f"  HA count distribution of excluded: {dict(sorted(Counter(ha_counts[~passes]).items()))}")

    # ── Per-cluster stats ─────────────────────────────────────────────────────
    stats = []
    for ci, cl in enumerate(clusters):
        excl = sum(1 for idx in cl if not passes[idx])
        stats.append({"rank": ci + 1, "total": len(cl),
                      "excluded": excl, "kept": len(cl) - excl})

    # ── Plot: stacked bar of kept vs excluded per cluster ─────────────────────
    cmap = plt.get_cmap("tab10")

    big  = [s for s in stats if s["total"] >= SMALL_THRESHOLD]
    small_agg = {
        "rank": "small",
        "total":    sum(s["total"]    for s in stats if s["total"] < SMALL_THRESHOLD),
        "excluded": sum(s["excluded"] for s in stats if s["total"] < SMALL_THRESHOLD),
        "kept":     sum(s["kept"]     for s in stats if s["total"] < SMALL_THRESHOLD),
    }
    plot_data = big + [small_agg]
    xlabels   = [f"C{s['rank']}" for s in big] + ["Small\nclusters"]
    kept_vals = [s["kept"]     for s in plot_data]
    excl_vals = [s["excluded"] for s in plot_data]

    fig, ax = plt.subplots(figsize=(max(8, len(xlabels) * 0.85), 5))
    fig.patch.set_facecolor("#f9f9f9")
    ax.set_facecolor("#f9f9f9")

    x = list(range(len(xlabels)))
    ax.bar(x, kept_vals, label="Complete", color="#4878d0",
           edgecolor="white", linewidth=0.8, width=0.65)
    ax.bar(x, excl_vals, bottom=kept_vals, label=f"Excluded (≠{EXPECTED_HA} heavy atoms)",
           color="#ee854a", edgecolor="white", linewidth=0.8, width=0.65)

    for xi, (kv, ev) in enumerate(zip(kept_vals, excl_vals)):
        if ev > 0:
            ax.text(xi, kv + ev * 1.05, str(ev),
                    ha="center", va="bottom", fontsize=8,
                    color="#cc5500", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.set_ylabel("Number of conformers")
    ax.set_title(
        f"Incomplete ATP conformers per cluster\n"
        f"{n_excl:,} of {N:,} excluded (expected {EXPECTED_HA} heavy atoms)  ·  cutoff {cutoff} Å",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    plot_path = f"{args.plot_dir}/cluster_filter_exclusions.png"
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved {plot_path}")

    # ── Write filtered + aligned representative SDFs ──────────────────────────
    out_dir = Path(d) / "cluster_representatives_filtered"
    out_dir.mkdir(exist_ok=True)
    print(f"\nWriting filtered + aligned SDFs to {out_dir}/")

    total_aligned, total_failed = 0, 0

    for s in stats:
        ci = s["rank"] - 1
        cl = clusters[ci]
        valid = [idx for idx in cl if passes[idx]]

        if not valid:
            print(f"  cluster_{s['rank']:02d}: SKIPPED — no complete members")
            continue

        n_sample = min(args.n_reps, len(valid))
        centroid = cl[0] if passes[cl[0]] else valid[0]
        rest = [idx for idx in valid if idx != centroid]
        random.shuffle(rest)
        sampled = [centroid] + rest[:n_sample - 1]

        n_ok, n_fail = 0, 0
        sdf_path = out_dir / f"cluster_{s['rank']:02d}_representatives.sdf"
        with Chem.SDWriter(str(sdf_path)) as w:
            w.SetKekulize(False)
            for pos, idx in enumerate(sampled):
                mol = mols[idx]
                aligned = align_mol(mol, ref_mol, ref_pattern, ref_match)
                if aligned:
                    n_ok += 1
                else:
                    n_fail += 1
                mol.SetProp("Cluster_Rank",          str(s["rank"]))
                mol.SetProp("Cluster_Size",          str(s["total"]))
                mol.SetProp("Cluster_Size_Filtered", str(len(valid)))
                mol.SetProp("Is_Centroid",           "1" if pos == 0 else "0")
                mol.SetProp("Aligned",               "1" if aligned else "0")
                mol.SetProp("Conformer_Index",       str(idx))
                w.write(mol)

        total_aligned += n_ok
        total_failed  += n_fail
        excl_note  = f", {s['excluded']} incomplete" if s["excluded"] else ""
        align_note = f", {n_fail} unaligned" if n_fail else ""
        print(f"  cluster_{s['rank']:02d}: {n_sample} reps "
              f"({len(valid)} complete{excl_note}{align_note})")

    print(f"\nAlignment summary: {total_aligned} aligned, {total_failed} unaligned")
    print(f"Done. Filtered + aligned SDFs in {out_dir}/")


if __name__ == "__main__":
    main()
