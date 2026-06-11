"""
Download all ATP-containing PDB structures from RCSB and extract/filter/align
ATP conformers. Saves:
  pdb_cache/<pdb_id>.pdb          — raw PDB files
  output/atp_conformers_raw.sdf   — filtered + aligned ATP molecules (31 heavy atoms)
  output/pdb_ids.json             — list of fetched PDB IDs

Run via submit/01_download.sh
"""

import os, json, time, argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdFMCS
from rdkit.ML.Cluster import Butina

# ── Config ───────────────────────────────────────────────────────────────────
PDB_CACHE   = Path('pdb_cache')
OUTPUT_DIR  = Path('output')
RAW_SDF     = OUTPUT_DIR / 'atp_conformers_raw.sdf'
IDS_JSON    = OUTPUT_DIR / 'pdb_ids.json'

MAX_STRUCTS = 9999
N_WORKERS   = 10      # download threads — keep modest to be polite to RCSB
SLEEP_S     = 0.05    # seconds between requests per worker

ATP_SMILES  = ('c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H]'
               '(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N')

# ── Helpers ──────────────────────────────────────────────────────────────────

def get_atp_pdb_ids(max_results=9999):
    url  = 'https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb/ATP'
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    atp_data = resp.json().get('ATP', [])
    ids = []
    for item in atp_data:
        if isinstance(item, str):
            ids.append(item.upper())
        elif isinstance(item, dict):
            for key in ('pdb_id', 'entry_id', 'id', 'pdbId'):
                if key in item:
                    ids.append(item[key].upper())
                    break
    return ids[:max_results]


def download_pdb(pdb_id):
    local = PDB_CACHE / f'{pdb_id.lower()}.pdb'
    if local.exists():
        return local.read_text()
    try:
        resp = requests.get(
            f'https://files.rcsb.org/download/{pdb_id}.pdb', timeout=30
        )
        if resp.status_code == 200:
            local.write_text(resp.text)
            return resp.text
    except Exception:
        pass
    return None


def extract_atp_from_pdb_text(pdb_text, pdb_id):
    atp_lines = [l for l in pdb_text.split('\n')
                 if (l.startswith('HETATM') and l[17:20].strip() == 'ATP')
                 or l.startswith('CONECT') or l.startswith('END')]
    if not atp_lines:
        return []
    groups = defaultdict(list)
    for l in atp_lines:
        if l.startswith('HETATM'):
            groups[f'{l[21]}_{l[22:26].strip()}'].append(l)
    mols = []
    for key, lines in groups.items():
        mol = Chem.MolFromPDBBlock('\n'.join(lines + ['END']),
                                   sanitize=False, removeHs=False)
        if mol and mol.GetNumAtoms() > 5:
            try:
                Chem.SanitizeMol(mol)
                mol.SetProp('_Name', f'{pdb_id}_{key}')
                mols.append(mol)
            except Exception:
                pass
    return mols


def build_reference_atp():
    mol = Chem.MolFromSmiles(ATP_SMILES)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.RemoveHs(mol)


def make_alignment_pattern(ref_mol, probe_mol):
    mcs = rdFMCS.FindMCS(
        [probe_mol, ref_mol],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=10,
    )
    pattern   = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref_mol.GetSubstructMatch(pattern)
    return pattern, ref_match


def align_mol(mol, ref_mol, pattern, ref_match):
    probe_match = mol.GetSubstructMatch(pattern)
    if not probe_match or not ref_match:
        return False
    try:
        rdMolAlign.AlignMol(mol, ref_mol, atomMap=list(zip(probe_match, ref_match)))
        return True
    except Exception:
        return False


def _download_worker(pdb_id):
    time.sleep(SLEEP_S)
    text = download_pdb(pdb_id)
    if text is None:
        return pdb_id, []
    return pdb_id, extract_atp_from_pdb_text(text, pdb_id)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    PDB_CACHE.mkdir(exist_ok=True)
    OUTPUT_DIR.mkdir(exist_ok=True)

    if RAW_SDF.exists() and IDS_JSON.exists():
        print(f'Checkpoint found at {RAW_SDF} — delete to rerun.')
        return

    print('Fetching ATP PDB IDs from PDBe...')
    pdb_ids = get_atp_pdb_ids(MAX_STRUCTS)
    print(f'  {len(pdb_ids)} structures found.')

    IDS_JSON.write_text(json.dumps(pdb_ids, indent=2))
    print(f'  IDs saved to {IDS_JSON}')

    print(f'\nDownloading & extracting ATP ({N_WORKERS} workers)...')
    mols = []
    done = 0
    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(_download_worker, pid): pid for pid in pdb_ids}
        for fut in futures:
            pid, mol_list = fut.result()
            mols.extend(mol_list)
            done += 1
            if done % 200 == 0:
                print(f'  {done}/{len(pdb_ids)} processed | {len(mols)} molecules')

    print(f'\nExtracted {len(mols):,} raw ATP conformers.')

    expected_ha = Chem.MolFromSmiles(ATP_SMILES).GetNumHeavyAtoms()
    mols_filtered = [m for m in mols if m.GetNumHeavyAtoms() == expected_ha]
    print(f'After 31-heavy-atom filter: {len(mols_filtered):,} kept '
          f'({len(mols)-len(mols_filtered):,} excluded)')

    print('Aligning to MMFF94 reference ATP...')
    ref_mol = build_reference_atp()
    pattern, ref_match = make_alignment_pattern(ref_mol, mols_filtered[0])
    print(f'  MCS: {pattern.GetNumAtoms()} atoms, {pattern.GetNumBonds()} bonds')
    n_aligned = sum(align_mol(m, ref_mol, pattern, ref_match) for m in mols_filtered)
    print(f'  Aligned: {n_aligned}/{len(mols_filtered)}')

    print(f'Saving to {RAW_SDF}...')
    with Chem.SDWriter(str(RAW_SDF)) as w:
        w.SetKekulize(False)
        for mol in mols_filtered:
            w.write(mol)
    print(f'Done. {len(mols_filtered):,} conformers written.')


if __name__ == '__main__':
    main()
