"""
Run PLIP on all cached PDB structures and parse ATP–protein interactions.
Results (XML per structure) are cached — re-runs skip completed structures.

Saves:
  output/binding/plip_results/<PDB_ID>/*_report.xml
  output/binding/atp_interactions.csv

Run via submit/03_plip.sh
"""

import os, sys, json, shutil, subprocess, argparse
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from tqdm import tqdm

PDB_CACHE_DIR = Path('pdb_cache')
OUTPUT_DIR    = Path('output/binding')
PLIP_DIR      = OUTPUT_DIR / 'plip_results'
CSV_PATH      = OUTPUT_DIR / 'atp_interactions.csv'

AA3 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

ITYPE_LABEL = {
    'hydrogen_bonds':          'HBond',
    'hydrophobic_interactions': 'HPhob',
    'salt_bridges':            'Salt',
    'pi_stacks':               'PiStack',
    'pi_cation_interactions':  'PiCat',
    'metal_complexes':         'Metal',
    'water_bridges':           'WaterBr',
    'halogen_bonds':           'Halogen',
}

ATP_MOIETY = {
    'PG':  'gamma-PO4', 'O1G': 'gamma-PO4', 'O2G': 'gamma-PO4', 'O3G': 'gamma-PO4',
    'O3B': 'gamma-PO4',
    'PB':  'beta-PO4',  'O1B': 'beta-PO4',  'O2B': 'beta-PO4',
    'O3A': 'beta-PO4',
    'PA':  'alpha-PO4', 'O1A': 'alpha-PO4', 'O2A': 'alpha-PO4',
    "C1'": 'ribose', "C2'": 'ribose', "C3'": 'ribose', "C4'": 'ribose',
    "C5'": 'ribose', "O4'": 'ribose', "O2'": 'ribose', "O3'": 'ribose',
    "O5'": 'ribose',
    'N1': 'adenine', 'C2': 'adenine', 'N3': 'adenine', 'C4': 'adenine',
    'C5': 'adenine', 'C6': 'adenine', 'N6': 'adenine', 'N7': 'adenine',
    'C8': 'adenine', 'N9': 'adenine',
}


def atp_moiety(atom_name: str) -> str:
    a = atom_name.strip().upper().replace('*', "'")
    return ATP_MOIETY.get(a, 'other')


def _extract_atp_atom(item, igroup_tag: str) -> str:
    atom = item.findtext('ligatom', '').strip()
    if atom:
        return atom
    if igroup_tag == 'salt_bridges':
        lig_atoms = item.find('lig_atoms')
        if lig_atoms is not None:
            first = lig_atoms.find('atom')
            if first is not None and first.text:
                return first.text.strip()
        return item.findtext('lig_group', '').strip()
    if igroup_tag in ('pi_stacks', 'pi_cation_interactions'):
        return item.findtext('lig_group', 'ring').strip() or 'ring'
    return ''


def run_plip(plip_exe, pdb_path: str, plip_dir: str):
    Path(plip_dir).mkdir(parents=True, exist_ok=True)
    existing = list(Path(plip_dir).glob('*_report.xml'))
    if existing:
        return str(existing[0])
    try:
        subprocess.run(
            [plip_exe, '-f', str(pdb_path), '-x', '-o', str(plip_dir)],
            capture_output=True, text=True, timeout=120, check=True,
        )
    except Exception:
        return None
    found = list(Path(plip_dir).glob('*_report.xml'))
    return str(found[0]) if found else None


def parse_atp_xml(xml_path: str, pdb_id: str) -> pd.DataFrame:
    try:
        root = ET.parse(xml_path).getroot()
    except Exception:
        return pd.DataFrame()
    records = []
    for bs in root.findall('.//bindingsite'):
        if bs.findtext('./identifiers/hetid', '') != 'ATP':
            continue
        atp_id = (f"{bs.findtext('./identifiers/chain','')}"
                  f":{bs.findtext('./identifiers/position','')}")
        interacts = bs.find('interactions')
        if interacts is None:
            continue
        for igroup in interacts:
            raw   = igroup.tag
            label = ITYPE_LABEL.get(raw, raw)
            for item in igroup:
                restype  = item.findtext('restype',  '')
                resnr    = item.findtext('resnr',    '')
                reschain = item.findtext('reschain', '')
                if raw == 'hydrogen_bonds':
                    pdon  = item.findtext('protisdon', 'True')
                    label = 'HBond_Don' if pdon == 'True' else 'HBond_Acc'
                atom_name = _extract_atp_atom(item, raw)
                moiety    = atp_moiety(atom_name)
                records.append({
                    'pdb_id':    pdb_id,
                    'atp_id':    atp_id,
                    'res_nr':    resnr,
                    'res_type':  restype,
                    'res_1':     AA3.get(restype, 'X'),
                    'itype':     label,
                    'itype_raw': raw,
                    'atp_atom':  atom_name,
                    'atp_moiety': moiety,
                })
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n-workers', type=int, default=16)
    parser.add_argument('--n-structures', type=int, default=None,
                        help='Cap number of structures (default: all)')
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PLIP_DIR.mkdir(parents=True, exist_ok=True)

    plip_exe = shutil.which('plip') or os.path.join(os.path.dirname(sys.executable), 'plip')
    if not Path(plip_exe).exists():
        raise EnvironmentError(f'PLIP not found at {plip_exe!r}')
    print(f'PLIP: {plip_exe}')

    pdb_files = sorted(PDB_CACHE_DIR.glob('*.pdb'))
    if args.n_structures:
        pdb_files = pdb_files[:args.n_structures]
    pdb_ids = [f.stem.upper() for f in pdb_files]
    print(f'Structures to process: {len(pdb_ids)}')

    n_cached = sum(1 for pid in pdb_ids
                   if list((PLIP_DIR / pid).glob('*_report.xml')))
    print(f'Already cached: {n_cached}  |  To run: {len(pdb_ids)-n_cached}')

    def _worker(pdb_id):
        pdb_path = PDB_CACHE_DIR / f'{pdb_id.lower()}.pdb'
        if not pdb_path.exists():
            return pdb_id, None
        xml_path = run_plip(plip_exe, str(pdb_path), str(PLIP_DIR / pdb_id))
        if xml_path is None:
            return pdb_id, None
        df = parse_atp_xml(xml_path, pdb_id)
        return pdb_id, df if not df.empty else None

    results, failed = {}, []
    with ThreadPoolExecutor(max_workers=args.n_workers) as pool:
        futures = {pool.submit(_worker, pid): pid for pid in pdb_ids}
        for fut in tqdm(as_completed(futures), total=len(pdb_ids), desc='PLIP'):
            pid, df = fut.result()
            if df is not None:
                results[pid] = df
            else:
                failed.append(pid)

    print(f'\nSucceeded: {len(results):,}  |  Failed/no ATP: {len(failed):,}')

    if results:
        df_all = pd.concat(results.values(), ignore_index=True)
        df_all.to_csv(CSV_PATH, index=False)
        print(f'Saved {len(df_all):,} interaction records to {CSV_PATH}')
    else:
        print('No results to save.')


if __name__ == '__main__':
    main()
