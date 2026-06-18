import gzip, re, os, sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "pipeline"))
import config

CMDS_FILE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(config.CMDS_DIR, "mpnn.cmds")

def get_motif_residues(cifgz_path):
    """Return sorted list of 'ChainResnum' strings for motif residues in the CIF."""
    seen = set()
    order = []
    with gzip.open(cifgz_path, 'rt') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            # Columns (0-indexed): 0=group_PDB 1=type_symbol 2=label_atom_id 3=label_alt_id
            # 4=label_comp_id 5=label_asym_id 6=label_entity_id 7=label_seq_id
            # 8=pdbx_PDB_ins_code 9=auth_seq_id 10=auth_comp_id 11=auth_asym_id
            # 12=auth_atom_id 13=id 14=B_iso 15=occupancy 16=formal_charge 17=is_motif
            is_motif = parts[17]
            if is_motif != 'True':
                continue
            auth_chain = parts[11]
            auth_res = parts[9]
            key = (auth_chain, int(auth_res))
            if key not in seen:
                seen.add(key)
                order.append(key)
    # Sort by chain then residue number
    order.sort()
    return " ".join(f"{c}{r}" for c, r in order)

with open(CMDS_FILE) as f:
    lines = [l.rstrip() for l in f if l.strip()]

new_lines = []
for line in lines:
    pdb_match = re.search(r'--pdb_path "([^"]+)"', line)
    if not pdb_match:
        new_lines.append(line)
        continue

    cif_path = pdb_match.group(1)
    if not os.path.exists(cif_path):
        print(f"WARNING: missing CIF: {os.path.basename(cif_path)}")
        new_lines.append(line)
        continue

    fixed_res = get_motif_residues(cif_path)
    new_line = re.sub(r'--fixed_residues "[^"]*"', f'--fixed_residues "{fixed_res}"', line)
    new_lines.append(new_line)
    print(f"{os.path.basename(cif_path)}: {fixed_res}")

with open(CMDS_FILE, "w") as f:
    f.write("\n".join(new_lines) + "\n")

print(f"\nDone. Updated {len(new_lines)} commands in mpnn.cmds")
