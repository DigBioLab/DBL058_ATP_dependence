#!/usr/bin/env python3
"""
Prepare OKT3 Fv inputs for the ATP-gated anti-CD3 campaign.

Fetches PDB 1SY6 (OKT3 Fab + CD3 gamma/epsilon heterodimer), strips the constant
domains, renames chains for the RFD3/LigandMPNN pipeline, and emits the residue-set
annotation (paratope / VH-VL interface core / redesignable) that the campaign config
and `fill_fixed_residues.py` both read.

Chain mapping
    1SY6 H 1-113   (VH)   -> A
    1SY6 L 1-107   (VL)   -> B
    1SY6 A         (CD3)  -> C      (renamed to avoid colliding with VH)

Outputs (default --out-dir inputs/, which is gitignored):
    okt3_fv.cif                 VH + VL only, apo design input
    okt3_fv_cd3.cif             VH + VL + CD3, for the paratope-check RF3 state
    okt3_1sy6_residue_sets.json residue annotation + per-residue distances

Stdlib only, so it runs on the cluster without the conda env.

Usage
    python3 scripts/prepare_okt3_fv.py --out-dir inputs/
    python3 scripts/prepare_okt3_fv.py --cif /path/to/1sy6.cif --out-dir inputs/
"""
import argparse
import json
import math
import os
import sys
import urllib.request
from collections import defaultdict

PDB_ID = "1SY6"
RCSB_URL = f"https://files.rcsb.org/download/{PDB_ID}.cif"

# Variable-domain boundaries in 1SY6 auth numbering, and target chain IDs.
VH_LAST, VL_LAST = 113, 107
CHAIN_MAP = {"H": ("A", 1, VH_LAST), "L": ("B", 1, VL_LAST), "A": ("C", None, None)}

# Distance thresholds (Angstrom)
PARATOPE_CUT = 4.5      # VH/VL residue -> CD3
PARATOPE_WIDE = 5.0     # wider set, used for MPNN fixed_residues
INTERFACE_CUT = 4.5     # VH residue -> VL and vice versa
REDESIGN_MIN_D_CD3 = 12.0   # interface residue is safe to rebuild beyond this
BUFFER_MIN_D_CD3 = 8.0      # 8-12 A: touch only with justification

CIF_COLUMN_ORDER = [
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
    "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code",
    "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
    "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num",
]


def fetch_cif(path):
    """Download 1SY6 unless already present."""
    if os.path.exists(path):
        print(f"using cached {path}")
        return path
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    print(f"downloading {RCSB_URL}")
    with urllib.request.urlopen(RCSB_URL, timeout=60) as r, open(path, "wb") as f:
        f.write(r.read())
    return path


def parse_atom_site(path):
    """Return (column_names, list_of_atom_dicts) from the _atom_site loop."""
    cols, rows = [], []
    in_loop = reading = False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop and not reading:
                if s.startswith(("ATOM", "HETATM")):
                    reading = True
                else:
                    continue
            if reading:
                if s.startswith("#") or not s:
                    break
                parts = s.split()
                if len(parts) == len(cols):
                    rows.append(dict(zip(cols, parts)))
    if not rows:
        sys.exit(f"no _atom_site records parsed from {path}")
    return cols, rows


def xyz(a):
    return float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"])


def min_dist(group_a, group_b):
    """Minimum interatomic distance between two atom lists."""
    best = float("inf")
    for a in group_a:
        ax, ay, az = xyz(a)
        for b in group_b:
            bx, by, bz = xyz(b)
            d = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
            if d < best:
                best = d
    return math.sqrt(best)


def remap(rows):
    """Keep VH/VL/CD3, rename chains, drop constant domains."""
    out = []
    for d in rows:
        if d.get("group_PDB") != "ATOM":
            continue
        src = d["auth_asym_id"]
        if src not in CHAIN_MAP:
            continue
        new, lo, hi = CHAIN_MAP[src]
        n = int(d["auth_seq_id"])
        if lo is not None and not (lo <= n <= hi):
            continue
        d = dict(d)
        d["auth_asym_id"] = d["label_asym_id"] = new
        out.append(d)
    return out


def write_cif(path, atoms, cols, data_name):
    order = [c for c in CIF_COLUMN_ORDER if c in cols]
    with open(path, "w") as f:
        f.write(f"data_{data_name}\n#\nloop_\n")
        for c in order:
            f.write(f"_atom_site.{c}\n")
        for i, d in enumerate(atoms, 1):
            d = dict(d)
            d["id"] = str(i)
            f.write(" ".join(d[c] for c in order) + "\n")
        f.write("#\n")
    print(f"wrote {path}  ({len(atoms)} atoms)")


def residue_sets(atoms):
    """Annotate VH/VL residues by distance to CD3 and to the partner domain."""
    by_chain = defaultdict(list)
    for a in atoms:
        by_chain[a["auth_asym_id"]].append(a)
    cd3 = by_chain.get("C", [])
    if not cd3:
        sys.exit("CD3 chain missing; cannot annotate paratope")

    res = defaultdict(list)
    for a in atoms:
        if a["auth_asym_id"] in ("A", "B"):
            res[(a["auth_asym_id"], int(a["auth_seq_id"]), a["auth_comp_id"])].append(a)

    per_res = {}
    for key, group in res.items():
        chain, num, aa = key
        partner = by_chain["B" if chain == "A" else "A"]
        label = ("H" if chain == "A" else "L") + str(num)
        per_res[label] = {
            "aa": aa,
            "design_chain": chain,
            "d_cd3": round(min_dist(group, cd3), 2),
            "d_partner": round(min_dist(group, partner), 2),
        }

    def srt(keys):
        return sorted(keys, key=lambda s: (s[0], int(s[1:])))

    paratope = srt([k for k, v in per_res.items() if v["d_cd3"] <= PARATOPE_CUT])
    paratope_wide = srt([k for k, v in per_res.items() if v["d_cd3"] <= PARATOPE_WIDE])
    core = srt([k for k, v in per_res.items() if v["d_partner"] <= INTERFACE_CUT])
    redesign = srt([k for k in core if per_res[k]["d_cd3"] >= REDESIGN_MIN_D_CD3])
    buffer_ = srt([k for k in core
                   if BUFFER_MIN_D_CD3 <= per_res[k]["d_cd3"] < REDESIGN_MIN_D_CD3])

    return {
        "source": f"PDB {PDB_ID} (OKT3 Fab + CD3 gamma/epsilon), 1SY6 auth numbering",
        "chain_map": {"VH": "H 1-113 -> A", "VL": "L 1-107 -> B", "CD3ge": "A -> C"},
        "thresholds": {
            "paratope": PARATOPE_CUT, "paratope_wide": PARATOPE_WIDE,
            "interface": INTERFACE_CUT, "redesign_min_d_cd3": REDESIGN_MIN_D_CD3,
        },
        "PARATOPE": paratope,
        "PARATOPE_WIDE": paratope_wide,
        "INTERFACE_CORE": core,
        "REDESIGNABLE": redesign,
        "BUFFER": buffer_,
        "per_residue": per_res,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="inputs",
                    help="output directory (default: inputs/, gitignored)")
    ap.add_argument("--cif", default=None,
                    help="existing 1SY6 mmCIF; downloaded to <out-dir>/1sy6.cif if omitted")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    cif = args.cif or fetch_cif(os.path.join(args.out_dir, "1sy6.cif"))

    cols, rows = parse_atom_site(cif)
    atoms = remap(rows)
    if not atoms:
        sys.exit("chain remapping produced no atoms; check VH/VL boundaries")

    fv_cd3 = atoms
    fv_only = [a for a in atoms if a["auth_asym_id"] != "C"]

    write_cif(os.path.join(args.out_dir, "okt3_fv_cd3.cif"), fv_cd3, cols, "OKT3_FV_CD3")
    write_cif(os.path.join(args.out_dir, "okt3_fv.cif"), fv_only, cols, "OKT3_FV_APO")

    sets = residue_sets(fv_cd3)
    json_path = os.path.join(args.out_dir, "okt3_1sy6_residue_sets.json")
    with open(json_path, "w") as f:
        json.dump(sets, f, indent=1)
    print(f"wrote {json_path}")

    for name in ("PARATOPE", "REDESIGNABLE", "BUFFER"):
        print(f"{name:<14}: {' '.join(sets[name])}")


if __name__ == "__main__":
    main()
