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


# ── ATP seeding (Route B) ────────────────────────────────────────────────────
# RFD3 extracts a ligand from the *input structure*; it does not fetch it from
# the CCD (rfd3/utils/inference.py, fetch_mask_from_component). Partial
# diffusion therefore needs ATP already placed in the input, even though
# section 1.3 of the plan measured that no pocket exists yet to place it in.
#
# The seed is deliberately crude: ATP is parked over the shallow groove section
# 1.3 identified, just clear of the surface, and partial diffusion is expected
# to excavate around it. It is a starting point, not a docked pose.
CCD_ATP_URL = "https://files.rcsb.org/ligands/download/ATP.cif"

# Groove lining from plan section 1.3, in design-chain numbering (H->A, L->B).
GROOVE_RESIDUES = ["A39", "A43", "A93", "A95", "B37", "B38", "B41", "B42", "B86"]

ATP_CLASH_MIN = 3.2     # Angstrom; closest allowed ATP-protein heavy-atom contact
ATP_PROBE_MAX = 14.0    # how far out along the groove normal to search
ATP_N_ROT = 60          # random orientations tried per offset


def fetch_atp_template(path):
    """Ideal-geometry ATP heavy atoms from the RCSB chemical component dictionary."""
    if not os.path.exists(path):
        print(f"downloading {CCD_ATP_URL}")
        with urllib.request.urlopen(CCD_ATP_URL, timeout=60) as r, open(path, "wb") as f:
            f.write(r.read())

    cols, rows, in_loop, reading = [], [], False, False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_chem_comp_atom."):
                cols.append(s.split(".", 1)[1]); in_loop = True; continue
            if in_loop:
                if s.startswith("#") or not s:
                    if reading:
                        break
                    continue
                parts = s.split()
                if len(parts) == len(cols):
                    rows.append(dict(zip(cols, parts))); reading = True

    def pick(d, *names):
        for n in names:
            v = d.get(n)
            if v not in (None, "?", "."):
                return v
        return None

    atoms = []
    for d in rows:
        if d.get("type_symbol") == "H":
            continue
        x = pick(d, "model_Cartn_x", "pdbx_model_Cartn_x_ideal")
        y = pick(d, "model_Cartn_y", "pdbx_model_Cartn_y_ideal")
        z = pick(d, "model_Cartn_z", "pdbx_model_Cartn_z_ideal")
        if None in (x, y, z):
            continue
        atoms.append({"name": d["atom_id"], "element": d["type_symbol"],
                      "xyz": (float(x), float(y), float(z))})
    if len(atoms) != 31:
        print(f"WARNING: ATP template has {len(atoms)} heavy atoms, expected 31")
    return atoms


def _centroid(points):
    n = len(points)
    return tuple(sum(p[i] for p in points) / n for i in range(3))


def _rotation(a, b, c):
    """Rotation matrix from three Euler angles."""
    ca, sa, cb, sb, cc, sc = (math.cos(a), math.sin(a), math.cos(b),
                              math.sin(b), math.cos(c), math.sin(c))
    return (
        (cb * cc, -cb * sc, sb),
        (sa * sb * cc + ca * sc, -sa * sb * sc + ca * cc, -sa * cb),
        (-ca * sb * cc + sa * sc, ca * sb * sc + sa * cc, ca * cb),
    )


def place_atp(fv_atoms, template, seed=0, clash_min=None, force_offset=None):
    """Park ATP over the groove: outward along the groove normal, closest
    non-clashing offset, best of ATP_N_ROT orientations by contact count.

    clash_min   override ATP_CLASH_MIN. Lowering it seeds ATP deeper, at the
                cost of starting from a clash that diffusion must resolve.
    force_offset  pin the offset along the groove normal instead of scanning
                for the first clear one; use with a low clash_min to seed
                inside the groove.

    Returns (placed_atoms, report_dict).
    """
    clash_min = ATP_CLASH_MIN if clash_min is None else clash_min
    import random
    rng = random.Random(seed)

    prot = [xyz(a) for a in fv_atoms]
    groove = [xyz(a) for a in fv_atoms
              if a["auth_asym_id"] + a["auth_seq_id"] in GROOVE_RESIDUES]
    if not groove:
        sys.exit("no groove-lining atoms matched; check GROOVE_RESIDUES")

    G, P = _centroid(groove), _centroid(prot)
    u = tuple(G[i] - P[i] for i in range(3))
    norm = math.sqrt(sum(v * v for v in u)) or 1.0
    u = tuple(v / norm for v in u)

    tpl = [a["xyz"] for a in template]
    tpl_c = _centroid(tpl)
    tpl = [tuple(p[i] - tpl_c[i] for i in range(3)) for p in tpl]

    rots = [_rotation(rng.uniform(0, 2 * math.pi),
                      rng.uniform(0, 2 * math.pi),
                      rng.uniform(0, 2 * math.pi)) for _ in range(ATP_N_ROT)]

    best = None
    fallback = None                    # best pose ignoring the clash floor
    offset = 0.0 if force_offset is None else force_offset
    while offset <= ATP_PROBE_MAX:
        centre = tuple(G[i] + offset * u[i] for i in range(3))
        for R in rots:
            pose = [tuple(centre[i] + sum(R[i][j] * p[j] for j in range(3))
                          for i in range(3)) for p in tpl]
            dmin, contacts = 1e9, 0
            for q in pose:
                for r in prot:
                    d2 = ((q[0]-r[0])**2 + (q[1]-r[1])**2 + (q[2]-r[2])**2)
                    if d2 < dmin:
                        dmin = d2
                    if d2 < 25.0:      # 5 A contact shell
                        contacts += 1
            dmin = math.sqrt(dmin)
            if fallback is None or dmin > fallback[1]:
                fallback = (contacts, dmin, offset, pose)
            if dmin >= clash_min and (best is None or contacts > best[0]):
                best = (contacts, dmin, offset, pose)
        if best is not None or force_offset is not None:
            break                      # first offset that clears; most buried
        offset += 0.5

    if best is None:
        if force_offset is None:
            sys.exit(f"could not place ATP without a clash below {ATP_PROBE_MAX} A")
        # Pinned offset with no clash-free orientation: take the least-clashing
        # pose and let diffusion resolve it. This is the point of a deep seed.
        best = fallback
        print(f"  note: no orientation clears {clash_min} A at offset "
              f"{force_offset}; seeding with a deliberate clash")

    contacts, dmin, offset, pose = best
    placed = [{"name": t["name"], "element": t["element"], "xyz": p}
              for t, p in zip(template, pose)]
    return placed, {"groove_centroid": [round(v, 2) for v in G],
                    "offset_along_normal": round(offset, 2),
                    "min_protein_distance": round(dmin, 2),
                    "contacts_within_5A": contacts}


def atp_cif_rows(placed, cols, chain="L", entity="9"):
    """ATP heavy atoms as _atom_site dicts matching the Fv CIF columns."""
    rows = []
    for i, a in enumerate(placed, 1):
        x, y, z = a["xyz"]
        d = {c: "?" for c in cols}
        d.update({
            "group_PDB": "HETATM", "id": str(i), "type_symbol": a["element"],
            "label_atom_id": a["name"], "label_alt_id": ".", "label_comp_id": "ATP",
            "label_asym_id": chain, "label_entity_id": entity, "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?",
            "Cartn_x": f"{x:.3f}", "Cartn_y": f"{y:.3f}", "Cartn_z": f"{z:.3f}",
            "occupancy": "1.00", "B_iso_or_equiv": "30.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "ATP",
            "auth_asym_id": chain, "auth_atom_id": a["name"],
            "pdbx_PDB_model_num": "1",
        })
        rows.append(d)
    return rows


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
    ap.add_argument("--with-atp", action="store_true",
                    help="also write okt3_fv_atp.cif with ATP seeded over the "
                         "section-1.3 groove (required for Route B partial diffusion)")
    ap.add_argument("--atp-offset", type=float, default=None,
                    help="pin the ATP offset along the groove normal (A) instead "
                         "of taking the first clash-free one; use with "
                         "--atp-clash-min to seed deep inside the groove")
    ap.add_argument("--atp-clash-min", type=float, default=None,
                    help=f"minimum ATP-protein contact (default {ATP_CLASH_MIN} A); "
                         "lower it to allow a deliberately clashing deep seed")
    ap.add_argument("--atp-out", default="okt3_fv_atp.cif",
                    help="filename for the ATP-seeded structure")
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

    if args.with_atp:
        template = fetch_atp_template(os.path.join(args.out_dir, "ATP_ccd.cif"))
        placed, report = place_atp(fv_only, template,
                                   clash_min=args.atp_clash_min,
                                   force_offset=args.atp_offset)
        write_cif(os.path.join(args.out_dir, args.atp_out),
                  fv_only + atp_cif_rows(placed, cols), cols, "OKT3_FV_ATP")
        print("ATP seed placement:")
        for k, v in report.items():
            print(f"  {k:<22} {v}")
        stem = os.path.splitext(args.atp_out)[0]
        with open(os.path.join(args.out_dir, f"{stem}_placement.json"), "w") as f:
            json.dump(report, f, indent=1)

    sets = residue_sets(fv_cd3)
    json_path = os.path.join(args.out_dir, "okt3_1sy6_residue_sets.json")
    with open(json_path, "w") as f:
        json.dump(sets, f, indent=1)
    print(f"wrote {json_path}")

    for name in ("PARATOPE", "REDESIGNABLE", "BUFFER"):
        print(f"{name:<14}: {' '.join(sets[name])}")


if __name__ == "__main__":
    main()
