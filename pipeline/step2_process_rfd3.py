#!/usr/bin/env python3
"""
Step 2 (LSF CPU job): Process RFD3 outputs, generate MPNN cmds, submit MPNN + step 3.
Runs automatically after the RFD3 job completes via LSF dependency.
"""
import sys, os, glob, json, csv, gzip, re, subprocess
sys.path.insert(0, "/work3/s212001/ATP_dependent/pipeline")
import config

sys.path.insert(0, config.ENZDES_SCRIPTS)
import jupyter_utils

sys.path.append("/dtu/projects/dbl/proteusai/latest/src/")
from pai.metrics.structure import compute_structure_metrics
from atomworks.io import parse


def _bsub_jid(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(r.stdout.strip(), flush=True)
    if r.returncode != 0:
        raise RuntimeError(f"bsub failed:\nSTDOUT: {r.stdout}\nSTDERR: {r.stderr}")
    m = re.search(r"Job <(\d+)>", r.stdout)
    if not m:
        raise RuntimeError(f"No job ID in bsub output: {r.stdout}")
    return m.group(1)


def _motif_residues_from_cif(cif_gz_path):
    """Return space-separated 'ChainRes' labels for atoms with is_motif_atom_with_fixed_seq=True."""
    seen, order = set(), []
    with gzip.open(cif_gz_path, "rt") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 18 or parts[17] != "True":
                continue
            key = (parts[11], int(parts[9]))   # (auth_asym_id, auth_seq_id)
            if key not in seen:
                seen.add(key)
                order.append(key)
    order.sort()
    return " ".join(f"{c}{r}" for c, r in order)


# ── 1. Collect RFD3 metrics and save CSV ─────────────────────────────────────
print("=== Step 2: Processing RFD3 outputs ===", flush=True)
rfd3_csv = f"{config.SCORES_DIR}/rfd3_metrics_with_json_path.csv"
json_paths = sorted(glob.glob(f"{config.DIFFUSION_OUT_DIR}/*/*json"))
if not json_paths:
    sys.exit("ERROR: No RFD3 JSON files found — did the RFD3 job complete successfully?")

rows, all_keys = [], {"json_path"}
for jp in json_paths:
    with open(jp) as f:
        data = json.load(f)
    row = {"json_path": os.path.abspath(jp)}
    for k, v in data.get("metrics", {}).items():
        if k in ("diffused_com", "fixed_com") and isinstance(v, (list, tuple)) and len(v) == 3:
            row[f"{k}_x"], row[f"{k}_y"], row[f"{k}_z"] = v
            all_keys.update({f"{k}_x", f"{k}_y", f"{k}_z"})
        else:
            row[k] = v
            all_keys.add(k)
    rows.append(row)

fieldnames = ["json_path"] + sorted(k for k in all_keys if k != "json_path")
with open(rfd3_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for r in rows:
        writer.writerow(r)
print(f"Saved RFD3 metrics: {rfd3_csv}  ({len(rows)} designs)", flush=True)

# ── 2. Generate MPNN cmds (filtering bad structures) ─────────────────────────
mpnn_cmds = f"{config.CMDS_DIR}/mpnn.cmds"
bad = 0
with open(mpnn_cmds, "w") as f:
    for structure in sorted(glob.glob(f"{config.DIFFUSION_OUT_DIR}/*/*cif.gz")):
        if config.MPNN_REMOVE_BAD:
            struct = parse(structure)
            m = compute_structure_metrics(struct["asym_unit"][0])
            if m["n_chainbreaks"] > 1 or m["n_clashes_sidechain"] > 0 or m["n_clashes_backbone"] > 0:
                bad += 1
                continue

        bn = os.path.basename(structure).replace(".cif.gz", "").replace(".pdb", "")
        this_outdir = f"{config.MPNN_OUT_DIR}/{bn}"
        os.makedirs(this_outdir, exist_ok=True)
        fixed_resis = _motif_residues_from_cif(structure)

        cmd = (
            f"python /dtu/projects/dbl/LigandMPNN/25.3.1-0/run.py"
            f" --seed {config.MPNN_SEED}"
            f" --pdb_path \"{structure}\""
            f" --out_folder \"{this_outdir}\""
            f" --fixed_residues \"{fixed_resis}\""
            f" --batch_size {config.MPNN_BATCH_SIZE}"
            f" --number_of_batches {config.MPNN_N_BATCHES}"
            f" --model_type \"ligand_mpnn\""
            f" --checkpoint_ligand_mpnn \"{config.MPNN_CKPT}\""
            f" --chains_to_design \"{config.MPNN_CHAINS_TO_DESIGN}\""
        )
        f.write(cmd + "\n")

n_tasks = sum(1 for _ in open(mpnn_cmds))
print(f"Filtered {bad} bad structures. Generated {n_tasks} MPNN commands.", flush=True)

if n_tasks == 0:
    sys.exit("ERROR: 0 MPNN commands generated — all structures were filtered out.")

jupyter_utils.make_sub_script(
    mpnn_cmds, n_task=n_tasks, group_size=10, mem="10G", queue="hpc",
    job_name="mpnn", cores=4, time_limit="2:00",
    python_path="/dtu/projects/dbl/LigandMPNN/25.3.1-0/",
    env="/dtu/projects/dbl/LigandMPNN/25.3.1-0//miniforge3/envs/ligandmpnn_env",
    logs_dir=config.LOGS_DIR,
)
mpnn_sh = mpnn_cmds.replace(".cmds", ".sh")

# ── 3. Submit MPNN ────────────────────────────────────────────────────────────
mpnn_jid = _bsub_jid(f"bsub < {mpnn_sh}")
print(f"Submitted MPNN job: {mpnn_jid}", flush=True)

# ── 4. Submit step 3 (waits for MPNN) ────────────────────────────────────────
step3_script = f"{config.PIPELINE_DIR}/step3_process_mpnn.py"
step3_jid = _bsub_jid(
    f'bsub -w "done({mpnn_jid})" -J step3_process_mpnn -q hpc -n 4 '
    f'-W 02:00 -R "rusage[mem=10GB]" -R "span[hosts=1]" '
    f'-o {config.LOGS_DIR}/step3_%J.out -e {config.LOGS_DIR}/step3_%J.err '
    f'"{config.PYTHON_EXEC} {step3_script}"'
)
print(f"Submitted step3 job: {step3_jid}  (waits for MPNN {mpnn_jid})", flush=True)
print("=== Step 2 complete ===", flush=True)
