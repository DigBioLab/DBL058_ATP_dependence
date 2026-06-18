#!/usr/bin/env python3
"""
Step 4 (LSF CPU job): Collect RF3 metrics, shortlist designs, generate no-ligand
RF3 configs, submit RF3-nol + step 5.
Runs automatically after the RF3 job completes via LSF dependency.
"""
import sys, os, glob, json, gzip, re, subprocess, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for HPC
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

sys.path.insert(0, config.ENZDES_SCRIPTS)
sys.path.insert(0, f"{config.ENZDES_SCRIPTS}/scripts")
import jupyter_utils
from rf3_sm_metrics import gather_rf3_small_molecule_metrics

from Bio import SeqIO


def _bsub_jid(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(r.stdout.strip(), flush=True)
    if r.returncode != 0:
        raise RuntimeError(f"bsub failed:\nSTDOUT: {r.stdout}\nSTDERR: {r.stderr}")
    m = re.search(r"Job <(\d+)>", r.stdout)
    if not m:
        raise RuntimeError(f"No job ID in bsub output: {r.stdout}")
    return m.group(1)


def get_chain_lengths_from_cif(cif_gz_path):
    chains, order = {}, []
    with gzip.open(cif_gz_path, "rt") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            chain, res = parts[11], parts[9]
            if chain not in chains:
                chains[chain] = set()
                order.append(chain)
            chains[chain].add(res)
    return [(c, len(chains[c])) for c in order]


def _seqs_for_design(example_id):
    """Return {chain_id: seq} for a shortlisted design from its MPNN FASTA."""
    mpnn_model, seq_idx = example_id.rsplit("_", 1)
    seq_idx = int(seq_idx)
    fastas = glob.glob(f"{config.MPNN_OUT_DIR}/{mpnn_model}/seqs/*.cif.gz.fa")
    if not fastas:
        return None
    cif_matches = glob.glob(f"{config.DIFFUSION_OUT_DIR}/*/{mpnn_model}.cif.gz")
    chain_lens  = get_chain_lengths_from_cif(cif_matches[0]) if cif_matches else []
    for i, r in enumerate(SeqIO.parse(fastas[0], "fasta")):
        if i == seq_idx:
            seq, pos, seqs = str(r.seq), 0, {}
            for chain_id, n_res in chain_lens:
                seqs[chain_id] = seq[pos:pos + n_res]
                pos += n_res
            return seqs
    return None


# ── 1. Collect RF3 metrics ────────────────────────────────────────────────────
print("=== Step 4: Processing RF3 outputs ===", flush=True)
N = config.RF3_N_MODELS

df = gather_rf3_small_molecule_metrics(
    parent=config.RF3_OUT_DIR,
    pattern="*_metrics.csv",
    n_models=N,
    add_stats=False,
)

pl_cols    = [f"iptm.iptm_protein_ligand_{i}" for i in range(N)]
pp_cols    = [f"iptm.iptm_protein_protein_{i}" for i in range(N)]
clash_cols = [f"count_clashing_chains.has_clash_{i}" for i in range(N)]

pl_arr    = df[[c for c in pl_cols    if c in df.columns]].to_numpy(float)
pp_arr    = df[[c for c in pp_cols    if c in df.columns]].to_numpy(float)
clash_arr = df[[c for c in clash_cols if c in df.columns]].to_numpy(float)

rows = np.arange(len(df))
best_idx  = np.nanargmax(pl_arr, axis=1)
worst_idx = np.nanargmin(pl_arr, axis=1)

df["best_model_idx"] = best_idx
df["best_iptm_pl"]   = pl_arr[rows, best_idx]
df["best_iptm_pp"]   = pp_arr[rows, best_idx]
df["best_clash"]     = clash_arr[rows, best_idx]
df["worst_iptm_pl"]  = pl_arr[rows, worst_idx]
df["worst_iptm_pp"]  = pp_arr[rows, worst_idx]
df["worst_clash"]    = clash_arr[rows, worst_idx]
df["any_clash"]      = (clash_arr > 0).any(axis=1)

summary_cols = [
    "example_id",
    "best_iptm_pl", "best_iptm_pp", "best_clash",
    "worst_iptm_pl", "worst_iptm_pp", "worst_clash",
    "any_clash",
]
out_csv = f"{config.SCORES_DIR}/rf3_iptm_summary.csv"
df[summary_cols].to_csv(out_csv, index=False)
print(f"Saved RF3 summary: {out_csv}  ({len(df)} designs)", flush=True)
print(f"Designs with no clash in any model: {(~df['any_clash']).sum()} / {len(df)}", flush=True)

# ── 2. Scatter plots ──────────────────────────────────────────────────────────
for label, x_col, y_col, clash_col, title in [
    ("best",  "best_iptm_pl",  "best_iptm_pp",  "best_clash",
     "RF3 — best model per design (highest ipTM protein–ligand)"),
    ("worst", "worst_iptm_pl", "worst_iptm_pp", "worst_clash",
     "RF3 — worst model per design (lowest ipTM protein–ligand)"),
]:
    no_clash  = df[df[clash_col] == 0]
    has_clash = df[df[clash_col] >  0]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(no_clash[x_col],  no_clash[y_col],  s=25, alpha=0.75, color="steelblue",
               label=f"no clash (n={len(no_clash)})")
    ax.scatter(has_clash[x_col], has_clash[y_col], s=25, alpha=0.5,  color="tomato",
               marker="x", label=f"clash (n={len(has_clash)})")
    ax.set_xlabel("ipTM protein–ligand", fontsize=12)
    ax.set_ylabel("ipTM protein–protein", fontsize=12)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out_png = f"{config.PICS_DIR}/rf3_iptm_{label}.png"
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"Saved: {out_png}", flush=True)

# ── 3. Shortlist ──────────────────────────────────────────────────────────────
df_short = df[
    (df["best_iptm_pl"] >= config.IPTM_PL_CUTOFF) &
    (df["best_iptm_pp"] >= config.IPTM_PP_CUTOFF) &
    (~df["any_clash"])
].copy()
print(f"Shortlisted: {len(df_short)} / {len(df)} designs "
      f"(pl≥{config.IPTM_PL_CUTOFF}, pp≥{config.IPTM_PP_CUTOFF}, no clash)", flush=True)

# Copy best-model CIFs
shortlisted_structs = f"{config.SHORTLISTED_DIR}/structures"
os.makedirs(shortlisted_structs, exist_ok=True)

def _best_cif_path(example_id, best_model_idx):
    design_folder = example_id.rsplit("_", 1)[0]
    fname = f"{example_id}_model_{int(best_model_idx)}.cif.gz"
    return os.path.join(config.RF3_OUT_DIR, design_folder, fname)

df_short["best_cif"] = df_short.apply(
    lambda r: _best_cif_path(r["example_id"], r["best_model_idx"]), axis=1)

copied = 0
for _, row in df_short.iterrows():
    if os.path.exists(row["best_cif"]):
        shutil.copy2(row["best_cif"], os.path.join(shortlisted_structs, os.path.basename(row["best_cif"])))
        copied += 1
print(f"Copied {copied} CIFs → {shortlisted_structs}", flush=True)

shortlisted_csv = f"{config.SHORTLISTED_DIR}/shortlisted_scores.csv"
df_short[summary_cols + ["best_model_idx", "best_cif"]].to_csv(shortlisted_csv, index=False)
print(f"Saved: {shortlisted_csv}", flush=True)

if len(df_short) == 0:
    sys.exit("No designs passed shortlisting thresholds — adjust IPTM_PL_CUTOFF / IPTM_PP_CUTOFF in config.py.")

# ── 4. Generate no-ligand RF3 configs and cmds ───────────────────────────────
nol_configs_dir = f"{config.SHORTLISTED_DIR}/rf3_nol_configs"
nol_out_dir     = f"{config.SHORTLISTED_DIR}/rf3_nol_out"
nol_cmds        = f"{config.CMDS_DIR}/rf3_nol.cmds"
os.makedirs(nol_configs_dir, exist_ok=True)
os.makedirs(nol_out_dir, exist_ok=True)

skipped = 0
with open(nol_cmds, "w") as f_out:
    for _, row in df_short.iterrows():
        eid  = row["example_id"]
        seqs = _seqs_for_design(eid)
        if not seqs:
            print(f"WARNING: no sequence found for {eid}", flush=True)
            skipped += 1
            continue

        json_f  = f"{nol_configs_dir}/{eid}.json"
        out_dir = f"{nol_out_dir}/{eid}"
        os.makedirs(out_dir, exist_ok=True)

        config_data = [{"name": eid, "components": [
            {"seq": seqs[c], "chain_id": c} for c in sorted(seqs)
            # ligand intentionally omitted
        ]}]
        with open(json_f, "w") as jf:
            json.dump(config_data, jf, indent=4)

        cmd = (
            f"rf3 fold inference_engine=rf3"
            f" inputs={json_f}"
            f" out_dir={out_dir}"
            f" ckpt_path={config.RF3_CKPT}"
            f" num_steps={config.RF3_NUM_STEPS}"
            f" annotate_b_factor_with_plddt=True"
            f" early_stopping_plddt_threshold=0"
        )
        f_out.write(cmd + "\n")

n_nol = sum(1 for _ in open(nol_cmds))
print(f"Generated {n_nol} no-ligand RF3 commands (skipped {skipped}).", flush=True)

jupyter_utils.make_sub_script(
    nol_cmds, n_task=n_nol, group_size=10, mem="10GB", time_limit="4:00",
    cores=4, gpu="num=1:mode=exclusive_process",
    queue="cabl40s", job_name="rf3_nol",
    env="/dtu/projects/dbl/rf3/release/conda_envs/rf3-env",
    python_path="/dtu/projects/dbl/rf3/release",
    logs_dir=config.LOGS_DIR,
)
nol_sh = nol_cmds.replace(".cmds", ".sh")

# ── 5. Submit RF3-nol ─────────────────────────────────────────────────────────
nol_jid = _bsub_jid(f"bsub < {nol_sh}")
print(f"Submitted RF3-nol job: {nol_jid}", flush=True)

# ── 6. Submit step 5 (waits for RF3-nol) ─────────────────────────────────────
step5_script = f"{config.PIPELINE_DIR}/step5_analyze.py"
step5_jid = _bsub_jid(
    f'bsub -w "done({nol_jid})" -J step5_analyze -q hpc -n 2 '
    f'-W 01:00 -R "rusage[mem=8GB]" -R "span[hosts=1]" '
    f'-o {config.LOGS_DIR}/step5_%J.out -e {config.LOGS_DIR}/step5_%J.err '
    f'"{config.PYTHON_EXEC} {step5_script}"'
)
print(f"Submitted step5 job: {step5_jid}  (waits for RF3-nol {nol_jid})", flush=True)
print("=== Step 4 complete ===", flush=True)
