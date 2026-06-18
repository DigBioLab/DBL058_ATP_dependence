#!/usr/bin/env python3
"""
Step 3 (LSF CPU job): Process MPNN outputs, generate RF3 cmds, submit RF3 + step 4.
Runs automatically after the MPNN job completes via LSF dependency.
"""
import sys, os, glob, json, gzip, re, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

sys.path.insert(0, config.ENZDES_SCRIPTS)
import jupyter_utils

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
    """Return ordered list of (chain_id, n_residues) for protein chains (ATOM records)."""
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


# ── Ligand entry ──────────────────────────────────────────────────────────────
if config.RF3_CCD_LIGAND:
    ligand_entry = config.LIGAND_NAME
elif config.RF3_SDF_PATH:
    ligand_entry = {"path": config.RF3_SDF_PATH}
else:
    ligand_entry = {"smiles": config.RF3_SMILES}

# ── 1. Generate RF3 configs and cmds ─────────────────────────────────────────
print("=== Step 3: Generating RF3 inputs ===", flush=True)
rf3_configs_dir = f"{config.CONFIGS_DIR}/rf3"
os.makedirs(rf3_configs_dir, exist_ok=True)

rf3_cmds = f"{config.CMDS_DIR}/rf3.cmds"
fasta_files = sorted(glob.glob(f"{config.MPNN_OUT_DIR}/*/seqs/*cif.gz.fa"))
if not fasta_files:
    sys.exit("ERROR: No MPNN FASTA files found — did the MPNN job complete successfully?")

with open(rf3_cmds, "w") as f:
    for fasta_f in fasta_files:
        name    = os.path.basename(fasta_f).split(".")[0]
        json_f  = f"{rf3_configs_dir}/{name}.json"
        out_dir = f"{config.RF3_OUT_DIR}/{name}"
        os.makedirs(out_dir, exist_ok=True)

        cif_matches = glob.glob(f"{config.DIFFUSION_OUT_DIR}/*/{name}.cif.gz")
        chain_lens  = get_chain_lengths_from_cif(cif_matches[0]) if cif_matches else []

        results = []
        for i, r in enumerate(SeqIO.parse(fasta_f, "fasta")):
            if i == 0:
                continue   # skip original RFD3 alanine-filled backbone

            seq, pos, components = str(r.seq), 0, []
            for chain_id, n_res in chain_lens:
                components.append({"seq": seq[pos:pos + n_res], "chain_id": chain_id})
                pos += n_res
            components.append(ligand_entry)
            results.append({"name": f"{name}_{i}", "components": components})

        with open(json_f, "w") as jf:
            json.dump(results, jf, indent=4)

        cmd = (
            f"rf3 fold inference_engine=rf3"
            f" inputs={json_f}"
            f" out_dir={out_dir}"
            f" ckpt_path={config.RF3_CKPT}"
            f" num_steps={config.RF3_NUM_STEPS}"
            f" annotate_b_factor_with_plddt=True"
            f" early_stopping_plddt_threshold=0"
        )
        f.write(cmd + "\n")

n_task = sum(1 for _ in open(rf3_cmds))
print(f"Generated {n_task} RF3 commands across {len(fasta_files)} FASTA files.", flush=True)

if n_task == 0:
    sys.exit("ERROR: 0 RF3 commands generated.")

jupyter_utils.make_sub_script(
    cmds=rf3_cmds, n_task=n_task, mem="10GB", time_limit="4:00",
    group_size=10, cores=4, gpu="num=1:mode=exclusive_process",
    queue="cabl40s", job_name="rf3",
    env="/dtu/projects/dbl/rf3/release/conda_envs/rf3-env",
    python_path="/dtu/projects/dbl/rf3/release",
    logs_dir=config.LOGS_DIR,
)
rf3_sh = rf3_cmds.replace(".cmds", ".sh")

# ── 2. Submit RF3 ─────────────────────────────────────────────────────────────
rf3_jid = _bsub_jid(f"bsub < {rf3_sh}")
print(f"Submitted RF3 job: {rf3_jid}", flush=True)

# ── 3. Submit step 4 (waits for RF3) ─────────────────────────────────────────
step4_script = f"{config.PIPELINE_DIR}/step4_process_rf3.py"
step4_jid = _bsub_jid(
    f'bsub -w "done({rf3_jid})" -J step4_process_rf3 -q hpc -n 4 '
    f'-W 02:00 -R "rusage[mem=10GB]" -R "span[hosts=1]" '
    f'-o {config.LOGS_DIR}/step4_%J.out -e {config.LOGS_DIR}/step4_%J.err '
    f'"{config.PYTHON_EXEC} {step4_script}"'
)
print(f"Submitted step4 job: {step4_jid}  (waits for RF3 {rf3_jid})", flush=True)
print("=== Step 3 complete ===", flush=True)
