#!/usr/bin/env python3
"""
Step 0: Generate RFD3 config JSON and submit script.
Run locally (not as an LSF job) before launching run_pipeline.sh.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

sys.path.insert(0, config.ENZDES_SCRIPTS)
import jupyter_utils

# ── Create all working directories ────────────────────────────────────────────
for d in [config.CMDS_DIR, config.LOGS_DIR, config.PICS_DIR, config.CONFIGS_DIR,
          config.SCORES_DIR, config.INPUTS_DIR, config.DIFFUSION_OUT_DIR,
          config.RF3_OUT_DIR, config.MPNN_OUT_DIR, config.SHORTLISTED_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Write RFD3 config JSON ────────────────────────────────────────────────────
bn     = os.path.basename(config.INPUT_PDB)[:-4]
json_f = f"{config.CONFIGS_DIR}/{bn}.json"

jupyter_utils.write_rfd3_json(
    json_f, config.INPUT_PDB, config.CONTIG, config.LENGTH,
    ligand=config.LIGAND_NAME,
    fixed_atoms=config.FIXED_ATOMS,
    dialect=config.DIALECT,
    redesign_motif_sc=config.REDESIGN_MOTIF_SIDECHAINS,
    hotspots=config.HOTSPOTS,
    infer_ori_strategy=config.INFER_ORI_STRATEGY,
)

# ── Write RFD3 cmds file ──────────────────────────────────────────────────────
outdir = f"{config.DIFFUSION_OUT_DIR}/{bn}"
os.makedirs(outdir, exist_ok=True)

cmd = (
    f"rfd3 design"
    f" out_dir={outdir}"
    f" inputs={json_f}"
    f" ckpt_path={config.RFD3_CKPT}"
    f" diffusion_batch_size={config.RFD3_DIFFUSION_BATCH}"
    f" n_batches={config.RFD3_N_BATCHES}"
    f" inference_sampler.step_scale=3"
    f" inference_sampler.gamma_0=0.2"
)

cmds_file = f"{config.CMDS_DIR}/rfd3.cmds"
with open(cmds_file, "w") as f:
    f.write(cmd + "\n")

# ── Generate LSF submit script ────────────────────────────────────────────────
jupyter_utils.make_sub_script(
    cmds=cmds_file, n_task=1, mem="10GB", time_limit="12:30",
    cores=4, gpu="num=1:mode=exclusive_process",
    queue="cabl40s", job_name="rfd3",
    env="/dtu/projects/dbl/foundry/miniforge3/envs/rfd3",
    python_path="/dtu/projects/dbl/foundry",
    logs_dir=config.LOGS_DIR,
)

sh_file = cmds_file.replace(".cmds", ".sh")
print(f"Config: {json_f}")
print(f"Cmds:   {cmds_file}")
print(f"Script: {sh_file}")
print("\nSetup complete. Launch the full pipeline with:")
print(f"  bash {os.path.join(os.path.dirname(config.PIPELINE_DIR), 'run_pipeline.sh')}")
