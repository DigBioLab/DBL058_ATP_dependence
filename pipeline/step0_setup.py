#!/usr/bin/env python3
"""
Step 0: Generate RFD3 config JSON and submit script.
Run locally (not as an LSF job) before launching run_pipeline.sh.

Two routes, selected by config.PARTIAL_DIFFUSION:

  False (Route A, default)  inpainting — CONTIG/LENGTH replace the redesignable
                            face with unindexed segments. One job.
  True  (Route B)           partial diffusion — the whole Fv is noised by
                            partial_t Angstroms with the paratope pinned, and the
                            redesignable face unindexed so it can be rebuilt.
                            One job per value in PARTIAL_T_SWEEP.

Route B writes its own JSON rather than calling jupyter_utils.write_rfd3_json,
which has no `partial_t` or `unindex` parameter. That helper lives in the shared
DBL_secrets checkout used by other campaigns, so it is left untouched.
"""
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import submit_utils

sys.path.insert(0, config.ENZDES_SCRIPTS)
import jupyter_utils


def collapse_residues(residues):
    """['A37','A39','A43','A44','A45'] -> 'A37,A39,A43-45'.

    RFD3 reads consecutive residues written as a range (A43-45) as tied in
    sequence, which is the correct hint here: these are contiguous stretches of
    one framework face, not independent single residues.
    """
    parsed = []
    for r in residues:
        m = re.match(r"^([A-Za-z]+)(\d+)$", r)
        if not m:
            raise SystemExit(f"cannot parse residue label {r!r}")
        parsed.append((m.group(1), int(m.group(2))))

    out = []
    for chain in sorted({c for c, _ in parsed}):
        nums = sorted(n for c, n in parsed if c == chain)
        start = prev = nums[0]
        for n in nums[1:] + [None]:
            if n is not None and n == prev + 1:
                prev = n
                continue
            out.append(f"{chain}{start}" if start == prev else f"{chain}{start}-{prev}")
            if n is not None:
                start = prev = n
    return ",".join(out)


def write_partial_diffusion_json(path, name, partial_t, unindex, fixed_atoms):
    """RFD3 input for one partial-diffusion arm.

    No contig/length: under partial diffusion RFD3 skips contig sampling
    entirely (input_parsing.py, `if not self.is_partial_diffusion`), so a contig
    would be silently ignored.

    ATP is rigid. Under partial diffusion the input array is subset to protein
    before `select_fixed_atoms` is applied, and the ligand is then re-appended
    with its fixed-coordinate flag hard-coded on. There is no input-spec knob to
    free it — see the note in configs/okt3_atp_fv.py.
    """
    payload = {name: {
        "dialect":            config.DIALECT,
        "input":              config.INPUT_PDB,
        "ligand":             config.LIGAND_NAME,
        "partial_t":          float(partial_t),
        "unindex":            unindex,
        "select_fixed_atoms": fixed_atoms,
        "is_non_loopy":       config.IS_NON_LOOPY,
    }}
    # Only understood by the patched RFD3 at /work3/s212001/rfd3_patched.
    # Stock RFD3 rejects unknown keys, so this is omitted unless configured.
    lig_fixed = getattr(config, "LIGAND_FIXED_ATOMS", None)
    if lig_fixed is not None:
        payload[name]["ligand_fixed_atoms"] = lig_fixed
    with open(path, "w") as f:
        json.dump(payload, f, indent=4)
    return path

# ── Create all working directories ────────────────────────────────────────────
for d in [config.CMDS_DIR, config.LOGS_DIR, config.PICS_DIR, config.CONFIGS_DIR,
          config.SCORES_DIR, config.INPUTS_DIR, config.DIFFUSION_OUT_DIR,
          config.RF3_OUT_DIR, config.MPNN_OUT_DIR, config.SHORTLISTED_DIR]:
    os.makedirs(d, exist_ok=True)

bn = os.path.basename(config.INPUT_PDB)[:-4]


def rfd3_cmd(json_f, outdir):
    return (
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
cmds, arms = [], []

if getattr(config, "PARTIAL_DIFFUSION", False):
    # ── Route B: partial diffusion, one arm per noise level ──────────────────
    unindex = collapse_residues(config.REDESIGNABLE_RESIDUES)
    # The paratope is the only thing pinned. Everything else — framework,
    # CDR backbones, the VH/VL packing — is free to move by partial_t, which is
    # the deformation section 1.3 of the plan says is needed to open a cavity.
    fixed_atoms = {r: "ALL" for r in config.PARATOPE_RESIDUES}

    print(f"Route B: partial diffusion, {len(config.PARTIAL_T_SWEEP)} arms")
    print(f"  unindexed ({len(config.REDESIGNABLE_RESIDUES)} residues): {unindex}")
    print(f"  pinned    ({len(config.PARATOPE_RESIDUES)} paratope residues, all atoms)")
    print(f"  ligand    : {config.LIGAND_NAME} (rigid — not configurable, see config)")

    for t in config.PARTIAL_T_SWEEP:
        tag    = f"{bn}_pt{t:g}"
        json_f = write_partial_diffusion_json(
            f"{config.CONFIGS_DIR}/{tag}.json", tag, t, unindex, fixed_atoms)
        outdir = f"{config.DIFFUSION_OUT_DIR}/{tag}"
        os.makedirs(outdir, exist_ok=True)
        cmds.append(rfd3_cmd(json_f, outdir))
        arms.append((tag, t, json_f, outdir))
        print(f"  arm partial_t={t:g} -> {json_f}")
else:
    # ── Route A: inpainting from CONTIG/LENGTH ───────────────────────────────
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
    outdir = f"{config.DIFFUSION_OUT_DIR}/{bn}"
    os.makedirs(outdir, exist_ok=True)
    cmds.append(rfd3_cmd(json_f, outdir))
    arms.append((bn, None, json_f, outdir))
    print("Route A: inpainting (CONTIG/LENGTH), 1 arm")

with open(cmds_file, "w") as f:
    f.write("\n".join(cmds) + "\n")

# Record which arm is which, so the sweep comparison does not have to re-derive
# it from directory names later.
with open(f"{config.CONFIGS_DIR}/rfd3_arms.json", "w") as f:
    json.dump([{"tag": t, "partial_t": pt, "config": j, "out_dir": o}
               for t, pt, j, o in arms], f, indent=4)

# ── Generate LSF submit script ────────────────────────────────────────────────
submit_utils.make_sub_script(
    cmds=cmds_file, n_task=len(cmds), mem="10GB", time_limit="12:30",
    cores=4, gpu="num=1:mode=exclusive_process",
    queue=getattr(config, "RFD3_QUEUE", "cabl40s"), job_name="rfd3",
    env="/dtu/projects/dbl/foundry/miniforge3/envs/rfd3",
    python_path="/dtu/projects/dbl/foundry",
    logs_dir=config.LOGS_DIR,
)

sh_file = cmds_file.replace(".cmds", ".sh")
print(f"\nArms:   {len(cmds)}")
print(f"Cmds:   {cmds_file}")
print(f"Script: {sh_file}")
print("\nSetup complete. Launch the full pipeline with:")
print(f"  bash {os.path.join(os.path.dirname(config.PIPELINE_DIR), 'run_pipeline.sh')}")
