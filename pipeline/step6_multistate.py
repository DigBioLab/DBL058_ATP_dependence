#!/usr/bin/env python3
"""
Step 6 (LSF): the third gating state — holo + CD3.

Steps 3-5 score two states: holo (VH+VL+ATP) and apo (VH+VL). Those answer
"does ATP glue the Fv together?" but say nothing about whether the rebuilt Fv
still binds its antigen. The S2 gate in docs/ATP_dependent_Fv_plan_v2.md section 3
is "CDR RMSD < 1.0 A, CD3 ipTM retained", and no state in the pipeline so far
contains CD3. This step adds it.

Two phases, chained through LSF:

  --phase generate   read step5's survivors, build VH+VL+CD3+ATP configs,
                     submit RF3, then submit this script again as `analyze`
  --phase analyze    collect per-chain-pair ipSAE, apply the CD3-retention
                     gate, write the final CSV and plot

Why ipSAE and not iptm_protein_protein: in a three-chain complex the aggregate
protein-protein ipTM mixes the VH/VL interface with the Fv/CD3 interface, so a
design that assembles well but has lost CD3 scores the same as one that kept it.
The RF3 .score files carry per-chain-pair rows (chain_i_interface /
chain_j_interface / ipsae_interface), which separate the two.

Run after step5. Can also be run by hand:
    python3 pipeline/step6_multistate.py --phase generate
"""
import sys, os, glob, json, re, argparse, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for HPC
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import submit_utils

sys.path.insert(0, config.ENZDES_SCRIPTS)
sys.path.insert(0, f"{config.ENZDES_SCRIPTS}/scripts")
import jupyter_utils

CD3_OUT_DIR     = f"{config.SHORTLISTED_DIR}/rf3_cd3_out"
CD3_CONFIGS_DIR = f"{config.SHORTLISTED_DIR}/rf3_cd3_configs"
NOL_CONFIGS_DIR = f"{config.SHORTLISTED_DIR}/rf3_nol_configs"
FINAL_CSV       = f"{config.SHORTLISTED_DIR}/final_ligand_dependent.csv"
REF_NAME        = "_wt_okt3_fv_cd3"   # leading underscore keeps it out of design globs


def _bsub_jid(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(r.stdout.strip(), flush=True)
    if r.returncode != 0:
        raise RuntimeError(f"bsub failed:\nSTDOUT: {r.stdout}\nSTDERR: {r.stderr}")
    m = re.search(r"Job <(\d+)>", r.stdout)
    if not m:
        raise RuntimeError(f"No job ID in bsub output: {r.stdout}")
    return m.group(1)


def _seqs_from_nol_config(example_id):
    """{chain_id: seq} for a design, reused from the no-ligand config step4 wrote.

    Cheaper and less fragile than re-deriving from the MPNN FASTA, and guaranteed
    to be the same sequence that was folded in the apo state.
    """
    f = f"{NOL_CONFIGS_DIR}/{example_id}.json"
    if not os.path.exists(f):
        return None
    with open(f) as fh:
        data = json.load(fh)
    return {c["chain_id"]: c["seq"] for c in data[0]["components"] if "seq" in c}


def _seq_from_cif(cif_path, chain_id):
    """One-letter sequence of one chain, read from a reference mmCIF."""
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    structure = MMCIFParser(QUIET=True).get_structure("ref", cif_path)
    for chain in structure[0]:
        if chain.id != chain_id:
            continue
        seq = ""
        for res in chain:
            try:
                seq += index_to_one(three_to_index(res.get_resname()))
            except Exception:
                seq += "X"
        return seq
    raise SystemExit(f"chain {chain_id} not found in {cif_path}")


def _rf3_cmd(json_f, out_dir):
    return (
        f"rf3 fold inference_engine=rf3"
        f" inputs={json_f}"
        f" out_dir={out_dir}"
        f" ckpt_path={config.RF3_CKPT}"
        f" num_steps={config.RF3_NUM_STEPS}"
        f" annotate_b_factor_with_plddt=True"
        f" early_stopping_plddt_threshold=0"
    )


def _write_state(name, components):
    """Write one RF3 input JSON and return (json_path, out_dir)."""
    json_f  = f"{CD3_CONFIGS_DIR}/{name}.json"
    out_dir = f"{CD3_OUT_DIR}/{name}"
    os.makedirs(out_dir, exist_ok=True)
    with open(json_f, "w") as jf:
        json.dump([{"name": name, "components": components}], jf, indent=4)
    return json_f, out_dir


# ══════════════════════════════════════════════════════════════════════════════
# Phase: generate
# ══════════════════════════════════════════════════════════════════════════════
def phase_generate():
    print("=== Step 6 (generate): holo + CD3 state ===", flush=True)

    if not os.path.exists(FINAL_CSV):
        sys.exit(f"missing {FINAL_CSV} — run step5 first.")
    df = pd.read_csv(FINAL_CSV)
    if len(df) == 0:
        sys.exit("step5 retained no designs — nothing to score against CD3.")
    print(f"Designs carried over from step5: {len(df)}", flush=True)

    os.makedirs(CD3_CONFIGS_DIR, exist_ok=True)
    os.makedirs(CD3_OUT_DIR, exist_ok=True)

    cd3_seq = _seq_from_cif(config.CD3_REF_CIF, config.CD3_SOURCE_CHAIN)
    print(f"CD3 antigen: {len(cd3_seq)} residues from "
          f"{config.CD3_REF_CIF} chain {config.CD3_SOURCE_CHAIN} "
          f"-> chain {config.CD3_CHAIN_ID}", flush=True)

    ligand = (config.RF3_GATING_LIGAND if config.RF3_CCD_LIGAND
              else {"path": config.RF3_SDF_PATH})

    cmds_file = f"{config.CMDS_DIR}/rf3_cd3.cmds"
    skipped = 0
    with open(cmds_file, "w") as f_out:
        # ── Designs: VH + VL + CD3 + ATP ──────────────────────────────────────
        for _, row in df.iterrows():
            eid  = row["example_id"]
            seqs = _seqs_from_nol_config(eid)
            if not seqs:
                print(f"WARNING: no apo config for {eid} — skipped", flush=True)
                skipped += 1
                continue
            components  = [{"seq": seqs[c], "chain_id": c} for c in sorted(seqs)]
            components += [{"seq": cd3_seq, "chain_id": config.CD3_CHAIN_ID}, ligand]
            json_f, out_dir = _write_state(eid, components)
            f_out.write(_rf3_cmd(json_f, out_dir) + "\n")

        # ── Reference: wild-type OKT3 Fv + CD3, no ligand ─────────────────────
        # Without this, CD3_IPSAE_MIN is an arbitrary absolute number. With it,
        # "retained" is measurable as a fraction of what the parent antibody
        # achieves under the same model and settings. One extra job.
        if config.CD3_WT_REFERENCE:
            vh = _seq_from_cif(config.CD3_REF_CIF, "A")
            vl = _seq_from_cif(config.CD3_REF_CIF, "B")
            ref_components = [
                {"seq": vh,      "chain_id": "A"},
                {"seq": vl,      "chain_id": "B"},
                {"seq": cd3_seq, "chain_id": config.CD3_CHAIN_ID},
            ]
            json_f, out_dir = _write_state(REF_NAME, ref_components)
            f_out.write(_rf3_cmd(json_f, out_dir) + "\n")
            print("Added wild-type OKT3 Fv + CD3 reference run.", flush=True)

    n_cmds = sum(1 for _ in open(cmds_file))
    if n_cmds == 0:
        sys.exit("No commands generated — every design was skipped.")
    print(f"Generated {n_cmds} holo+CD3 RF3 commands (skipped {skipped}).", flush=True)

    submit_utils.make_sub_script(
        cmds_file, n_task=n_cmds, group_size=10, mem="16GB", time_limit="4:00",
        cores=4, gpu="num=1:mode=exclusive_process",
        queue=config.RF3_QUEUE, job_name="rf3_cd3",
        env="/dtu/projects/dbl/rf3/release/conda_envs/rf3-env",
        python_path="/dtu/projects/dbl/rf3/release",
        logs_dir=config.LOGS_DIR,
    )
    cd3_sh  = cmds_file.replace(".cmds", ".sh")
    cd3_jid = _bsub_jid(f"bsub < {cd3_sh}")
    print(f"Submitted RF3 holo+CD3 job: {cd3_jid}", flush=True)

    this = os.path.abspath(__file__)
    analyze_jid = _bsub_jid(
        f'bsub -w "done({cd3_jid})" -J step6_analyze -q hpc -n 2 '
        f'-W 01:00 -R "rusage[mem=8GB]" -R "span[hosts=1]" '
        f'-o {config.LOGS_DIR}/step6_analyze_%J.out '
        f'-e {config.LOGS_DIR}/step6_analyze_%J.err '
        f'"{config.PYTHON_EXEC} {this} --phase analyze"'
    )
    print(f"Submitted step6 analyze job: {analyze_jid}  (waits for {cd3_jid})", flush=True)
    print("=== Step 6 (generate) complete ===", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# Phase: analyze
# ══════════════════════════════════════════════════════════════════════════════
def _cd3_engagement(example_id):
    """Best-case CD3 engagement for one design.

    Reads the per-chain-pair rows of the RF3 .score files. For each model, CD3
    engagement is min(ipSAE[VH,CD3], ipSAE[VL,CD3]) — the *weaker* of the two
    arms, because section 1.1 of the plan measures the paratope as 59% VH / 41% VL,
    so a design that engages CD3 with one domain only has not retained the
    epitope. Across models the best such value is taken.

    Returns (cd3_ipsae, vh_cd3, vl_cd3) or (nan, nan, nan).
    """
    score_files = glob.glob(f"{CD3_OUT_DIR}/{example_id}/*.score")
    if not score_files:
        return np.nan, np.nan, np.nan

    A, B, C = "A", "B", config.CD3_CHAIN_ID
    best = (np.nan, np.nan, np.nan)
    for sf in score_files:
        try:
            d = pd.read_csv(sf)
        except Exception:
            continue
        need = {"chain_i_interface", "chain_j_interface", "ipsae_interface"}
        if not need.issubset(d.columns):
            continue
        d = d.dropna(subset=["chain_i_interface", "chain_j_interface"]).copy()
        # Chain labels carry a copy suffix, e.g. "A_1".
        d["ci"] = d["chain_i_interface"].astype(str).str.split("_").str[0]
        d["cj"] = d["chain_j_interface"].astype(str).str.split("_").str[0]

        def pair(frame, x, y):
            m = frame[((frame.ci == x) & (frame.cj == y)) |
                      ((frame.ci == y) & (frame.cj == x))]
            return float(m["ipsae_interface"].max()) if len(m) else np.nan

        group_key = "batch_idx" if "batch_idx" in d.columns else None
        groups = d.groupby(group_key) if group_key else [(None, d)]
        for _, grp in groups:
            vh_cd3, vl_cd3 = pair(grp, A, C), pair(grp, B, C)
            if np.isnan(vh_cd3) or np.isnan(vl_cd3):
                continue
            cd3 = min(vh_cd3, vl_cd3)
            if np.isnan(best[0]) or cd3 > best[0]:
                best = (cd3, vh_cd3, vl_cd3)
    return best


def phase_analyze():
    print("=== Step 6 (analyze): CD3-retention gate ===", flush=True)

    df = pd.read_csv(FINAL_CSV)
    print(f"Designs from step5: {len(df)}", flush=True)

    vals = [_cd3_engagement(eid) for eid in df["example_id"]]
    df["cd3_ipsae"]    = [v[0] for v in vals]
    df["cd3_ipsae_vh"] = [v[1] for v in vals]
    df["cd3_ipsae_vl"] = [v[2] for v in vals]

    n_missing = int(df["cd3_ipsae"].isna().sum())
    if n_missing:
        print(f"WARNING: {n_missing} / {len(df)} designs have no usable CD3 "
              f".score output — treated as failures, not dropped.", flush=True)

    # ── Reference ─────────────────────────────────────────────────────────────
    ref = _cd3_engagement(REF_NAME)[0] if config.CD3_WT_REFERENCE else np.nan
    if not np.isnan(ref):
        df["cd3_retention"] = df["cd3_ipsae"] / ref
        print(f"Wild-type OKT3 Fv + CD3 reference ipSAE: {ref:.3f}", flush=True)
    else:
        df["cd3_retention"] = np.nan
        if config.CD3_WT_REFERENCE:
            print("WARNING: reference run produced no score — retention is "
                  "absolute-threshold only, and the threshold is a guess.", flush=True)

    # ── Gate ──────────────────────────────────────────────────────────────────
    # Absolute floor always applies. The relative test applies only when the
    # reference actually ran, so a missing reference cannot silently pass designs.
    passes = df["cd3_ipsae"] >= config.CD3_IPSAE_MIN
    if not np.isnan(ref):
        passes &= df["cd3_retention"] >= config.CD3_RETENTION_MIN
    df["passes_cd3"] = passes.fillna(False)

    df_final = df[df["passes_cd3"]].copy()
    print(f"\nThree-state survivors: {len(df_final)} / {len(df)}", flush=True)
    print(f"  gates: holo-apo gain >= {config.IPTM_PP_DIFF_MIN} (step5), "
          f"holo iptm_pp >= {config.IPTM_PP_CUTOFF} (step5), "
          f"CD3 ipSAE >= {config.CD3_IPSAE_MIN}"
          + (f", CD3 retention >= {config.CD3_RETENTION_MIN}"
             if not np.isnan(ref) else " (no retention test — reference missing)"),
          flush=True)

    all_csv = f"{config.SHORTLISTED_DIR}/multistate_all.csv"
    out_csv = f"{config.SHORTLISTED_DIR}/final_multistate.csv"
    df.to_csv(all_csv, index=False)
    df_final.to_csv(out_csv, index=False)
    print(f"Saved: {all_csv}", flush=True)
    print(f"Saved: {out_csv}", flush=True)

    # ── Plot: switch strength vs antigen retention ────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(df["iptm_pp_gain"], df["cd3_ipsae"],
               s=25, alpha=0.6, color="steelblue", label="step5 survivors")
    ax.scatter(df_final["iptm_pp_gain"], df_final["cd3_ipsae"],
               s=70, alpha=0.9, color="tomato", marker="*",
               label=f"retained CD3 (n={len(df_final)})")
    ax.axvline(config.IPTM_PP_DIFF_MIN, color="gray", linestyle=":",
               label=f"switch gate ({config.IPTM_PP_DIFF_MIN})")
    ax.axhline(config.CD3_IPSAE_MIN, color="gray", linestyle="--",
               label=f"CD3 gate ({config.CD3_IPSAE_MIN})")
    if not np.isnan(ref):
        ax.axhline(ref, color="green", linestyle="-", alpha=0.5,
                   label=f"wild-type OKT3 ({ref:.2f})")
    ax.set_xlabel("ipTM protein-protein gain  (holo - apo)", fontsize=12)
    ax.set_ylabel("CD3 engagement  (min ipSAE over VH/CD3, VL/CD3)", fontsize=12)
    ax.set_title("Three-state gate: ATP switch vs antigen retention", fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()
    out_png = f"{config.PICS_DIR}/rf3_multistate.png"
    plt.savefig(out_png, dpi=150)
    plt.close()
    print(f"Saved: {out_png}", flush=True)

    # ── Top designs ───────────────────────────────────────────────────────────
    cols = ["example_id", "best_iptm_pl", "best_iptm_pp", "best_iptm_pp_nol",
            "iptm_pp_gain", "cd3_ipsae", "cd3_retention"]
    top = df_final[[c for c in cols if c in df_final.columns]].sort_values(
        "iptm_pp_gain", ascending=False)
    print(f"\nTop {min(20, len(top))} three-state designs:", flush=True)
    print(top.head(20).to_string(index=False), flush=True)
    print("\n=== Step 6 complete ===", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", choices=["generate", "analyze"], default="generate")
    args = ap.parse_args()
    (phase_generate if args.phase == "generate" else phase_analyze)()
