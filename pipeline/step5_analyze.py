#!/usr/bin/env python3
"""
Step 5 (LSF CPU job): Collect no-ligand RF3 metrics, apply ligand-dependence
filter, save final CSV and scatter plot. Terminal step — no further jobs submitted.
Runs automatically after the RF3-nol job completes via LSF dependency.
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for HPC
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

sys.path.insert(0, f"{config.ENZDES_SCRIPTS}/scripts")
from rf3_sm_metrics import gather_rf3_small_molecule_metrics

NOL_OUT_DIR = f"{config.SHORTLISTED_DIR}/rf3_nol_out"
N = config.RF3_N_MODELS

# ── 1. Collect no-ligand metrics ──────────────────────────────────────────────
print("=== Step 5: Ligand-dependence analysis ===", flush=True)
df_nol = gather_rf3_small_molecule_metrics(
    parent=NOL_OUT_DIR,
    pattern="*_metrics.csv",
    n_models=N,
    add_stats=False,
)

pp_cols_nol = [f"iptm.iptm_protein_protein_{i}" for i in range(N)]
pp_arr_nol  = df_nol[[c for c in pp_cols_nol if c in df_nol.columns]].to_numpy(float)
df_nol["best_iptm_pp_nol"] = np.nanmax(pp_arr_nol, axis=1)

# ── 2. Merge with shortlisted results ─────────────────────────────────────────
shortlisted_csv = f"{config.SHORTLISTED_DIR}/shortlisted_scores.csv"
df_short = pd.read_csv(shortlisted_csv)

df_merged = df_short.merge(
    df_nol[["example_id", "best_iptm_pp_nol"]],
    on="example_id", how="left",
)

# ── 3. Apply ligand-dependence filter ────────────────────────────────────────
# A differential, not two independent thresholds. The original pair
# (nol < 0.50, pp >= 0.20) admitted designs whose switch runs backwards: apo 0.45
# with holo 0.25 satisfies both, yet that design assembles BETTER without ligand.
# Configs that do not set IPTM_PP_DIFF_MIN keep the legacy behaviour.
diff_min = getattr(config, "IPTM_PP_DIFF_MIN", None)

df_merged["iptm_pp_gain"] = df_merged["best_iptm_pp"] - df_merged["best_iptm_pp_nol"]

if diff_min is None:
    df_final = df_merged[
        (df_merged["best_iptm_pp_nol"] < config.IPTM_PP_NOL_CUTOFF) &
        (df_merged["best_iptm_pp"]     >= config.IPTM_PP_CUTOFF)
    ].copy()
    print(f"Ligand-dependent designs: {len(df_final)} / {len(df_merged)} "
          f"(nol iptm_pp < {config.IPTM_PP_NOL_CUTOFF}, "
          f"iptm_pp ≥ {config.IPTM_PP_CUTOFF})  [legacy filter]", flush=True)
else:
    df_final = df_merged[
        (df_merged["iptm_pp_gain"] >= diff_min) &
        (df_merged["best_iptm_pp"] >= config.IPTM_PP_CUTOFF)
    ].copy()
    print(f"Ligand-dependent designs: {len(df_final)} / {len(df_merged)} "
          f"(holo − apo ≥ {diff_min}, holo iptm_pp ≥ "
          f"{config.IPTM_PP_CUTOFF})", flush=True)

n_inverted = int((df_merged["iptm_pp_gain"] < 0).sum())
if n_inverted:
    print(f"WARNING: {n_inverted} / {len(df_merged)} designs assemble better WITHOUT "
          f"ligand (negative iptm_pp gain)", flush=True)

# ── 4. Save final CSV ─────────────────────────────────────────────────────────
final_csv = f"{config.SHORTLISTED_DIR}/final_ligand_dependent.csv"
df_final.to_csv(final_csv, index=False)
print(f"Saved: {final_csv}", flush=True)

# ── 5. Scatter: with-ligand vs without-ligand iptm_pp ────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(df_merged["best_iptm_pp"], df_merged["best_iptm_pp_nol"],
           s=25, alpha=0.6, color="steelblue", label="shortlisted")
ax.scatter(df_final["best_iptm_pp"], df_final["best_iptm_pp_nol"],
           s=60, alpha=0.9, color="tomato", marker="*",
           label=f"retained — ligand-dependent (n={len(df_final)})")
ax.axhline(config.IPTM_PP_NOL_CUTOFF, color="gray", linestyle="--",
           label=f"no-ligand cutoff ({config.IPTM_PP_NOL_CUTOFF})")
ax.axvline(config.IPTM_PP_CUTOFF, color="gray", linestyle=":",
           label=f"with-ligand pp cutoff ({config.IPTM_PP_CUTOFF})")
ax.plot([0, 1], [0, 1], "k--", alpha=0.25, label="equal (no ligand effect)")
ax.set_xlabel("ipTM protein–protein  (with ligand)",    fontsize=12)
ax.set_ylabel("ipTM protein–protein  (without ligand)", fontsize=12)
ax.set_title("Ligand-dependence filter", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
out_png = f"{config.PICS_DIR}/rf3_ligand_dependence.png"
plt.savefig(out_png, dpi=150)
plt.close()
print(f"Saved: {out_png}", flush=True)

# ── 6. Print top designs ──────────────────────────────────────────────────────
top_cols = ["example_id", "best_iptm_pl", "best_iptm_pp", "best_iptm_pp_nol", "best_cif"]
top = df_final[[c for c in top_cols if c in df_final.columns]].sort_values(
    "best_iptm_pl", ascending=False)
print(f"\nTop {min(20, len(top))} ligand-dependent designs:", flush=True)
print(top.head(20).to_string(index=False), flush=True)
print("\n=== Pipeline complete ===", flush=True)
print(f"Final results: {final_csv}", flush=True)
