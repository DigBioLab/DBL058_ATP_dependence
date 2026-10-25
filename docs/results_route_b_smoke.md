# Route B smoke test — partial diffusion at partial_t 2–5 Å

**Date:** 2026-10-25 · **Branch:** `feature/darwol-atp-gated-cd3-fv`
**Campaign:** `okt3_atp_fv_smoke` · **Plan:** [`ATP_dependent_Fv_plan_v2.md`](ATP_dependent_Fv_plan_v2.md) §4 Route B

**Headline: Route B did not excavate a pocket. 0 of 80 designs reached ATP scale,
and mean enclosed volume fell below the starting structure.** The approach is not
yet ruled out — three untested variables are listed in §6 — but as configured it
does not work.

---

## 1. What was run

RFD3 partial diffusion of the whole OKT3 Fv with ATP present, four noise levels,
20 designs each (80 total). Scope deliberately small: a feasibility check before
committing the 10³–10⁴ budget.

| Setting | Value |
|---|---|
| Input | `okt3_fv_atp.cif` — VH (A, 113 res) + VL (B, 107 res) + ATP (L, 31 atoms) |
| Unindexed | all 18 REDESIGNABLE residues, `A37,A39,A43-46,A95,A109-112,B37,B41-43,B45,B86,B97` |
| Pinned | 19 PARATOPE_WIDE residues, all atoms (328 atoms) |
| `partial_t` | 2.0, 3.0, 4.0, 5.0 Å |
| ATP | rigid — not configurable, see §5 |
| Designs | 20 per arm (`RFD3_N_BATCHES=2` × batch 10) |
| Compute | 1 GPU (cabgpu), ~3 min wall for all 80 |

All 80 completed. No tracebacks, no CUDA errors, no ligand clashes.
`diffused_index_map` in every output confirms exactly the intended 18 residues
were diffused.

## 2. Structural metrics — the design mechanism works

From RFD3's own per-model metrics.

| partial_t | insertion RMSD | max Cα dev | ligand clashes | bb clashes | chainbreaks | sheet frac |
|---|---|---|---|---|---|---|
| 2.0 | 0.60 | 0.92 | 0 | 0 | 1.95 | 0.29 |
| 3.0 | 0.69 | 0.93 | 0 | 0 | 2.15 | 0.29 |
| 4.0 | 0.69 | 0.93 | 0 | 0 | 2.05 | 0.28 |
| 5.0 | 0.78 | 1.08 | 0 | 0 | 2.20 | 0.27 |

Deformation scales with noise as expected, and the fold survives: sheet content
is stable and there are no backbone clashes. Chainbreaks ≈ 2 is consistent with
the Fv being two chains, though **this has not been verified to be the A/B break
rather than new breaks inside the diffused segments** — worth checking before
any of these backbones are reused.

## 3. ATP contact counts — misleading, recorded for honesty

| | seed | pt 2 | pt 3 | pt 4 | pt 5 |
|---|---|---|---|---|---|
| ATP heavy-atom contacts <5 Å (mean) | 55 | 68.9 | 67.9 | 65.6 | 65.6 |
| (best single model) | — | 76 | 79 | 79 | 82 |

Contacts rose in every arm, which initially read as "the protein closed around
ATP". The cavity scan in §4 shows that reading was wrong: the protein draped
over the ligand surface without hollowing anything out underneath. **Contact
count rewards exactly the failure mode we care about avoiding** and should not
be used as a pocket metric again.

## 4. Cavity scan — the decisive measurement

Method: [`scripts/cavity_scan.py`](../scripts/cavity_scan.py), reimplementing the
§1.3 protocol (1 Å voxels, 1.4 Å probe, ray-occlusion buriedness over 42
directions). ATP is excluded from the atom set so it cannot fill its own cavity.
Volumes are the **largest single connected cavity**, which is what matters —
ATP is ~350–400 Å³ as one contiguous molecule.

| buriedness | seed | pt 2 | pt 3 | pt 4 | pt 5 |
|---|---|---|---|---|---|
| ≥70 % (mean) | 160 | 155 | 154 | 153 | 149 |
| ≥80 % (mean) | **33** | 16.3 | 14.4 | 13.8 | 15.7 |
| ≥80 % (best) | 33 | 28 | 32 | 23 | **44** |
| ≥90 % (mean) | 1 | 1.1 | 1.1 | 1.1 | 1.1 |

- Designs beating the seed's 33 Å³ enclosed pocket: **1 of 80** (44 Å³).
- Designs reaching ATP scale (≥300 Å³ at ≥80 %): **0 of 80**.
- At ≥90 % buriedness every arm sits at ~1 Å³ — identical to the seed and to apo
  1SY6. Nothing enclosed formed anywhere.

Mean enclosed volume **decreased** relative to the seed (33 → 14–16 Å³).
Partial diffusion tightened packing rather than opening it. The trend across
`partial_t` is flat and within noise.

### 4.1 Caveat: these volumes are not comparable to §1.3's

`cavity_scan.py` reproduces §1.3's *conclusion* but not its *numbers*. On the
same apo Fv:

| buriedness | §1.3 (total/largest) | this script (total/largest) |
|---|---|---|
| ≥90 % | 0 / — | 9 / 2 |
| ≥80 % | 49 / 17 | 156 / 30 |
| ≥70 % | 208 / 75 | 526 / 154 |
| ≥60 % | 578 / 339 | 1044 / 596 |

Systematically 2–3× larger, because §1.3's script was never committed and its
region definition, ray count, ray length and vdW radii are unrecorded. The
qualitative finding is identical — essentially nothing enclosed at ≥90 %.

**Use these volumes for relative comparison only** (seed vs design, arm vs arm,
all scanned identically). Do not quote them against §1.3 or against literature
pocket sizes.

## 5. Constraint discovered: ATP cannot be made mobile

Under partial diffusion in this RFD3 build the ligand is **always
coordinate-fixed**, and no input-spec setting changes it:

- [`input_parsing.py:672`](/dtu/projects/dbl/foundry/models/rfd3/src/rfd3/inference/input_parsing.py) calls
  `extract_ligand_array(..., fixed_atoms={}, ...)` with a hardcoded empty dict,
  so the unfix path in `utils/inference.py:119-134` can never fire.
- The flag is then forced on: `is_motif_atom_with_fixed_coord = np.ones(...)`.

Confirmed empirically — the built pipeline input has 31/31 ATP atoms fixed.

The rigid-vs-free comparison requested as a design variable therefore cannot be
run through configuration. It needs either Route A or a patched RFD3.

Related: RFD3 extracts the ligand **from the input structure** rather than from
the CCD, so the input CIF must already contain ATP. `prepare_okt3_fv.py --with-atp`
now seeds it; an apo input fails with `Invalid contig format: 'ATP'`.

## 6. Interpretation

**Confirmed by measurement.** Route B at `partial_t` 2–5 Å, paratope pinned and
ATP rigid, does not excavate an ATP-scale pocket in 80 attempts. Mean enclosed
volume is worse than the starting structure.

**Inference, not established.** Three candidate causes, none tested:

1. **Rigid ATP may be the core obstacle.** The ligand cannot move aside while the
   backbone reorganises, so draping the protein over it is the only accessible
   solution. This promotes ligand mobility from a nice-to-have to a central
   variable (§5).
2. **The seed may be too exposed.** ATP sits 10.5 Å out along the groove normal —
   the first clash-free offset, because the §1.3 groove centroid is ~10 Å inside
   the protein. 2–5 Å of noise plausibly cannot translate a rigid ligand into a
   pocket that would have to form 10 Å away.
3. **2–5 Å may be too timid.** RFD3 permits up to 15 Å and the docs note the
   response is nonlinear. Nothing here tests 8–15 Å.

**What this does not say.** Nothing about whether the switch works — no MPNN, no
RF3, no apo/holo comparison was run. This tests only whether RFD3 will open a
cavity.

## 7. Reproducing

```bash
export ENZDES_SCRIPTS=/work3/s212001/DBL_secrets/enzdes_pipeline-main
python3 scripts/prepare_okt3_fv.py --out-dir inputs/ --with-atp
# set PARTIAL_DIFFUSION = True in configs/okt3_atp_fv.py, then
ATP_WORKING_DIR=<dir> python3 pipeline/step0_setup.py
bsub < <dir>/cmds/rfd3.sh
python3 scripts/cavity_scan.py --file-list <files> --around-ligand ATP \
    --half-extent 10 --out cavity.csv
```

Data: `okt3_atp_fv_smoke/partial_t_sweep.csv` (per-model metrics),
`okt3_atp_fv_smoke/cavity/cavity_all.csv` (per-model cavity volumes),
`inputs/atp_seed_placement.json` (seed geometry).

## 8. Open items

- [x] `partial_t` 8–15 Å sweep (§6.3) — **no effect**, see
      [`results_excavation_options.md`](results_excavation_options.md)
- [x] Deeper ATP seed, accepting initial clashes (§6.2) — **marginal only**, ibid.
- [x] Patched RFD3 with mobile ligand (§5) — **no effect**; ligand rigidity was
      not the obstacle, ibid. §3.1
- [x] Route A comparison (§6.1) — **substantially worse than Route B**, ibid. §3.2
- [ ] Verify the ~2 chainbreaks are the A/B chain break, not new breaks (§2).

**All three hypotheses in §6 were tested and none explains the failure.**
The follow-up concludes that the constant across every arm — only 18 residues
rebuildable, paratope pinned — is the remaining untested factor. See
[`results_excavation_options.md`](results_excavation_options.md) §5 for the
recommended next steps.
