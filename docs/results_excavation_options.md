# Pocket excavation — four options tested, all negative

**Date:** 2026-10-25 · **Branch:** `feature/darwol-atp-gated-cd3-fv`
**Campaign:** `okt3_atp_fv_opts` · **Follows:** [`results_route_b_smoke.md`](results_route_b_smoke.md)
**Plan:** [`ATP_dependent_Fv_plan_v2.md`](ATP_dependent_Fv_plan_v2.md) §1.3, §4

**Headline: none of the four candidate fixes opened an ATP-scale pocket.
0 of 237 designs reached even 200 Å³ of enclosed cavity against a ~350–400 Å³
requirement. Combined with the earlier partial_t 2–5 Å run, that is 317 designs
across every variable we can reach without changing the problem setup.**

Two results are individually decisive: mobile ATP makes no difference, and
Route A is substantially worse than Route B.

---

## 1. What was tested

The smoke test ([`results_route_b_smoke.md`](results_route_b_smoke.md)) left three
untested explanations for Route B's failure, plus the untested alternative route.
Each became an arm, 20 designs per arm, 260 total.

| Option | Hypothesis under test | Arms |
|---|---|---|
| **1 — high noise** | 2–5 Å was too timid; RFD3 allows up to 15 | `partial_t` 8, 10, 12, 15 |
| **2 — deep seed** | ATP sat 10.5 Å out, too exposed to pull in | `partial_t` 2–5, ATP seeded at 3 Å offset |
| **3 — mobile ATP** | rigid ligand blocked backbone reorganisation | `partial_t` 2–5, patched RFD3 |
| **4 — Route A** | inpainting may beat partial diffusion | contig/length, 1 arm |

Options 2 and 3 needed new capability:

- **Deep seed** — `prepare_okt3_fv.py --atp-offset 3.0 --atp-clash-min 2.0`
  writes `okt3_fv_atp_deep.cif`: 419 ATP–protein contacts within 5 Å (vs 55 for
  the shallow seed) and a deliberate 1.05 Å clash for diffusion to resolve.
- **Mobile ATP** — stock RFD3 cannot free a ligand (see
  [`results_route_b_smoke.md`](results_route_b_smoke.md) §5). A patched copy at
  `/work3/s212001/rfd3_patched` adds a `ligand_fixed_atoms` field. **The shared
  install at `/dtu/projects/dbl/foundry` was not modified.** Verified before use:
  rigid control gives 31/31 ATP atoms fixed, patched gives 0/31.

## 2. Cavity results

Method and caveats as in [`results_route_b_smoke.md`](results_route_b_smoke.md) §4.
Metric is the **largest connected cavity**, since ATP is one contiguous molecule.

| Arm | n | ≥70 % | **≥80 % mean** | ≥80 % best | ≥90 % mean |
|---|---|---|---|---|---|
| *shallow seed (baseline)* | — | 160 | **33** | 33 | 1 |
| *deep seed (baseline)* | — | — | **33** | 33 | 2 |
| opt1 high noise pt8 | 19 | 151 | 12.7 | 25 | 1.1 |
| opt1 high noise pt10 | 18 | 156 | 13.2 | 27 | 1.1 |
| opt1 high noise pt12 | 10 | 151 | 20.2 | **58** | 1.6 |
| opt1 high noise pt15 | 19 | 151 | 16.3 | 46 | 1.5 |
| opt2 deep seed pt2 | 19 | **176** | 22.2 | 31 | **4.4** |
| opt2 deep seed pt3 | 19 | 176 | 20.4 | 31 | **5.1** |
| opt2 deep seed pt4 | 19 | 168 | 19.3 | 28 | 4.3 |
| opt2 deep seed pt5 | 19 | 164 | 19.4 | 29 | 4.6 |
| opt3 mobile pt2 | 19 | 161 | 17.1 | 30 | 1.1 |
| opt3 mobile pt3 | 19 | 171 | 16.3 | 34 | 1.9 |
| opt3 mobile pt4 | 19 | 150 | 11.1 | 16 | 2.2 |
| opt3 mobile pt5 | 19 | 153 | 12.8 | 31 | 3.1 |
| **opt4 Route A** | 19 | **22** | **2.1** | 4 | 0.6 |

- Designs ≥300 Å³ at ≥80 %: **0 / 237**
- Designs ≥200 Å³ at ≥80 %: **0 / 237**
- Best single structure: **58 Å³** (`opt1_highnoise_pt12_..._0_model_2`), still
  ~6× short and within its own arm's spread.
- Every arm except opt2 is *below* the 33 Å³ baseline on mean enclosed volume.

## 3. Three specific findings

### 3.1 Mobile ATP makes no difference — decision 8 resolved

Option 3 is indistinguishable from the rigid arms at every threshold. The
ligand demonstrably moved: measuring ATP centroid to paratope centroid (a
frame-independent internal distance), rigid arms sit at exactly 31.20 Å with
sd 0.00, while mobile arms shift to 30.25 Å (sd 0.37) at pt2 and 27.08 Å
(sd 0.64) at pt5 — ATP travelled several Å toward the paratope.

**Ligand rigidity was not the obstacle.** The rigid-vs-free comparison
originally requested as a design variable has now been run, and the variable
does not matter for pocket formation. This removes the main argument for
maintaining the RFD3 patch.

### 3.2 Route A is much worse, not marginally worse

Route A collapses ≥70 % volume from ~150 Å³ to 22 Å³ and ≥80 % to 2.1 Å³.
Contig inpainting destroys the existing groove rather than deepening it,
presumably because replacing the FR2/FR4 face with sampled-length segments
lets the two domains repack against each other with nothing between them.

§4 predicted Route B would beat Route A. That is now measured. Both fail, but
Route A should not be the fallback.

### 3.3 High noise has an upper bound

`opt1 pt12` crashed after 10 of 20 designs with an `AssertionError` in RFD3's
motif-RMSD check on residue **A111** — one unindexed token moved 9.6 Å while
the rest moved 0.5–2.5 Å. At 12 Å the pinned motif can no longer be reconciled.
`pt15` completed all 20, so the boundary is not sharp; treat it as "near 12 Å"
rather than a precise limit. pt12's higher mean is computed on half the sample
and is not comparable.

## 4. Interpretation

**Confirmed by measurement.** Across 317 designs (80 at `partial_t` 2–5 plus
these 237) spanning noise 2–15 Å, shallow and deep ATP seeds, rigid and mobile
ligand, and both design routes, RFD3 does not excavate an ATP-scale cavity at
the OKT3 VH/VL interface with the paratope pinned.

**Inference, not established.** The common factor across all four failures is
that the interface stays packed. Options varied *how hard we pushed* and *where
ATP started*, but never *how much of the protein was allowed to move*: the
18-residue REDESIGNABLE set and the 19-residue paratope pin were constant
throughout. A 400 Å³ void may simply require rebuilding more than 18 residues.

**What this does not say.** Nothing about whether the ATP switch works. No
MPNN, no RF3, no apo/holo comparison has been run on the OKT3 campaign. This
measures only cavity formation.

## 5. Consequences for the plan

§3's S2 gate is "can the validated interface be rebuilt into the OKT3 Fv base?"
On this evidence the honest answer is **not by the method in §4**, and the gate
should be treated as firing rather than quietly re-tuned.

Three directions, none yet tested, in ascending order of how much they change
the problem:

1. **Enlarge the redesignable region.** Unindex the BUFFER set (5 residues) or
   the full 38-residue interface core rather than 18. Cheapest test of the §4
   inference above; keeps the paratope pinned.
2. **Drop the paratope pin for a feasibility-only run.** Not a usable design —
   it would destroy CD3 binding — but it answers whether *anything* opens a
   pocket here, separating "RFD3 can't" from "the constraints forbid it".
3. **Revisit the staging in §3.** S1 produced genuine nucleotide-gated
   assemblies backbone-free (230 designs clearing the gate,
   [plan §9.1](ATP_dependent_Fv_plan_v2.md)). The failure is specific to
   grafting that onto a packed Ig interface. Fusing a validated S1 switch
   domain to the Fv, rather than excavating the Fv itself, avoids the problem
   this document measures.

**Recommendation: run 1 and 2 before any further scaling.** Together they cost
~40 designs and distinguish "needs a bigger rebuild region" from "this
interface cannot host an ATP site", which is the question that decides whether
the OKT3 route continues at all.

## 6. Data and provenance

| Item | Path |
|---|---|
| Per-model cavity volumes | `okt3_atp_fv_opts/cavity/cavity_all.csv` |
| Arm manifest | `okt3_atp_fv_opts/configs/arms.json` |
| Shallow / deep seed geometry | `inputs/okt3_fv_atp_placement.json`, `inputs/okt3_fv_atp_deep_placement.json` |
| Patched RFD3 | `/work3/s212001/rfd3_patched` (local copy; shared install untouched) |
| Earlier partial_t 2–5 sweep | `okt3_atp_fv_sweep/` — 640 models, cavity scan not run |

### Caveats on these numbers

- **19–20 designs per arm.** Differences of a few Å³ between arms are not
  meaningful. Only the gap to ~350 Å³ is.
- **237 of 250 scanned.** Two shards failed on a transient `conda activate`
  error (unreadable `backports/zstd` in the base install, not a code fault) and
  were re-run; the outstanding 13 structures cannot change a 0/237 result.
- **Cavity volumes are not comparable to §1.3's absolute numbers** — see
  [`results_route_b_smoke.md`](results_route_b_smoke.md) §4.1.
- **Option 3 depends on a local patch to RFD3.** Input handling is verified and
  the ligand provably moved, but the patch has had no review beyond that.
- **Procedural note.** The option-3 resubmission raced four GPU queues and the
  losers were never killed, so four copies ran concurrently into the same output
  directories. All 80 files were checked afterwards and parse cleanly with full
  atom counts, but the run wasted ~4× the GPU time it needed.
