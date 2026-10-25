> **SUPERSEDED on 2026-10-15 by [`ATP_dependent_Fv_plan_v2.md`](ATP_dependent_Fv_plan_v2.md).**
>
> Kept for provenance: v2 cites this document, and several of its conclusions were
> overturned by measurements on PDB 1SY6. In particular §5 Phase 2 below proposes
> siting the ATP pocket by docking into the V<sub>H</sub>/V<sub>L</sub> cleft. A cavity
> scan later showed no enclosed ATP-scale pocket exists there, so that approach was
> replaced by generative excavation. §5 Phase 5's ADP/adenosine counter-screen was also
> de-scoped. Read v2 for current plans; read this only for history.

# ATP-Dependent Anti-CD3ε Fv — Conversation Record & Project Plan

**Date:** 2026-10-15
**Participants:** Darian (DigBioLab) · Claude
**Reference repo:** `DigBioLab/DBL058_ATP_dependence`

---

## 1. Background question — unindexed scaffolding in RFdiffusion 3

**Q:** What is "unindexed scaffolding"?

"Unindexed" scaffolding is the part of a designed chain that has **no correspondence to any input structure** — the model invents it, and it is given no residue numbering from a template.

### Indexed regions

Residues pinned to something real. You give a chain and a residue range from the input structure (e.g. `A30-45`), and those residues carry their original numbering and geometry into the output. The model is told *"these atoms exist, keep them."* Used for motifs, binding sites, and targets you design against.

### Unindexed regions

Specified only by a **length** or length range (e.g. `10-25`). There is no source index because there is no source. The model samples a length from that range, then generates backbone coordinates with no template constraint. Numbering is assigned after the fact when the output chain is written. These are the linkers, loops and helical arms — the scaffold that holds the indexed motif in place.

### Why the framing matters

- **Length is a sampled variable, not a fixed input.** Because the segment has no index to anchor it, residue indices of everything downstream shift between designs. You cannot write a constraint like "residue 57 must contact the ligand" across a batch unless 57 sits in an indexed block.
- **It is what makes motif scaffolding work.** Indexed residues are seen as having real sequence separation; unindexed stretches get indices that are internally consistent but carry no claim about input geometry. That is how two motif fragments from unrelated proteins can be joined, with the model deciding how much backbone to put between them and in what conformation.
- **Chain breaks vs. linkers.** A large index jump (or an explicit chain-break token, `/0` in the DBL058 contigs) tells the model two segments need not be covalently continuous. An unindexed segment of length *n* tells it they **must** be joined by *n* residues. Confusing the two is the usual cause of *"why did it build a linker I didn't want."*

> **Version caveat:** RFdiffusion3's input format differs from the original contig-string syntax. Flag names and schema keys below are taken from the DBL058 repo's actual working configs rather than from memory.

---

## 2. The proposed project

**Goal:** integrate one or two ATP binding sites into an anti-CD3ε scFv so that the paratope only forms when heavy and light chain assemble — giving tumour-microenvironment specificity, since extracellular ATP is elevated in the TME relative to healthy tissue.

**Originally proposed workflow:**

1. Download all PDB structures containing ATP.
2. Process with PLACER to extract binding residues + ATP molecules.
3. Fit them onto the scFv surface, focusing on the V<sub>H</sub>/V<sub>L</sub> interface, via a custom geometric-matching script.
4. Visualise best matches in PyMOL and make a selection.
5. For selected matches, define the ATP motif as unindexed and redesign that part of the scFv without altering the CD3ε binding site.
6. Compute ipTM (H↔L, and each towards the ligand) and RMSD.

---

## 3. What the existing repo already provides

The `DBL058_ATP_dependence` repo contains roughly 80% of the needed machinery. The `conditioned_assembly` campaign is the right precedent.

### Repo structure

```
run_pipeline.sh              # orchestrator: RFD3 → MPNN → RF3 → RF3-nol → analysis
pipeline/config.py           # campaign parameters
pipeline/step0_setup.py      # writes RFD3 config JSON + LSF submit script
pipeline/step2_process_rfd3.py
pipeline/step3_process_mpnn.py
pipeline/step4_process_rf3.py   # metrics, shortlisting, no-ligand config generation
pipeline/step5_analyze.py       # ligand-dependence filter
scripts/fill_fixed_residues.py  # populates MPNN --fixed_residues from motif atoms
scripts/03_run_plip.py          # PLIP interaction profiling
scripts/atp_conformer_analysis.py
submit/01_download.sh           # ~3,600 ATP-containing PDB entries
submit/02_rmsd_array.sh, 02_merge.sh, 03_plip.sh
```

### Key existing config (campaign `8btg_D501`)

```python
CONTIG = "20-30,D151-D159,3-8,D120-123,20-30,/0,20-30,D287-D289,20-30,D312-D314,10-20,C259-C265,20-30"
LENGTH = "136-203"

FIXED_ATOMS = {            # pin functional side-chain atoms, not whole backbone
    "D121": "N,O",   "D157": "N,NZ",   "D313": "NH1,NH2",
    "C260": "NH1,NH2", "C264": "NH1,NH2", ...
}

MPNN_CHAINS_TO_DESIGN = "A,B"
IPTM_PL_CUTOFF     = 0.80   # protein–ligand, with ligand
IPTM_PP_CUTOFF     = 0.20   # protein–protein, with ligand
IPTM_PP_NOL_CUTOFF = 0.50   # reject if it assembles without ligand
```

Already in place and directly reusable:

- Two-chain architecture via the `/0` chain break, with ATP contacts split across chains.
- Functional-atom pinning rather than rigid backbone grafting (`FIXED_ATOMS`).
- LigandMPNN for sequence design (ligand-aware, unlike plain ProteinMPNN).
- The **apo/holo rejection filter** — designs that assemble fine without ligand are discarded (`step5_analyze.py`).
- PLIP-based interaction profiling and conformer clustering over the whole ATP-containing PDB.

---

## 4. Four corrections to the premise

### 4.1 The indexed/unindexed assignment is inverted

Step 5 of the proposal says "define the ATP binding motif as unindexed motif." It is the opposite. The ATP-contacting residues are the **indexed** motif, pinned to input coordinates (as `D151-D159` etc. already are). The **unindexed** segments are the bare length ranges (`20-30`) that RFD3 invents to connect them.

For the scFv case: CDRs indexed, ATP-contact residues indexed, short unindexed stretches replacing selected framework segments at the V<sub>H</sub>/V<sub>L</sub> interface base.

### 4.2 An scFv cannot be a true assembly switch

V<sub>H</sub> and V<sub>L</sub> joined by a (G₄S)₃ linker associate *intramolecularly*. Effective concentration is in the millimolar range — roughly +4 to +6 kcal/mol of pairing bias that an ATP site would have to overcome. The apo form will still pair.

Also, `iptm_protein_protein` — the entire discriminating metric in `step5_analyze.py` — is undefined for a single chain.

**Therefore:** design and screen in **two-chain Fv format** (separate V<sub>H</sub> and V<sub>L</sub> chains, as the `/0` break already does), then decide the final molecule — a dsFv/Fab-like two-chain arm, or a split engager where each half carries one V domain. If a single-chain product is required, treat the switch as conformational/affinity-modulating rather than assembly-gating and expect a much smaller dynamic range.

### 4.3 PLIP ≠ PLACER

The proposal says PLACER; the repo uses **PLIP** (`scripts/03_run_plip.py`), which is the right tool for *extracting* contacts. **PLACER** is a different tool — a ligand-environment/pose plausibility network — and belongs in **evaluation**, as a cheap orthogonal filter catching designed pockets that cannot hold ATP in a sane pose. Use both, at different stages.

### 4.4 Do not write a bespoke geometric matcher

Step 3 as specified is the weakest link and largely unnecessary. Rigid multi-residue motif grafting onto an Ig β-sandwich will almost always fail — the fold will not accommodate a P-loop or Rossmann-style site. Two better framings, to run in parallel:

- **Track A (preferred) — ligand placement, not motif matching.** Reduce the problem from "superimpose a 6-residue motif" to "sample ATP's 6 rigid-body DOF in the interface cleft," then let RFD3 build the pocket. This is what the 8BTG campaign does, minus the fixed scaffold.
- **Track B — pharmacophore grafting.** Keep only *functional-atom positions* from natural sites, not backbone. `FIXED_ATOMS` already does this. Push further: discard backbone motif entirely, keep ATP plus a constellation of required functional-group positions, and let RFD3 + LigandMPNN find residues to satisfy them.

---

## 5. Proposed plan

### Phase 0 — Inputs and decisions (days)

- Pick the anti-CD3ε clone by whether a **co-crystal with CD3ε exists**. UCHT1, SP34 and OKT3 all have deposited structures — verify IDs rather than trusting memory. SP34 is preferable if cyno cross-reactivity matters for downstream tox work.
- Build the two-chain Fv input: strip the linker, assign V<sub>H</sub> = chain A, V<sub>L</sub> = chain B.
- Annotate three residue sets in a JSON the whole pipeline reads:
  - `PARATOPE` — CDRs by IMGT/Kabat, plus any framework residue within 5 Å of CD3ε
  - `INTERFACE_CORE` — V<sub>H</sub>/V<sub>L</sub> contacts
  - `REDESIGNABLE` — framework, interface-proximal, ≥8 Å from CD3ε
- **Decide where the pocket goes.** The paratope sits at the top of the V<sub>H</sub>/V<sub>L</sub> interface; CDR-H3 and L3 pack *into* that interface.
  - Pocket in the interface core → maximum coupling to assembly, maximum risk to CDR geometry.
  - Pocket at the base (former C<sub>H</sub>1/C<sub>L</sub> face) → safer paratope, weaker coupling.
  - Sample both; let the CDR-RMSD metric adjudicate.

### Phase 1 — ATP site library (reuse, ~1 day)

Nothing new to compute. Take from DBL058: `pdb_cache/`, `output/atp_conformers_raw.sdf`, the Butina conformer clusters, and the PLIP fingerprint clusters.

Add one filter the repo lacks:

- **Exclude catalytic sites.** Most PDB ATP sites are Mg·ATP with turnover machinery; grafting one risks building a slow ATPase that destroys its own ligand. Prefer non-catalytic, allosteric or transport sites. Flag Mg-dependence per cluster from the PLIP metal-coordination records.
- Reduce each retained cluster to a pharmacophore: functional-atom positions + ATP conformer, exported as a CST-like spec.

### Phase 2 — Pocket siting (replaces proposed step 3)

- Enumerate candidate ATP poses in the chosen Fv region by docking (GNINA/Vina for pose sampling, or RifGen if Rosetta is available). Only plausible *starting* placements are needed, not accurate affinity.
- Score candidates on: buriedness, clash with paratope, distance to nearest CDR atom, and how many contacts come from each chain.
- **Hard filter: contacts must be split across both chains** (e.g. ≥30% of contact surface from each). A pocket satisfied entirely by V<sub>H</sub> gates nothing.
- Emit one input CIF per surviving pose: Fv + ATP in the same frame, motif flags set.

### Phase 3 — Visual triage

Keep this, but gate it behind Phase 2 scores so ~20–50 poses are reviewed, not thousands. The repo's PyMOL session export already does the heavy lifting.

### Phase 4 — Design (RFD3 → LigandMPNN)

Reuse `run_pipeline.sh` wholesale, with these config changes:

```python
# Scaffold-conditioned inpainting: CDRs + ATP contacts indexed,
# short unindexed stretches replace interface-base framework.
CONTIG = "A1-37,6-10,A44-97,/0,B1-33,5-9,B39-88,8-12,B95-107"
REDESIGN_MOTIF_SIDECHAINS = True   # non-CDR motif positions may change identity
```

Two notes:

- Also try **partial diffusion** — whole Fv diffused at low noise with ATP present — as an alternative to hard inpainting. Letting the backbone relax around the ligand usually succeeds more often than strict segment replacement, at the cost of more CDR drift, which the RMSD metric will catch.
- **`fill_fixed_residues.py` needs extending.** It currently *overwrites* `--fixed_residues` with motif residues only:

  ```python
  new_line = re.sub(r'--fixed_residues "[^"]*"', f'--fixed_residues "{fixed_res}"', line)
  ```

  With `MPNN_CHAINS_TO_DESIGN = "A,B"`, LigandMPNN will happily redesign the CDRs away. It must emit `motif_residues ∪ PARATOPE`.

### Phase 5 — Multi-state evaluation (the real work)

The proposal's step 6 has two states; this needs five. For each design, run RF3 on:

| State | Purpose | Want |
|---|---|---|
| V<sub>H</sub> + V<sub>L</sub> + ATP | assembly, holo | high `iptm_pp`, high `iptm_pl` |
| V<sub>H</sub> + V<sub>L</sub> (apo) | assembly, apo | **low** `iptm_pp` |
| V<sub>H</sub> + V<sub>L</sub> + ATP + CD3ε | paratope intact | high ipTM to CD3ε, low CDR-RMSD vs parent |
| V<sub>H</sub> + V<sub>L</sub> + ADP / AMP / adenosine | TME selectivity | low `iptm_pl` |
| V<sub>H</sub> + V<sub>L</sub> + GTP | nucleotide specificity | low `iptm_pl` |

Three substantive additions:

#### Fix the ligand-dependence filter

With `IPTM_PP_CUTOFF = 0.20` and `IPTM_PP_NOL_CUTOFF = 0.50`, a design with apo `iptm_pp` = 0.45 and holo = 0.25 **passes** — the switch runs backwards. Replace the two independent thresholds with a differential:

```python
df_final = df_merged[
    (df_merged["best_iptm_pp"] - df_merged["best_iptm_pp_nol"] >= 0.30) &
    (df_merged["best_iptm_pp"] >= 0.60)
]
```

#### The ADP/adenosine counter-screen is not optional

The TME is CD39/CD73-rich; extracellular ATP is hydrolysed to ADP → AMP → adenosine, all of which accumulate there. A site that binds ADP equally well gives zero tumour selectivity — and most natural ATP sites do exactly that. This is the likeliest failure mode of the whole concept. Test it early, before the CD3ε state, as a cheap kill filter.

#### ipTM is a confidence metric, not a ΔG

The quantity that matters is ΔΔG_assembly between apo and holo; ipTM is a weak proxy. Add:

- **PLACER** on final pockets (pose plausibility)
- A **physics-based ΔΔG** on the top ~50 (Rosetta, or MM-GBSA / short MD)
- Target an ATP K_d around **1–20 µM**, so the switch is off at normal interstitial ATP and on at reported TME concentrations.

Also track two developability metrics, since germline framework at a V<sub>H</sub>/V<sub>L</sub> interface is being mutated: **germline deviation** (immunogenicity risk) and **mutation count** (cap it).

### Phase 6 — Experimental readout

- Display or express top designs.
- Titrate ATP, with ADP/AMP/adenosine counter-titrations.
- Use **ATPγS or AMP-PNP** alongside ATP to separate binding from hydrolysis; apyrase as a negative control.
- Functional readout: Jurkat-NFAT reporter ± ATP.
- Note ATP's instability in culture medium — this bites people.

---

## 6. Two honest caveats

**Negative design is the hard part.** RFD3 is trained for positive design; nothing in it destabilises the apo state. The repo's rejection-based approach (generate, then discard designs that assemble without ligand) is sound but inefficient — expect low single-digit percent pass rates, so budget 10³–10⁴ designs, not the current 100 × 10. Odds can be biased by deliberately designing an unsatisfied polar core or cavity into the apo interface.

**ATP is not tumour-exclusive.** It is released by inflammation, trauma, hypoxia and exercise. This narrows on-target off-tumour toxicity rather than eliminating it — worth stating plainly in the project framing.

---

## 7. Open decisions

1. **Final molecular format** — two-chain Fv / dsFv, Fab-like arm, split engager, or single-chain with a weakened switch?
2. **Anti-CD3ε clone** — UCHT1, SP34 or OKT3 (driven by available co-crystal and cyno cross-reactivity needs).
3. **Pocket location** — interface core (strong coupling, CDR risk) vs. interface base (safe paratope, weak coupling), or both.
4. **One ATP site or two** — two sites multiply the switch cooperativity but roughly square the design difficulty.

---

## 8. Next deliverables discussed

- Patched `config.py` for the scFv/Fv campaign
- Extended `fill_fixed_residues.py` (motif ∪ paratope)
- New `step6_multistate.py` running the five-state evaluation
