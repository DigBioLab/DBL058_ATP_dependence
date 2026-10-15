# ATP-Gated Anti-CD3ε Fv — Plan v2

**Date:** 2026-10-15 · **Supersedes:** `ATP_dependent_scFv_plan.md`
**Rev. 2026-10-15:** ATP-vs-ADP selectivity de-scoped by Darian. See §5.
**Grounding:** PDB 1SY6 (measured, this session) · repo `DigBioLab/DBL058_ATP_dependence` @ HEAD (read, this session)

Everything in §1 is measured from the deposited coordinates. §2 onwards separates
measurement from inference and recommendation.

---

## 1. What 1SY6 actually is, and what it measures out to

**Confirmed.** 1SY6 = *Crystal Structure of CD3γε Heterodimer in Complex with OKT3 Fab
Fragment* (2004). Three chains:

| Chain | Content | Residues (auth) |
|---|---|---|
| H | OKT3 Fab heavy — V<sub>H</sub> 1–113, C<sub>H</sub>1 114–219 | 219 |
| L | OKT3 Fab light — V<sub>L</sub> 1–107, C<sub>L</sub> 108–213 | 213 |
| A | CD3γ–ε ectodomain, single chain via a 26-residue linker | 168 modelled, numbered 0–203 |

Two consequences before anything else:

- **This settles clone choice: OKT3.** It is murine. No cyno cross-reactivity, and
  the framework is non-human at exactly the V<sub>H</sub>/V<sub>L</sub> positions this
  project wants to mutate. See §6.
- **The antigen is CD3γε, not εδ or εε.** Epitope and geometry are specific to the
  γε pair.

### 1.1 Paratope — both chains contribute, V<sub>H</sub> leads

Residues within 4.5 Å of chain A:

```
H:  T33  Y50  N52  S54  R55  Y57  N59  K74  Y99  D101 Y104      (11 residues, 38 atoms)
L:  S30  Y31  D49  W90  S91                                      ( 5 residues, 26 atoms)
```

At 5.0 Å, H31, L92 and L93 join. Buried-atom split is 59 % V<sub>H</sub> / 41 % V<sub>L</sub>.

**This is good news for the concept.** L31, L49 and W90 sit at 2.66, 2.81 and 3.10 Å
from CD3 — deep, not peripheral. Neither domain alone should bind CD3 appreciably, so
the apo (unassembled) state has a genuine off-state. The v1 plan assumed this; it is
now measured.

### 1.2 V<sub>H</sub>/V<sub>L</sub> interface, and how much of it is safe to touch

38 interface residues (21 V<sub>H</sub> + 17 V<sub>L</sub>, 4.5 Å). Classified by distance
to CD3:

```
REDESIGNABLE (>=12 A from CD3, 18 residues)
  H37 H39 H43 H44 H45 H46 H95 H109 H110 H111 H112
  L37 L41 L42 L43 L45 L86 L97

BUFFER (8-12 A, 5 residues)   H47 H106 H107 L35 L94

PARATOPE-PROXIMAL (<8 A, 15 residues)  — do not touch
  H35 H50 H59 H99 H103 H104 H105 L31 L33 L48 L49 L88 L90 L93 L95
```

Exported machine-readable to [inputs/okt3_1sy6_residue_sets.json](../inputs/okt3_1sy6_residue_sets.json).
This is the `PARATOPE` / `INTERFACE_CORE` / `REDESIGNABLE` annotation v1 asked for, done.

**Inference:** the redesignable set is the canonical FR2 ledge (H37/H39/H45, L37/L41/L43)
plus the FR4/J-segment face (H109–H112, L97). H39 and L37 are the conserved interface
glutamine pair. Rebuilding exactly this face is the standard lever on V<sub>H</sub>/V<sub>L</sub>
affinity, which is what a switch needs.

### 1.3 The decisive measurement: **there is no pocket to dock into**

I ran a grid cavity scan over the V<sub>H</sub>/V<sub>L</sub> interface base (1 Å voxels,
probe-accessible, ray-occlusion buriedness). Method is crude, so the threshold sweep
matters more than any single number:

| Buriedness threshold | Interpretation | Total volume | Largest single cavity |
|---|---|---|---|
| ≥ 90 % | fully enclosed | 0 Å³ | — |
| ≥ 80 % | enclosed pocket | 49 Å³ | 17 Å³ |
| ≥ 70 % | semi-enclosed | 208 Å³ | 75 Å³ |
| ≥ 60 % | surface groove | 578 Å³ | 339 Å³ |

ATP occupies roughly 350–400 Å³ and a µM-affinity site is typically 400–700 Å³ enclosed.

**Confirmed:** the Fv interface is tightly packed. There is a shallow groove of
about the right footprint at the interface base, 23 Å from CD3, lined by H39, H43,
H93, H95, L37, L38, L41, L42, L86 — i.e. precisely the redesignable set. But there is
no enclosed cavity anywhere near ATP scale.

**Therefore v1's Phase 2 does not work as written.** "Enumerate ATP poses by docking
into the chosen Fv region" will return only shallow, mostly solvent-exposed surface
poses, because that is all the input geometry offers. Docking cannot site a pocket
that does not exist. The pocket has to be *excavated* — backbone removed and rebuilt —
which is a generative task, not a sampling task.

### 1.4 The Fab format reinforces the wrong thing

C<sub>H</sub>1/C<sub>L</sub> bury 25 residues of contact and are covalently locked by an
interchain disulfide **L213–H134** (2.03 Å, confirmed in the coordinates). Internal
disulfides: H22–H96, H146–H201, L23–L87, L133–L193.

An ATP switch at the V<sub>H</sub>/V<sub>L</sub> interface cannot gate assembly of a Fab:
the constant domains hold the two chains together regardless. **Design and screen in
Fv format — H1–113 and L1–107 only, constant domains deleted, no engineered disulfide.**
This confirms v1 §4.2 and now quantifies why.

---

## 2. Two confirmed bugs in the repo, as applied to this campaign

Both are real; I read the current code.

### 2.1 `step5_analyze.py` — the ligand-dependence filter admits inverted switches

```python
df_final = df_merged[
    (df_merged["best_iptm_pp_nol"] < config.IPTM_PP_NOL_CUTOFF) &   # 0.50
    (df_merged["best_iptm_pp"]     >= config.IPTM_PP_CUTOFF)        # 0.20
]
```

With apo = 0.45 and holo = 0.25, both conditions hold and the design passes, despite
assembling *better without* ATP. Replace the independent thresholds with a differential:

```python
df_final = df_merged[
    (df_merged["best_iptm_pp"] - df_merged["best_iptm_pp_nol"] >= 0.30) &
    (df_merged["best_iptm_pp"] >= 0.60)
]
```

### 2.2 `fill_fixed_residues.py` — will delete the paratope

It rewrites the MPNN flag with motif residues only:

```python
new_line = re.sub(r'--fixed_residues "[^"]*"', f'--fixed_residues "{fixed_res}"', line)
```

With `MPNN_CHAINS_TO_DESIGN = "A,B"`, every non-motif position is designable, CDRs
included. It must emit `motif_residues ∪ PARATOPE_5A` from the JSON in §1.2.

### 2.3 Chain-ID collision, minor but will bite

The RFD3 campaign uses chains A/B. In 1SY6, **A is CD3**. Rename on input: V<sub>H</sub>→A,
V<sub>L</sub>→B, CD3→C. Do it in the input-prep step, not by hand.

---

## 3. The honest strategic problem

The project stacks three independently hard requirements:

1. Excavate a 400 Å³ pocket in a tightly packed Ig interface (§1.3).
2. Make assembly nucleotide-dependent — **negative design**, which RFD3 is not trained for.
3. Keep the CD3 paratope intact through that rebuild.

A fourth, discriminating ATP from ADP, has been de-scoped (§5). Each of the remaining
three is a plausible-but-uncertain research project. Attempting all three in one campaign
means a failure tells you nothing about which of the three broke.

**Recommendation: stage it, with a kill gate after each stage.** This is the main
structural change from v1.

| Stage | Question | Antigen involved? | Gate to proceed |
|---|---|---|---|
| **S1** | Can this pipeline produce *any* nucleotide-gated two-chain assembly, backbone free? | No | ≥10 designs with holo−apo ipTM ≥0.3, and ≥1 confirmed wet-lab |
| **S2** | Can the validated interface be rebuilt into the OKT3 Fv base? | Yes | CDR RMSD <1.0 Å, CD3 ipTM retained |
| **S3** | Format, developability, function | Yes | Jurkat-NFAT ATP dependence |

Three stages, not four: the former selectivity stage collapses into S1 (§5).
S1 is largely the existing `8btg_D501` campaign — `top9fixed_structures/` already exists.
**Before writing any new code, pull the realised holo/apo ipTM distribution from that
run.** If it never separated, S2 is premature and the whole plan should change. That is
the cheapest informative action available and it should happen first.

---

## 4. Revised design strategy for S2 (replaces v1 Phase 2)

Drop docking-based pocket siting. Two generative routes, run in parallel.

**Route A — interface-base rebuild by inpainting.** Keep CDRs and paratope indexed,
replace the redesignable FR2/FR4 face with unindexed segments, place ATP as ligand,
let RFD3 build a composite site. Starting contig, 1SY6 auth numbering after the
A/B rename of §2.3:

```python
# A = VH (1-113), B = VL (1-107). Unindexed lengths replace the redesignable face.
CONTIG = "A1-36,4-8,A47-94,6-10,A105-113,/0,B1-36,6-10,B46-85,4-8,B87-96,4-8,B98-107"
LENGTH = "210-240"
REDESIGN_MOTIF_SIDECHAINS = False      # paratope identities frozen
LIGAND_NAME = "ATP"
```

Note the `/0` break is doing real work here: V<sub>H</sub> and V<sub>L</sub> are genuinely
separate chains in Fv format, which is also what makes `iptm_protein_protein` meaningful.

**Route B — partial diffusion of the whole Fv with ATP present.** Low-noise diffusion
over both domains, paratope restrained. Given §1.3, letting the backbone breathe around
the ligand is more likely to open a real cavity than strict segment replacement. Costs
more CDR drift; the RMSD metric catches it.

**Expect Route B to outperform Route A**, because the required deformation is distributed
(opening a packed interface) rather than local. That is an inference, not a result; run
both and let the numbers decide.

**Deliberate apo destabilisation.** RFD3 will not do negative design, but you can bias
it: require the designed pocket to bury ≥2 unsatisfied polar groups or a charged residue
that is only neutralised by the ATP phosphates. Score this post hoc and enrich, rather
than trying to condition on it.

---

## 5. Nucleotide selectivity — de-scoped

**Decision (Darian, 2026-10-15): ATP-vs-ADP discrimination is not a requirement.**
Recorded as a scope decision, not a finding.

It is defensible on the biology. CD39 and CD73 make the TME rich in ADP, AMP and
adenosine *as well as* ATP, so the whole series is elevated together. What the switch
needs to read is total extracellular adenine nucleotide against healthy interstitium,
not the identity of one member. Accepting ADP cross-reactivity therefore matches the
biology rather than conceding to it.

Three consequences, all favourable:

- **The hardest sub-problem disappears.** Reading the γ-phosphate specifically would
  have meant building a charged subsite that is penalised when empty — negative design
  nested inside the negative design of §3. That requirement is gone.
- **The usable pharmacophore library grows substantially.** Most natural ATP sites bind
  adenine and ribose and tolerate phosphate count. The repo's PLIP fingerprint clusters
  can now be used broadly instead of filtered down to γ-phosphate readers.
- **A stage and a screening round disappear** (§3).

The affinity target loosens too. If ATP, ADP and AMP all contribute, the effective
ligand pool in the TME exceeds the ATP concentration alone, so a weaker site suffices.
**Suggest K_d 50–200 µM against ATP**, against v1's 1–20 µM.

### 5.1 Still worth measuring, as tracked metrics rather than gates

- **Adenosine.** The one member of the series worth watching. It circulates in plasma
  and normal tissue at appreciable levels and it is the smallest fragment, so a design
  that binds adenosine alone tightly has the least tumour contrast. Record
  `iptm_protein_ligand` for ADN; do not gate on it.
- **GTP.** Cheap, and tells you whether the site reads the base at all. Informative
  about design quality rather than tumour selectivity. Record it.

Run these with CCD codes directly (`RF3_CCD_LIGAND = True`): ATP, ADP, AMP, ADN, GTP.
The current config uses an SDF path; CCD codes are simpler for a panel.

Excluding catalytic Mg·ATP sites from the library still stands, and is now the *only*
remaining filter on it (v1 Phase 1). Grafting turnover machinery risks building a slow
ATPase that consumes its own ligand. Flag Mg-dependence per cluster from the PLIP
metal-coordination records.

### 5.2 Concentration targets — weak sources, flag in any grant text

Reported TME extracellular ATP is in the hundreds of µM against nanomolar-to-low-µM in
healthy interstitium. These values come largely from pmeLUC-type reporter probes, are
contested, and are method-dependent. ATP is also pericellular rather than bulk. Any
grant text should say plainly that this narrows rather than eliminates on-target
off-tumour risk, since adenine nucleotides are released by inflammation, trauma,
hypoxia and exercise.

One off-target is worth naming now that ADP is in scope: **platelet activation releases
high local ADP** at sites of clotting and vascular injury. Stating it in the project
framing is the right handling; it is not a reason to reinstate the selectivity
requirement.

## 6. OKT3-specific issues v1 did not have

- **Murine framework.** The redesignable set in §1.2 *is* framework, and it is mouse
  framework. Mutating it compounds an existing immunogenicity liability.
  **Recommendation:** build on a humanised OKT3 V-region (teplizumab lineage) rather
  than the 1SY6 sequence, using 1SY6 only for geometry. Verify the humanised sequence
  from a primary source — do not trust recalled residue numbers. If humanisation changes
  FR2/FR4 positions, those are the exact positions being redesigned, so do it *before*
  designing, not after.
- **Track germline deviation and mutation count** as developability metrics, capped.
- **CD3 agonism.** OKT3 is a potent agonist; an ATP-gated version still carries
  cytokine-release risk when it does engage. Relevant to how S4 is framed, not to design.

---

## 7. Carried over from v1, unchanged and still correct

- Indexed = pinned motif, unindexed = sampled-length linker. v1 §4.1's correction stands.
- PLIP for *extracting* contacts (repo already does this); PLACER for *evaluating*
  designed pockets. Different tools, different stages.
- Functional-atom pinning via `FIXED_ATOMS` rather than rigid backbone grafting.
- The multi-state RF3 evaluation, now three gating states — holo, apo, holo+CD3ε —
  plus a nucleotide panel that is recorded but not gated (§5.1).
- ipTM is a confidence metric, not a ΔG. Add a physics-based ΔΔG on the top ~50.
- Budget 10³–10⁴ designs, not 100 × 10. Rejection-based negative design is inefficient.
- Experimental: ATPγS / AMP-PNP alongside ATP to separate binding from hydrolysis,
  apyrase as negative control, and watch ATP instability in culture medium.

---

## 8. Open decisions, updated

| # | Decision | Status after this re-evaluation |
|---|---|---|
| 1 | Anti-CD3ε clone | **Closed — OKT3**, by availability of 1SY6. Open sub-question: murine vs humanised V-region (§6). |
| 2 | Molecular format | **Narrowed — Fv for design and screening** (§1.4). Final therapeutic format still open. |
| 3 | Pocket location | **Narrowed — interface base**, the 18-residue redesignable set (§1.2). Interface core is ruled out: it is paratope-proximal. |
| 4 | One ATP site or two | **Recommend one.** §1.3 shows accommodating even one requires real excavation. |
| 5 | *(new)* Stage or go direct | **Recommend staged**, S1→S3 (§3). |
| 6 | *(new)* ATP vs ADP selectivity | **Closed — not required** (§5). Collapses a stage and loosens the affinity target. |

## 9. Immediate next actions, in order

1. **Pull holo/apo ipTM from the existing `8btg_D501` run.** Cheapest decision-relevant
   data that exists. Everything downstream is conditional on it.
2. Patch `step5_analyze.py` (§2.1) and re-score that run with the differential filter —
   the current numbers may be inflated by the inverted-switch admission.
3. Patch `fill_fixed_residues.py` to emit `motif ∪ paratope` (§2.2).
4. Build the Fv input: strip C<sub>H</sub>1/C<sub>L</sub>, rename chains, emit `okt3_fv.cif`.
5. Write `step6_multistate.py` for the five-state panel, CCD ligands.
6. Only then run S2 design.
