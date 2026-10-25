# Branch `feature/darwol-atp-gated-cd3-fv`

A different design task from the one this repo was built for.

**Original campaign (`8btg_D501`, on `main`):** design a two-chain assembly around ATP
with a free backbone. Nothing is held fixed except the grafted ATP motif.

**This branch:** excavate an ATP site into the base of an *existing* antibody VH/VL
interface, so that the OKT3 anti-CD3 paratope only forms when the two chains assemble
in the presence of extracellular nucleotide. The paratope is held fixed throughout.

The machinery is shared; the constraints are not.

## What is here

| Path | Purpose |
|---|---|
| `docs/ATP_dependent_Fv_plan_v2.md` | The plan, grounded in measurements from PDB 1SY6 |
| `docs/ATP_dependent_scFv_plan.md` | The superseded v1 plan, kept for provenance |
| `scripts/prepare_okt3_fv.py` | Builds the Fv inputs and residue annotation from 1SY6 |
| `configs/okt3_atp_fv.py` | Campaign config. Unvalidated starting point |

## Start here

```bash
python3 scripts/prepare_okt3_fv.py --out-dir inputs/
cp configs/okt3_atp_fv.py pipeline/config.py
python3 pipeline/step0_setup.py
```

`inputs/` is gitignored. It is regenerated from 1SY6 by the prep script, so nothing
structural needs to be tracked.

## Three findings that drive the design

Measured from 1SY6, detail in the plan.

1. **No pre-formed pocket.** The VH/VL interface holds at most ~17 Å³ of enclosed
   cavity against the ~400 Å³ ATP needs. The site must be generated, not docked into.
   This is why `FIXED_ATOMS` is empty and no natural ATP motif is grafted.
2. **The constant domains must go.** CH1/CL bury 25 contact residues and are locked by
   the L213–H134 interchain disulfide, so a Fab cannot be assembly-gated. Design runs
   in Fv format.
3. **Both domains contribute to the paratope**, 59 % VH / 41 % VL by buried atoms, with
   three light-chain contacts under 3.2 Å. Neither half alone should bind CD3, so the
   off-state is genuine.

## Scope decision

ATP-versus-ADP discrimination is **not** a requirement. CD39 and CD73 make the tumour
microenvironment rich in ADP, AMP and adenosine alongside ATP, so the switch reads total
extracellular adenine nucleotide. `RF3_TRACKED_LIGANDS` are recorded, not gated.

## Caveats

- `CONTIG`, `LENGTH` and every threshold in the config are unrun first guesses.
- OKT3 is murine, and the redesign target is framework. Consider rebasing on a
  humanised V-region before committing to a sequence.
- Per DBL convention a genuinely separate project would get its own `DBLxxx` number and
  repo. This lives as a branch by request; worth revisiting if it outgrows the shared
  pipeline.
