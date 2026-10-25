"""
Campaign config: ATP-gated anti-CD3 OKT3 Fv.

Different task from the repo's original `8btg_D501` campaign. There, a two-chain
assembly was designed around ATP with a free backbone. Here, an ATP site is excavated
into the base of an existing antibody VH/VL interface while the CD3 paratope is held
fixed.

To activate:
    cp configs/okt3_atp_fv.py pipeline/config.py      # pipeline does `import config`
or point PYTHONPATH at configs/ before running step0_setup.py.

Prerequisite:
    python3 scripts/prepare_okt3_fv.py --out-dir inputs/

STATUS: unvalidated starting point. CONTIG, LENGTH and the thresholds below are
first guesses from the 1SY6 geometry analysis in docs/ATP_dependent_Fv_plan_v2.md.
They have not been run. Treat every value marked TUNE as a parameter to sweep.
"""
import json
import os

# ── Infrastructure ────────────────────────────────────────────────────────────
PIPELINE_DIR   = os.path.dirname(os.path.abspath(__file__))
REPO_DIR       = os.path.dirname(PIPELINE_DIR)
SCRIPTS_DIR    = os.path.join(REPO_DIR, "scripts")
PYTHON_EXEC    = os.environ.get("ATP_PYTHON_EXEC", "python3")
# jupyter_utils / rf3_sm_metrics live here, imported via sys.path by steps 0,3,4,5,6.
# An empty default made `import jupyter_utils` fail with a confusing traceback.
ENZDES_SCRIPTS = os.environ.get(
    "ENZDES_SCRIPTS", "/work3/s212001/DBL_secrets/enzdes_pipeline-main")

# ── Campaign ──────────────────────────────────────────────────────────────────
user     = os.environ.get("DBL_USER", "s212001")
family   = "atp_binding"
target   = "conditioned_assembly"
campaign = "okt3_atp_fv"

WORKING_DIR = os.environ.get(
    "ATP_WORKING_DIR",
    f"/work3/{user}/projects/enzyme_design/{family}/{target}/{campaign}",
)

# ── Derived directories (do not edit) ────────────────────────────────────────
CMDS_DIR          = f"{WORKING_DIR}/cmds"
LOGS_DIR          = f"{WORKING_DIR}/logs"
PICS_DIR          = f"{WORKING_DIR}/pics"
CONFIGS_DIR       = f"{WORKING_DIR}/configs"
SCORES_DIR        = f"{WORKING_DIR}/scores"
INPUTS_DIR        = f"{WORKING_DIR}/inputs"
DIFFUSION_OUT_DIR = f"{WORKING_DIR}/diffusion_out"
RF3_OUT_DIR       = f"{WORKING_DIR}/rf3_out"
MPNN_OUT_DIR      = f"{WORKING_DIR}/mpnn_out"
SHORTLISTED_DIR   = f"{WORKING_DIR}/shortlisted"

# ── Residue annotation (from scripts/prepare_okt3_fv.py) ─────────────────────
RESIDUE_SETS_JSON = os.path.join(REPO_DIR, "inputs", "okt3_1sy6_residue_sets.json")


def _load_sets():
    if not os.path.exists(RESIDUE_SETS_JSON):
        raise SystemExit(
            f"missing {RESIDUE_SETS_JSON}\n"
            "run: python3 scripts/prepare_okt3_fv.py --out-dir inputs/"
        )
    with open(RESIDUE_SETS_JSON) as f:
        return json.load(f)


def _to_design_chain(labels):
    """'H37' -> 'A37', 'L41' -> 'B41' (pipeline chain IDs)."""
    return [("A" if s[0] == "H" else "B") + s[1:] for s in labels]


_SETS = _load_sets()

# Held fixed through LigandMPNN: the CD3 paratope at the wider 5.0 A cut.
PARATOPE_RESIDUES = _to_design_chain(_SETS["PARATOPE_WIDE"])
# The VH/VL interface face that may be rebuilt: >=12 A from CD3.
REDESIGNABLE_RESIDUES = _to_design_chain(_SETS["REDESIGNABLE"])

# ── RFD3 parameters ───────────────────────────────────────────────────────────
# Built by scripts/prepare_okt3_fv.py into the repo's inputs/, not WORKING_DIR —
# same convention as RESIDUE_SETS_JSON above. inputs/ is gitignored and regenerates.
#
# MUST be the ATP-seeded file, for BOTH routes. RFD3 extracts the ligand from the
# input structure rather than fetching it from the CCD, so an apo input fails with
# "Invalid contig format: 'ATP'". Regenerate with:
#     python3 scripts/prepare_okt3_fv.py --out-dir inputs/ --with-atp
INPUT_PDB = os.path.join(REPO_DIR, "inputs", "okt3_fv_atp.cif")
DIALECT   = 2

# A = VH (1-113), B = VL (1-107). Unindexed length ranges replace the FR2 ledge and
# the FR4/J face, i.e. the REDESIGNABLE set. CDRs and paratope stay indexed.
# TUNE: segment boundaries and lengths are the main lever on pocket volume.
CONTIG = "A1-36,4-8,A47-94,6-10,A105-113,/0,B1-36,6-10,B46-85,4-8,B87-96,4-8,B98-107"
LENGTH = "210-240"

LIGAND_NAME               = "ATP"
REDESIGN_MOTIF_SIDECHAINS = False   # paratope identities frozen
IS_NON_LOOPY              = True
INFER_ORI_STRATEGY        = None
HOTSPOTS                  = None

# No natural ATP motif is grafted: there is no pre-formed pocket to pin to
# (docs/ATP_dependent_Fv_plan_v2.md section 1.3). The site is generated de novo, so
# FIXED_ATOMS stays empty until a pharmacophore is chosen in a later round.
FIXED_ATOMS = {}

RFD3_CKPT            = "/dtu/projects/dbl/foundry/ckpt/rfd3_latest.ckpt"
RFD3_DIFFUSION_BATCH = 10
RFD3_N_BATCHES       = 100   # TUNE: negative design has a low pass rate; budget 1e3-1e4

# ── Route B: partial diffusion ───────────────────────────────────────────────
# Low-noise partial diffusion of the whole Fv with ATP present. Expected to open
# a cavity more reliably than segment replacement, because the deformation needed
# is distributed rather than local (plan section 4).
#
# When True, step0_setup.py ignores CONTIG/LENGTH — RFD3 skips contig sampling
# under partial diffusion — and instead unindexes REDESIGNABLE_RESIDUES while
# pinning PARATOPE_RESIDUES. One RFD3 job per value in PARTIAL_T_SWEEP.
PARTIAL_DIFFUSION = False

# partial_t is the standard deviation of re-added noise in ANGSTROMS, not a 0-1
# fraction. RFD3 recommends <= 15 and says to start small and raise, because the
# response is nonlinear. The previous value here (0.3) was a near no-op.
PARTIAL_T_SWEEP = [2.0, 3.0, 4.0, 5.0]

# ATP MOBILITY IS NOT CONFIGURABLE under partial diffusion in this RFD3 build.
# input_parsing.py subsets the input to protein (line ~614) before applying
# select_fixed_atoms, then re-appends the ligand with
# is_motif_atom_with_fixed_coord hard-coded to 1 (line ~678). ATP is therefore
# always rigid here. A rigid-vs-free comparison has to come from Route A, or
# from patching RFD3.

# ── MPNN parameters ───────────────────────────────────────────────────────────
MPNN_SEED             = 42
MPNN_BATCH_SIZE       = 10
MPNN_N_BATCHES        = 1
MPNN_CKPT             = "/dtu/projects/dbl/LigandMPNN/25.3.1-0//model_params/ligandmpnn_v_32_010_25.pt"
MPNN_CHAINS_TO_DESIGN = "A,B"
MPNN_REMOVE_BAD       = True
# Consumed by scripts/fill_fixed_residues.py, which unions this with the motif
# residues. Without it LigandMPNN redesigns the CDRs away.
MPNN_EXTRA_FIXED_RESIDUES = PARATOPE_RESIDUES

# ── RF3 parameters ────────────────────────────────────────────────────────────
RF3_CKPT      = "/dtu/projects/dbl/rf3/release/ckpt/rf3_latest.pt"
RF3_NUM_STEPS = 50
RF3_N_MODELS  = 5

RF3_CCD_LIGAND = True      # CCD codes are simpler than an SDF for a nucleotide panel
RF3_SDF_PATH   = None
RF3_SMILES     = None

RF3_QUEUE = "cabl40s"   # GPU queue for RF3 array jobs

# Multi-state panel. Only ATP gates; the rest are recorded for context.
# ADP cross-reactivity is acceptable by design decision (plan section 5).
RF3_GATING_LIGAND   = "ATP"
# DEFERRED 2026-10-21: the tracked nucleotide panel is out of scope for now.
# step6_multistate.py scores three states only — holo, apo, holo+CD3. Left
# declared so re-enabling it is a code change in one place, not a config hunt.
RF3_TRACKED_LIGANDS = []   # was ["ADP", "AMP", "ADN", "GTP"]

# ── Third gating state: holo + CD3 (step6_multistate.py) ─────────────────────
# NOTE ON NAMING: the plan calls this state "holo+CD3e", but 1SY6 and therefore
# inputs/okt3_fv_cd3.cif contain the CD3 gamma-epsilon ectodomain as a single
# chain joined by a 26-residue linker (plan section 1). What is folded here is
# CD3ge, not CD3e alone.
CD3_REF_CIF       = os.path.join(REPO_DIR, "inputs", "okt3_fv_cd3.cif")
CD3_SOURCE_CHAIN  = "C"    # chain holding CD3 in the reference CIF
CD3_CHAIN_ID      = "C"    # chain ID given to CD3 in the RF3 input
CD3_WT_REFERENCE  = True   # also fold wild-type OKT3 Fv + CD3 as a baseline

# TUNE — both are guesses until the reference run reports a wild-type value.
# CD3_IPSAE_MIN is an absolute floor; CD3_RETENTION_MIN is the fraction of the
# wild-type value a design must keep, and is only applied if the reference ran.
CD3_IPSAE_MIN     = 0.30
CD3_RETENTION_MIN = 0.70

RUN_STEP6_MULTISTATE = True   # step5 chains step6 instead of being terminal

# ── Shortlisting thresholds ───────────────────────────────────────────────────
IPTM_PL_CUTOFF = 0.80   # protein-ligand, holo
IPTM_PP_CUTOFF = 0.60   # protein-protein, holo. TUNE

# ── Ligand-dependence filter ──────────────────────────────────────────────────
# Differential, not two independent thresholds. The original pair
# (nol < 0.50, pp >= 0.20) admits designs that assemble BETTER without ligand.
IPTM_PP_DIFF_MIN = 0.30   # require holo - apo >= this
IPTM_PP_NOL_CUTOFF = 0.50  # retained for backwards compatibility with step5
