"""
Enzyme design pipeline configuration — edit this file to match your campaign.
All other pipeline scripts import from here; you should not need to change them.
"""

# ── Infrastructure ────────────────────────────────────────────────────────────
PIPELINE_DIR   = "/work3/s212001/ATP_dependent/pipeline"
PYTHON_EXEC    = "/work3/s212001/miniforge3/envs/atp_project_env/bin/python3"
ENZDES_SCRIPTS = "/work3/s212001/DBL_secrets/enzdes_pipeline-main"

# ── Campaign ──────────────────────────────────────────────────────────────────
user     = "s212001"
family   = "atp_binding"
target   = "conditioned_assembly"
campaign = "8btg_D501"

WORKING_DIR = f"/work3/{user}/projects/enzyme_design/{family}/{target}/{campaign}"

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

# ── RFD3 parameters ───────────────────────────────────────────────────────────
INPUT_PDB = f"{WORKING_DIR}/inputs/ATP_binding_8btg.cif"
DIALECT   = 2

# Scaffold layout: free ranges (e.g. 30-50) interspersed with fixed motif residues (e.g. A58-58)
CONTIG = "20-30,D151-D159,3-8,D120-123,20-30,/0,20-30,D287-D289,20-30,D312-D314,10-20,C259-C265,20-30"
LENGTH = "136-203"

LIGAND_NAME               = "ATP"
REDESIGN_MOTIF_SIDECHAINS = False
IS_NON_LOOPY              = True
INFER_ORI_STRATEGY        = None
HOTSPOTS                  = None

# Sidechain atoms to pin during diffusion (preserves catalytic geometry)
FIXED_ATOMS = {
    "D121": "N,O",
    "D152": "O",
    "D154": "N",
    "D155": "N",
    "D156": "CA,C,N",
    "D157": "N,NZ",
    "D158": "N,OG1",
    "D159": "N,ND1",
    "D288": "NZ",
    "D313": "NH1,NH2",
    "C260": "NH1,NH2",
    "C264": "NH1,NH2",
}

RFD3_CKPT             = "/dtu/projects/dbl/foundry/ckpt/rfd3_latest.ckpt"
RFD3_DIFFUSION_BATCH  = 10
RFD3_N_BATCHES        = 10

# ── MPNN parameters ───────────────────────────────────────────────────────────
MPNN_SEED             = 42
MPNN_BATCH_SIZE       = 10   # sequences per structure
MPNN_N_BATCHES        = 1    # total = MPNN_BATCH_SIZE × MPNN_N_BATCHES
MPNN_CKPT             = "/dtu/projects/dbl/LigandMPNN/25.3.1-0//model_params/ligandmpnn_v_32_010_25.pt"
MPNN_CHAINS_TO_DESIGN = "A,B"
MPNN_REMOVE_BAD       = True   # discard structures with chain breaks or clashes before design

# ── RF3 parameters ────────────────────────────────────────────────────────────
RF3_CKPT      = "/dtu/projects/dbl/rf3/release/ckpt/rf3_latest.pt"
RF3_NUM_STEPS = 50   # 50 for speed; 200 for production
RF3_N_MODELS  = 5

# Ligand spec: set RF3_CCD_LIGAND=True if ligand has a CCD code (e.g. "ATP"),
# else provide RF3_SDF_PATH (CIF/SDF file) or RF3_SMILES as fallback.
RF3_CCD_LIGAND = False
RF3_SDF_PATH   = f"{WORKING_DIR}/inputs/cluster_01_top1_8BTG_D_501.cif"
RF3_SMILES     = None

# ── Shortlisting thresholds ───────────────────────────────────────────────────
IPTM_PL_CUTOFF = 0.80   # minimum ipTM protein–ligand (best model)
IPTM_PP_CUTOFF = 0.20   # minimum ipTM protein–protein (best model)

# ── No-ligand filter threshold ────────────────────────────────────────────────
# Designs scoring ABOVE this without ligand are rejected (they fold fine without it)
IPTM_PP_NOL_CUTOFF = 0.50
