#!/bin/bash
# Enzyme design pipeline orchestrator: RFD3 → fill_fixed_residues → MPNN → RF3 → RF3-nol → analysis
#
# Usage:
#   bash /path/to/ATP_dependent/run_pipeline.sh
#
# What this script does:
#   1. Runs step0_setup.py locally to generate RFD3 configs and submit script
#   2. Submits the RFD3 GPU job
#   3. Submits step2 as a CPU job that waits for RFD3
#      — step2 generates MPNN cmds, calls scripts/fill_fixed_residues.py to populate
#        fixed residues from motif atoms, then submits MPNN and chains step3 → step4 → step5
#      — step3 generates RF3 cmds, submits RF3, then submits step4
#      — step4 shortlists, generates RF3-nol cmds, submits RF3-nol, then submits step5
#      — step5 runs the final ligand-dependence analysis (terminal step)
#
# Edit pipeline/config.py to change campaign parameters (user, PYTHON_EXEC, thresholds).

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_DIR="${SCRIPT_DIR}/pipeline"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

# Read campaign-specific paths from config.py — edit pipeline/config.py to change these.
_cfg() { python3 -c "import sys; sys.path.insert(0, '${PIPELINE_DIR}'); import config; print($1)"; }
PYTHON_EXEC="$(_cfg 'config.PYTHON_EXEC')"
WORKING_DIR="$(_cfg 'config.WORKING_DIR')"
CMDS_DIR="$(_cfg 'config.CMDS_DIR')"
LOGS_DIR="$(_cfg 'config.LOGS_DIR')"

mkdir -p "${LOGS_DIR}"

get_jid() {
    grep -oP '(?<=Job <)\d+(?=>)'
}

# ── Step 0: Generate RFD3 config and submit script ────────────────────────────
echo "Running setup (generating RFD3 config and submit script)..."
"${PYTHON_EXEC}" "${PIPELINE_DIR}/step0_setup.py"

RFD3_SH="${CMDS_DIR}/rfd3.sh"
if [[ ! -f "${RFD3_SH}" ]]; then
    echo "ERROR: ${RFD3_SH} was not created by step0_setup.py" >&2
    exit 1
fi

# ── Step 1: Submit RFD3 (GPU) ─────────────────────────────────────────────────
RFD3_JID=$(bsub < "${RFD3_SH}" | get_jid)
echo "Submitted RFD3 job:    ${RFD3_JID}"

# ── Step 2: Submit post-RFD3 processing (CPU, waits for RFD3) ────────────────
# step2 will internally submit MPNN and then chain step3 → step4 → step5
STEP2_JID=$(bsub \
    -w "done(${RFD3_JID})" \
    -J step2_process_rfd3 \
    -q hpc \
    -n 4 \
    -W 02:00 \
    -R "rusage[mem=10GB]" \
    -R "span[hosts=1]" \
    -o "${LOGS_DIR}/step2_%J.out" \
    -e "${LOGS_DIR}/step2_%J.err" \
    "${PYTHON_EXEC} ${PIPELINE_DIR}/step2_process_rfd3.py" \
    | get_jid)
echo "Submitted step2 job:   ${STEP2_JID}  (waits for RFD3 ${RFD3_JID})"

echo ""
echo "Pipeline queued. The full chain will run automatically:"
echo "  RFD3 (${RFD3_JID}) → step2 (${STEP2_JID}) → fill_fixed_residues → MPNN → step3 → RF3 → step4 → RF3-nol → step5"
echo ""
echo "Monitor with:"
echo "  bjobs ${RFD3_JID} ${STEP2_JID}"
echo "  bqueues"
echo ""
echo "Logs:"
echo "  ${LOGS_DIR}/"
