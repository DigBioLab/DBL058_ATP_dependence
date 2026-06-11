# ATP Binding Analysis

Characterises ATP binding across PDB structures to identify conserved interaction archetypes and conformational diversity. Two notebooks form a pipeline: one builds an ATP conformer library, the other maps protein–ATP contacts across those structures.

---

## Notebooks

### 1. `atp_conformer_analysis.ipynb`
**What it does:** Queries the PDBe API for all PDB structures containing ATP (~3,600 entries), downloads and caches the PDB files, extracts every ATP molecule, and clusters them by structural similarity.

**Key steps:**
- Fetches all ATP-containing PDB IDs from the PDBe REST API
- Extracts and filters ATP conformers to those with exactly 31 heavy atoms (complete structures only)
- Aligns all conformers to an MMFF94-optimised reference via maximum common substructure (MCS)
- Computes an all-vs-all pairwise RMSD matrix (expensive — checkpointed to `output/rmsd_matrix_filtered.npy`)
- Clusters by Butina algorithm (default cutoff 1.5 Å) and embeds in 2D via MDS for visualisation
- Exports up to 10 representative conformers per cluster, plus MMFF94 energy-minimised canonical representatives

**Outputs:** `output/atp_conformers_raw.sdf`, `output/rmsd_matrix_filtered.npy`, `output/cluster_representatives/`, `output/cluster_minimised/`, plots in `pictures/`

**Run this notebook first.** It populates `pdb_cache/` which the binding analysis notebook reads.

---

### 2. `atp_binding_analysis.ipynb`
**What it does:** Runs PLIP (protein–ligand interaction profiler) on the cached PDB structures to detect and classify all ATP–protein contacts, then clusters structures by their interaction fingerprints to identify binding archetypes.

**Key steps:**
- Runs PLIP on each PDB file (parallelised; results cached as XML — re-runs are instant)
- Detects hydrogen bonds, hydrophobic contacts, π-stacking, π-cation interactions, salt bridges, metal coordination, and water bridges
- Encodes each structure as a binary fingerprint of (residue type, interaction type) pairs
- Clusters by Ward linkage on Jaccard distances (default 7 clusters)
- Produces heatmaps, per-cluster profiles, ATP moiety contact maps, and enrichment analysis
- Exports annotated PyMOL sessions for representative structures per cluster

**Outputs:** `output/binding/atp_interactions.csv`, `output/binding/cluster_assignments.csv`, `pymol_sessions/`, plots in `pictures/`

---

## Running on HPC (Recommended)

The heavy computational steps are extracted into scripts and LSF submission scripts so the notebooks only handle fast analysis and visualisation.

### Pipeline overview

```
submit/01_download.sh          →  pdb_cache/*.pdb
                                   output/atp_conformers_raw.sdf
                                   output/pdb_ids.json

submit/02_rmsd_array.sh        →  output/rmsd_bands/band_XXXX.npy  (40 parallel tasks)
submit/02_merge.sh             →  output/rmsd_matrix_filtered.npy

submit/03_plip.sh              →  output/binding/plip_results/*/
                                   output/binding/atp_interactions.csv
```

The notebooks load these checkpoints and run clustering, MDS, and visualisation interactively.

### Step-by-step

```bash
cd /work3/s212001/ATP_dependent

# 1. Download ~3,600 PDB files and extract/filter/align ATP conformers
#    Runtime: ~1–2 h  |  Queue: hpc  |  4 cores, 16 GB
bsub < submit/01_download.sh

# 2a. Compute the N×N RMSD matrix in 40 parallel bands
#     Runtime: ~2–6 h per band  |  Queue: hpc  |  1 core × 40 tasks, 16 GB each
#     Run only after step 1 completes
bsub < submit/02_rmsd_array.sh

# 2b. Merge bands into a single matrix (run after all 40 array tasks finish)
#     Runtime: < 5 min  |  Queue: hpc  |  1 core, 32 GB
bsub < submit/02_merge.sh

# 3. Run PLIP interaction profiling on all structures
#     Runtime: ~2–4 h  |  Queue: hpc  |  16 cores, 32 GB
#     Can run in parallel with steps 2a/2b
bsub < submit/03_plip.sh

# 4. Open notebooks — checkpoints are detected automatically, heavy steps are skipped
jupyter notebook
```

### Monitoring jobs

```bash
bjobs          # list running/pending jobs
bpeek <jobid>  # tail stdout of a running job
bhist -l <jobid>  # job history and exit status
```

---

## Environment Setup

### Requirements
- [Miniforge](https://github.com/conda-forge/miniforge) or Miniconda with `mamba`
- HPC note: the setup below is needed on systems with an older system `libstdc++` (RHEL/CentOS 7–8)

### Install

```bash
# 1. Create the conda environment
mamba env create -f environment.yml

# 2. Activate
mamba activate atp_project_env

# 3. Install plip (pip-only; must be after conda step so openbabel is already present)
pip install plip==3.0.0 --no-build-isolation

# 4. Fix libstdc++ on HPC (system library too old for conda-forge binaries)
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"' \
  > $CONDA_PREFIX/etc/conda/activate.d/ld_library_path.sh

# 5. Reactivate to apply the hook
mamba activate atp_project_env
```

### Verify

```bash
python -c "import openbabel; import plip; import rdkit; import pymol2; print('all OK')"
```

### Why the extra steps?

| Step | Reason |
|------|--------|
| `pip install --no-build-isolation` | plip's build script calls `pip install openbabel` internally; `--no-build-isolation` lets it see conda's already-installed openbabel instead of trying to build it from source |
| `ld_library_path.sh` hook | Conda-forge binaries require `GLIBCXX_3.4.31+`; HPC system `/lib64/libstdc++.so.6` is typically older. The hook makes the conda env's newer `libstdc++` take precedence on activation |

---

## Dependencies

| Package | Source | Purpose |
|---------|--------|---------|
| `rdkit` | conda-forge | Molecule handling, RMSD, clustering, MMFF94 minimisation |
| `openbabel` | conda-forge | Structure format conversion (plip dependency) |
| `pymol-open-source` | conda-forge | PyMOL sessions for binding site visualisation |
| `biopython` | conda-forge | Structure parsing (plip dependency) |
| `plip` | pip | Protein–ligand interaction profiling |
| `numpy`, `scipy`, `pandas` | conda-forge | Numerics, clustering, data handling |
| `matplotlib`, `seaborn` | conda-forge | Plotting |
| `requests` | conda-forge | PDBe API queries |
| `tqdm` | conda-forge | Progress bars |
| `lxml` | conda-forge | PLIP XML output parsing |
| `ipykernel` | conda-forge | Jupyter kernel |
