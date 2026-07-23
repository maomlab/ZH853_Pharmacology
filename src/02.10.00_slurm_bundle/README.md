# ZH853–MOR MD SLURM bundle

Self-contained membrane-MD workflow for the MOR–ZH853 (+Gi) complex, to run on an academic
SLURM/GPU cluster. Prepared locally in Phase 2; **cluster-specific values are marked `# TODO(OQ-3)`**
(GPU type, partition, wall-time, account, module/conda paths) and must be filled from the cluster
specs before submission.

## Inputs (copy from the repo)
- `intermediate/02.03.00_receptor/receptorR_fixed_heavy.pdb` — rebuilt receptor (Phase 2).
- `intermediate/02.04.00_ligand/ZH853_prepared.sdf` — protonated ligand (+1), bond orders assigned.
- Protonation/His/D2.50 decisions: `product/02.02.00_protonation_*.md`.

## Systems to build (SPECIFICATION D-10/D-11)
- **System A — active-state complex:** MOR + ZH853 + Gi(αβγ) in POPC:chol (9:1), scFv16 removed.
  For Objective 1–2 dynamics/occupancy.
- **System B — binding/FEP:** MOR + ZH853 in POPC:chol, intracellular half Cα-restrained (or Gα
  α5-helix retained) to hold the active state. For Objective 4 throughput.
- Each built **twice**: D2.50 (Asp116) charged vs protonated (ASH).

## Workflow (run in order)
| Step | Script | What | Where |
|------|--------|------|-------|
| 0 | `00_install.sh` | build the `zh853mor-sim` conda env | login node |
| 1 | `ligand_resp/run_resp.sh` | GAFF2+AM1-BCC (quick) or RESP (rigorous) ligand params | GPU/CPU |
| 2 | `01_build_system.sh` | PACKMOL-Memgen bilayer + tleap assembly → prmtop/rst7 | CPU |
| 3 | `submit_equilibrate.sbatch` → `02_equilibrate.py` | 6-stage restrained equilibration | GPU |
| 4 | `submit_production.sbatch` → `03_production.py` | production (≥3 replicas) | GPU |
| 5 | `04_analyze.py` | QC: RMSD/RMSF, ligand RMSD, contact occupancy, membrane | CPU |

## Force field (SPECIFICATION D-12)
ff19SB (protein) + Lipid21 (membrane) + OPC water + GAFF2/RESP ligand; 0.15 M NaCl. HMR → 4 fs.
(CHARMM36m + CHARMM-GUI is the documented alternative — do not mismatch the water model to the FF.)

## Reproducibility
Pin exact versions in `environment-cluster.yml` once the cluster CUDA/driver stack is known.
Record `openmm.version`, GPU, and CUDA in each run log. Trajectories stay on the cluster; copy back
QC summaries + representative frames to `product/`.
