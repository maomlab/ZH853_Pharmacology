# ZH853–MOR MD SLURM bundle

Self-contained membrane-MD workflow for the MOR–ZH853 (+Gi) complex, to run on an academic
SLURM/GPU cluster. Prepared locally in Phase 2; **cluster-specific values are marked `# TODO(OQ-3)`**
(GPU type, partition, wall-time, account, module/conda paths) and must be filled from the cluster
specs before submission.

## Inputs (copy from the repo)
- `intermediate/02.05.00_oriented/receptorR_oriented.pdb` — rebuilt receptor (Phase 2), superposed
  onto the OPM reference so the membrane normal is z and the OPM midplane is z = 0. `01_build_system.sh`
  stages this automatically; the un-oriented `02.03.00_receptor/receptorR_fixed_heavy.pdb` must **not**
  be used with `--preoriented`.
- `intermediate/02.04.00_ligand/ZH853_prepared.sdf` — protonated ligand (+1), bond orders assigned.
- Protonation/His/D2.50 decisions: `product/02.02.00_protonation_*.md`.

## Systems to build (SPECIFICATION D-10/D-11)
- **System A — active-state complex:** MOR + ZH853 + Gi(αβγ) in POPC:chol (9:1), scFv16 removed.
  For Objective 1–2 dynamics/occupancy.
- **System B — binding/FEP:** MOR + ZH853 in POPC:chol, intracellular half Cα-restrained (or Gα
  α5-helix retained) to hold the active state. For Objective 4 throughput.
- Each built **twice**: D2.50 (Asp116) charged vs protonated (ASH).

> **OPEN — the build is currently apo, and neither System A nor System B.** `01_build_system.sh`
> stages `receptorR_oriented.pdb`, which is the **receptor alone**: no ZH853 and no Gi. `tleap.in`
> does `LIG = loadmol2 ZH853.mol2` but never combines `LIG` into `sys`, so that unit is dropped and
> `saveamberparm` writes a ligand-free system — silently, since loading parameters for an absent
> residue is not an error. `make_tleap.py` now warns when no ligand residue is present in the packed
> box. Fixing this is a build-design decision, not a one-line edit: the ligand (and Gi, for System A)
> must be packed *with* the protein, so `02.05.00` needs to emit a complex whose ligand residue name
> matches `ZH853.mol2`, and PACKMOL-Memgen's `reduce` preprocessing has to be kept from mangling the
> HETATM records. `intermediate/02.05.00_oriented/complex_oriented.pdb` (receptor + ligand, same
> OPM transform) is the natural starting point. Decide this before spending GPU time.

## Conda environments (three, task-specific)
They are split because their `openmm` pins are mutually incompatible (`openmm-plumed` lags the
newest `openmm` that `openmmforcefields` requires):
| Env | Spec | Used by | Notes |
|-----|------|---------|-------|
| `zh853mor-prep` | `environment-prep.yml` | steps 1–2 | CPU; AmberTools/PACKMOL-Memgen |
| `zh853mor-sim` | `environment-cluster.yml` | steps 3–5 + FEP | GPU; openmm + openmmforcefields |
| `zh853mor-plumed` | `environment-plumed.yml` | metadynamics (3.9) | GPU; older openmm + openmm-plumed; optional |

## Workflow (run in order)
| Step | Script | Env | Where |
|------|--------|-----|-------|
| 0 | `00_install.sh` | — | login node (builds the envs) |
| 0.5 | `check_gpu_env.sh` | `zh853mor-sim` | **GPU (pre-flight)** |
| 1 | `ligand_resp/run_resp.sh` | `zh853mor-prep` | CPU |
| 2 | `01_build_system.sh` | `zh853mor-prep` | CPU |
| 2a | `check_placement.py` | `zh853mor-prep` | CPU (called by step 2) |
| 2b | `make_tleap.py` | `zh853mor-prep` | CPU (called by step 2) |

### Building the system (step 2)

`01_build_system.sh` creates a **fresh timestamped `build_*/` directory** and works there. This is
not cosmetic: PACKMOL-Memgen writes to the CWD *and silently reuses anything it finds there* — the
component PDBs (`POPC.pdb`, `CHL1.pdb`, `WAT.pdb`, `Na+.pdb`, `Cl-.pdb`; watch for its "Using
WAT.pdb in the folder" message) and the preprocessed protein files (`receptor_Trim_H.pdb`,
`*.grid.pdb`, `receptorin_EMBED*.pdb`). Once a stale component PDB is present whose atom order no
longer matches the `atoms 1 20` / `atoms 88 131` head/tail constraints in the generated
`packmol.inp`, PACKMOL fails with

```
ERROR: Packmol was unable to put the molecules in the desired regions
       even without considering distance tolerances.
Maximum violation of the restraints:  26.12
```

and **no combination of `--lipids` / `--ratio` / `--dist` / `--preoriented` will fix it**, because
those flags do not affect the cached files. The 26 Å residual is one lipid length and the report says
100 % of that type violates — the signature of contradictory constraints, not of an overpacked box.
Overpacking instead shows up as a GENCAN convergence failure. Set `BUILD_DIR=` to override the
location; the script refuses a non-empty directory.

Two post-build steps then run automatically:

- **`check_placement.py`** — guards the OPM placement from 02.05.00. PACKMOL-Memgen re-centres the
  solute on its own *z* bounding box when it orients the protein itself, and our receptor's bbox
  centre is ~5 Å below the OPM midplane (the intracellular face — H8, ICL3, C-term — protrudes
  further than the extracellular face), so a re-centred build would sit ~5 Å too high in the membrane
  with nothing downstream complaining. Measured 2026-07-26 with `--preoriented` on packmol-memgen
  2025.1.29: translation `(+0.02, −0.18, +0.00) Å`, misregistration **+0.07 Å** — honoured, so this
  is a regression guard rather than a workaround. It measures the receptor's rigid-body shift against
  the lipid phosphate planes, **fails the build** past 1.5 Å, and reports the Trp/Tyr girdle for
  comparison with the OPM ±15.7 Å slab (built: −16.5/+13.7 Å, P–P thickness 38.6 Å, matching the
  37–40 Å predicted in Methods §3.7). Override with `SKIP_PLACEMENT_CHECK=1` only to inspect a
  known-bad build.
- **`make_tleap.py`** — fills the two `@PLACEHOLDERS@` in `tleap.in` that cannot be known until the
  system is packed, writing `tleap_run.in` plus a `bilayer_system_ff.pdb`:
  - the **disulfide**, because `loadpdb` renumbers every residue sequentially from 1 across the whole
    system, so the OPRM1 numbering in `bond sys.142.SG sys.219.SG` no longer exists. C142–C219 is
    found geometrically (SG–SG 2.02 Å) and emitted in tleap's numbering — 74/151 for the current
    69–349 construct, but that shifts with ACE/NME caps, so it is never hardcoded. The two cysteines
    are renamed `CYS`→`CYX` so ff19SB does not build an HG onto a bonded sulfur.
  - the **periodic box**, from PACKMOL's own cell (CRYST1, else the `inside box` bounds in
    `packmol.inp`) — 91.38 × 91.38 × 108.22 Å for the current setup. The previous `setBox sys vdw`
    derived a box from van der Waals extents, which does not reproduce the packing cell.
| 3 | `submit_equilibrate.sbatch` → `02_equilibrate.py` | `zh853mor-sim` | GPU |
| 4 | `submit_production.sbatch` → `03_production.py` | `zh853mor-sim` | GPU |
| 5 | `04_analyze.py` | `zh853mor-sim` | CPU |

## CUDA / modules
OpenMM (conda-forge) bundles its own CUDA runtime, so it needs only the node's NVIDIA **driver** —
usually **no `module load cuda` required**. **cuDNN is never needed** for OpenMM MD (only for ML
potentials). Before the first real run, submit the pre-flight **`sbatch check_gpu_env.sh`** (step 0.5):
it prints the modules/`nvidia-smi`/env, runs `openmm.testInstallation`, and does a real 200-step CUDA
run, ending in a clear **PASS/FAIL** with a fix checklist. If CUDA is missing, load a CUDA module
≤ the `nvidia-smi` "CUDA Version"; if the JIT compiler is not found,
`export OPENMM_CUDA_COMPILER=$(which nvcc)`.

## Force field (SPECIFICATION D-12)
ff19SB (protein) + Lipid21 (membrane) + OPC water + GAFF2/RESP ligand; 0.15 M NaCl. HMR → 4 fs.
(CHARMM36m + CHARMM-GUI is the documented alternative — do not mismatch the water model to the FF.)

## Reproducibility
Pin exact versions in `environment-cluster.yml` once the cluster CUDA/driver stack is known.
Record `openmm.version`, GPU, and CUDA in each run log. Trajectories stay on the cluster; copy back
QC summaries + representative frames to `product/`.
