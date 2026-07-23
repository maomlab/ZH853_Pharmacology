# ZH853 · MOR Docking & Structure Analysis

Computational structural-biology study of how the cyclic peptide agonist **ZH853** binds and modulates
the **μ-opioid receptor (MOR / OPRM1)**, starting from a 3.5 Å cryo-EM structure of the
**MOR–ZH853–Gi–scFv16** complex.

## Goals
1. Characterize the **unique structural interactions** of ZH853 vs other MOR agonist complexes.
2. Identify **mutations** that selectively abrogate ZH853 but spare other full agonists.
3. Propose **drug-likeness improvements** (peptide PK/PD, GLP-1-style strategies).
4. Predict relative affinities of ZH853 and analogs via **ABFE / relative FEP / metadynamics**.

## Where things live
| Path | Contents |
|------|----------|
| `data/` | Raw inputs — the cryo-EM PDB, comparator structures, ligand SMILES |
| `src/` | Workflow + code, numbered `src/##.##.##_name` (sequential + hierarchical) |
| `intermediate/` | Cached/temporary results (git-ignored) |
| `product/` | Outputs — figures, tables, reports (`_YYYYMMDD.ext`), manuscript |
| `docs/` | [`PLAN.md`](docs/PLAN.md) (roadmap), [`references.md`](docs/references.md) (comparators + methods) |
| `OBJECTIVES.md` | Original project brief |
| `SPECIFICATION.md` | Decisions log + open questions |

## Status
Phases 0–1 and static interaction analysis (Phase 3) complete — see [`docs/PLAN.md`](docs/PLAN.md)
and [`docs/RESULTS_interactions.md`](docs/RESULTS_interactions.md). Key findings: ZH853 keeps the
conserved D149 (D3.32) anchor but uniquely engages **E231 (ECL2, salt bridge)** and **H321 (7.36)** —
yielding a ranked mutation panel led by **E231Q/E231A** (Objective 2). Reproduce with `make analysis`.
Analog design (Objective 3) also complete — see [`docs/RESULTS_analog_design.md`](docs/RESULTS_analog_design.md):
all analogs are beyond-Ro5; two design series proposed (N-methylation for permeability, lipidation for
half-life). Open items: cluster specs (OQ-3) for MD/prep; MD occupancy validation; Phase 6 free energy.

## Environment
Analysis (local, macOS): Python 3.10 + numpy, biotite, rdkit, MDAnalysis (present); add ProLIF, PLIP,
pdbfixer, mdtraj. Simulation (SLURM cluster): OpenMM, OpenFF, openmmforcefields, PLUMED, OpenFE/BAT2.
Env specs and a `Makefile` (`make lint typecheck test figures`) land in Phase 0.
