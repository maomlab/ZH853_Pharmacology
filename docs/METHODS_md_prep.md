# Methods — MD system preparation (Phase 2)

_Living methods section. Pipeline: `src/02.*` and `src/02.10.00_slurm_bundle/`
(`make prep`). Last updated 2026-07-22._

## Overview
Preparation of the cryo-EM MOR–Gi–scFv16–ZH853 model for membrane MD, done locally as far as the
tooling allows (assessment, protonation, receptor rebuild, ligand prep) with the compute-heavy steps
(membrane packing, parameterization, simulation) staged as a self-contained SLURM bundle.

## 1. Assessment (`02.01.00`)
Chain inventory, missing atoms, termini, His sites, disulfides, and the membrane frame.
- MOR chain R spans 69–349 with **no internal gaps** (only sidechain truncations to rebuild).
- **45 incomplete residues** (truncated sidechains, mostly surface Lys/Glu/Arg — expected at 3.5 Å).
- Conserved disulfide **C142–C219** present (enforce explicitly).
- Membrane normal from the TM-bundle principal axis; **3 modeled cholesterols** (84 atoms) span ~28 Å along it
  (bilayer-thickness cue); receptor spans ~69 Å along the normal. Use to pre-orient (PPM/OPM).

## 2. Protonation (`02.02.00`, PROPKA @ pH 7.4)
- **D2.50 = Asp116, pKa 7.61** — essentially at physiological pH ⇒ genuinely ambiguous. **Build
  parallel systems** (charged vs protonated ASH); constant-pH MD is the rigorous fallback (D-11).
- D3.32 (Asp149, pKa 6.02) and E231 (ECL2, pKa 4.93) are **charged** — consistent with their salt
  bridges to ZH853 (Objective 1). D3.49 (Asp166, DRY) charged.
- All four chain-R histidines are **neutral** at pH 7.4 (H299 pKa 3.97, H321 pKa 5.15); assign
  HID/HIE tautomers from the local H-bond network (H299/H321 line the pocket).
- Only one non-standard state overall (D2.50) — the system is otherwise standard.

## 3. Receptor rebuild (`02.03.00`, PDBFixer)
Rebuilt **69 missing sidechain atoms** across the 45 truncated residues → **0 incomplete residues**
remain. Outputs a heavy-atom receptor for the membrane builder and a pH-7.4 protonated copy.
Termini (T69/F349) are **capped ACE/NME** in the tleap step (they are internal fragments of
full-length OPRM1, so neutral caps, not charged termini). Cholesterol is dropped here; the membrane
builder places lipids.

## 4. Ligand preparation (`02.04.00`, RDKit)
Bond orders/aromaticity transferred from the reference SMILES onto the deposited coordinates;
protonated to **net +1** (Tyr1 α-amine — the D149 salt-bridge partner; no free carboxylate, the Glu is
in the lactam). Explicit H added on the 3D pose → SDF/PDB for parameterization. Two FF routes:
GAFF2 + AM1-BCC (quick, for MD) and multi-conformer RESP (rigorous, for FEP; D-3), plus an OpenFF
cross-check.

## 5. System assembly & simulation (SLURM bundle `02.10.00_slurm_bundle/`)
- **Membrane placement (D-14):** production orientation via the **OPM/PPM transfer-energy method**
  (community standard; hydrophobic thickness **~32 Å for MOR**, OPM 4DKL 32.0±1.0 Å; class-A GPCRs
  31–35 Å). The local `02.05.00` cholesterol-centred orientation is a quick first-pass proxy only —
  3 site-specific cholesterols fix the midplane to ~2 Å and their ~28 Å span (≈ POPC hydrocarbon core
  2Dc=28.8 Å, Kučerka 2011) is thin. `03.04.00` validates the placement against the **Trp/Tyr aromatic
  girdle** (~30 Å, agrees with OPM) and experiment; build to the ~31–32 Å OPM slab (P-P ~37–40 Å).
- **Membrane build:** PACKMOL-Memgen, POPC:cholesterol 9:1 (D-4), 0.15 M NaCl, OPC water, ≥15 Å pad,
  extra intracellular Z for the Gi domain; tleap assembly with **ff19SB + Lipid21 + OPC + GAFF2/RESP
  ligand** (D-12).
- **Systems** (D-10): A = MOR + ZH853 + Gi in bilayer (scFv16 removed) for Objective 1–2 dynamics;
  B = MOR + ZH853 (intracellular half restrained / α5 retained) for Objective 4 FEP throughput. Each
  built twice for the D2.50 variants.
- **Equilibration** (`02_equilibrate.py`): CHARMM-GUI-style **6-stage restraint release**
  (backbone 4000→10, sidechain 2000→0, lipid 1000→0 kJ/mol/nm²), NVT→NPT, MonteCarloMembraneBarostat
  (semi-isotropic, γ=0), 2 fs.
- **Production** (`03_production.py`): LangevinMiddle 310 K, MonteCarloMembraneBarostat 1 bar,
  **HMR → 4 fs**, ≥3 replicas × ~500 ns.
- **Per-replica QC** (`04_analyze.py`): backbone RMSD, pre-aligned RMSF, **receptor-aligned ligand
  RMSD**, key-contact occupancy, membrane APL — a smoke test that one run produced something sane.

## 6. Production analysis (`02.11.00_analyze_simulations/`)
Two machines by design (D-18): each replica is reduced **on the cluster** (one SLURM array task;
the panel is ~100 GB of trajectory) to a few hundred kB of observables, and aggregation, tables and
figures run **locally** on those.

- **Numbering (D-19):** tleap renumbers the system from 1, so human OPRM1 numbers are transferred
  positionally from the staged `receptor.pdb` (69–349) and the transfer is refused if the residue
  sequences disagree.
- **Selections:** by mass, not by name or chain — an Amber prmtop has no chains and no elements
  here, and Lipid21 has no POPC residue (PC + PA + OL), so lipids are counted one per phosphorus.
- **Structural QC:** Cα RMSD vs the staged OPM-oriented receptor (whole and TM-only, superposed),
  Cα RMSF, receptor-aligned ligand RMSD against the **deposited pose**, ligand internal RMSD and
  R_g, C142–C219, bilayer thickness (P–P), area per lipid (gross and convex-hull-corrected), OPM
  z-registration, and T/density/volume from the OpenMM state log.
- **Interactions (D-21):** per-residue minimum heavy-atom distance to the ligand every frame, with
  the **same cutoffs as the static fingerprint** (4.5 Å contact, 3.5 Å polar), so MD occupancy and
  the cryo-EM fingerprint are directly comparable; plus bridging-water counts, since the canonical
  H6.52 contact is water-mediated and would otherwise read as absent.
- **State observables:** R3.50–T6.34 and R3.50–Y7.53 Cα rulers (Cα rather than sidechain tips: at
  3.5 Å the rotamers are the least reliable coordinates), and Na⁺ occupancy of the D2.50 site —
  the direct test of the ASP/ASH pair (D-11).
- **Convergence (D-20):** equilibration detected per replica by maximising the effective sample
  size (Chodera 2016) on the TM Cα RMSD, then statistical inefficiency g, n_eff, blocking curves,
  R-hat across replicas and the Hess (2002) cosine content of the leading principal components.
  Means are quoted with correlation-corrected standard errors and replicate spread; occupancies
  with SEM over replicas.

## Open items (cluster-dependent — OQ-3)
GPU type / partition / wall-time / account and QM engine (Psi4 vs Gaussian for RESP) are marked
`TODO(OQ-3)` in the sbatch and RESP scripts; fill from the cluster specs and pin
`environment_zh853mor-sim.yml` before submission.
