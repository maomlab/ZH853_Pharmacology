# Specification of Deliverables

This file records decisions, clarifications, and conventions made as the ZH853–MOR analysis proceeds
(per OBJECTIVES.md). Decisions (D-n) are resolved; Open Questions (OQ-n) await user input.

## Conventions
- **Residue numbering:** human OPRM1 (UniProt **P35372**), which matches the deposited construct
  (chain R). Mouse-numbered comparators (4DKL, 5C1M, 7T2G, 7SBF) are **−2** relative to ours.
- **Ballesteros-Weinstein** generic numbers used alongside construct numbers throughout.
- **Directory/versioning:** raw inputs → `data/`; code → `src/##.##.##_name`; caches →
  `intermediate/`; outputs → `product/…_YYYYMMDD.ext`; plan/decisions → `docs/`.

## Decisions (resolved)
- **D-1 (numbering):** Use human OPRM1 P35372 numbering as the project standard; verified against the
  construct's canonical orthosteric residues. (2026-07-22)
- **D-2 (comparator set):** Primary comparator = **8F7R (endomorphin-1)** as the closest chemical analog;
  plus DAMGO (8EFQ/6DDE), β-endorphin (8F7Q), and small-molecule/biased set (8EF5/6/B/L/O, 7T2G, 5C1M) and
  inactive 4DKL. (2026-07-22)
- **D-3 (ligand FF strategy):** Parameterize ZH853 via the **Amber residue-library route** (QM-RESP capped
  fragments + explicit D-chirality + LEaP ring closure), with an OpenFF/openmmforcefields build as an
  independent cross-check. (2026-07-22)
- **D-4 (membrane baseline):** POPC:cholesterol **9:1** baseline; 7:3 / asymmetric variant as sensitivity
  check. Embed receptor only; Gi in cytoplasmic slab with extra intracellular Z-padding. (2026-07-22)
- **D-5 (free-energy framing):** Treat ABFE/FEP/CTMD results as **rank-ordering** unless convergence is
  explicitly demonstrated, given documented macrocyclic-peptide non-convergence. (2026-07-22)
- **D-6 (ZH853 identity — resolves OQ-1):** ZH853 is the **monomeric macrocycle per the SMILES/structure**
  (Tyr–cyclo[D-Lys–Trp–Phe–Glu]–Gly-NH₂, C₄₂H₅₁N₉O₈, MW 810). The "(…)₂" dimer name in OBJECTIVES is a
  typo. All parameterization/analysis uses this species. (user, 2026-07-22)
- **D-7 (lead priority — resolves OQ-2):** **Objectives 1–2 lead** (comparative interaction analysis +
  mutation panel — highest-confidence, structure-driven). Free energy (Obj 4) proceeds later as
  **rank-ordering only**. (user, 2026-07-22)
- **D-8 (signaling bias — resolves OQ-4):** **Bias analysis is in scope.** Include Gz (9WST) and
  β-arrestin-1 (9WSV) comparators and transducer-state determinant analysis. (user, 2026-07-22)
- **D-9 (interaction-fingerprint method):** Primary static analysis uses a **transparent heavy-atom
  geometric fingerprint** (`zh853mor.interactions`) — appropriate at 3.5 Å where H-bond angle criteria
  are over-precise, and uniform across ligand classes (no per-ligand hydrogens/bond orders needed;
  aromatic rings found by planar-ring geometry). **PLIP + ProLIF cross-validation deferred** to when the
  analysis conda env is built (PLIP's pip build fails without conda OpenBabel). (2026-07-22)
- **D-10 (MD systems):** Build two systems. **A (active-state complex):** MOR + ZH853 + Gi(αβγ) in
  POPC:chol, **scFv16 removed** (crystallization aid, non-physiological) — for Objective 1–2 dynamics/
  occupancy. **B (binding/FEP):** MOR + ZH853 in bilayer with the intracellular half Cα-restrained (or
  Gα α5-helix retained) to hold the active state — for Objective 4 throughput. (2026-07-22)
- **D-11 (D2.50 protonation):** PROPKA gives **Asp116 (D2.50) pKa 7.61** — at physiological pH. Build
  **parallel systems (charged vs protonated ASH)** and compare; constant-pH MD is the rigorous fallback.
  All chain-R His are neutral at pH 7.4 (assign HID/HIE by H-bonding). (2026-07-22)
- **D-12 (force field):** **ff19SB + Lipid21 + OPC water + GAFF2/RESP ligand** (Amber route; matches the
  residue-library ligand plan D-3). CHARMM36m + CHARMM-GUI is the documented alternative. Water model is
  not mismatched to the protein FF. (2026-07-22)

## Open questions (need user input)
- **OQ-3 (compute environment):** SLURM cluster specs (GPU types/count, wall-time limits, queue), and which
  software is preinstalled vs must be built (OpenMM, PLUMED, OpenFE, phenix/MolProbity, Gaussian/Psi4 for RESP)?
  Needed before Phase 2 SLURM bundles; pins versions in `environment-cluster.yml`.
- **OQ-5 (analog set for FEP):** Confirm the analog panel for Objective 4 = ZH850, ZH831, ZH809 (from
  OBJECTIVES), plus any Phase-5 designed analogs; the OBJECTIVES analog names/SMILES have internal Trp/Phe
  inconsistencies (same class of typo as D-6) to reconcile — the **SMILES should be treated as authoritative**.

## Data provenance
- `data/mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb` — 3.5 Å cryo-EM real-space-refined
  model (PHENIX 2.0), MOR–Gi–scFv16–ZH853; provided by the user. Refinement metadata in REMARK 3.
