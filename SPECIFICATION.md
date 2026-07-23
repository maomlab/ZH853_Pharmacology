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

## Open questions (need user input)
- **OQ-1 (ZH853 identity — load-bearing):** OBJECTIVES names ZH853 as "(Tyr-[D-Lys-Phe-Phe-Asp]₂-NH₂)₂"
  (a dimer, Phe-Phe), but the accompanying SMILES and the modeled ligand `L01` (59 heavy atoms, C₄₂H₅₁N₉O₈)
  are a **monomeric macrocycle containing Trp** (Tyr–cyclo[D-Lys–Trp–Phe–Glu]–Gly-NH₂). The SMILES matches
  the structure; the name does not. **Please confirm the intended chemical identity and stereochemistry.**
- **OQ-2 (prioritization & depth):** Which of Objectives 1–4 should lead, and how deep on the high-risk
  free-energy work (rank-ordering only vs full ABFE campaign)?
- **OQ-3 (compute environment):** SLURM cluster specs (GPU types/count, wall-time limits, queue), and which
  software is preinstalled vs must be built (OpenMM, PLUMED, OpenFE, phenix/MolProbity, Gaussian/Psi4 for RESP)?
- **OQ-4 (signaling bias scope):** Is ZH853's reported signaling bias / reduced-tolerance profile in scope
  (would pull in Gz/β-arrestin comparators 9WST/9WSV and bias-focused analysis)?
- **OQ-5 (analog set for FEP):** Confirm the analog panel for Objective 4 = ZH850, ZH831, ZH809 (from
  OBJECTIVES), plus any Phase-5 designed analogs; note the OBJECTIVES analog names/SMILES also have internal
  Trp/Phe inconsistencies to reconcile alongside OQ-1.

## Data provenance
- `data/mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb` — 3.5 Å cryo-EM real-space-refined
  model (PHENIX 2.0), MOR–Gi–scFv16–ZH853; provided by the user. Refinement metadata in REMARK 3.
