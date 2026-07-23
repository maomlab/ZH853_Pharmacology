# References — comparator structures & computational methods

Compiled 2026-07-22 from cross-verified literature/PDB research. PDB IDs were confirmed against RCSB
structure pages; entries that could not be verified are flagged.

---

## 1. Comparator MOR (OPRM1) structures

**Adversarial catches to remember:** 6DG7 is *not* MOR (it is 5-HT3A); the DAMGO/Gi companion to 6DDF is
**6DDE**. TRV130 = oliceridine (single structure 8EFB). The stabilizer in MOR–Gi complexes is **scFv16**,
not Nb33/Nb39 (Nb39 appears only in 5C1M, G-protein-free).

### Active-state peptide-agonist complexes (cryo-EM; Gi + scFv16 unless noted)
| PDB | Ligand | Species | Transducer | Res. Å | Citation |
|-----|--------|---------|-----------|--------|----------|
| 6DDE / 6DDF | DAMGO | mouse | Gi + scFv16 | 3.5 | Koehl et al., Nature 558:547 (2018) |
| 8EFQ | DAMGO | human | Gi + scFv16 | 3.30 | Zhuang et al., Cell 185:4361 (2022) |
| **8F7R** | **Endomorphin-1** | human | Gi + scFv16 | 3.28 | Wang et al., Cell 186:413 (2023) |
| 8F7Q | β-endorphin | human | Gi + scFv16 | 3.22 | Wang et al., Cell 186:413 (2023) |
| 8K9L | DAMGO + BMS-986122 (PAM) | human | Gi + scFv16 | 3.05 | Kaneko et al., Nat Commun 15:3544 (2024) |
| 9WST / 9WSW | DAMGO / Endomorphin-1 | human | **Gz** + scFv16 | 2.80 | Zhang et al., Cell Res 35:1021 (2025) |
| 9WSV / 9WSX | DAMGO / Endomorphin-1 | human | **β-arrestin-1** + Fab30 | 2.80 | Zhang et al., Cell Res 35:1021 (2025) |

*8F7R (endomorphin-1) is the chemically closest deposited comparator: ZH853 is an endomorphin-1 analog.*

### Active-state small-molecule agonist complexes
| PDB | Ligand | Species | Transducer | Method Å | Citation |
|-----|--------|---------|-----------|----------|----------|
| 5C1M | BU72 (morphinan) | mouse | Nb39, no G protein | X-ray 2.07 | Huang et al., Nature 524:315 (2015) |
| 8EF5 | Fentanyl | human | Gi + scFv16 | cryo-EM 3.30 | Zhuang et al., Cell 2022 |
| 8EF6 | Morphine | human | Gi + scFv16 | cryo-EM 3.20 | Zhuang et al., Cell 2022 |
| 8EFB | Oliceridine / TRV130 (biased) | human | Gi + scFv16 | cryo-EM 3.20 | Zhuang et al., Cell 2022 |
| 8EFL | SR-17018 (biased) | human | Gi + scFv16 | cryo-EM 3.20 | Zhuang et al., Cell 2022 |
| 8EFO | PZM21 (biased) | human | Gi + scFv16 | cryo-EM 2.80 | Zhuang et al., Cell 2022 |
| 7SBF | PZM21 | mouse | Gi + scFv16 | cryo-EM 2.90 | Wang, Hetzer et al., Angew Chem 61:e202200269 (2022) |
| 7SCG | FH210 | mouse | Gi + scFv16 | cryo-EM 3.00 | Wang, Hetzer et al., 2022 |
| 7T2G | Mitragynine pseudoindoxyl (biased) | mouse | Gi | cryo-EM 2.50 | Qu et al., Nat Chem Biol 19:423 (2023) |
| 7T2H | Lofentanil | mouse | Gi + scFv16 | cryo-EM 3.20 | Qu et al., 2023 |
| 9ODE–9ODI | Lofentanil (activation series) | — | Gi | 2.4–4.3 | Robertson et al., Nature 652:794 (2026) |
| 9PY2/9PY3/9PY4 | Loperamide (GDP intermediates) | — | Gi(GDPβS) | ~3.2 | Khan et al. (Gati lab), Nature 648:755 (2025) |

### Inactive / antagonist (state comparison)
| PDB | Ligand / partner | Species | Method Å | Citation |
|-----|------------------|---------|----------|----------|
| 4DKL | β-funaltrexamine (β-FNA, covalent to K5.39) | mouse (T4L) | X-ray 2.80 | Manglik et al., Nature 485:321 (2012) |
| 8QOT | antagonist nanobody NbE | mouse | cryo-EM 3.20 | Yu et al., Nat Commun 15:8687 (2024) |
| 7UL4 | alvimopan | — | verify | *construct/res. not re-fetched — verify before citing* |
| 9BJK | naloxone + NAM | — | verify | *not re-fetched — verify before citing* |

**Verified ABSENT (no deposited MOR complex, as of 2026-07):** carfentanil (8TFP/9AXN are anti-fentanyl
Fabs, not MOR), sufentanil, remifentanil, ohmefentanyl, hydromorphone, oxymorphone, endomorphin-2, DALDA/
[Dmt¹]DALDA, biphalin, and any cyclic/stapled-peptide MOR agonist. **ZH853 itself is not deposited.**

## 2. Orthosteric pocket residues (BW → human OPRM1 P35372; verified matches our construct)
D149 (3.32, salt-bridge anchor) · Y150 (3.33) · M153 (3.36) · K235 (5.39, β-FNA covalent site) ·
W295 (6.48, toggle switch) · I298 (6.51) · **H299 (6.52, water-mediated phenol H-bond)** · V302 (6.55) ·
W320 (7.35, μ/δ selectivity) · I324 (7.39) · G327 (7.42) · Y328 (7.43, 3–7 lock) ·
Q126 (2.60) · N129 (2.63). Mouse numbers (4DKL/5C1M/7T2G/7SBF) are **−2** relative to these.

---

## 3–8. Computational methods (2024–2026 best practice)

Full detail in the plan; condensed pointers with maturity flags:

- **Membrane build (mature):** CHARMM-GUI Membrane Builder → OpenMM inputs, or PACKMOL-Memgen; pre-orient
  PPM/OPM; POPC:chol 9:1 baseline (7:3 sensitivity); Gi in cytoplasmic slab with extra Z-padding. OpenMM
  `Modeller.addMembrane()` cannot do cholesterol. Refs: JCTC 10.1021/acs.jctc.2c01246; JCIM 10.1021/acs.jcim.9b00269.
- **Force fields (mature):** ff19SB+OPC+Lipid21 **or** CHARMM36m+C36+mTIP3P — never mismatch water to FF.
  Refs: 10.1021/acs.jctc.9b00591; 10.1002/pro.4413; 10.1021/acs.jctc.1c01217.
- **Macrocyclic-peptide ligand FF (moderate/hard):** Amber residue-library route — capped-dipeptide
  fragments, multi-conformer RESP HF/6-31G(d), ff19SB backbone + GAFF2 side chains, explicit D-chirality,
  LEaP ring closure. openmmforcefields `SystemGenerator` (≥ v0.16.0 for multi-residue ligands) to combine
  in OpenMM. OpenFF-Sage is a small-molecule (not peptide) FF. cis/trans + ring conformers poorly handled →
  enhanced sampling + NMR validation. Refs: PMC12713365; 10.1021/acs.jcim.4c01120; openmmforcefields GitHub.
- **Equilibration (mature):** CHARMM-GUI 6-step restraint release; LangevinMiddle 310 K +
  MonteCarloMembraneBarostat (semi-isotropic, surface tension 0); HMR 4 fs. Refs: choderalab b2ar_membrane
  tutorial; 10.1021/acs.jctc.9b00160 (HMR).
- **ABFE (immature/high-risk for macrocycles):** OpenFE v1.7 `AbsoluteBindingProtocol`, BAT2; Boresch
  restraints ill-defined for floppy peptide; free-ligand non-convergence. Mitigate: ATM, Lambda-ABF-OPES,
  restrained ABFE; alchemlyb/MBAR. Refs: openfree.energy v1.7; 10.1021/acs.jctc.4c00205; Chem Sci 2020
  macrocycle FEP convergence.
- **Relative FEP (emerging/fragile for peptides):** OpenFE (hybrid-topology/SepTop), PMX, Perses; Kartograf
  mapping (LOMAP can't break rings — matters for Glu↔Asp edits). Refs: OpenFE tutorials; 10.1021/acs.jctc.3c01206.
- **CTMD / metadynamics (precedented, peptide-hard):** PLUMED + openmm-plumed; funnel metadynamics for
  qualitative ΔG/pose; c(t) time-independent estimator (10.1021/jp504920s); CTMD-scoring is single-paper
  (PMC12889651). Watch PLUMED 1-indexed vs OpenMM 0-indexed atoms.
- **Interaction analysis (mature):** ProLIF (RDKit+MDAnalysis) for trajectory fingerprints/occupancy
  (`fp.to_dataframe().mean()`); getcontacts (GPCR standard); PLIP for static snapshots. Refs: ProLIF &
  getcontacts GitHub; PLIP PMC4489249.
- **QC (mature; convergence criteria contested):** Cα RMSD plateau, RMSF (pre-align!), receptor-aligned
  ligand RMSD, contact occupancy, membrane APL (POPC≈0.68 nm²)/thickness/S_CD; multiple replicas + block
  averaging; RMSD-plateau alone is insufficient. Tools: MDAnalysis, LiPyphilic.
