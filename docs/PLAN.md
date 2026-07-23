# ZH853–MOR Structure Analysis — Project Plan

**Status:** living document. Last updated 2026-07-22.
Tracks what is planned, in progress, done, and considered-but-deferred.
Decisions and clarifications are logged in [`../SPECIFICATION.md`](../SPECIFICATION.md).

> **User-set direction (2026-07-22):** ZH853 = monomeric macrocycle (D-6). **Objectives 1–2 lead**
> (interactions + mutations); free energy is later rank-ordering only (D-7). **Signaling bias is in
> scope** — Gz/β-arrestin comparators included (D-8).
>
> **Progress:** Phase 0 **done** (scaffold, envs, Makefile, CI gate green). Phase 1 **done**
> (13 comparator PDBs fetched; QC report; BW map). Phase 3 static analysis **done** — see
> [RESULTS_interactions.md](RESULTS_interactions.md): **Objective 1** identified E231 (ECL2, unique
> ionic) + H321 (7.36) as ZH853-distinctive vs the conserved D149 anchor; **Objective 2** yields a
> ranked mutation panel led by **E231Q/E231A**. Phase 5 analog design **done** — see
> [RESULTS_analog_design.md](RESULTS_analog_design.md): all analogs beyond-Ro5 (TPSA 235–280, HBD 8–10);
> structure-based map (46 buried / 13 exposed atoms) yields two series — N-methylation (permeability) from
> ZH831, C-terminal lipidation (half-life) from ZH853. **Phase 2 local prep done** — see
> [METHODS_md_prep.md](METHODS_md_prep.md): receptor rebuilt (69 atoms, 0 incomplete), protonation
> resolved (D2.50/Asp116 pKa 7.61 → parallel systems), ligand prepped (+1), and a SLURM bundle
> (`src/02.10.00_slurm_bundle/`) with OpenMM equilibration/production/QC scripts. **Phase 7 manuscript
> draft started** — `product/manuscript/manuscript_20260722.md` synthesizes Objectives 1–3 + methods
> ([prospective] items flagged pending MD/FEP). Remaining: cluster submission (OQ-3), MD occupancy
> validation, PLIP/ProLIF cross-check, Phase 6 FEP, then fold results into the manuscript.

---

## 1. Scientific framing

We have a **3.5 Å cryo-EM structure of the μ-opioid receptor (MOR / OPRM1) bound to the
cyclic peptide agonist ZH853, in complex with the Gi heterotrimer and the scFv16 stabilizing
fragment** (`data/mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb`).

ZH853 is a **macrocyclic, D-amino-acid endomorphin-1 analog** from the Zadina lab
(Zadina et al., 2016, DOI:10.1016/j.neuropharm.2015.12.024). Two facts shape the entire project:

1. **This is (to our knowledge) the first cyclic-peptide MOR complex.** No macrocyclic-peptide
   MOR structure is deposited in the PDB; the nearest cyclic-peptide opioid structure (8FEG) is
   κ-opioid receptor. There is therefore **no direct structural precedent for the ligand class** —
   the comparators are linear peptides (DAMGO, endomorphin-1, β-endorphin) and small molecules.
2. **The free-energy component carries most of the technical risk.** ZH853 is a large, flexible
   macrocyclic peptide; the shared failure mode across ABFE, relative FEP, and metadynamics is that
   the ligand's conformational ensemble (bound *and especially free*) does not converge on routine MD
   timescales. Free-energy numbers should be treated as **rank-ordering** unless convergence is
   explicitly demonstrated (see [references](references.md) §4–6).

### 1.1 What the structure actually contains (verified locally)

| Chain | Contents | Modeled span |
|-------|----------|--------------|
| A | Gαi | 3–354 |
| B | Gβ1 | 5–340 |
| C | Gγ2 | 10–63 |
| D | scFv16 | 1–248 |
| R | **MOR (OPRM1)** + 84 cholesterol (CLR) atoms | **69–349, no internal gaps** |
| E | **ZH853** (ligand `L01`) | 1 macrocyclic residue, 59 heavy atoms |

- **Numbering is human OPRM1 (UniProt P35372).** Verified: every canonical orthosteric position
  matches (D149/D3.32, Y150/3.33, M153/3.36, K235/5.39, W295/6.48, I298/6.51, **H299/6.52**,
  V302/6.55, W320/7.35, I324/7.39, G327/7.42, Y328/7.43). **No +2 mouse offset applies here** — but
  the offset *does* apply when importing residue numbers from mouse structures (4DKL, 5C1M, 7T2G, 7SBF).
- **Model quality (from REMARK 3):** CC_mask 0.80, CC_volume 0.76, clashscore 12.3, rotamer outliers
  2.8 %, Ramachandran outliers 0.0 %, bond RMSD 0.022 Å. Adequate but 3.5 Å — sidechain rotamers and
  the ligand pose have real uncertainty; MD is partly a refinement/validation exercise.
- **Conserved ECL2–TM3 disulfide C142–C219 present** (2.03 Å).
- **ZH853 = C₄₂H₅₁N₉O₈, MW 810** (RDKit, from the OBJECTIVES SMILES) — heavy-atom count matches the
  modeled ligand exactly (59), confirming the deposited species is the **monomeric macrocycle**
  Tyr–cyclo[D-Lys–Trp–Phe–Glu]–Gly-NH₂, **not** the dimer implied by the "(…)₂" name in OBJECTIVES
  (see SPECIFICATION open question OQ-1).

### 1.2 ZH853 binding-pocket interaction map (verified locally, 4.5 Å shell)

Direct polar anchors (heavy–heavy distances):

| Receptor residue | BW | Ligand atom | Distance | Interpretation |
|---|---|---|---|---|
| **Asp149** | 3.32 | N09 | 3.16 Å | Salt bridge to ligand protonated amine — **universal opioid anchor** |
| **Tyr328** | 7.43 | N09 | 3.17 Å | 3–7 "lock"; H-bond to same amine |
| **Glu231** | ECL2 | N51 | 3.12 Å | Charged ECL2 contact — candidate ZH853-distinctive |
| **Asn129** | 2.63 | N31 | 3.02 Å | Upper-pocket polar network |
| His299 | 6.52 | (4.89 Å) | — | Canonical water-mediated phenol contact; just outside 4.5 Å shell |
| **His321** | ~7.36/ECL3 | contact | 3.04 Å | **Second His, direct contact — candidate ZH853-distinctive** |

Full hydrophobic/aromatic shell (4.5 Å): Y77, Q126, N129, Y130, W135, V145, I146, D149, Y150, M153,
C219, T220, L221, F223, E231, L234, K235, V238, W295, I298, V302, W320, H321, I324, Y328.

This map is the seed for Objective 2 (mutation design): **shared anchors** (D149, H299, Y328) likely
cannot discriminate ZH853 from other agonists, whereas **ZH853-enriched contacts** (E231, H321, N129,
the macrocycle-extended contacts around ECL2/TM7) are the candidate selectivity determinants.

### 1.3 Objectives → deliverables mapping

| # | Objective (from OBJECTIVES.md) | Primary deliverable |
|---|---|---|
| 1 | Unique structural interactions vs other MOR complexes | Comparative interaction-fingerprint analysis + figures (Phase 3) |
| 2 | Mutations that abrogate ZH853 but spare other full agonists | Ranked mutation panel + MD/FEP validation (Phases 3, 6) |
| 3 | Drug-likeness improvements (peptide PK/PD, GLP-1-style) | Analog design report + property/permeability predictions (Phase 5) |
| 4 | ABFE / FEP / CTMD affinity predictions for ZH853 + analogs | Free-energy task bundles + rank-ordering (Phase 6) |
| — | MD of receptor+ligand complex; QC of simulation quality | Equilibrated system + production + QC dashboard (Phases 2, 4) |
| — | Manuscript | `product/manuscript/` (Phase 7) |

---

## 2. Comparator structures (verified — full table in [references.md](references.md))

**Primary comparators (human, active, Gi + scFv16):**
- **8F7R — endomorphin-1** (3.28 Å): chemically the closest deposited analog (ZH853 is an
  endomorphin-1 analog). **Top-priority comparator.**
- **8EFQ — DAMGO** (human, 3.30 Å) and **6DDE/6DDF — DAMGO** (mouse, 3.5 Å, same resolution as ours).
- **8F7Q — β-endorphin** (3.22 Å).

**Small-molecule / biased-agonist comparators:** 8EF5 (fentanyl), 8EF6 (morphine), 8EFB (oliceridine/
TRV130), 8EFL (SR-17018), 8EFO/7SBF (PZM21), 7T2G (mitragynine pseudoindoxyl, 2.5 Å high-res),
7T2H/9ODE-series (lofentanil). **5C1M (BU72, 2.07 Å X-ray)** and **7T2G** are the high-resolution pocket
references.

**State / transducer comparators:** 4DKL (β-FNA, inactive), 8QOT (antagonist Nb); 9WST/9WSW (Gz),
9WSV/9WSX (β-arrestin-1) — relevant if ZH853 signaling bias is in scope.

---

## 3. Project structure & conventions

Per OBJECTIVES.md organizing principles:

```
data/            # raw/external inputs (the cryoEM PDB; downloaded comparator PDBs; SMILES)
src/             # workflow + supporting code, numbered src/##.##.##_name (sequential + hierarchical)
intermediate/    # cached/temporary results, named after the generating script (git-ignored)
product/         # outputs: figures/tables/reports, <name>_YYYYMMDD.<ext>; manuscript/
docs/            # this plan, references, decisions
```

**Planned `src/` layout (numbering is provisional; fill in as scripts land):**

| Index | Stage |
|-------|-------|
| `src/00.00.00_env/` | Conda/mamba env specs (analysis env; cluster env); Makefile targets |
| `src/01.*` | **Data acquisition & QC** — fetch comparator PDBs; parse/validate cryoEM model; MolProbity-style checks |
| `src/02.*` | **Structure preparation** — protonation, caps, disulfides, ligand parameterization, membrane build |
| `src/03.*` | **Static interaction analysis** — PLIP/ProLIF fingerprints; comparative pocket analysis (Obj 1, 2) |
| `src/04.*` | **MD equilibration & production** — SLURM task bundles; QC (Obj: MD + QC) |
| `src/05.*` | **Analog cheminformatics & design** — RDKit properties, permeability, GLP-1-style modifications (Obj 3) |
| `src/06.*` | **Free-energy calculations** — ABFE / relative FEP / metadynamics task bundles (Obj 4) |
| `src/07.*` | **Manuscript assembly** — figures, tables, text |

**Engineering standards:** ruff (lint) + mypy (type-check) + pytest, driven from a `Makefile`
(`make lint typecheck test figures`). Python 3.10 analysis env (local: numpy, biotite, rdkit,
MDAnalysis already present; add ProLIF, PLIP, pdbfixer, mdtraj). OpenMM/OpenFF/PLUMED live in the
**cluster** env, not required locally.

---

## 4. Phased plan

Each phase lists concrete tasks and an explicit **decision gate**. Cluster-bound work is packaged as
self-contained SLURM bundles (data + `submit.sbatch` + install/prepare/gather/post-process scripts) per
OBJECTIVES; the user runs these manually.

### Phase 0 — Repo scaffold & environment  *(local, ~0.5 day)*
- Initialize `src/`, `intermediate/`, `product/`, `docs/` with the numbering scheme; `Makefile`;
  `environment.yml` (analysis) and `environment-cluster.yml` (OpenMM/OpenFF/PLUMED); pre-commit
  (ruff+mypy); `README.md`.
- **Gate:** `make lint typecheck` green on a trivial module.

### Phase 1 — Data acquisition & structure QC  *(local, 1–2 days)*
- Fetch primary comparator PDBs (8F7R, 8EFQ, 6DDE, 8F7Q, 5C1M, 7T2G, 4DKL, 8EFB, 8EFL, 8EFO) into `data/`.
- Parse the cryoEM model: per-chain inventory, B-factor/occupancy profile, ligand sanity, cholesterol
  positions; reproduce the numbering-verification and pocket-contact analysis as versioned scripts.
- Independent geometry QC (MolProbity via `phenix`/`molprobity` if available, else Biotite-based clash/
  rotamer summary) to characterize where the 3.5 Å model is least reliable (esp. ligand rotamers, ECL2/3).
- Build a **Ballesteros-Weinstein ↔ construct-residue map** by aligning chain R to OPRM1 (P35372) so all
  downstream analyses can speak BW numbers.
- **Gate:** confirmed residue-numbering map; documented list of low-confidence regions feeding prep decisions.

### Phase 2 — Structure preparation  *(local build → cluster equilibration, 3–5 days)*
System prep is the mature, low-risk backbone of the project. Key decisions logged to SPECIFICATION.
- **Receptor:** cap truncated termini (R69 / R349); model any missing sidechain atoms; decide ICL3 /
  fusion handling; assign protonation with PDB2PQR/PROPKA at pH 7.4 **plus explicit parallel D2.50
  (Asp) protonated/deprotonated systems** (Na⁺-pocket ambiguity, see references §1); enforce C142–C219
  and other disulfides.
- **ZH853 ligand parameterization (the crux):** use the **Amber residue-library route** (not
  whole-macrocycle small-molecule typing): capped-dipeptide fragments for each non-canonical/D residue,
  multi-conformer **RESP at HF/6-31G(d)**, ff19SB backbone + GAFF2 novel side-chain terms, explicit
  D-chirality, head-to-tail ring closure in LEaP. Cross-check against an OpenFF-Sage/openmmforcefields
  build (≥ v0.16.0 for multi-residue ligands) as an independent sanity check.
- **Membrane:** build MOR in a **POPC:cholesterol bilayer** (9:1 control; a 7:3 / asymmetric
  plasma-membrane-mimetic variant as a sensitivity check) via **CHARMM-GUI Membrane Builder** (emits
  OpenMM inputs) or PACKMOL-Memgen; pre-orient with PPM/OPM; Gi heterotrimer in the cytoplasmic water
  slab with **extra intracellular Z-padding** to avoid PBC self-image.
- **FF / solvent:** ff19SB + OPC + Lipid21 **or** CHARMM36m + C36 lipids + mTIP3P (do not mismatch water
  to FF); 0.15 M NaCl; neutralize.
- **Deliverable:** `src/02.*` build scripts + a **prep SLURM bundle**; equilibrate with the CHARMM-GUI
  6-step restraint schedule; production settings LangevinMiddle 310 K + MonteCarloMembraneBarostat, HMR
  4 fs.
- **Gate:** equilibrated system passing basic QC (stable box/density, APL ≈ 0.65–0.68 nm² for POPC,
  intact pocket); the ligand pose stable under light restraints before free production.

### Phase 3 — Static & short-MD interaction analysis  *(local + short cluster runs, 3–5 days)* — **Objectives 1 & 2**
- **Static fingerprints** with **PLIP** (representative cryoEM pose) and cross-check with **ProLIF**.
- **Trajectory occupancy** with **ProLIF** (or getcontacts): interaction persistence over equilibration/
  short production, reported as fraction-of-frames (0–1).
- **Comparative analysis:** overlay ZH853 fingerprints against 8F7R (endomorphin-1), 8EFQ/6DDE (DAMGO),
  8F7Q (β-endorphin), and small-molecule agonists — identify **shared anchors vs ZH853-enriched contacts**
  (current candidates: E231/ECL2, H321, N129, macrocycle-extended TM7/ECL contacts).
- **Objective 2 mutation panel:** rank candidate mutations by (a) ZH853-specific contact strength/occupancy
  and (b) predicted sparing of DAMGO/endomorphin/morphine — e.g. E231A/Q, H321A/F, N129A, W320 variants,
  vs the "do-not-touch" shared anchors D149/H299/Y328. Produce a ranked, rationalized table with hypotheses.
- **Gate:** Objective-1 comparative interaction figure + Objective-2 candidate mutation table (feeds Phase 6).

### Phase 4 — Production MD & simulation QC  *(cluster, wall-clock weeks; analysis local)*
- Production: **≥ 3 independent replicas**, several hundred ns each (µs aspirational), for the ZH853
  complex; optionally a DAMGO or endomorphin-1 complex as a comparator ensemble.
- **QC dashboard** (`product/…_qc_YYYYMMDD`): backbone Cα RMSD (plateau), RMSF (pre-aligned),
  **receptor-aligned ligand RMSD**, key-contact occupancy, membrane APL / thickness / S_CD order
  parameters, box/density; error bars from replicate spread + block averaging. Report per references §8.
- **Mutant simulations** for the top Phase-3 candidates (in-silico mutagenesis + MD) to check pocket
  integrity and ZH853-contact disruption.
- **Gate:** documented convergence assessment (not RMSD-plateau alone) and a QC report establishing
  simulation quality before any quantitative claims.

### Phase 5 — Analog cheminformatics & drug-likeness design  *(local, 3–5 days)* — **Objective 3**
- Property panel with RDKit for ZH853 + analogs (already: TPSA 235–280, HBD 8–10, MW 714–810,
  cLogP −0.15…0.73 — all beyond-rule-of-5). Quantify the permeability liabilities.
- **Design strategies** (peptide PK/PD, GLP-1-inspired), each with a concrete enumerated set + rationale:
  backbone **N-methylation** to mask HBDs and rigidify; reduce HBD/TPSA; **lipidation** (C16–C18 fatty-acid
  + γGlu/linker for albumin binding & half-life extension, GLP-1 style); explore stapling/PEGylation;
  retain the D149/Tyr pharmacophore. Predict permeability/EPSA-like descriptors and flag synthetic
  feasibility.
- Feed the most promising analogs to Phase 6 relative FEP.
- **Gate:** analog design report with ranked, property-annotated candidates.

### Phase 6 — Free-energy calculations  *(cluster, high risk — see references §4–6)* — **Objective 4**
Treat as **rank-ordering**, not absolute affinities, unless convergence is demonstrated. Staged by risk:
- **6a Relative FEP (lowest risk):** ZH853 vs ZH850/ZH831/ZH809. These differ by **single-residue edits**
  (Trp↔Phe, Glu↔Asp ring-size, ±Gly cap). R-group-like edits are the tractable regime; **ring-size
  changes (Glu↔Asp) break LOMAP ring-mapping** — use **Kartograf** mapping and expect fragility. Tooling:
  OpenFE (hybrid-topology / SepTop) or PMX. Multiple replicas; enhanced sampling (REST2 on the ligand);
  explicit convergence diagnostics.
- **6b ABFE (higher risk):** OpenFE v1.7 `AbsoluteBindingProtocol` or BAT2; Boresch restraints are
  **ill-defined for a floppy macrocycle** (no rigid 3-atom anchor) — mitigate with restrained/enhanced-
  sampling ABFE (ATM, Lambda-ABF-OPES). Budget heavy free-ligand sampling; report convergence honestly.
- **6c CTMD / funnel metadynamics (exploratory):** PLUMED + openmm-plumed. Funnel metadynamics for
  qualitative ΔG_bind / pose exploration; multiple walkers + richer CVs than a single RMSD (peptide DOF).
  c(t)-scoring is a single-paper method — exploratory only.
- **Design task bundles per calculation:** each a self-contained SLURM directory (inputs, install/prepare/
  submit/gather/post-process, `alchemlyb`/MBAR analysis).
- **Gate:** at minimum, a defensible **rank-ordering** of ZH853 vs analogs with convergence diagnostics;
  quantitative ΔG only if demonstrably converged.

### Phase 7 — Synthesis & manuscript  *(local, ongoing)*
- Assemble into `product/manuscript/`: intro + related work (opioid structures, biased agonism,
  cyclic-peptide PK), methods (full protocols + FF + convergence), results (Obj 1–4), conclusions.
- High-quality summary figures + quantitative tables; QC appendix.
- **Gate:** internally consistent manuscript draft with every quantitative claim traceable to a script/QC.

---

## 5. Risk register

| Risk | Severity | Mitigation |
|------|----------|------------|
| Macrocyclic-peptide ligand FF parameterization wrong (D-AA, ring, cis/trans) | High | Amber residue-library + QM-RESP route; NMR validation of conformers if available; enhanced sampling |
| Free-energy non-convergence (bound & free ensembles) | High | Rank-ordering framing; REST2/enhanced sampling; replicas; explicit diagnostics; qualitative claims |
| 3.5 Å pose/rotamer uncertainty propagates to interaction claims | Medium | MD refinement; report occupancy w/ error bars; cross-check vs higher-res comparators (5C1M, 7T2G) |
| D2.50 / pocket protonation ambiguity | Medium | Parallel protonation-state systems; constant-pH as fallback |
| Relative FEP ring-size edits (Glu↔Asp) unmappable | Medium | Kartograf 3D mapping; ATM/SepTop; restrict network to R-group-like edits where possible |
| ABFE Boresch anchors ill-defined for floppy macrocycle | Medium | ATM / enhanced-sampling ABFE; careful anchor selection + standard-state correction checks |
| Cluster env drift / reproducibility | Low | Pinned `environment-cluster.yml`; self-contained bundles; logged software versions |

## 6. Considered but deferred
- **FEP+ (Schrödinger)** has the only mature macrocycle machinery (ring open/close moves) but is
  commercial — deferred unless open-source convergence proves inadequate and a license is available.
- **Espaloma-0.3** self-consistent GNN parameterization — promising but unvalidated for D/non-natural
  macrocycles; deferred to a sensitivity check.
- **Full asymmetric plasma-membrane bilayer (PIP2/PSM, physiological cholesterol)** — deferred to a
  sensitivity variant after the POPC:chol 9:1 baseline is established.
- *(Signaling-bias structural analysis is now **in scope** per D-8 — moved out of "deferred".)*

---

## 7. Open questions for the user
See SPECIFICATION.md (OQ-1 … OQ-n). The load-bearing ones: (OQ-1) ZH853 identity/stereochemistry
confirmation; (OQ-2) objective prioritization & depth; (OQ-3) cluster resources & software availability;
(OQ-4) whether signaling bias is in scope.
