# Structural basis of ZH853 recognition by the µ-opioid receptor: a distinctive extracellular-loop contact defines a selectivity handle and guides analog design

**Authors:** M. J. O'Meara¹ *et al.*
¹ Department of Computational Medicine and Bioinformatics, University of Michigan

**Status:** working draft, 2026-07-22. Results for Objectives 1–3 are complete; the molecular-dynamics
and free-energy validation (Methods §2.5) is prepared but not yet run — claims dependent on it are
flagged **[prospective]**.

---

## Abstract

The cyclic endomorphin-1 analog ZH853 is a potent µ-opioid receptor (MOR) agonist with an
attractive preclinical profile, but the structural basis of its binding — and how it differs from
other opioid agonists — has not been characterized. Using a 3.5 Å cryo-EM structure of the
MOR–ZH853–Gi–scFv16 complex, we computed heavy-atom interaction fingerprints for ZH853 and benchmarked
them against thirteen deposited MOR complexes spanning peptide agonists (DAMGO, endomorphin-1,
β-endorphin), morphinan and fentanyl-class small molecules, biased agonists, and an antagonist. ZH853
engages the canonical opioid anchor (a salt bridge to Asp3.32 satisfied by its Tyr1 α-amine, plus the
Tyr7.43 "3–7 lock") shared by all agonists, but additionally forms a **salt bridge to Glu231 in
extracellular loop 2 (ECL2) that no other agonist in the set makes**, together with an
aromatic/cation-π contact to His7.36 shared only with its endomorphin-1 parent. This extended ECL2/ECL3
engagement is the structural signature of the cyclic scaffold. We nominate **E231Q/E231A** as
mutations predicted to disrupt ZH853 selectively while sparing related full agonists, and identify the
molecule's beyond-rule-of-5 liabilities (TPSA 235–280 Å², 8–10 H-bond donors), proposing two orthogonal
optimization series — backbone N-methylation for passive permeability and semaglutide-style C-terminal
lipidation for half-life — targeted to solvent-exposed positions established from the bound pose. A
membrane molecular-dynamics and free-energy protocol to test these hypotheses is established and
released. This is, to our knowledge, the first structural analysis of a cyclic-peptide MOR complex.

---

## 1. Introduction

Opioid analgesics acting at the µ-opioid receptor (MOR, gene *OPRM1*) remain the mainstay of severe
pain management, but their utility is limited by respiratory depression, tolerance, and dependence.
A central goal of contemporary opioid pharmacology is to retain analgesic efficacy while separating it
from these liabilities — through signaling bias, partial agonism, or novel chemotypes. Endomorphins,
the endogenous MOR-selective tetrapeptides (Tyr-Pro-Trp-Phe-NH₂), are attractive scaffolds because of
their high selectivity, but their linear form is metabolically labile and poorly bioavailable. The
Zadina laboratory addressed this by engineering cyclic, D-amino-acid endomorphin analogs; **ZH853**
(Tyr-cyclo[D-Lys-Trp-Phe-Glu]-Gly-NH₂) is a lead from this series with potent antinociception and a
reduced side-effect profile in rodents [Zadina 2016].

Despite this promise, no structure of ZH853 — or of any cyclic-peptide MOR complex — has been reported,
leaving open three questions this work addresses computationally from a 3.5 Å cryo-EM structure of the
MOR–ZH853–Gi–scFv16 complex: (1) **what, if anything, is structurally distinctive about how ZH853
engages MOR** relative to other agonists; (2) **which receptor mutations would abrogate ZH853 while
sparing related full agonists**, providing a pharmacological fingerprint and a route to mechanistic
tests; and (3) **how ZH853's drug-like properties could be improved**, drawing on modern peptide-drug
strategies. We further establish and release a membrane molecular-dynamics (MD) and free-energy
protocol to validate these structure-derived hypotheses.

## 2. Related work

**MOR structural pharmacology.** The inactive MOR structure (β-FNA, 4DKL [Manglik 2012]) and the active
BU72 complex (5C1M [Huang 2015]) established the orthosteric pocket and the activation-associated
rearrangements of the W6.48 toggle and TM6. Cryo-EM of MOR–Gi complexes then resolved a series of
agonists — the peptide DAMGO (6DDE [Koehl 2018]; human 8EFQ) and endomorphin-1/β-endorphin
(8F7R/8F7Q [Wang 2023]), and the small molecules fentanyl, morphine, oliceridine, SR-17018, PZM21
(8EF5/8EF6/8EFB/8EFL/8EFO [Zhuang 2022]) and mitragynine pseudoindoxyl/lofentanil (7T2G/7T2H
[Qu 2023]). All share the D3.32 anchor and a conserved hydrophobic subpocket; biased agonists differ
mainly in TM6/ECL2 engagement. Notably, **every deposited peptide-agonist MOR structure carries a
*linear* peptide** — ZH853 would be the first cyclic-peptide complex.

**Cyclic peptides as drugs.** Macrocyclic peptides occupy "beyond-rule-of-5" (bRo5) chemical space
[Doak 2014]; their oral/passive permeability is governed less by static polarity than by
conformational shielding of backbone donors ("molecular chameleons" [Whitty 2016]), which
N-methylation and cyclization promote (as in cyclosporine). Half-life, by contrast, is engineered
through albumin-binding **lipidation** — the strategy that gives semaglutide a ~1-week half-life via a
γGlu-2×AEEA-C18-diacid conjugate [Knudsen 2019]. These two axes — permeability and duration — are
distinct and are treated separately here.

## 3. Methods

**3.1 Structure and numbering.** We used the deposited real-space-refined MOR–Gi–scFv16–ZH853 model
(3.5 Å; CC_mask 0.80; 0.0% Ramachandran outliers). Chains: Gαi (A), Gβ1 (B), Gγ2 (C), scFv16 (D),
MOR/*OPRM1* residues 69–349 with 84 modeled cholesterols (R), and ZH853 as HETATM ligand L01 (E). The
construct uses **human OPRM1 (UniProt P35372) numbering**, verified against all canonical orthosteric
positions; mouse-numbered comparators were offset by +2 to this frame. The ligand's 59 heavy atoms and
formula (C₄₂H₅₁N₉O₈, 810 Da) match the reference SMILES, confirming the monomeric macrocycle.

**3.2 Comparator set.** Thirteen structures were retrieved from the RCSB PDB: peptide agonists
6DDE/8EFQ (DAMGO), 8F7R (endomorphin-1), 8F7Q (β-endorphin), 9WST/9WSV (DAMGO with Gz/β-arrestin);
small molecules 5C1M (BU72), 8EF5 (fentanyl), 8EFB (oliceridine), 8EFL (SR-17018), 8EFO (PZM21),
7T2G (mitragynine pseudoindoxyl); and antagonist 4DKL (β-FNA). For each, the receptor chain and bound
agonist were isolated programmatically and mapped to human OPRM1 numbering.

**3.3 Interaction fingerprints.** Because the 3.5 Å model lacks hydrogens, we used a resolution-matched
**heavy-atom geometric fingerprint** (rather than angle-based H-bond criteria), classifying each
receptor–ligand residue contact as ionic (opposite-charge groups ≤4.0 Å), H-bond (polar N/O ≤3.5 Å),
hydrophobic (C–C ≤4.0 Å), π-stacking (aromatic centroids ≤5.5 Å), or cation-π (≤6.0 Å). Aromatic rings
of non-standard ligands were detected by planar-ring geometry, a template-free method that generalizes
across chemotypes. The same method was applied uniformly to all fourteen complexes.

**3.4 Mutation panel and cheminformatics.** Candidate selectivity mutations were scored by ZH853
interaction strength × rarity across the twelve agonists × a bonus for sparing the two closest full
agonists (endomorphin-1, DAMGO). Physicochemical descriptors (2D, plus 3D radius of gyration, polar
SASA fraction, and intramolecular H-bond count from macrocycle-aware conformers) were computed with
RDKit; modification variants (N-methylation, lipid conjugation, halogenation) were enumerated
programmatically and their properties predicted. Buried vs solvent-exposed ligand atoms were
determined from the bound pose (≤4.5 Å to receptor).

**3.5 MD/free-energy protocol [prospective].** A membrane MD system was prepared: receptor sidechains
rebuilt (PDBFixer), protonation assigned by PROPKA at pH 7.4, and ZH853 protonated to +1 and prepared
for GAFF2/RESP parameterization. PROPKA places **D2.50 (Asp116) at pKa 7.61**, so parallel
protonated/charged systems are built. The production protocol (POPC:cholesterol 9:1, ff19SB/Lipid21/OPC,
CHARMM-GUI-style six-stage equilibration, LangevinMiddle 310 K with a semi-isotropic membrane barostat,
hydrogen-mass repartitioning at 4 fs, ≥3 replicas) and free-energy calculations (relative FEP across
the ZH850/831/809 analogs; ABFE and funnel metadynamics as exploratory) are released as a SLURM bundle.
Given documented convergence difficulties for flexible macrocycles, free-energy results are framed as
rank-ordering.

*All analyses are scripted and version-controlled (Data & code availability).*

## 4. Results

### 4.1 ZH853 satisfies the conserved opioid anchor

ZH853 contacts 25 receptor residues. Its Tyr1 α-amine forms the **canonical Asp3.32 (Asp149) salt
bridge** (3.16 Å) and hydrogen-bonds Tyr7.43 (Tyr328, 3.17 Å) — the "3–7 lock" — while the Tyr1 ring
and the macrocyclic Trp/Phe pack into the conserved hydrophobic subpocket (Met3.36, Val3.28/Ile3.29,
Ile6.51, Val6.55, Val5.42). Every one of these contacts is shared by all twelve comparator agonists
(Figure 1), confirming that ZH853 presents the same "message" pharmacophore as morphine, fentanyl,
DAMGO, and endomorphin-1. Consistent with this, PROPKA assigns Asp149 a depressed pKa (6.0), i.e.
charged and salt-bridge-competent.

### 4.2 A distinctive ECL2 salt bridge distinguishes ZH853 (Objective 1)

Against this conserved background, three contacts set ZH853 apart (Figure 1; Table 1):

- **Glu231 (ECL2): an ionic + H-bond contact made by ZH853 alone (0/12 agonists).** A ZH853 basic
  group reaches ECL2 to form a salt bridge (3.12 Å) that no linear peptide or small molecule in the set
  reproduces. PROPKA confirms Glu231 is charged (pKa 4.9). This is the single strongest structural
  discriminator.
- **His7.36 (His321): aromatic/cation-π, shared only with endomorphin-1 and mitragynine (2/12).** An
  ECL3-proximal contact retained from the endomorphin-1 lineage.
- Weaker ZH853-unique contacts at Leu221/Phe223 (ECL2) and Leu234 (TM5), gained by the extended/cyclic
  scaffold.

The picture is coherent: **ZH853 keeps the deep conserved anchor but reaches "up and out" toward
ECL2/ECL3**, engaging a rim region that the compact linear agonists do not. This extended loop
engagement is the structural hallmark of the macrocycle.

### 4.3 A selectivity handle: E231 mutations (Objective 2)

Ranking pocket residues by ZH853 interaction strength, rarity across agonists, and sparing of the
endomorphin-1/DAMGO parents nominates **Glu231 as the top selective target** (Table 1). We predict
**E231Q** (charge-neutral isostere) and **E231A** (removal) will attenuate ZH853 binding while leaving
agonists that do not contact ECL2 Glu231 largely unaffected — a testable pharmacological signature.
His7.36 substitutions (H321A/F) are a secondary, less clean handle (they also perturb endomorphin-1).
By contrast, the universal anchors (Asp3.32, Tyr7.43, Tyr3.33, Met3.36, and the hydrophobic cage) are
**do-not-mutate controls**: disrupting them abrogates all agonists and cannot confer selectivity.
**[prospective]** magnitudes await MD occupancy and relative-FEP/in-silico-mutagenesis validation.

### 4.4 Drug-likeness liabilities and a two-axis design strategy (Objective 3)

All four analogs sit firmly beyond rule-of-5 (MW 714–810, TPSA 235–280 Å², 8–10 H-bond donors; Figure
2); the dominant liabilities are high polar surface area and donor count. Their 3D conformers show
4–6 intramolecular H-bonds — partial donor self-shielding consistent with chameleonic behavior. From
the bound pose, **46 of 59 ligand heavy atoms are buried** (the Tyr1/aromatic pharmacophore) and **13
are solvent-exposed**, defining where derivatization is structurally safe. We therefore propose two
orthogonal series (Figure 3):

1. **Permeability** — backbone **N-methylation** of solvent-facing amides (each removes one donor and
   ~9 Å² TPSA; a 2–3-site subset reaches TPSA ≈253 Å², crossing the oral-bRo5 threshold), optionally
   with 4-F-Phe for metabolic stability. Best pursued from **ZH831**, the least liable analog.
2. **Half-life** — **C-terminal semaglutide-style lipidation** on the exposed cap (γGlu-2×AEEA-C18-diacid),
   trading polarity for albumin-mediated duration; pursued from the most potent analog, **ZH853**.

Both series must be potency-checked against MOR by relative FEP, and N-methyl sites chosen to avoid the
receptor-facing backbone amides identified in §4.1–4.2.

## 5. Discussion

The central finding is that ZH853's distinctiveness lies not in the conserved orthosteric anchor — which
it shares with all opioid agonists — but in a **cyclic-scaffold-enabled ECL2 salt bridge (Glu231)** and
adjacent ECL3 aromatic contact (His7.36). ECL2 is a recognized determinant of ligand selectivity and
kinetics across GPCRs, and its engagement here provides a plausible structural correlate for ZH853's
distinct pharmacology and a concrete, low-ambiguity selectivity handle (E231Q). Because Glu231 is
peripheral to the conserved message, its mutation is well suited to a clean pharmacological separation
of ZH853 from morphinan/fentanyl chemotypes and even from its linear endomorphin-1 parent.

The design analysis reframes ZH853 optimization as two independent problems. The molecule already
banks the hard-won macrocyclic and D-amino-acid protease resistance; what remains is permeability
(a donor/PSA problem, addressable by N-methylation guided by which amides face solvent vs receptor) and
duration (a half-life problem, addressable by lipidation on the exposed cap). Coupling the design to the
bound pose avoids the common error of derivatizing the pharmacophore.

**Limitations.** All results derive from a single 3.5 Å static model; sidechain rotamers and the exact
ligand pose carry real uncertainty, and interaction assignments use heavy-atom geometry rather than
explicit hydrogens. The comparative analysis is a contact census, not an energetic ranking. The
mutation and design predictions are hypotheses; their magnitudes require the MD/FEP validation whose
protocol we establish here, and macrocycle free-energy convergence is itself a known challenge.

## 6. Conclusions and future work

From the first structural analysis of a cyclic-peptide MOR complex, we identify an ECL2 Glu231 salt
bridge as the defining, ZH853-specific interaction, nominate E231Q/E231A as selectivity-conferring
mutations, and propose bound-pose-guided N-methylation and lipidation series to address the molecule's
permeability and half-life liabilities. The immediate next step is to run the released membrane-MD and
free-energy protocol to (i) confirm contact occupancy and pose stability, (ii) quantify the E231
selectivity, and (iii) rank the analog and designed-variant affinities. Longer term, the Gz and
β-arrestin comparators included here support a follow-on analysis of whether ECL2/ECL3 engagement
correlates with ZH853's reported signaling bias.

## Figures

- **Figure 1.** Pocket interaction fingerprint of ZH853 vs 13 MOR complexes (`figures/fig1_interaction_heatmap.png`).
  Rows = receptor residues (human OPRM1/BW); columns = complexes; cells colored by interaction strength
  and labeled by type (I ionic, P cation-π, A π-stack, H H-bond, h hydrophobic). Asp3.32 is engaged
  across all columns; Glu231 (ECL2) is engaged in the ZH853 column alone.
- **Figure 2.** ZH853 analogs in Ro5/bRo5 property space (`figures/fig2_property_space.png`).
- **Figure 3.** Predicted property shifts of modification strategies by PK axis (`figures/fig3_design_shifts.png`).

## Table 1 — ZH853 contacts ranked by distinctiveness (excerpt)

| Residue | BW | ZH853 interaction | # of 12 agonists sharing | note |
|---|---|---|---|---|
| Glu231 | ECL2 | ionic + H-bond | **0** | ZH853-unique; primary selectivity handle (E231Q/E231A) |
| His321 | 7.36 | π-stack + cation-π | 2 | shared only with endomorphin-1, mitragynine |
| Tyr130 | 2.64 | hydrophobic | 1 | endomorphin-1 only |
| Leu221/Phe223 | 45.52/45.54 | contact | 0 | weak ECL2 contacts (cyclic scaffold) |
| Asp149 | 3.32 | ionic + H-bond | 12 | universal anchor — do-not-mutate control |
| Tyr328 | 7.43 | H-bond | 12 | universal (3–7 lock) — control |

*(Full matrix: `product/03.01.00_interaction_matrix_20260722.csv`; full panel:
`product/03.02.00_mutation_panel_20260722.md`.)*

## Data and code availability

All analyses are scripted and version-controlled in the project repository (`src/`, numbered
workflows; `zh853mor` package; `make analysis` reproduces every figure and table). Comparator
structures are re-fetchable via `make fetch`. Decisions and clarifications are logged in
`SPECIFICATION.md`; detailed methods in `docs/METHODS_md_prep.md`, `docs/RESULTS_interactions.md`,
and `docs/RESULTS_analog_design.md`.

## References

1. Zadina JE, et al. Endomorphin analog analgesics with reduced abuse liability. *Neuropharmacology*
   2016. doi:10.1016/j.neuropharm.2015.12.024
2. Manglik A, et al. Crystal structure of the µ-opioid receptor bound to a morphinan antagonist.
   *Nature* 2012;485:321. (4DKL)
3. Huang W, et al. Structural insights into µ-opioid receptor activation. *Nature* 2015;524:315. (5C1M)
4. Koehl A, et al. Structure of the µ-opioid receptor–Gi protein complex. *Nature* 2018;558:547. (6DDE/6DDF)
5. Zhuang Y, et al. Molecular recognition of morphine and fentanyl by the µ-opioid receptor. *Cell*
   2022;185:4361. (8EF5/8EF6/8EFB/8EFL/8EFO/8EFQ)
6. Qu Q, et al. Insights into distinct signaling profiles of the µOR. *Nat Chem Biol* 2023;19:423. (7T2G/7T2H)
7. Wang Y, et al. Structures of the entire human opioid receptor family (endomorphin-1/β-endorphin).
   *Cell* 2023;186:413. (8F7R/8F7Q)
8. Zhang Y, et al. Structures of MOR with Gz and β-arrestin. *Cell Res* 2025;35:1021. (9WST/9WSV)
9. Doak BC, et al. Oral druggable space beyond the rule of 5. *Chem Biol* 2014;21:1115.
10. Whitty A, et al. Quantifying the chameleonic properties of macrocycles. *Drug Discov Today* 2016.
11. Knudsen LB, Lau J. The discovery and development of semaglutide. *Front Endocrinol* 2019;10:155.
