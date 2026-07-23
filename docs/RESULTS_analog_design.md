# Results — ZH853 drug-likeness & analog design (Objective 3)

_Living results section. Products carry date codes in `product/`. Last updated 2026-07-22._
Pipeline: `src/05.01.00_analog_properties.py`, `src/05.02.00_design_modifications.py`
(`make analogs design`).

## The liability

All four cyclic-peptide analogs sit firmly in **beyond-rule-of-5** space:

| Analog | MW | TPSA | HBD | HBA | cLogP | intramol H-bonds (3D) |
|---|---|---|---|---|---|---|
| ZH853 | 810 | 280 | 10 | 9 | −0.15 | 6 |
| ZH850 | 753 | 251 | 9 | 8 | 0.73 | 5 |
| ZH809 | 739 | 251 | 9 | 8 | 0.34 | 5 |
| ZH831 | 714 | 235 | 8 | 8 | 0.25 | 4 |

The dominant liability is **high TPSA (235–280 Å²) and HBD count (8–10)** — far above the oral Ro5
limits (TPSA 140, HBD 5). Cyclization + D-amino acids already give protease resistance, so the two
remaining gaps are **(1) passive permeability** and **(2) plasma half-life**. Encouragingly, the 3D
conformers show 4–6 **intramolecular H-bonds**, i.e. some backbone donors self-shield ("molecular
chameleon" behavior) — the property that lets macrocycles of this size be permeable at all.
Figure: `product/05.01.00_analog_property_space_*.png`.

## Structure-based derivatization map

From the bound pose, **46 of 59 ZH853 heavy atoms are buried** (the pharmacophore: Tyr1 amine →
D149 salt bridge, aromatic pocket contacts) and **13 are solvent-exposed** — the safe handles for
conjugation. This directly couples Objective 3 to Objective 1: derivatize the exposed C-terminal-cap
region; do **not** touch the buried Tyr1 message or the backbone amides that H-bond D149/Y328/E231.

## Two orthogonal PK axes and the recommended edits

**Axis 1 — passive permeability / oral-CNS exposure** (lower TPSA & HBD):
- **N-methylation of solvent-facing backbone amides** is the highest-value lever. Each N-methyl removes
  one donor and ~9 Å² TPSA; a **2–3 site subset** brings TPSA to ~253 Å² (below the ~250 oral-bRo5
  threshold) and HBD to 7 while limiting the potency risk. hexa-N-methyl (TPSA 227, HBD 4) is the
  property-space extreme but likely perturbs the bound conformation — test subsets, rank by FEP.
- **4-F-Phe** (para-fluoro-Phe): near property-neutral; blocks Phe para-hydroxylation (metabolic
  stability) at low structural risk.

**Axis 2 — plasma half-life / duration** (albumin binding, GLP-1/semaglutide strategy):
- **C-terminal lipidation on the exposed cap**: a C16 palmitoyl (minimal) or the validated
  **γ-Glu-2×AEEA-C18-diacid** albumin-binder. This *raises* size/polarity (keeps the molecule
  injectable, not oral) but extends half-life from minutes toward days. Orthogonal to axis 1.

## Recommendation
- **Permeability series:** start from **ZH831** (least liable) + a 2–3-site N-methyl scan (solvent-facing
  amides only), optionally + 4-F-Phe.
- **Long-acting series:** keep **ZH853** (most potent) + C-terminal semaglutide-style lipidation.
- Both series feed **Objective 4 relative FEP** to confirm the edits preserve MOR affinity.

## Caveats
Descriptor-level predictions only. Macrocycle permeability is conformation-dependent; confirm with
**3D-PSA over conformer ensembles** and experimental PAMPA/Caco-2. The N-methyl site selection depends
on which backbone amides face solvent vs the receptor — resolve this precisely with **MD occupancy
(Phase 4)** before synthesis. Any pharmacophore-adjacent edit must be potency-checked by FEP.
