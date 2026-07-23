# Manuscript

Working draft synthesizing the completed analyses (Objectives 1–3) plus the established
MD/free-energy methods.

- `manuscript_20260722.md` — the draft in Markdown (intro, related work, methods, results,
  discussion, conclusions, references).
- `manuscript.tex` — LaTeX source (article class, `authblk`, `booktabs`, embedded figures,
  `thebibliography`). Compile with **`make manuscript`** (uses `tectonic`, which auto-fetches
  packages) → `manuscript.pdf` (~7 pp).
- `manuscript.pdf` — compiled PDF.
- `figures/` — figures copied from dated `product/` outputs:
  - `fig1_interaction_heatmap.png` ← `product/03.01.00_fingerprint_heatmap_*.png`
  - `fig2_property_space.png` ← `product/05.01.00_analog_property_space_*.png`
  - `fig3_design_shifts.png` ← `product/05.02.00_design_property_shifts_*.png`

## Provenance / regeneration
Every figure, table, and number traces to a scripted analysis. Regenerate the underlying
products with `make analysis` (interactions + mutations + analogs + design) and `make prep`
(Phase-2 methods numbers), then refresh the copies in `figures/`.

## Status of claims
Results for Objectives 1–3 are complete and computed. Claims marked **[prospective]** in the
text depend on the molecular-dynamics / free-energy validation (Methods §2.5 / §3.5), whose
protocol is prepared and released (`src/02.10.00_slurm_bundle/`) but **not yet run**. Do not
present [prospective] items as established results.

## To do before submission
- Author list, affiliations, and corresponding-author details (currently placeholder).
- Run the MD/FEP protocol and fold in occupancy / ΔΔG results, upgrading the [prospective] claims.
- PLIP/ProLIF cross-validation of the interaction fingerprints (SPECIFICATION D-9).
- Resolve analog-SMILES naming (SPECIFICATION OQ-5) in the Methods if analogs are foregrounded.
- Final reference formatting to target-journal style; verify every PDB ID and DOI.
