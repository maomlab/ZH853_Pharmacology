# Results — ZH853 interaction analysis (Objectives 1 & 2)

_Living results section. Generated products carry date codes in `product/`; this file
synthesizes the current interpretation. Last updated 2026-07-22._

Method: heavy-atom geometric interaction fingerprints (no modeled hydrogens; resolution-matched
to 3.5 Å) computed uniformly for the ZH853 complex and 13 comparator structures, all mapped onto
human OPRM1 (P35372) numbering. Pipeline: `src/03.01.00_interaction_fingerprints.py`,
`src/03.02.00_mutation_panel.py`. See [PLAN](PLAN.md) and [references](references.md).

## Objective 1 — what is structurally unique about ZH853 binding

ZH853 contacts **25 receptor residues**. Benchmarked against 12 agonist complexes (endomorphin-1,
DAMGO ×4 incl. Gz/arrestin, β-endorphin, BU72, fentanyl, oliceridine, SR-17018, PZM21,
mitragynine-pseudoindoxyl) plus the β-FNA antagonist, the contacts partition cleanly:

**Universal anchors (shared by all 12 agonists) — the conserved opioid pharmacophore:**
- **D149 (3.32)** — ionic salt bridge to the ligand protonated amine (dark across the entire heatmap row).
- **Y328 (7.43), Y150 (3.33)** — the 3–7 lock / amine H-bond network.
- Hydrophobic cage: **M153 (3.36), V145/I146 (3.28/3.29), I298 (6.51), V302 (6.55), V238 (5.42), Q126 (2.60)**.

These reproduce the canonical MOR orthosteric recognition seen in every deposited agonist complex —
i.e. ZH853 satisfies the same "message" anchor as morphine, DAMGO, and endomorphin-1.

**ZH853-distinctive contacts (the interesting part):**

| Contact | BW | ZH853 interaction | Shared with | Significance |
|---|---|---|---|---|
| **E231** | ECL2 | **ionic + H-bond** | **0/12 — unique** | A second salt bridge from a ZH853 basic group to ECL2 that *no other agonist forms*. The single strongest distinguishing feature. |
| **H321** | 7.36 | **π-stack + cation-π** | endomorphin-1, mitragynine (2/12) | An aromatic/ECL3 contact shared only with its direct parent and one biased alkaloid. |
| Y130 | 2.64 | hydrophobic | endomorphin-1 (1/12) | Upper-pocket contact retained from the endomorphin-1 scaffold. |
| L221, F223 | 45.52/45.54 | contact-only | 0/12 — unique | ECL2 contacts gained by the extended/cyclized ZH853 scaffold (weak). |
| L234 | 5.38 | contact-only | 0/12 — unique | TM5 contact absent in linear agonists (weak). |

**Interpretation:** ZH853 keeps the conserved D3.32 anchor but reaches "up and out" toward ECL2/ECL3
(E231, H321, L221/F223), a region the linear peptides and small molecules do not engage. This extended
ECL2 engagement is the structural signature of the cyclic scaffold and the most likely origin of any
distinct pharmacology. Figure: `product/03.01.00_fingerprint_heatmap_*.png`.

## Objective 2 — mutations to abrogate ZH853 while sparing related full agonists

Candidates ranked by ZH853 interaction strength × rarity among agonists × sparing of the two
closest full agonists (endomorphin-1, the parent, and DAMGO). Full table:
`product/03.02.00_mutation_panel_*.md`.

**Primary recommendation — E231 (ECL2):**
- **E231Q** (charge-neutral isostere) and **E231A** (removal). ZH853 forms an ionic + H-bond contact
  here that *none* of the 12 agonists share, so loss should be ZH853-selective. E231Q is the cleaner
  test (removes charge, preserves size/H-bonding geometry), E231A the harder knock-out.

**Secondary candidates (with explicit caveats):**
- **H321A / H321F (7.36):** removes the π/cation-π contact, but H321 is *also* engaged by endomorphin-1,
  so expect partial loss of the parent — a discriminating-but-not-clean test.
- **L221A / F223A / L234A:** ZH853-unique but only weak (contact-only) contacts; likely low-magnitude
  and less specific. Useful as a combination/secondary panel.

**Do-not-mutate controls (universal anchors):** D149, Y328, Y150, M153, Q126, I146, V145, I298, V302,
V238 — mutating these abrogates all agonists and should be used as negative (non-selective) controls.

## Caveats & next steps
- These are **static, single-model** inferences from a 3.5 Å structure; sidechain rotamers and the
  ligand pose carry real uncertainty. Every claim above should be confirmed by **MD contact-occupancy**
  (Phase 4) — a contact present in one model may be transient — and the selectivity of E231Q/H321F
  quantified by **relative FEP / in-silico mutagenesis** (Phase 6).
- The geometric fingerprint should be **cross-validated with PLIP and ProLIF** once the analysis conda
  env is built (PLIP needs OpenBabel, unavailable in the current base env; see SPECIFICATION D-9).
- β-endorphin (a 21-mer) makes many extra-orthosteric contacts (e.g. R213/D218 salt bridges) that are a
  property of its length, not directly comparable to the compact ZH853 macrocycle.
