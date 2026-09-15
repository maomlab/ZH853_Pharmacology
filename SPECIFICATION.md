# Specification of Deliverables

This file records decisions, clarifications, and conventions made as the ZH853–MOR analysis proceeds
(per OBJECTIVES.md). Decisions (D-n) are resolved; Open Questions (OQ-n) await user input.

## Conventions
- **Residue numbering:** human OPRM1 (UniProt **P35372**), which matches the deposited construct
  (chain R). Mouse-numbered comparators (4DKL, 5C1M, 7T2G, 7SBF) are **−2** relative to ours.
- **Ballesteros-Weinstein** generic numbers used alongside construct numbers throughout.
- **Directory/versioning:** raw inputs → `data/`; code → `src/##.##.##_name`; caches →
  `intermediate/`; outputs → `product/…_YYYYMMDD.ext`; plan/decisions → `docs/`.
- **`src/` holds only workflows; everything generated goes to `intermediate/`** (D-16). Scripts run
  from `src/` must never write beside themselves — a build directory under `src/` is both untracked
  clutter and, for PACKMOL-Memgen, an active hazard (it reuses component PDBs it finds in its CWD).
  This includes **SLURM job logs**: a stage submitted from `src/` directs `--output` into the
  intermediate directory of the step that produces the data (see D-22).
- **Naming:** a workflow step is named for what it *does* — `##.##.##_<verb>_<object>`
  (`02.03.00_prepare_receptor.py`, `02.11.00_analyze_simulations/`) — and the intermediate
  directory it writes to carries the **same** name (`intermediate/02.11.00_analyze_simulations/`),
  so a step, its cache and its products are found by one string.

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
- **D-13 (three cluster conda envs):** The cluster stack is split into **`zh853mor-prep`** (AmberTools/
  PACKMOL-Memgen, CPU), **`zh853mor-sim`** (openmm + openmmforcefields, GPU: equilibration/production/
  FEP), and **`zh853mor-plumed`** (openmm 8.4 + openmm-plumed, GPU: metadynamics only). Reason:
  `openmmforcefields>=0.16` needs openmm>=8.5.1 while `openmm-plumed` (latest) supports only openmm<=8.4,
  so they are unsatisfiable together; splitting by task resolves it and keeps every capability. Python is
  left unpinned so conda tracks the latest openmm (needs >=3.12). (2026-07-24)
- **D-14 (membrane placement — supersedes the cholesterol-primary approach):** Production orientation uses
  **OPM/PPM** (transfer-energy minimization, the community standard); hydrophobic thickness **~32 Å for MOR**
  (OPM 4DKL 32.0±1.0 Å; class-A GPCRs 31–35 Å), the value built to. The local `02.05.00` cholesterol-centred
  orientation is a first-pass proxy only (3 site-specific cholesterols fix the midplane to ~2 Å; their ~28 Å
  span ≈ the POPC hydrocarbon core 2Dc=28.8 Å [Kučerka 2011] and is thin). Placement is cross-checked in
  `03.04.00` against the **Trp/Tyr aromatic girdle** (~30 Å, agrees with OPM) and experiment. Reason: bound
  cholesterols are a weak, biased ruler; OPM/PPM + the aromatic belt are the recognized methods. (2026-07-25)

- **D-15 (receptor finalisation happens in prep, not assembly):** His tautomers and neutral ACE/NME
  termini are written into the receptor by `02.03.00`, not left to the membrane-builder/tleap step.
  Reason: both were silently overridden by downstream defaults in the first build. tleap maps a residue
  named `HIS` to **HIE** regardless of the tautomer determined upstream — the 2026-07-26 build got HIE
  at all four sites although OpenMM's H-bond-network assignment (the standard method, run inside
  PDBFixer) had chosen **HID for H225, H299/H6.52 and H321/H7.36**, i.e. three wrong, including both
  pocket histidines. Likewise a chain starting/ending on a standard residue is built with charged
  termini (`NTHR`/`CPHE`), adding two formal charges that full-length OPRM1 (69–349 of 400 aa) does not
  have. Neither raises an error, so both are now fixed in the receptor and verifiable in
  `product/02.03.00_receptor_prep_*.md`. Cap torsions are chosen by a clash scan (closest cap–protein
  contact 4.3 Å). (2026-07-26)
- **D-16 (generated artefacts live under `intermediate/`):** `src/` contains workflows only; build
  directories, packing intermediates and staged inputs go to `intermediate/##.##.##_name/`. See
  Conventions. (2026-07-26)
- **D-17 (D2.50 variant is a build-time rename):** The ASH116 system is produced by
  `D250=ASH ./01_build_system.sh`, which renames the residue in the staged receptor, rather than by
  duplicating the prep/orient chain. The two variants differ only in protonation, so sharing one
  prepared and OPM-oriented receptor guarantees the comparison is not confounded by a differing
  starting geometry (D-11). (2026-07-26)
- **D-18 (production analysis is a two-machine pipeline):** Trajectory reduction
  (`02.11.00/01_reduce_trajectory.py`) runs **on the cluster**, one SLURM array task per replica, and
  writes a few hundred kB per replica; aggregation and figures run **locally** on those files. Reason:
  the panel is ~100 GB of DCD (3 replicas x 500 ns x up to 10 systems at 100 ps sampling), which is
  not worth moving, while every table and figure must be re-derivable in seconds without a
  trajectory. Consequence: anything the figures may ever want must be computed in the single
  reduction pass, so the reduction stores the full per-residue contact matrix, not just the
  anchors. (2026-09-12)
- **D-19 (MD residue numbers come from `receptor.pdb`, not the prmtop):** tleap renumbers the system
  from 1, so a prmtop residue id is not the human OPRM1 number the rest of the project speaks
  (D149, E231, H299...). `zh853mor.md.construct_residue_map` transfers the numbering positionally
  from the staged `receptor.pdb`, which keeps 69-349, and **refuses** if the residue names disagree.
  Reason: analysing by raw prmtop resid reports the wrong residues while still producing plausible
  numbers -- there is no error to notice. This is also why the 4.5 A shell is recomputed for every
  residue rather than only for a hard-coded key set. (2026-09-12)
- **D-20 (convergence is judged on four diagnostics, not an RMSD plateau):** effective sample size
  after automatic equilibration detection (Chodera 2016), a blocking curve, R-hat across replicas,
  and the Hess cosine content of the leading PCs. Thresholds: n_eff >= 20 per replica,
  R-hat <= 1.2, PC1 cosine content <= 0.5. Reason: a plateaued, correlated series yields error bars
  that are too small by sqrt(g), with g routinely 10-100 for pocket observables at 100 ps sampling;
  and the first PC of a too-short run reproduces the half-cosine of free diffusion, which looks
  exactly like a slow collective motion. Consequence: with 3 replicas R-hat is a flag, not a
  measurement, and the QC report says so. (2026-09-12)
- **D-21 (MD contact criterion is the static one):** occupancy uses the same heavy-atom cutoffs as
  `zh853mor.interactions` (contact 4.5 A, polar/H-bond 3.5 A, ionic 4.0 A; D-9), imported rather
  than restated, so an MD occupancy and a cryo-EM fingerprint entry mean the same thing and
  Objective 1 can be read across the two. Water-mediated contacts are counted separately (a water O
  within 3.5 A of both partners) because the canonical H6.52 interaction is water-bridged and a
  direct-contact count alone scores it as absent. (2026-09-12)
- **D-22 (SLURM logs live with the data, not in `src/`):** `--output` is resolved relative to the
  submission directory, so the stages submitted from a `src/` directory (`./submit.sh params` and
  `build`, and `02.11.00`'s `submit_reduce.sh`) write their `.out` files into
  `intermediate/<the step's directory>/logs/`; the GPU stages, submitted from a build directory,
  already land beside their trajectory and are unchanged. The submit scripts create the directory
  first, because SLURM cannot start a job whose output file it cannot open, and each `.sbatch`
  header names the same path relative to its own directory so a hand-submitted job behaves the
  same. Reason: D-16 — `src/` is source, and a job array otherwise drops one `.out` per task into
  it. (2026-09-12)

- **D-23 (trajectory movies are a first-class analysis step, and a movie is an EXPORTED FILE
  before it is a video):** `02.12.00` writes a `<replica>_movie.{pdb,xtc,json}` triple on the
  cluster and renders locally, the same split as `02.11.00` and for the same reason (the DCDs stay
  where they are). The exported pair is the deliverable: it opens in VMD/PyMOL/ChimeraX/MolStar,
  so the two rendered movies are conveniences on top of it rather than the only way to look.
  Reason: exploratory inspection is how failures that nobody thought to measure get found -- a
  lipid tail entering the orthosteric site, ECL2 peeling off the lid, the receptor sliding along
  the membrane normal -- and all of those reach `02.11.00` only as "the RMSD went up". (2026-09-14)
- **D-24 (what the export fixes, and what it deliberately does not):** frames are unwrapped by
  fragment and re-imaged around the protein (OpenMM wraps molecule by molecule, so the receptor
  and its lipids are torn at the box faces and the ligand can sit a box length from its pocket);
  the annular-lipid selection is the UNION over frames sampled across the whole run, because a
  trajectory file has one atom count for every frame and "within 8 A" cannot be re-evaluated per
  frame; pocket waters are instead RANKED by occupancy and truncated, because thousands pass
  through in 500 ns and drawing them all hides the few that bridge a contact. Superposition is on
  the same membrane-embedded CA set `02.11.00` measures its TM RMSD on. `--no-align` skips the
  rotational fit *and* the z-centring, so drift, tilt and membrane registration stay visible,
  while x and y are still centred (the box origin is arbitrary in the membrane plane). Both
  renderers read the receptor/ligand chain assignment out of the JSON rather than re-deriving it,
  so they cannot disagree about which chain is the subject. (2026-09-14)
- **D-25 (the movie carries its own sampling caveat):** the export records `frame_spacing_ns` and
  the caption on every rendered frame carries the simulated time, not just a frame index. Reason:
  a 500 ns replica in 300 frames is one frame per 1.7 ns, and nothing faster than that -- a
  rotamer flip, a water exchange, a transient contact -- is present rather than merely blurred.
  `submit_export.sh` prints the ns/frame column before the array is submitted, so the decision to
  raise `--frames` is made before the compute rather than after the movie disappoints. (2026-09-14)

- **D-26 (the conformational-feature matrices are produced by the reduction, not by a second
  pass):** `02.11.00`'s reduction now also stores, per frame, the pairwise Ca-Ca distances among
  the 36 key/functional residues (`key_pairs`) alongside the existing per-residue ligand distances
  (`min_dist`). Reason: it is the only stage that reads the trajectories, and a separate
  featurisation stage would be a second ~100 GB read for arrays it already has the coordinates to
  compute. Cost: the reduced files grow from a few hundred kB to ~7 MB per replica (~200 MB for
  the panel), which is still an scp rather than a data-management problem. **A reduction produced
  before this change has no `key_pairs` and must be re-run with `--force`**; `02.13.00` says so by
  name rather than falling back to a thinner feature set. Ca rather than sidechain tips, for the
  reason the activation rulers already give: at 3.5 A the rotamers are the least reliable
  coordinates in the model, and featurising them would let tICA find slow modes in the model's
  guesses. Ca distances are also defined for the **apo** arm, which ligand-contact features are
  not — that is what lets apo and holo share one landscape. (2026-09-14)
- **D-27 (landscapes are drawn on ONE basis, fitted once, and the between-system variance is
  reported):** `02.13.00` fits its tICA basis on the pooled post-equilibration frames of every
  system being compared and projects each into it. Fitting per system would give each its own
  axes, and two landscapes drawn on different linear combinations of distances cannot be laid
  beside each other however alike they look. The consequence has to be stated rather than hidden:
  systems do not interconvert, so a direction separating two of them has an autocorrelation of ~1
  at any lag and tICA returns it first with an infinite implied timescale. That is a useful
  DISCRIMINATIVE axis, not a relaxation time. Every component therefore carries the **fraction of
  its variance lying between systems** (flagged above 0.5), and its timescale is re-estimated
  separately within each system, where the coordinate can actually relax. `--fit-on <system>` is
  available for the stricter reading. (2026-09-14)
- **D-28 (a landscape is a sampling density, and is labelled as one):** the heatmaps are
  -kT ln P at 310 K over the frames the trajectories actually visited — not a reweighted or
  converged free energy (D-5). Bins below a frame-count floor are left unshaded rather than drawn
  as a barrier, because -kT ln(1/N) for one stray frame is a confident-looking wall built from a
  single sample; difference maps are shown only where BOTH systems are sampled, with the
  sampled-by-one-only region hatched, so a difference in coverage is never read as a difference in
  energetics. tICA inherits D-20's trap unchanged — on an unconverged replica the slowest apparent
  process is the drift — so `convergence.cosine_content` is reported per replica for the tIC
  projections exactly as it is for the PCs. (2026-09-14)

## Open questions (need user input)
- **OQ-3 (compute environment):** SLURM cluster specs (GPU types/count, wall-time limits, queue), and which
  software is preinstalled vs must be built (OpenMM, PLUMED, OpenFE, phenix/MolProbity, Gaussian/Psi4 for RESP)?
  Needed before Phase 2 SLURM bundles; pins versions in `environment_zh853mor-sim.yml`.
- **OQ-5 (analog set for FEP):** Confirm the analog panel for Objective 4 = ZH850, ZH831, ZH809 (from
  OBJECTIVES), plus any Phase-5 designed analogs; the OBJECTIVES analog names/SMILES have internal Trp/Phe
  inconsistencies (same class of typo as D-6) to reconcile — the **SMILES should be treated as authoritative**.

## Data provenance
- `data/mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb` — 3.5 Å cryo-EM real-space-refined
  model (PHENIX 2.0), MOR–Gi–scFv16–ZH853; provided by the user. Refinement metadata in REMARK 3.
