# ZH853 · MOR Docking & Structure Analysis

Computational structural-biology study of how the cyclic peptide agonist **ZH853** binds and modulates
the **μ-opioid receptor (MOR / OPRM1)**, starting from a 3.5 Å cryo-EM structure of the
**MOR–ZH853–Gi–scFv16** complex.

## Goals
1. Characterize the **unique structural interactions** of ZH853 vs other MOR agonist complexes.
2. Identify **mutations** that selectively abrogate ZH853 but spare other full agonists.
3. Propose **drug-likeness improvements** (peptide PK/PD, GLP-1-style strategies).
4. Predict relative affinities of ZH853 and analogs via **ABFE / relative FEP / metadynamics**.

## Where things live
| Path | Contents |
|------|----------|
| `data/` | Raw inputs — the cryo-EM PDB, comparator structures, ligand SMILES |
| `src/` | Workflow + code, numbered `src/##.##.##_name` (sequential + hierarchical) |
| `intermediate/` | Cached/temporary results (git-ignored) |
| `product/` | Outputs — figures, tables, reports (`_YYYYMMDD.ext`), manuscript |
| `docs/` | [`PLAN.md`](docs/PLAN.md) (roadmap), [`references.md`](docs/references.md) (comparators + methods) |
| `OBJECTIVES.md` | Original project brief |
| `SPECIFICATION.md` | Decisions log + open questions |

## Workflow

Two machines, and the split is deliberate. **Local** carries `openmm` + `pdbfixer` + RDKit for
structure preparation; the **cluster** carries AmberTools (CPU) and the driver-pinned GPU OpenMM in
two separate envs, and deliberately does *not* carry `pdbfixer` — so the receptor rebuild and
membrane orientation are local steps. `intermediate/` is git-ignored, so anything produced on one
machine must be **copied** to the other; `git pull` will not bring it.

### Configure `cluster.env` first

Before any step that runs on the cluster, copy the template at the **repository root** and fill in
the two required fields:

```bash
cp cluster.env.example cluster.env
$EDITOR cluster.env          # ZH_ACCOUNT and ZH_PARTITION have no default
```

Every batch stage refuses to run without it rather than falling back to defaults, which would
activate a different conda env or run for a different length and still exit 0. It is gitignored:
per-cluster and per-user.

It holds **only site facts** — account, partitions, resources, wall-times, conda env names. What
changes the science (run lengths, replica count) lives per build in the `sampling.env` that
`01_build_system.sh` writes into each build directory, so re-running one system with different
sampling cannot silently inherit another build's choices, and the settings a trajectory was
produced under sit beside the trajectory.

Run `make help` for the grouped target list.

| # | Step | Command | Where |
|---|------|---------|-------|
| 0 | Create the environments | `make env-local` · `make env-cluster` | local · cluster |
| 1 | Fetch comparator structures | `make fetch` | local |
| 2 | Static analysis (Objectives 1–3) | `make analysis` | local |
| 3 | Receptor, ligand and analog preparation | `make prep` | local |
| 4 | Copy prep outputs to the cluster | `scp` — see below | → cluster |
| 5 | Ligand force-field parameters | `./submit.sh params` (4-task array) | cluster (CPU) |
| 6 | Build the membrane systems | `./submit.sh build` (10-task array) | cluster (CPU) |
| 7 | Equilibrate → pre-produce → produce | `./submit.sh all` (or `check` / `eq` / `preprod` / `prod`) | cluster (GPU) |
| 8 | Trajectory QC | `04_analyze.py` | cluster |

**Steps 5–8 are the SLURM bundle. Follow
[`src/02.10.00_slurm_bundle/README.md`](src/02.10.00_slurm_bundle/README.md)** — it is the
authority on cluster settings (`cluster.env`), the CUDA/driver pinning, what each build check
guards against, how equilibration is judged, and the wall-clock estimates. Everything below is only
the hand-off into it.

### Step 3 — preparation, and the one part that also runs on the cluster

`make prep` runs the whole Phase-2 chain: component split, protonation states, receptor rebuild,
OPM orientation, ZH853 preparation, and the analog poses. All of it is local, because the receptor
rebuild needs `pdbfixer`.

The exception is the last step, `make prep-analogs-pose`, whose only inputs are `complex_oriented.pdb`
and RDKit — both present in `zh853mor-prep`. So the three analog poses can be regenerated **on the
cluster** rather than copied, which is the reproducible option. It gives ZH850/ZH831/ZH809 ZH853's
binding mode by constrained embedding on the common scaffold, and prints the scaffold coverage,
strain relative to the ZH853 template, and closest receptor contact for each — read those before
building. The method and its caveats are in the bundle README under *Ligands and the apo system*.

### Step 4 — what has to cross to the cluster

`intermediate/` is git-ignored by design (large, regenerable), so these are copied, not pulled:

```bash
CLUSTER=<user>@<login-node>; REPO=<repo path on the cluster>
scp intermediate/02.05.00_oriented/receptorR_oriented.pdb \
    intermediate/02.05.00_oriented/complex_oriented.pdb   $CLUSTER:$REPO/intermediate/02.05.00_oriented/
scp intermediate/02.04.00_ligand/ZH853_prepared.sdf       $CLUSTER:$REPO/intermediate/02.04.00_ligand/
```

`receptorR_oriented.pdb` is the finalised receptor (ACE/NME caps, named His tautomers, OPM-oriented)
and `complex_oriented.pdb` additionally carries the deposited ZH853 pose. `01_build_system.sh`
refuses to build from a stale copy of the first, so a missed sync fails loudly rather than quietly
producing a differently-protonated system. Copy the analog files too, or regenerate them there with
`make prep-analogs-pose`.

### Steps 5–6 — the panel

Five systems (`apo`, `ZH853`, `ZH850`, `ZH831`, `ZH809`) × two D2.50 protonation states
(`ASP`, `ASH`) = **10 builds**, each landing in
`intermediate/02.10.00_build/<LIGAND>_<D250>_<timestamp>/` as a self-contained run directory:

Both CPU stages are SLURM job arrays, one task per unit of work, so they run concurrently instead
of serially on a login node — parameterization is four `sqm` runs, and each build is ~12–15 min of
PACKMOL-Memgen:

```bash
cd src/02.10.00_slurm_bundle
./submit.sh params      # array 1-4:  one ligand per task  -> intermediate/02.08.00_ligand_params/
./submit.sh build       # array 1-10: one (ligand, D2.50) per task
```

Add `-n` to either to print the `sbatch` command without submitting. They request
`ZH_CPU_PARTITION` and **no GPU** — holding one through PACKMOL-Memgen would waste the allocation.
Parallel builds are safe by construction: each writes to its own timestamped directory, which is
why `01_build_system.sh` insists on a pristine one.

To run either without SLURM, the underlying commands still work directly:
`make prep-ligand-parameterize`, and `LIGAND=… D250=… ./01_build_system.sh`.

Then, from each build directory, `./submit.sh all`.

## Status
Phases 0–1 and static interaction analysis (Phase 3) complete — see [`docs/PLAN.md`](docs/PLAN.md)
and [`docs/RESULTS_interactions.md`](docs/RESULTS_interactions.md). Key findings: ZH853 keeps the
conserved D149 (D3.32) anchor but uniquely engages **E231 (ECL2, salt bridge)** and **H321 (7.36)** —
yielding a ranked mutation panel led by **E231Q/E231A** (Objective 2). Reproduce with `make analysis`.
Analog design (Objective 3) also complete — see [`docs/RESULTS_analog_design.md`](docs/RESULTS_analog_design.md):
all analogs are beyond-Ro5; two design series proposed (N-methylation for permeability, lipidation for
half-life).

**Phase 2/4 (MD) in progress.** The SLURM bundle builds and runs all ten systems; cluster settings
that were OQ-3 now live in one `cluster.env`. The **apo/ASP** arm is furthest along — built,
equilibrated, and through an unrestrained pre-production leg (measured on the H200 nodes: 94.4 ns/day
at 2 fs, 587 ns/day at 4 fs, so ~25 h end to end per system). The ligand arms are not yet built.
Open items: RESP charges still use the AM1-BCC route pending the QM-engine question (the remaining
`TODO(OQ-3)` in `src/02.08.00_ligand_parameterize.sh`); MD occupancy validation; Phase 6 free energy.

## Environment
Four conda environments, created from the specs in the repo root:

| Env | Spec | Where | Used by |
|-----|------|-------|---------|
| `zh853mor-local` | `environment_zh853mor-local.yml` | local | every `make` target here (`make env-local`) |
| `zh853mor-prep` | `environment_zh853mor-prep.yml` | cluster | steps 5–7: AmberTools, PACKMOL-Memgen, RDKit |
| `zh853mor-sim` | `environment_zh853mor-sim.yml` | cluster | steps 8–9: GPU OpenMM (CUDA pin tracks the driver) |
| `zh853mor-plumed` | `environment_zh853mor-plumed.yml` | cluster | metadynamics; optional |

Specs are named `environment_<env name>.yml`. `make env-local` creates the first;
`make env-cluster` creates all three cluster envs. The prep and sim envs are separate
because their `openmm` pins are mutually incompatible, and `zh853mor-prep` omits `pdbfixer`/`openmm`
on purpose — see the bundle README for the CUDA/driver pinning, which is the fiddliest part.
