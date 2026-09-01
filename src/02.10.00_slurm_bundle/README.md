# ZH853–MOR MD SLURM bundle

Self-contained membrane-MD workflow for the MOR–ZH853 (+Gi) complex, to run on an academic
SLURM/GPU cluster. Prepared locally in Phase 2.

**All cluster-specific values live in one file: `cluster.env`** (account, partition, GPU type,
wall-times, replica count, modules, conda env names). Copy the template and fill in the two
required fields before anything is submitted:

```bash
cp cluster.env.example cluster.env
$EDITOR cluster.env          # ZH_ACCOUNT and ZH_PARTITION are the only ones with no default
```

`cluster.env` is gitignored — it is per-cluster and per-user, so it must not travel with the repo;
`cluster.env.example` is the tracked template. `submit.sh` reads it and passes the values to
`sbatch` as **command-line flags**, which take precedence over `#SBATCH` directives in a job
script. That is why the `.sbatch` files carry no account or partition of their own: there is one
source of truth per cluster rather than a `CHANGEME` in each job script, duplicated again into
every timestamped build directory.

## Inputs (copy from the repo)
- `intermediate/02.05.00_oriented/receptorR_oriented.pdb` — rebuilt receptor (Phase 2), superposed
  onto the OPM reference so the membrane normal is z and the OPM midplane is z = 0. `01_build_system.sh`
  stages this automatically; the un-oriented `02.03.00_receptor/receptorR_fixed_heavy.pdb` must **not**
  be used with `--preoriented`.
- `intermediate/02.04.00_ligand/ZH853_prepared.sdf` — protonated ligand (+1), bond orders assigned.
- Protonation/His/D2.50 decisions: `product/02.02.00_protonation_*.md`.

## Systems to build (SPECIFICATION D-10/D-11)
- **System A — active-state complex:** MOR + ZH853 + Gi(αβγ) in POPC:chol (9:1), scFv16 removed.
  For Objective 1–2 dynamics/occupancy.
- **System B — binding/FEP:** MOR + ZH853 in POPC:chol, intracellular half Cα-restrained (or Gα
  α5-helix retained) to hold the active state. For Objective 4 throughput.
- Each built **twice**: D2.50 (Asp116) charged vs protonated (ASH).

> **OPEN — the build is currently apo, and neither System A nor System B.** `01_build_system.sh`
> stages `receptorR_oriented.pdb`, which is the **receptor alone**: no ZH853 and no Gi. `tleap.in`
> does `LIG = loadmol2 ZH853.mol2` but never combines `LIG` into `sys`, so that unit is dropped and
> `saveamberparm` writes a ligand-free system — silently, since loading parameters for an absent
> residue is not an error. `make_tleap.py` now warns when no ligand residue is present in the packed
> box. Fixing this is a build-design decision, not a one-line edit: the ligand (and Gi, for System A)
> must be packed *with* the protein, so `02.05.00` needs to emit a complex whose ligand residue name
> matches `ZH853.mol2`, and PACKMOL-Memgen's `reduce` preprocessing has to be kept from mangling the
> HETATM records. `intermediate/02.05.00_oriented/complex_oriented.pdb` (receptor + ligand, same
> OPM transform) is the natural starting point. Decide this before spending GPU time.

## Conda environments (three, task-specific)
They are split because their `openmm` pins are mutually incompatible (`openmm-plumed` lags the
newest `openmm` that `openmmforcefields` requires):
| Env | Spec | Used by | Notes |
|-----|------|---------|-------|
| `zh853mor-prep` | `environment-prep.yml` | steps 1–2 | CPU; AmberTools/PACKMOL-Memgen |
| `zh853mor-sim` | `environment-cluster.yml` | steps 3–5 + FEP | GPU; openmm + openmmforcefields |
| `zh853mor-plumed` | `environment-plumed.yml` | metadynamics (3.9) | GPU; older openmm + openmm-plumed; optional |

## Workflow (run in order)

Steps 0–2 are run **from this directory** (`src/02.10.00_slurm_bundle/`) on the login/CPU node.
Steps 3–5 are run **from the build directory** that step 2 creates — `01_build_system.sh` copies
the run scripts and `.sbatch` files in there, so each build is a self-contained run directory that
you `cd` into and submit from. Nothing in steps 3–5 is invoked out of `src/`.

| Step | Script | Env | Where | Run from |
|------|--------|-----|-------|----------|
| 0 | `00_install.sh` | — | login node (builds the envs) | `src/02.10.00_slurm_bundle/` |
| 0.5 | `./submit.sh check` → `check_gpu_env.sh` | `zh853mor-sim` | **GPU (pre-flight)** | either |
| 1 | `ligand_resp/run_resp.sh` | `zh853mor-prep` | CPU | `src/02.10.00_slurm_bundle/` |
| 2 | `01_build_system.sh` | `zh853mor-prep` | CPU | `src/02.10.00_slurm_bundle/` |
| 2a | `check_placement.py` | `zh853mor-prep` | CPU (called by step 2) | — |
| 2b | `check_piercing.py` | `zh853mor-prep` | CPU (called by step 2) | — |
| 2c | `fix_caps.py` | `zh853mor-prep` | CPU (called by step 2) | — |
| 2d | `make_tleap.py` | `zh853mor-prep` | CPU (called by step 2) | — |
| 3 | `./submit.sh eq` → `submit_equilibrate.sbatch` → `02_equilibrate.py` | `zh853mor-sim` | GPU | the build directory |
| 4 | `./submit.sh prod` → `submit_production.sbatch` → `03_production.py` | `zh853mor-sim` | GPU | the build directory |
| 5 | `04_analyze.py` | `zh853mor-sim` | CPU | the build directory |

### Building the system (step 2)

```bash
./01_build_system.sh              # D2.50 charged (default)
D250=ASH ./01_build_system.sh     # D2.50 protonated variant (D-11/D-17)
```

Both systems come from the **same** prepared and OPM-oriented receptor — the variant is a residue
rename applied to the staged copy, so the comparison cannot be confounded by a different starting
geometry. Output goes to `intermediate/02.10.00_build/<D250>_<timestamp>/`, never inside `src/`
(D-16).

`01_build_system.sh` creates a **fresh timestamped build directory** and works there. This is
not cosmetic: PACKMOL-Memgen writes to the CWD *and silently reuses anything it finds there* — the
component PDBs (`POPC.pdb`, `CHL1.pdb`, `WAT.pdb`, `Na+.pdb`, `Cl-.pdb`; watch for its "Using
WAT.pdb in the folder" message) and the preprocessed protein files (`receptor_Trim_H.pdb`,
`*.grid.pdb`, `receptorin_EMBED*.pdb`). Once a stale component PDB is present whose atom order no
longer matches the `atoms 1 20` / `atoms 88 131` head/tail constraints in the generated
`packmol.inp`, PACKMOL fails with

```
ERROR: Packmol was unable to put the molecules in the desired regions
       even without considering distance tolerances.
Maximum violation of the restraints:  26.12
```

and **no combination of `--lipids` / `--ratio` / `--dist` / `--preoriented` will fix it**, because
those flags do not affect the cached files. The 26 Å residual is one lipid length and the report says
100 % of that type violates — the signature of contradictory constraints, not of an overpacked box.
Overpacking instead shows up as a GENCAN convergence failure. Set `BUILD_DIR=` to override the
location; the script refuses a non-empty directory.

Two post-build steps then run automatically:

- **`check_placement.py`** — guards the OPM placement from 02.05.00. PACKMOL-Memgen re-centres the
  solute on its own *z* bounding box when it orients the protein itself, and our receptor's bbox
  centre is ~5 Å below the OPM midplane (the intracellular face — H8, ICL3, C-term — protrudes
  further than the extracellular face), so a re-centred build would sit ~5 Å too high in the membrane
  with nothing downstream complaining. Measured 2026-07-26 with `--preoriented` on packmol-memgen
  2025.1.29: translation `(+0.02, −0.18, +0.00) Å`, misregistration **+0.07 Å** — honoured, so this
  is a regression guard rather than a workaround. It measures the receptor's rigid-body shift against
  the lipid phosphate planes, **fails the build** past 1.5 Å, and reports the Trp/Tyr girdle for
  comparison with the OPM ±15.7 Å slab (built: −16.5/+13.7 Å, P–P thickness 38.6 Å, matching the
  37–40 Å predicted in Methods §3.7). Override with `SKIP_PLACEMENT_CHECK=1` only to inspect a
  known-bad build.
- **`check_piercing.py`** — finds lipid tails threaded through rings. PACKMOL enforces a minimum
  pairwise distance, which does not stop a tail passing *through* a Phe/Tyr/Trp/His/Pro ring or one of
  cholesterol's fused rings: every atom can be >2 Å from every other while the tail is threaded. This
  is what packmol-memgen's own finder tries to catch and could not here ("Lipid piercing finder
  failed"). It matters more than a clash because it is topological — minimisation separates atoms
  along straight lines and no such path unthreads a ring, so the lipid stays trapped for the whole
  trajectory. Rings are found as 5-/6-cycles in each residue's distance-inferred bond graph (so
  cholesterol and the ligand need no hardcoded atom names) and every bond from another residue is
  tested for segment/disc intersection. Reports without failing the build: the fix (re-pack with a
  different seed) is a human decision.
- **`fix_caps.py`** — normalises the ACE/NME caps after packing, before tleap. packmol-memgen's
  preprocessing mangles them on this stack: reduce protonates the capping amides and a stray H
  (`HN2`) survives the H-strip on the NME nitrogen, and the ff19SB libraries name the NME methyl
  carbon `C` (H1/H2/H3) while 02.03.00 writes the older pdbfixer convention `CH3` — `loadpdb`
  matches atoms by name, so the file `CH3` becomes an unknown typeless atom and tleap dies with
  "Atom .R<NME …>.A<CH3> does not have a type" (2026-08-25 build). The script reads the expected
  heavy atoms from `$AMBERHOME/dat/leap/lib/aminoct12.lib` itself (hardcoded ff19SB fallback),
  renames unambiguous mismatches in place — `CH3`→`C` moves nothing — drops stray cap hydrogens
  (leap rebuilds every H anyway), and can restore a genuinely missing heavy atom from the staged
  receptor rigidly translated by the median CA offset. Idempotent; ACE is audited too.
- **`make_tleap.py`** — fills the two `@PLACEHOLDERS@` in `tleap.in` that cannot be known until the
  system is packed, writing `tleap_run.in` plus a `bilayer_system_ff.pdb`:
  - the **disulfide**, because `loadpdb` renumbers every residue sequentially from 1 across the whole
    system, so the OPRM1 numbering in `bond sys.142.SG sys.219.SG` no longer exists. C142–C219 is
    found geometrically (SG–SG 2.02 Å) and emitted in tleap's numbering — 74/151 for the current
    69–349 construct, but that shifts with ACE/NME caps, so it is never hardcoded. The two cysteines
    are renamed `CYS`→`CYX` so ff19SB does not build an HG onto a bonded sulfur.
  - the **periodic box**. Candidate cells are collected from CRYST1, `packmol-memgen.json` and the
    `inside box` bounds in `packmol.inp`, and each is checked against the packed coordinates — the
    first that actually contains them wins. This check exists because the 2026-07-26 build produced a
    system spanning 87.0 × 87.5 × 113.2 Å inside a nominal 81.40 × 81.40 × 111.72 Å cell: a cell
    smaller than its own contents wraps atoms onto their periodic images, and no amount of
    minimisation repairs that. The original `setBox sys vdw` was worse still — it derives a box from
    van der Waals extents and does not reproduce the packing cell at all.
    In practice none of the recorded sources ever fires with packmol-memgen 2025.1.29: it writes no
    CRYST1 record at all, `packmol-memgen.json` is only produced for `--dry_run`, and its own box
    summary goes to the DEBUG log. Worse, PACKMOL ends both all-together runs at maxiter with rim
    lipids and water pushed partway past the `inside box` walls (which are per-atom constraints —
    memgen writes `inside box` at `structure` level — so overflow atoms are real restraint
    violations), making the parsed union systematically ~6 Å narrower (x+y) than what was packed.
    The offset is constant in `--dist`: 82.62 → 87.85 and 92.62 → 97.87 Å for `--dist 15` vs `20`
    (+5.24 both), i.e. systematic maxiter squeeze, not random failure — repacking does not remove it.
    When every candidate fails, `make_tleap.py` classifies the overflow by how deep atoms reach past
    the absolute `inside box` walls: ≤ 6 Å is **rim squeeze** (edge molecules' tails/waters poking
    out, most of each molecule still inside — absorbed with the packed extent + 2 × 1.25 Å and a
    loud WARNING; the content already fills that cell so no low-density gap for the membrane
    barostat to collapse is created), anything deeper means displaced molecules (holes inside,
    clashes across the periodic seam) and aborts. `python make_tleap.py --diagnose
    bilayer_system.pdb` prints every candidate and this verdict without building anything.
    Convergence context: `grep -nE "STOP|SUCCESS|Maximum violation" packmol.log`.

### Running the simulation (steps 3–5)

> **Check what you built before spending GPU time.** `tleap.in` never combines `LIG` into `sys`
> (see the OPEN note at the top), so a build made from the receptor-only `receptorR_oriented.pdb`
> is **apo** — no ZH853, no Gi — and step 5's `resname LIG` selection will be empty. `make_tleap.py`
> warns about this during step 2. Confirm before submitting:
> ```bash
> python -c "import parmed; print(sorted({r.name for r in parmed.load_file('system.prmtop').residues}))"
> ```
> If no ligand residue name (`ZH8`/`LIG`/…) appears, the trajectory cannot answer Objectives 1–2 and
> only the receptor dynamics are usable.

**Everything below runs from the build directory**, e.g.
`intermediate/02.10.00_build/ASP_20260825_154710/`. Step 2 stages `02_equilibrate.py`,
`03_production.py`, `04_analyze.py` and both `.sbatch` files there next to `system.prmtop`, so the
paths inside the sbatch scripts are all bare filenames and the job's working directory is the run
directory. Do not submit from `src/` and do not point `--prmtop` across directories: the analysis
and checkpoint outputs are written relative to the CWD and belong with the build they came from.

```bash
cd intermediate/02.10.00_build/ASP_20260825_154710

./submit.sh check     # step 0.5  short GPU pre-flight (openmm on CUDA) -> gpucheck_<jid>.out
./submit.sh all       # steps 3+4 equilibration, then production chained with afterok
```

or one stage at a time, to inspect `system_eq.pdb` in between:

```bash
./submit.sh eq        # step 3 -> system_eq.xml, system_eq.pdb
./submit.sh prod      # step 4 -> prod_r{1,2,3}.dcd  (only after eq has finished)
```

Add `-n` to any of these to print the `sbatch` command without submitting. `submit.sh all` chains
the two with `--dependency=afterok:` — not `afterany` — because production loads the state
equilibration writes, so a failed or cancelled step 3 must not launch replicas that would die on a
missing `.xml`. The wrapper also checks that `system.prmtop`/`system.rst7` are actually in the
current directory and that `ZH_ACCOUNT`/`ZH_PARTITION` are set, so a wrong directory or an unfilled
`cluster.env` fails immediately with a pointer, rather than as a queued job that dies minutes later
or an `Invalid account` rejection with no indication of which file to edit.

Then, once production finishes (CPU is fine, one call per replica):

```bash
python 04_analyze.py --top system.prmtop --traj prod_r1.dcd --lig LIG --out qc_r1
```

**The file naming is load-bearing.** Both job scripts key off `ZH_SYS` (default `system`, matching
what tleap writes) and equilibration writes `${ZH_SYS}_eq.xml`, which is exactly what
`submit_production.sbatch` reads back. Passing `--out eq` when calling `02_equilibrate.py` by hand
produces `eq.xml` and the production job will then fail to find `system_eq.xml`. Prefer the wrapper;
if you do run a stage manually, keep the `system_eq` basename:

```bash
python 02_equilibrate.py --prmtop system.prmtop --inpcrd system.rst7 --out system_eq
```

Wall-time: equilibration is 1.125 M steps at 2 fs (2.25 ns over six restrained stages); production
is `ZH_PROD_NS` (500) ns/replica at 4 fs with HMR. Size `ZH_EQ_TIME`/`ZH_PROD_TIME` from the ns/day
your GPU reports in `prod_r*.log` — the 12 h / 48 h in the template are only a starting guess, and
production will need restarting from `prod_r<N>.chk` if it does not fit.

### cluster.env reference

| Variable | Default | Notes |
|----------|---------|-------|
| `ZH_ACCOUNT` | *(required)* | allocation to charge — `sacctmgr show assoc user=$USER` |
| `ZH_PARTITION` | *(required)* | GPU queue — `sinfo -s` |
| `ZH_GRES` | `gpu:1` | pin the model (`gpu:a100:1`) if the partition is mixed |
| `ZH_CPUS` / `ZH_MEM` | `8` / `32G` | |
| `ZH_EQ_TIME` / `ZH_PROD_TIME` / `ZH_CHECK_TIME` | `12:00:00` / `48:00:00` / `00:05:00` | |
| `ZH_QOS` / `ZH_CONSTRAINT` / `ZH_EXTRA_SBATCH` | empty | passed through to `sbatch` when set |
| `ZH_REPLICAS` / `ZH_PROD_NS` | `3` / `500` | array size and ns per replica |
| `ZH_MODULES` | empty | normally stays empty — see CUDA / modules below |
| `ZH_CONDA_SH` | `$(conda info --base)/…` | set only if conda is not on `PATH` in the job |
| `ZH_SIM_ENV` / `ZH_PREP_ENV` | `zh853mor-sim` / `zh853mor-prep` | |
| `ZH_SYS` | `system` | basename of the tleap products |

`submit.sh` looks for `cluster.env` in `$ZH_CLUSTER_ENV`, then the current (build) directory, then
beside itself in `src/02.10.00_slurm_bundle/` — so filling in the bundle copy once covers every
build, and a build directory may override it with a local copy for a one-off (a longer wall-time,
a different partition) without touching the shared file. The resolved path is exported into the
job, so the job body reads the same settings the submission used.

### Applying this to an existing build

Builds made before the run scripts were staged automatically will not have them. Copy them in:

```bash
cd intermediate/02.10.00_build/ASP_20260825_154710
cp ../../../src/02.10.00_slurm_bundle/{02_equilibrate.py,03_production.py,04_analyze.py} .
cp ../../../src/02.10.00_slurm_bundle/{submit_equilibrate,submit_production}.sbatch .
cp ../../../src/02.10.00_slurm_bundle/{submit.sh,check_gpu_env.sh} .
```

`cluster.env` need not be copied — `submit.sh` falls back to the bundle's copy.

## CUDA / modules
OpenMM (conda-forge) bundles its own CUDA runtime, so it needs only the node's NVIDIA **driver** —
`module load cuda` does nothing for it and cuDNN is never needed. The env's `cuda-version` pin
must be **≤ the driver's ceiling** (the "CUDA Version" in the `nvidia-smi` header): conda cannot
see the driver and otherwise grabs the newest toolkit build, whose PTX the driver then refuses to
JIT — `CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`; this bit the H200 nodes (driver ceiling 13.2)
with an env solved at 13.3 on 2026-08-25, fixed by `cuda-version=13.2` in
`environment-cluster.yml`. Bump that pin when the cluster driver moves.
Subtler and also hit on 2026-08-25: conda-forge publishes **two openmm builds per python version,
built against different CUDA toolchains**, and the toolchain baked into the build decides which
PTX ISA it emits — not the `cuda-nvrtc` installed beside it. The driver accepted PTX ≤ ISA 9.2 and
rejected 9.3+ (probed via raw `cuModuleLoadData`), so the working build is pinned exactly in
`environment-cluster.yml`; if CUDA breaks again after a solver upgrade, re-probe the boundary and
re-pick the hash. Before the first real
run, submit the pre-flight **`./submit.sh check`** (step 0.5): it prints the
modules/`nvidia-smi`/env, runs `openmm.testInstallation`, and does a real 200-step CUDA run,
ending in a clear **PASS/FAIL** with a fix checklist. If the JIT compiler itself is not found,
`export OPENMM_CUDA_COMPILER=$(which nvcc)`.

## Force field (SPECIFICATION D-12)
ff19SB (protein) + Lipid21 (membrane) + OPC water + GAFF2/RESP ligand; 0.15 M NaCl. HMR → 4 fs.
(CHARMM36m + CHARMM-GUI is the documented alternative — do not mismatch the water model to the FF.)

## Reproducibility
Pin exact versions in `environment-cluster.yml` once the cluster CUDA/driver stack is known.
Record `openmm.version`, GPU, and CUDA in each run log. Trajectories stay on the cluster; copy back
QC summaries + representative frames to `product/`.
