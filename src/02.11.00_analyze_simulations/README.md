# 02.11.00 — Production simulation analysis

Turns the production trajectories written by
[`src/02.10.00_slurm_bundle/`](../02.10.00_slurm_bundle/README.md) into (a) a defensible
statement about **how reliable the sampling is** and (b) the **scientific comparisons the
objectives ask for**. Phase 4 of [`docs/PLAN.md`](../../docs/PLAN.md).

Where `04_analyze.py` in the simulation bundle is a per-replica smoke test ("did this run produce
something sane?"), this stage is the analysis: it pools replicas, propagates uncertainty, and
judges convergence rather than assuming it.

## Why it is split in two

| Step | Runs | Env | Reads | Writes |
|---|---|---|---|---|
| `01_reduce_trajectory.py` | **cluster** (CPU, job array) | `zh853mor-prep` | `intermediate/02.10.00_build/<system>/prod_r*.dcd` | `intermediate/02.11.00_analyze_simulations/<system>/<replica>.{npz,json}` |
| `02_aggregate.py` | **local** | `zh853mor-local` | those `.npz`/`.json` | `product/02.11.00_*.{csv,md}` |
| `03_figures.py` | **local** | `zh853mor-local` | those `.npz`/`.json` | `product/02.11.00_*.png` |

A 500 ns replica sampled every 100 ps is ~5 GB, and the panel is 3 replicas x up to 10 systems —
on the order of 100 GB, which stays on the cluster. The reduction writes a few hundred kB per
replica, so **only the reduced files need to be copied back**, and every table and figure can then
be regenerated on a laptop in seconds without touching a trajectory.

## Running it

On the cluster, after production has finished:

```bash
cd src/02.11.00_analyze_simulations
./submit_reduce.sh -n          # dry run: shows the inventory and the sbatch command, submits nothing
./submit_reduce.sh             # one array task per replica; resources from the root cluster.env
```

Both print what they found before doing anything, with **how long each replica actually is** —
frames and ns from the DCD header, against the `ZH_PROD_NS` its build was configured for:

```
  #   system       replica     frames        ns    target  status
  1   apo_ASH      prod_r1       5000     500.0       500  complete
  2   apo_ASH      prod_r2       1431     143.1       500  partial-29%
  3   apo_ASH      prod_r3       4880     488.0       500  writing
note: ZH831_ASH has 2 of 3 configured replicas (missing prod_r3)
```

This is the check a file listing cannot make. A replica killed at its wall-time, one still being
written, and a finished one are the same `prod_r*.dcd` to a `glob`; reducing the first two
silently analyses a partial run, and a build that produced two replicas instead of three rests
its replicate spread on two points without saying so. `writing` means the DCD was touched in the
last 15 minutes, `partial-NN%` that it is short of its target, and `unreadable` that the file
could not be opened at all. Everything is read from the DCD header and the state log, so listing
the whole panel costs milliseconds however large the trajectories are.

The submit script warns but does not refuse: reducing a partial run is often what you want. To
take only the finished ones, call `01_reduce_trajectory.py --build <dir>` (repeatable) or
`--replica <name>` directly instead of submitting the array.

or, for a single replica without SLURM (a few minutes):

```bash
python 01_reduce_trajectory.py --build ../../intermediate/02.10.00_build/apo_ASH_20260901_133002
```

Then bring `intermediate/02.11.00_analyze_simulations/` to the machine with the local env and:

```bash
make sim-aggregate sim-figures      # or run 02_aggregate.py / 03_figures.py directly
```

Job logs go to `intermediate/02.11.00_analyze_simulations/logs/reduce_<jobid>_<task>.out` — the
stage is submitted from `src/`, where generated files do not belong (D-16/D-22). `submit_reduce.sh`
creates that directory; submitting `submit_reduce.sbatch` by hand needs the `mkdir -p` first.

`01_reduce_trajectory.py` skips a replica whose `.json` already exists — rerun with `--force` after
changing what is measured. `--stride N` subsamples frames for a quick look.

## What is measured, and why

**Reliability.** Thermodynamics (T, density and its drift) from the OpenMM state log; membrane
area per lipid, thickness and OPM registration; Ca RMSD against the staged OPM-oriented receptor,
whole and TM-only; the C142–C219 disulfide; receptor-aligned ligand RMSD against the **deposited
pose**; and — the part an RMSD plateau cannot answer — statistical inefficiency, effective sample
size, residual drift with a standard error, blocking curves, R-hat across replicas, and the cosine
content of the leading principal components. Thresholds and their justification are constants at
the top of `02_aggregate.py`; the physical ranges match `check_equilibration.py`, the structural
tolerances are deliberately looser because 500 ns of free dynamics is *supposed* to move.

**Science.** Contact occupancy for every receptor residue at the same 4.5 Å heavy-atom criterion
as the static cryo-EM fingerprint in `03.01.00`, so MD and structure are directly comparable
(Objective 1); direct **and water-mediated** polar contacts for the anchors, because the canonical
H6.52 interaction is water-bridged and a direct-contact count alone scores it as absent;
ZH853-minus-analog occupancy differences, which is the part of Objective 2 a single structure
cannot supply; ligand pose retention; the R3.50–T6.34 and R3.50–Y7.53 activation rulers, which ask
whether the ligand holds the active state in a system built without a transducer (D-10 system B);
and Na+ occupancy of the D2.50 site, the direct test of the ASP/ASH pair built under D-11.

## Reading the numbers

* **Occupancy** is the fraction of *post-equilibration* frames in contact, reported as
  mean ± SEM **over replicas** — three replicas give a coarse error bar, and it is the replicate
  spread, not the frame count, that sets it.
* **`t0`** is chosen per replica by maximising the effective sample size (Chodera 2016) on the TM
  Ca RMSD, and the same window is then used for every observable of that replica.
* **`n_eff`**, not the number of frames, is the sample size. A 5000-frame replica with g = 100
  carries 50 independent samples.
* **R-hat > 1.2** or **PC1 cosine content > 0.5** means the result is not converged, whatever the
  RMSD trace looks like.

## Known limits

* Area per lipid is reported gross and hull-corrected. The convex hull overestimates the
  receptor's cross-section, so the corrected value is a **lower** bound; neither is a clean
  comparison against a pure-POPC literature area.
* The activation rulers are Ca–Ca distances. At 3.5 Å the sidechain rotamers are the least
  reliable coordinates in the starting model, so a sidechain-tip definition would measure the
  model's guesses rather than the receptor's state.
* Occupancy is geometric and hydrogen-free, matching `zh853mor.interactions` (D-9). PLIP/ProLIF
  cross-validation is still deferred.
* With three replicas, R-hat is a flag rather than a measurement.

## Files

| File | Role |
|---|---|
| `01_reduce_trajectory.py` | one pass per replica over the DCD -> compact `.npz` + `.json` |
| `submit_reduce.sh` / `submit_reduce.sbatch` | SLURM array, one replica per task (CPU) |
| `02_aggregate.py` | QC verdicts, occupancy and activation tables, the markdown report |
| `03_figures.py` | QC dashboard, convergence, contact occupancy, pocket dynamics |

The analysis primitives live in the package, where lint, type checking and the unit tests reach
them: [`zh853mor/md.py`](../zh853mor/md.py) (selections, the construct-numbering map, per-residue
contacts, membrane geometry) and [`zh853mor/convergence.py`](../zh853mor/convergence.py)
(statistical inefficiency, equilibration detection, blocking, R-hat, cosine content), tested in
`tests/test_md.py` and `tests/test_convergence.py`.
