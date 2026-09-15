# 02.13.00 — Conformational landscapes

Reduces every replica to **two collective variables**, and draws each system's density on that
plane as a free-energy heatmap — so the panel of simulations becomes a set of pictures that can be
laid side by side and subtracted, instead of ten tables of means.

Where [`02.11.00`](../02.11.00_analyze_simulations/README.md) answers *is each observable
converged and what is its mean*, and [`02.12.00`](../02.12.00_render_trajectory_movies/README.md)
answers *what is actually happening in there*, this stage answers the one in between: **which
regions of conformational space does each system occupy, and where do they differ?** That is the
Objective-1/2 question in its most direct form.

Both steps **run locally**. Nothing here touches a trajectory: the per-frame feature matrices come
from `02.11.00`'s reduction, which is the only stage that reads the DCDs.

| Step | Reads | Writes |
|---|---|---|
| `01_fit_landscape.py` | `intermediate/02.11.00_analyze_simulations/*/*.npz` | `intermediate/02.13.00_map_conformational_landscapes/<space>/{projections.npz,model.json}` + `product/02.13.00_<space>_model_*.md` |
| `02_figures.py` | those | `product/02.13.00_<space>_{landscapes,differences,diagnostics}_*.png` |

```bash
make landscape                                  # receptor space, both steps
python 01_fit_landscape.py --space ligand --scan && python 02_figures.py --space ligand
```

## The three coordinate spaces

| `--space` | Features | Covers | What it is for |
|---|---|---|---|
| `receptor` | pairwise Cα–Cα distances among the 36 key + functional residues (`key_pairs`) | **every** system, apo included | receptor conformational state — pocket shape, DRY/NPxxY, TM6 |
| `ligand` | every receptor residue's minimum heavy-atom distance to the ligand (`min_dist`) | holo only | the binding-mode landscape: Objectives 1 and 2 |
| `pair` | two named features, no fitting | whatever defines them | an interpretable plane you choose |

`receptor` and `ligand` use **tICA** — time-lagged independent component analysis, implemented in
[`zh853mor.landscape`](../zh853mor/landscape.py). PCA finds the directions of largest *variance*;
tICA finds the directions of slowest *decorrelation*. For a receptor those differ: the
largest-amplitude motion in 500 ns is usually a floppy loop, while the coordinate worth a figure —
the ligand leaving a subsite, TM6 opening — is slower and smaller.

`pair` skips the fitting and takes two features you name, which is the honest option when the
question is already specific. `--x` and `--y` accept any series the reduction stored, plus three
prefixed forms:

```bash
--x activation:tm3_tm6_ca --y activation:npxxy_tm3_ca    # the classic GPCR activation plane
--x lig:149 --y lig:299                                  # D3.32 anchor vs the H6.52 contact
--x ca:167-281 --y na_d250_dist                          # TM3–TM6 opening vs the D-11 Na+ site
```

## Three things that decide whether the picture means anything

**The basis is fitted once and every system is projected into it.** Fitting per system gives each
its own axes, and two landscapes drawn on different linear combinations of distances cannot be
compared however alike they look. `--fit-on <system>` instead fits on one system and projects the
rest — the stricter reading, where the axes are one system's kinetic modes and the others are
being viewed in coordinates chosen without them.

**A pooled fit's slowest "process" is partly not a process.** Systems do not interconvert: no
trajectory turns apo into holo. A direction that separates two systems therefore has an
autocorrelation of ~1 at any lag, and tICA will return it first, with an infinite implied
timescale. That is a useful *discriminative* axis — it is what makes the comparison legible — but
it is not a relaxation time. So the report gives, per component, the **fraction of its variance
that lies between systems**; above ~0.5 the component is mostly telling you which simulation you
are looking at. Timescales are additionally re-estimated *within* each system, where the
coordinate can actually relax.

**tICA inherits the D-20 trap.** On a trajectory that has not converged, the slowest apparent
process *is* the drift of an unequilibrated coordinate, and tIC1 reproduces the half-cosine of
free diffusion exactly as PC1 does. `convergence.cosine_content` applies to a tIC projection
unchanged, and is reported per replica; above 0.5 that replica's density is a picture of its own
drift. Everything is fitted on **post-equilibration frames only**, using each replica's own `t0`
from `02.11.00`.

## Reading the figures

* **`landscapes`** — one −kT ln P heatmap per system, all on the same bins and the same colour
  scale. Dark is low free energy, which inverts the usual "lightest means near zero" on purpose:
  here near-zero is the most-populated and best-determined region, and the pale rim is where the
  estimate rests on a few frames. Bins below `--min-count` frames are left unshaded rather than
  drawn as a barrier — −kT ln(1/N) for one stray frame is a confident-looking wall built from one
  sample.
* **`differences`** — each system minus a reference, on a diverging scale through a neutral grey.
  Shown only where **both** are sampled; hatching marks where this system goes and the reference
  does not, which is a difference in sampling, not a measured one.
* **`diagnostics`** — the implied timescale against lag (it must stop depending on the lag before
  the lag is defensible, and it must stay above the *t* = lag line to be resolved at all); the
  region holding 50 % of each replica's frames, as a replicate-agreement check; cosine content
  against the 0.5 threshold; and what the axes are made of, as the features most correlated with
  each component.

The binning adapts to the sample size unless you fix `--bins`: a grid fine enough for 15 000
frames leaves a short replica as a scatter of single-sample specks.

## Known limits

* **These are sampling densities in energy units, not converged free energies** (D-5). The depth
  of a basin visited twice is a statement about those two visits. Nothing here is reweighted,
  because production is unbiased MD; a barrier between two basins is a barrier the trajectory
  failed to cross, which is not the same as a barrier that is high.
* Features are Cα-based (D-19's reasoning): at 3.5 Å the sidechain rotamers are the least reliable
  coordinates in the starting model, and featurising them would let tICA find slow modes in the
  model's guesses. The cost is that a pure rotamer rearrangement is invisible in `receptor` space.
* The `ligand` space is undefined for apo, and the `receptor` space is only comparable between
  builds whose key-residue set is identical — a replica whose feature axis differs is dropped with
  a warning rather than silently projected onto someone else's coordinates.
* tICA is linear. A collective variable that is a nonlinear function of the input distances is not
  in the reachable set, whatever the lag.
* With three replicas the replicate-agreement panel is a flag, not a measurement — as R-hat is in
  `02.11.00` (D-20).
