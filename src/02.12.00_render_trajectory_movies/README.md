# 02.12.00 — Trajectory movies

Turns a production replica into something you can **watch**. Where
[`02.11.00`](../02.11.00_analyze_simulations/README.md) answers *is this trajectory trustworthy,
and what does it say*, this stage answers the question that comes before both: **what is actually
happening in there?** It is deliberately exploratory — its job is to surface the things nobody
thought to measure.

Things a scalar time series will not tell you, and a movie will in ten seconds:

* the ligand leaves the pocket, or rotates in place, or its macrocycle flips a peptide bond;
* a lipid tail slides into the orthosteric site through the TM5/TM6 gap;
* the receptor tilts in the bilayer, or slides along the membrane normal;
* a helix end frays, or ECL2 peels off the pocket lid;
* the bilayer buckles, or one leaflet thins under the receptor;
* the whole thing is fine and the RMSD step at 230 ns was a rigid-body wobble of the Gi lobe.

Every one of those shows up as "the RMSD went up" in a plot.

## The three steps

| Step | Runs | Env | Reads | Writes |
|---|---|---|---|---|
| `01_export_movie_trajectory.py` | **cluster** (CPU, job array) | `zh853mor-prep` | `intermediate/02.10.00_build/<system>/prod_r*.dcd` | `intermediate/02.12.00_render_trajectory_movies/<system>/<replica>_movie.{pdb,xtc,json}` |
| `02_render_dashboard.py` | local **or cluster** (CPU, job array) | `zh853mor-local` / `zh853mor-prep` | those movie files (+ `02.11.00`'s `.npz` if present) | `product/02.12.00_*_dashboard_*.mp4` |
| `03_render_molstar.js` | local **or cluster** (see below) | Node.js ≥ 18 | those movie files | `product/02.12.00_*_molstar_*.mp4` |

The split is the same one `02.11.00` makes, for the same reason: a 500 ns replica is ~5 GB and the
panel is ~100 GB, so the trajectory stays on the cluster. **Only the exported movie files come
back** — a few tens of MB per replica.

## Running it

On the cluster, after production has finished:

```bash
cd src/02.12.00_render_trajectory_movies
./submit_export.sh -n          # dry run: shows the inventory and the sbatch command
./submit_export.sh             # one array task per replica; resources from the root cluster.env
```

The inventory printed before anything is submitted includes the two columns specific to this
stage — how many frames the movie will have, and **how much simulated time each frame spans**:

```
  #   system       replica     frames        ns    movie   ns/frame  stride  status
  1   apo_ASH      prod_r1       5000     500.0      300       1.67   1/16  complete
  2   apo_ASH      prod_r2       1431     143.1      300       0.48   1/4   partial-29%
```

Read `ns/frame` before submitting. At 1.67 ns per frame nothing faster than that is in the
result: a sidechain rotamer flip, a water exchange, a transient contact — all gone. Raise
`--frames` (`./submit_export.sh -- --frames 1000`) when the question is a fast event.

Then render. Both passes are **single-core and serial within a replica** — matplotlib draws
every frame in turn and pipes it to ffmpeg, and MolStar draws its frames one at a time in software
WebGL — so a panel of fifteen replicas takes fifteen times one replica. There is no threading to
be had inside one; the way to make it fast is to run the replicas concurrently.

**On the cluster** (one array task per replica; fifteen then take about as long as the slowest
one):

```bash
cd src/02.12.00_render_trajectory_movies
./submit_render.sh -n                    # dry run: inventory + the sbatch command
./submit_render.sh                       # both passes
./submit_render.sh --only dashboard      # matplotlib only — needs no node
```

**Locally**, after copying `intermediate/02.12.00_render_trajectory_movies/` back:

```bash
make movies            # both renders; or movie-dashboard / movie-molstar separately
```

Either way **`ffmpeg` must be on PATH**, or `02_render_dashboard.py` degrades to an animated GIF —
several times larger and not seekable — with a warning that is easy to miss until the files exist.
It is declared in both environment files; an env created before that needs
`conda env update -f environment_zh853mor-prep.yml` (or `…-local.yml`).

The MolStar pass additionally needs `node` **and** an installed `node_modules`. `node` ships
inside `zh853mor-prep`; the install has to happen **from a login node**, because compute nodes
usually have no outbound network and puppeteer downloads its own Chromium:

```bash
make env-cluster                                    # creates the envs and does this
conda activate zh853mor-prep && make env-node       # or just this, if the envs exist
```

(The root [README](../../README.md#nodejs-and-why-it-has-to-be-installed-on-a-login-node) says the
same thing next to the rest of the installation.) Even with it installed, that Chromium links
against X/NSS shared libraries that not every cluster image carries. `submit_render.sh` checks all
of this before submitting and each task checks again, skipping the MolStar pass with a message
rather than failing the array — so the dashboard movies are still produced. If it will not run
there, render that pass locally; it is the presentable one, not the diagnostic one.

`01_export_movie_trajectory.py` skips a replica whose `.json` already exists — rerun with
`--force` after changing what is exported.

## What the export actually does, and why

A DCD written by OpenMM is not something a viewer can show. Four things are fixed here.

**Molecules are made whole.** OpenMM wraps coordinates into the periodic box *molecule by
molecule*. Rendered as written, the receptor and its lipids are torn in half at the box faces, and
the ligand can appear on the far side of the cell while still sitting in the pocket. Every frame
is unwrapped by fragment, centred on the protein and re-wrapped around it.

**The selection is fixed once.** A trajectory file has one atom count for every frame, so
"lipids within 8 Å" cannot be re-evaluated per frame. The annular shell is taken as the **union**
over frames sampled across the whole run — so a lipid that only arrives at 400 ns is present from
the start and can be watched arriving. Waters are handled the other way: over 500 ns thousands
pass within 5 Å of the pocket, so they are **ranked by occupancy** and only the most persistent
`--waters` (default 40) are kept, which is what leaves the water-bridged H6.52 contact visible
instead of buried in bulk solvent.

**Frames are superposed** on the membrane-embedded Cα — the same set `02.11.00` measures its TM
RMSD on. Rotational and translational diffusion of the whole system is real, uninformative, and
completely dominates an unsuperposed movie.

**…unless you ask for the raw box.** `--no-align` skips the rotational fit *and* the centring
along z, so whole-box drift, receptor tilt and membrane registration are left in. x and y are
still centred, because the box origin is arbitrary in the membrane plane and the receptor would
otherwise wander out of frame. Export both when a QC verdict is in doubt: the aligned movie shows
internal motion, the raw one shows whether the system as a whole is going somewhere.

The exported PDB also carries three things the raw topology does not: **chain IDs** (tleap writes
every protein chain into one unbroken numbering, which draws a cartoon bond straight from the
receptor's C-terminus into Gαi's N-terminus), **HETATM records** for the ligand, lipids, waters and
ions (left as ATOM, MolStar reads the macrocycle as part of the polymer), and the **human OPRM1
numbering**, so a residue clicked in the movie is the residue the manuscript names. Reserved chain
letters — `L` ligand, `M` membrane, `W` water, `I` ion — are excluded from the protein alphabet,
and the JSON records which chain is the receptor so both renderers agree on the subject.

## The two renders

**`02_render_dashboard.py` — the diagnostic one.** The structure beside its own time series on one
clock. A membrane view (x–z) with the Cα trace coloured by per-residue displacement from frame 0,
which says *which part* of the bundle is moving; a pocket view (x–y) with the key residues
coloured red while they are within 4.5 Å of the ligand; a contact strip that is the occupancy
analysis accumulating as the movie plays; and the RMSD, ligand-RMSD and bilayer-thickness traces
with a cursor on the current frame. Pure matplotlib — no Node, no GPU.

Where `02.11.00`'s reduced `.npz` is present its traces are preferred, because they were measured
on *every* frame of the full-resolution trajectory and the ligand RMSD there is against the
deposited pose rather than against frame 0. Where it is absent everything is recomputed from the
movie itself, so a replica that has not been through `02.11.00` still renders. **The panel titles
say which source was used.**

**`03_render_molstar.js` — the presentable one.** Headless MolStar in software WebGL, driven the
same way as [`03.10.00`](../03.10.00_molstar_render/README.md): receptor cartoon, ligand in
ball-and-stick, annular lipids as lines, pocket waters and ions. MolStar loads the PDB + XTC pair
natively (its own topology/coordinates route), the camera is set **once** before the first frame
and never touched again — a per-frame camera reset would follow the receptor's own drift and hide
exactly the motion the movie is being made for.

```bash
npm install                                     # once
node 03_render_molstar.js --width 1920 --height 1080 --fps 25
```

## Reading the result

* The caption carries the **simulated time**, not the frame number alone, and says whether the
  frames were superposed. A movie without that line is unreadable a week later.
* A motion that takes fewer than ~3 frames is an artefact of the stride, not an event. Check
  `frame_spacing_ns` in the JSON before believing anything fast.
* The dashboard's pocket view is a **slab**, ±9 Å around the ligand in z. A residue that leaves
  the slab disappears; it has not left the protein.
* Both renders are **orthographic-ish projections of one fixed camera**. Depth is not cued, so
  two things that overlap on screen are not necessarily in contact — that is what the contact
  strip and `02.11.00`'s occupancy table are for.

## Known limits

* The movie is a *subsample*. The default 300 frames of a 500 ns replica is one frame per 1.7 ns;
  fast events are absent rather than slow.
* Bulk water is not exported at all, so the movie cannot show hydration of anything except the
  pocket waters that were ranked in.
* The `--waters` ranking is computed over `--scan-frames` (default 40) sampled frames, not every
  frame: it finds persistent waters reliably and short-lived ones only by luck.
* The dashboard's 3D panels are projections, not renders — no depth sorting, no perspective. For
  anything that has to be looked at rather than measured, use the MolStar render or open the
  PDB/XTC pair in ChimeraX.
* MolStar renders here in **software** WebGL (swiftshader), around 0.3 s per frame. A 300-frame
  replica is ~2 minutes; the whole panel of 30 replicas is an hour.

## The files are the point

`<replica>_movie.pdb` + `<replica>_movie.xtc` is a standard topology/trajectory pair. Both movies
are conveniences on top of it — the pair itself opens directly in **VMD, PyMOL, ChimeraX and the
MolStar web viewer**, where you can rotate it, and that is usually the fastest way to chase down
something a rendered movie only hinted at.

```bash
chimerax intermediate/02.12.00_render_trajectory_movies/apo_ASH/prod_r1_movie.pdb \
         intermediate/02.12.00_render_trajectory_movies/apo_ASH/prod_r1_movie.xtc
```
