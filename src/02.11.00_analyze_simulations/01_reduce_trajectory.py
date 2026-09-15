#!/usr/bin/env python3
"""Reduce one production replica to a compact set of QC and interaction observables.

**Runs on the cluster**, in `zh853mor-prep` (MDAnalysis), where the trajectories are: a 500 ns
replica sampled every 100 ps is ~5 GB, and the panel is 3 replicas x up to 10 systems. Everything
downstream works on what this writes -- ~7 MB per replica, nine tenths of it the two per-frame
feature matrices -- so the aggregation, the figures and the conformational landscapes run locally
on a laptop without moving a single DCD. (~200 MB for the whole panel; `--stride` shrinks it
proportionally if that ever matters more than the time resolution of the tICA lag.)

One pass over the trajectory computes, per frame:

  * Ca RMSD against the staged OPM-oriented receptor, whole and TM-only, superposed
  * receptor-aligned ligand RMSD (pose retention) and the ligand's own internal RMSD and Rg
  * minimum heavy-atom distance from EVERY receptor residue to the ligand (the contact matrix)
  * pairwise CA-CA distances among the key/functional residues -- the shared conformational
    feature set `02.13.00` fits its tICA basis on, and the only one apo and holo have in common
  * polar (N/O-N/O) distances and bridging-water counts for the anchor residues
  * activation rulers (R3.50-T6.34, R3.50-Y7.53), the C142-C219 disulfide, the D2.50 Na+ site
  * bilayer thickness, gross area per lipid, box, and the TM z-registration

and then, after locating the end of the relaxation with `convergence.detect_equilibration`,
per-residue RMSF, a Ca PCA with its cosine content, and correlation-corrected means.

Residue numbers in the outputs are HUMAN OPRM1 (D149, E231, H299 ...), not the prmtop's, which
tleap renumbers from 1 -- see `zh853mor.md.construct_residue_map`.

Usage (from anywhere; paths are resolved against the repository):

    python 01_reduce_trajectory.py --all                       # every replica of every build
    python 01_reduce_trajectory.py --build intermediate/02.10.00_build/apo_ASH_20260901_133002
    python 01_reduce_trajectory.py --list                      # (build, replica) pairs, numbered
    python 01_reduce_trajectory.py --index 4                   # the 4th pair (SLURM array task)

Writes intermediate/02.11.00_analyze_simulations/<ligand>_<D250>/<replica>.{npz,json}.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402
from MDAnalysis.analysis.align import rotation_matrix  # noqa: E402
from MDAnalysis.lib.distances import self_distance_array  # noqa: E402

from zh853mor import convergence as cv  # noqa: E402
from zh853mor import md, paths  # noqa: E402

BUILD_ROOT = paths.INTERMEDIATE / "02.10.00_build"
OUT_ROOT = paths.INTERMEDIATE / "02.11.00_analyze_simulations"

# Equilibration is detected on ONE observable and applied to all of them, so that every mean for a
# replica describes the same window. The TM-only Ca RMSD is the right choice: it is the slowest
# thing that must settle before any pocket measurement means anything, and unlike the whole-
# receptor RMSD it is not dominated by the floppy termini and ICL3.
EQUILIBRATION_OBSERVABLE = "rmsd_ca_tm"

# A ligand that starts this far from its reference pose has not moved -- the reference is wrong.
# 15 A is well past any real first-frame deviation (the deposited pose IS frame 0's starting
# point, up to minimisation and the equilibration restraints) and well short of the ~90 A that a
# reference in the wrong coordinate frame produces.
LIG_POSE_SANITY_A = 15.0


@dataclass
class Job:
    build: md.BuildDir
    replica: str
    dcd: Path


def jobs(builds: list[md.BuildDir]) -> list[Job]:
    return [Job(b, name, dcd) for b in builds for name, dcd in b.replicas()]


def listing(todo: list[Job]) -> list[str]:
    """One tab-separated line per job: index, system, replica, frames, ns, target, status, path.

    Read from the DCD header and the state log, so listing the whole panel is instant. The length
    columns are the point: a replica that died at its wall-time or is still being written looks
    exactly like a finished one to a `glob`, and reducing it silently analyses a partial run.
    """
    targets: dict[str, float | None] = {}
    rows = []
    for n, job in enumerate(todo, 1):
        name = job.build.name
        if name not in targets:
            targets[name] = job.build.sampling().get("ZH_PROD_NS")
        try:
            span = md.replica_span(job.dcd, targets[name])
            frames, ns, status = f"{span.n_frames}", f"{span.ns:.1f}", span.status
        except (OSError, ValueError, KeyError, IndexError) as exc:
            frames, ns, status = "?", "?", f"unreadable ({type(exc).__name__})"
        target = targets[name]
        rows.append("\t".join([str(n), name, job.replica, frames, ns,
                                "?" if target is None else f"{target:.0f}", status,
                                str(job.dcd)]))
    return rows


def inventory_notes(builds: list[md.BuildDir]) -> list[str]:
    """Per-build remarks that are not about any single replica -- chiefly MISSING ones.

    A build configured for 3 replicas that produced 2 is invisible in a per-replica listing:
    everything present looks fine, and the replicate spread quietly rests on two points.
    """
    notes = []
    for b in builds:
        present = {name for name, _ in b.replicas()}
        want = int(b.sampling().get("ZH_REPLICAS", 0))
        if want and len(present) < want:
            missing = [f"prod_r{i}" for i in range(1, want + 1) if f"prod_r{i}" not in present]
            notes.append(f"note: {b.name} has {len(present)} of {want} configured replicas"
                         + (f" (missing {', '.join(missing)})" if missing else ""))
    return notes


def _minimum_image(positions: np.ndarray, anchor: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Shift a whole molecule by box vectors so its centroid is nearest `anchor`.

    OpenMM writes coordinates wrapped into the periodic box MOLECULE BY MOLECULE, so the ligand
    -- a separate molecule from the receptor -- can be imaged to the far side of the box while
    still sitting in the pocket. Distances computed with `box=` are immune, but superposition is
    not: without this the ligand RMSD jumps by a box length and the pose looks lost.
    """
    lengths = np.asarray(box[:3], dtype=float)
    disp = positions.mean(axis=0) - np.asarray(anchor, dtype=float)
    return positions - np.round(disp / lengths) * lengths


def summarize(series: np.ndarray, t0: int, frames_per_ns: float) -> dict[str, float]:
    """Post-equilibration mean, correlation-corrected SEM, and residual drift (per ns)."""
    tail = np.asarray(series, dtype=float)[t0:]
    tail = tail[np.isfinite(tail)]
    if tail.size < 3:
        return {"mean": float("nan"), "sem": float("nan"), "g": float("nan"),
                "drift_per_ns": float("nan"), "drift_se": float("nan")}
    g = cv.statistical_inefficiency(tail)
    mean, sem = cv.mean_sem(tail, g)
    # linear_drift works in frames; the thresholds and the report are per ns.
    slope, slope_se = cv.linear_drift(tail)
    return {"mean": mean, "sem": sem, "g": g,
            "drift_per_ns": slope * frames_per_ns, "drift_se": slope_se * frames_per_ns}


def reduce_replica(job: Job, out_dir: Path, stride: int = 1,
                   topology: Path | None = None) -> dict[str, object]:
    """Compute every observable for one replica and write the .npz/.json pair."""
    build = job.build
    warnings: list[str] = []
    u = mda.Universe(str(topology or build.prmtop), str(job.dcd))
    ref_u = mda.Universe(str(build.receptor_pdb))

    # Human OPRM1 numbering, transferred from the staged receptor; `warnings` collects the
    # protonation/disulfide form differences tleap introduces (CYX at the disulfide, always).
    resids = md.construct_residue_map(u, build.receptor_pdb, notes=warnings)
    prot = u.select_atoms("protein")
    heavy = prot.select_atoms("prop mass > 2.0")               # H is invisible at 3.5 A anyway
    ca = prot.select_atoms("name CA")
    polar_heavy = u.select_atoms(
        f"protein and ((prop mass > {md.MASS_N[0]} and prop mass < {md.MASS_N[1]}) or "
        f"(prop mass > {md.MASS_O[0]} and prop mass < {md.MASS_O[1]}))")

    ref_ca = ref_u.select_atoms("protein and name CA").positions.copy()
    if len(ref_ca) != len(ca):
        raise SystemExit(f"ERROR: {build.receptor_pdb.name} has {len(ref_ca)} CA, the topology "
                         f"{len(ca)}. The build directory is inconsistent.")

    # --- ligand -------------------------------------------------------------------------------
    resname = build.ligand_resname()
    lig = u.select_atoms(f"resname {resname} and prop mass > 2.0") if resname else u.atoms[[]]
    has_ligand = len(lig) > 0
    if resname and not has_ligand:
        warnings.append(f"system.json declares ligand resname {resname} but no atoms match it")
    lig_ref: np.ndarray | None = None
    lig_ref_source = "none"
    if has_ligand:
        # The deposited pose, transferred through the build: receptor.pdb holds the grafted
        # ligand under its INPUT resname (L01), so match on atom names rather than on order.
        cand = ref_u.select_atoms("not protein")
        ref_names = [a.name.strip() for a in cand]
        lig_names = [a.name.strip() for a in lig]
        # The name -> position map is only a valid transfer if each name it is asked for appears
        # EXACTLY ONCE on each side. `len(by_name) >= len(lig)` was too weak: a duplicated name
        # silently keeps whichever atom came last, and the reference pose is then a scrambled
        # molecule -- a wrong RMSD that still looks like an RMSD. The reference is allowed to be a
        # superset (it may carry hydrogens the mass-filtered topology selection dropped).
        dupes = {n for n in ref_names if ref_names.count(n) > 1}
        missing = sorted({n for n in lig_names if n not in set(ref_names)})
        ambiguous = sorted({n for n in lig_names if n in dupes})
        repeated = len(set(lig_names)) != len(lig_names)
        if not (missing or ambiguous or repeated):
            by_name = dict(zip(ref_names, cand.positions, strict=True))
            lig_ref = np.array([by_name[n] for n in lig_names], dtype=float)
            lig_ref_source = "receptor.pdb (deposited pose)"
        else:
            why = (f"{len(missing)} topology ligand atom name(s) absent from receptor.pdb "
                   f"({', '.join(missing[:5])}{'...' if len(missing) > 5 else ''})" if missing
                   else f"ambiguous name(s) {', '.join(ambiguous[:5])}" if ambiguous
                   else "the topology's ligand has repeated atom names")
            warnings.append(
                f"cannot transfer the deposited pose: {why}. lig_rmsd_pose therefore measures "
                "displacement from production frame 0, NOT from the deposited pose -- it says "
                "whether the ligand stayed where it started, not whether it kept the cryo-EM "
                "binding mode. `ligand_reference` in this replica's JSON records which was used")

    # --- fixed selections ---------------------------------------------------------------------
    phos = md.select_by_mass(u, *md.MASS_P)
    water_o = md.select_by_mass(u, *md.MASS_O, extra=f"resname {' '.join(sorted(md.WATER_RESN))}")
    sodium = md.select_by_mass(u, *md.MASS_NA)
    n_phos, n_sterol = md.count_lipids(u)
    if not n_phos:
        warnings.append("no phosphorus found: membrane observables are undefined")

    polar_sel = (f"(prop mass > {md.MASS_N[0]} and prop mass < {md.MASS_N[1]}) or "
                 f"(prop mass > {md.MASS_O[0]} and prop mass < {md.MASS_O[1]})")
    lig_polar = lig.select_atoms(polar_sel) if has_ligand else u.atoms[[]]

    contact = md.ResidueMinDistance(heavy)
    if contact.offsets.size != len(resids):
        raise SystemExit(f"ERROR: {contact.offsets.size} contact rows for {len(resids)} protein "
                         "residues; the per-residue occupancy would be shifted.")
    polar_contact = md.ResidueMinDistance(polar_heavy)
    # Which construct residue each row of `polar_contact` belongs to: not every residue has a
    # polar sidechain atom, so the polar rows are a subset of the contact rows. Likewise the Ca
    # set, which the RMSF and the PCA are computed on: the ACE and NME caps have no Ca, so it is
    # TWO residues shorter than the receptor and needs its own residue axis.
    polar_resids = md.construct_ids_for(polar_heavy, prot, resids)
    ca_resids = md.construct_ids_for(ca, prot, resids)

    index_of = {int(r): i for i, r in enumerate(resids)}
    anchor_rows = [index_of[r] for r in md.ANCHORS if r in index_of]
    anchor_ids = [r for r in md.ANCHORS if r in index_of]
    if len(anchor_ids) != len(md.ANCHORS):
        warnings.append("some anchor residues are absent from the receptor construct")
    anchor_polar_groups = [prot.residues[index_of[r]].atoms.select_atoms(polar_sel)
                           for r in anchor_ids]

    def ca_of(resid: int):
        i = index_of.get(resid)
        return prot.residues[i].atoms.select_atoms("name CA") if i is not None else u.atoms[[]]

    activation_pairs = [(label, ca_of(a), ca_of(b)) for label, a, b in md.ACTIVATION_DISTANCES]
    for label, a, b in activation_pairs:
        if not (len(a) and len(b)):
            warnings.append(f"activation ruler {label} is missing a Ca atom")
    # The sodium site is the D2.50 SIDECHAIN carboxylate. Selecting every oxygen of the residue
    # would include the backbone O, which points along the helix and would report a shorter Na+
    # distance than the site actually has.
    d250_carboxylate = u.atoms[[]]
    if 116 in index_of:
        d250_res = prot.residues[index_of[116]].atoms
        d250_carboxylate = d250_res.select_atoms("name OD1 OD2")
        if not len(d250_carboxylate):
            d250_carboxylate = d250_res.select_atoms(
                f"({polar_sel}) and not backbone and not name N O")
            warnings.append("D2.50 has no OD1/OD2; the Na+ distance uses its non-backbone "
                            f"polar atoms instead ({len(d250_carboxylate)} atoms)")
    # --- the shared conformational feature set ------------------------------------------------
    # Pairwise CA-CA distances among the key/functional residues, kept per frame. This is the
    # feature axis `02.13.00` fits its tICA basis on, and it has to be produced HERE because it is
    # the only place the trajectories are read: a separate featurisation pass would be a second
    # ~100 GB read for arrays this step already has the coordinates to compute.
    #
    # CA rather than sidechain tips, for the reason the activation rulers give (D-19): at 3.5 A
    # the rotamers are the least reliable coordinates in the starting model, so a sidechain
    # featurisation would let tICA find slow modes in the model's guesses. CA distances are also
    # defined for EVERY system including apo, which is what lets the apo arm share a landscape
    # with the holo ones -- the ligand-contact features (`min_dist`) cannot.
    key_ids = [r for r in sorted(md.KEY_RESIDUES) if r in index_of]
    key_ca = u.atoms[[]]
    for r in key_ids:
        key_ca = key_ca + ca_of(r)
    if len(key_ca) != len(key_ids):
        missing = len(key_ids) - len(key_ca)
        warnings.append(f"{missing} key residue(s) have no CA; the conformational feature set is "
                        "that much smaller and is NOT comparable with another build's")
        key_ids = [r for r, g in zip(key_ids, (ca_of(r) for r in key_ids), strict=True) if len(g)]
    key_pair_ids = np.array([(a, b) for i, a in enumerate(key_ids) for b in key_ids[i + 1:]],
                            dtype=int)

    ss_pair = [prot.residues[index_of[r]].atoms.select_atoms("name SG")
               for r in md.C142_C219 if r in index_of]
    # Na+ to the D2.50 carboxylate: one row (the single residue), minimised over every Na+.
    na_to_d250 = md.ResidueMinDistance(d250_carboxylate)

    # TM subset from the first frame: |z - midplane| < TM_HALF_A.
    u.trajectory[0]
    mid0 = md.membrane_midplane(phos.positions) if n_phos else 0.0
    tm_mask = np.abs(ca.positions[:, 2] - mid0) < md.TM_HALF_A
    if tm_mask.sum() < 100:
        warnings.append(f"only {int(tm_mask.sum())} Ca inside the membrane slab; "
                        "TM RMSD is computed over an unusually small set")

    dt_ps = float(u.trajectory.dt) or md.REPORT_PS
    if not np.isfinite(dt_ps) or dt_ps <= 0:
        dt_ps = md.REPORT_PS
        warnings.append(f"trajectory reports no timestep; assuming {md.REPORT_PS} ps per frame")

    # --- the frame loop -------------------------------------------------------------------
    n_res = len(resids)
    rows: dict[str, list[float]] = {k: [] for k in (
        "time_ns", "rmsd_ca_all", "rmsd_ca_tm", "registration_z", "thickness", "apl", "apl_net",
        "box_x", "box_y", "box_z", "lig_rmsd_pose", "lig_rmsd_frame0", "lig_internal_rmsd",
        "lig_rgyr", "na_d250_dist", "ss_dist")}
    min_dist = np.empty((0, n_res), dtype=np.float16)
    dist_rows, polar_rows, bridge_rows, activation_rows, ca_frames = [], [], [], [], []
    key_pair_rows: list[np.ndarray] = []
    lig_frame0: np.ndarray | None = None
    started = time.time()

    for i, ts in enumerate(u.trajectory[::stride]):
        box = ts.dimensions
        # OpenMM's DCDReporter writes the first frame AFTER one reporting interval, so frame 0
        # is at t = dt, not t = 0. Using (i+1) keeps this axis aligned with the state log's own
        # Time column, which the thermodynamic series are read from.
        rows["time_ns"].append((i * stride + 1) * dt_ps / 1000.0)
        rows["box_x"].append(float(box[0]))
        rows["box_y"].append(float(box[1]))
        rows["box_z"].append(float(box[2]))

        mob = ca.positions
        mob_c = mob - mob.mean(axis=0)
        ref_c = ref_ca - ref_ca.mean(axis=0)
        rot, rms_all = rotation_matrix(mob_c, ref_c)
        rows["rmsd_ca_all"].append(float(rms_all))
        # The TM subset is re-centred on ITSELF: superposing a subset about the whole-receptor
        # centroid leaves an uncorrected translation and reports an RMSD that is too large.
        tm_mob = mob[tm_mask] - mob[tm_mask].mean(axis=0)
        tm_ref = ref_ca[tm_mask] - ref_ca[tm_mask].mean(axis=0)
        rows["rmsd_ca_tm"].append(float(rotation_matrix(tm_mob, tm_ref)[1]))
        ca_frames.append((mob_c @ rot.T).astype(np.float32))   # superposed, for RMSF and PCA

        if n_phos:
            p = phos.positions
            mid = md.membrane_midplane(p)
            rows["thickness"].append(md.bilayer_thickness(p))
            # Gross AND protein-corrected: the receptor occupies a good fraction of the box
            # cross-section, so the gross value is not comparable with a pure-bilayer literature
            # area, while the hull correction overestimates the protein and makes the corrected
            # value a lower bound. Reporting both is the honest option (as at equilibration).
            xy = (float(box[0]), float(box[1]))
            rows["apl"].append(md.area_per_lipid(xy, n_phos + n_sterol))
            hull = md.protein_cross_section(heavy.positions, mid)
            rows["apl_net"].append(md.area_per_lipid(xy, n_phos + n_sterol, hull)
                                   if hull is not None else np.nan)
            rows["registration_z"].append(float(mob[tm_mask][:, 2].mean() - mid))
        else:
            rows["thickness"].append(np.nan)
            rows["apl"].append(np.nan)
            rows["apl_net"].append(np.nan)
            rows["registration_z"].append(np.nan)

        if has_ligand:
            lig_pos = _minimum_image(lig.positions, mob.mean(axis=0), box)
            # Into the REFERENCE's frame, by the same superposition the CA RMSD just used.
            aligned = (lig_pos - mob.mean(axis=0)) @ rot.T + ref_ca.mean(axis=0)
            if lig_frame0 is None:
                lig_frame0 = lig_pos.copy()          # box frame, for lig_rmsd_frame0 below
                if lig_ref is None:
                    # The fallback reference has to be the ALIGNED frame 0, not the raw one.
                    # `aligned` lives in receptor.pdb's OPM frame, which is centred near the
                    # origin, while `lig_pos` lives in the simulation box, which runs 0..L -- so
                    # a raw frame-0 fallback measured the ~90 A offset between those two origins
                    # and reported it as the ligand's pose RMSD, near-constant for the whole run.
                    lig_ref = aligned.copy()
                    lig_ref_source = "production frame 0 (receptor-aligned)"
                # Whatever the reference, frame 0 must sit close to it: a pose RMSD that starts
                # tens of angstroms out is a reference in the wrong frame or matched to the wrong
                # atoms, never a ligand that has moved. Cheap, and it fails loudly at frame 0
                # rather than producing a plausible-looking series nobody questions.
                start = float(np.sqrt(((aligned - lig_ref) ** 2).sum(axis=1).mean()))
                if start > LIG_POSE_SANITY_A:
                    warnings.append(
                        f"ligand pose RMSD is {start:.1f} A at the FIRST frame, against "
                        f"'{lig_ref_source}'. A reference in the right frame starts within a "
                        f"couple of angstroms; {LIG_POSE_SANITY_A} A or more means the reference "
                        "is in the wrong coordinate frame or matched to the wrong atoms, and "
                        "lig_rmsd_pose is not a pose RMSD for this replica")
            assert lig_ref is not None
            rows["lig_rmsd_pose"].append(
                float(np.sqrt(((aligned - lig_ref) ** 2).sum(axis=1).mean())))
            rows["lig_rmsd_frame0"].append(
                float(np.sqrt(((lig_pos - lig_frame0) ** 2).sum(axis=1).mean())))
            rows["lig_internal_rmsd"].append(float(rotation_matrix(
                lig_pos - lig_pos.mean(axis=0), lig_ref - lig_ref.mean(axis=0))[1]))
            rows["lig_rgyr"].append(float(lig.radius_of_gyration()))
            dist_rows.append(contact(lig_pos, box).astype(np.float16))
            polar_rows.append(polar_contact(lig_pos, box))
            bridge_rows.append(md.water_bridges(
                water_o.positions, [g.positions for g in anchor_polar_groups],
                lig_polar.positions, box))
        else:
            for key in ("lig_rmsd_pose", "lig_rmsd_frame0", "lig_internal_rmsd", "lig_rgyr"):
                rows[key].append(np.nan)

        activation_rows.append([
            float(np.linalg.norm(a.positions[0] - b.positions[0])) if len(a) and len(b) else np.nan
            for _, a, b in activation_pairs])
        if len(key_ca) > 1:
            # No `box=`: these are intra-protein distances and the receptor spans more than half
            # the box, where the minimum-image convention would fold a genuine 60 A separation
            # back to 31 A. The protein arrives whole because OpenMM wraps molecule by molecule,
            # which is the same assumption the CA RMSD above already rests on.
            key_pair_rows.append(self_distance_array(key_ca.positions).astype(np.float16))
        rows["na_d250_dist"].append(
            float(na_to_d250(sodium.positions, box)[0])
            if len(d250_carboxylate) and len(sodium) else np.nan)
        rows["ss_dist"].append(
            float(np.linalg.norm(ss_pair[0].positions[0] - ss_pair[1].positions[0]))
            if len(ss_pair) == 2 and all(len(g) for g in ss_pair) else np.nan)

    n_frames = len(rows["time_ns"])
    if n_frames < 10:
        raise SystemExit(f"ERROR: {job.dcd.name} has only {n_frames} frames; nothing to analyse.")
    series = {k: np.array(v, dtype=float) for k, v in rows.items()}
    if dist_rows:
        min_dist = np.array(dist_rows, dtype=np.float16)
    key_pairs = (np.array(key_pair_rows, dtype=np.float16) if key_pair_rows
                 else np.empty((0, len(key_pair_ids)), dtype=np.float16))
    activation = np.array(activation_rows, dtype=float)
    frames_per_ns = 1000.0 / (dt_ps * stride)

    # --- equilibration, then everything that depends on the window --------------------------
    t0, g_eq, n_eff = cv.detect_equilibration(series[EQUILIBRATION_OBSERVABLE])
    ca_arr = np.array(ca_frames, dtype=np.float32)
    prod = ca_arr[t0:]
    rmsf = np.sqrt(((prod - prod.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))
    if rmsf.size != ca_resids.size:  # the two must stay paired or every residue label shifts
        raise SystemExit(f"ERROR: {rmsf.size} RMSF values for {ca_resids.size} Ca atoms.")
    proj, explained = md.principal_components(prod, n_components=3)
    cosine = [cv.cosine_content(proj[:, k], index=k + 1) for k in range(proj.shape[1])]

    means = {k: summarize(v, t0, frames_per_ns) for k, v in series.items() if k != "time_ns"}
    for j, (label, _, _) in enumerate(activation_pairs):
        means[label] = summarize(activation[:, j], t0, frames_per_ns)

    occ: dict[str, float] = {}
    polar_occ: dict[str, float] = {}
    bridges: dict[str, float] = {}
    per_residue_occupancy = np.zeros(n_res)
    if len(min_dist):
        per_residue_occupancy = md.occupancy(min_dist[t0:].astype(float))
        polar_arr = np.array(polar_rows, dtype=float)[t0:]
        bridge_arr = np.array(bridge_rows, dtype=float)[t0:]
        for k, resid in enumerate(anchor_ids):
            occ[str(resid)] = float(per_residue_occupancy[anchor_rows[k]])
            bridges[str(resid)] = float(bridge_arr[:, k].mean())
        for k, resid in enumerate(polar_resids):
            if int(resid) in md.ANCHORS:
                polar_occ[str(int(resid))] = float((polar_arr[:, k] <= md.HBOND_CUT).mean())

    # --- thermodynamics from the state log ----------------------------------------------------
    log_series: dict[str, np.ndarray] = {}
    log_path = job.dcd.with_suffix(".log")
    if log_path.exists():
        try:
            log_series = md.read_state_log(log_path)
        except (OSError, ValueError) as exc:
            warnings.append(f"{log_path.name}: {exc}")
    else:
        warnings.append(f"{log_path.name} not found; thermodynamic QC is unavailable")
    log_means = {f"log_{k}": summarize(v, int(len(v) * t0 / max(n_frames, 1)), frames_per_ns)
                 for k, v in log_series.items() if k not in ("step", "time_ps")}

    summary: dict[str, object] = {
        "build": build.name, "build_dir": str(build.path), "ligand": build.ligand,
        "d250": build.d250, "stamp": build.stamp, "replica": job.replica,
        "n_frames": n_frames, "stride": stride, "dt_ps": dt_ps,
        "length_ns": float(series["time_ns"][-1]),
        # What the build was CONFIGURED to run, so the report can say 143 of 500 ns rather than
        # just 143 ns -- a replica killed at its wall-time is otherwise indistinguishable.
        "target_ns": build.sampling().get("ZH_PROD_NS"),
        "n_atoms": int(len(u.atoms)),
        "n_protein_residues": n_res, "n_phospholipids": n_phos, "n_sterols": n_sterol,
        "n_waters": int(len(water_o)), "n_sodium": int(len(sodium)),
        "n_key_residues": len(key_ids), "n_key_pairs": int(len(key_pair_ids)),
        "ligand_resname": resname, "n_ligand_atoms": int(len(lig)),
        "ligand_reference": lig_ref_source,
        "equilibration": {"observable": EQUILIBRATION_OBSERVABLE, "t0_frames": int(t0),
                          "t0_ns": float(series["time_ns"][t0]), "g": g_eq, "n_eff": n_eff},
        "means": {**means, **log_means},
        "occupancy": occ, "polar_occupancy": polar_occ, "water_bridges": bridges,
        "pca": {"explained": [float(x) for x in explained],
                "cosine_content": [float(c) for c in cosine]},
        "warnings": warnings,
        "wall_seconds": round(time.time() - started, 1),
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_dir / f"{job.replica}.npz",
        resids=resids, ca_resids=ca_resids, rmsf=rmsf, occupancy=per_residue_occupancy,
        min_dist=min_dist, anchor_resids=np.array(anchor_ids, dtype=int),
        key_pairs=key_pairs, key_pair_ids=key_pair_ids,
        activation=activation, activation_labels=np.array(
            [label for label, _, _ in activation_pairs]),
        pc_projection=proj.astype(np.float32), t0=t0,
        **{k: v.astype(np.float32) for k, v in series.items()},
        **{f"log_{k}": v.astype(np.float32) for k, v in log_series.items()})
    (out_dir / f"{job.replica}.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="append", default=[],
                    help="build directory (repeatable); default is every build found")
    ap.add_argument("--replica", help="analyse only this replica, e.g. prod_r1")
    ap.add_argument("--all", action="store_true", help="every replica of every build")
    ap.add_argument("--index", type=int, help="analyse the Nth (build, replica) pair, 1-based "
                                              "-- this is what the SLURM array task passes")
    ap.add_argument("--list", action="store_true", help="print the numbered pairs and exit")
    ap.add_argument("--stride", type=int, default=1, help="use every Nth frame")
    ap.add_argument("--topology", type=Path,
                    help="topology to use instead of the build's system.prmtop")
    ap.add_argument("--out", type=Path, default=OUT_ROOT)
    ap.add_argument("--force", action="store_true", help="recompute even if outputs exist")
    ap.add_argument("--all-builds", action="store_true",
                    help="include superseded builds, not just the newest per system")
    args = ap.parse_args()

    if args.build:
        builds = []
        for path in (Path(b).resolve() for b in args.build):
            found = md.discover_builds(path.parent, newest_only=False)
            match = [b for b in found if b.path == path]
            if not match:
                raise SystemExit(f"ERROR: {path} is not a build directory "
                                 "(expected <ligand>_<ASP|ASH>_<timestamp>/system.prmtop).")
            builds += match
    else:
        if not BUILD_ROOT.is_dir():
            raise SystemExit(f"ERROR: {BUILD_ROOT} does not exist. Build the systems first "
                             "(src/02.10.00_slurm_bundle/README.md).")
        builds = md.discover_builds(BUILD_ROOT, newest_only=not args.all_builds)

    todo = jobs(builds)
    if args.replica:
        todo = [j for j in todo if j.replica == args.replica]
    if not todo:
        raise SystemExit("ERROR: no production trajectories (prod_r*.dcd) found.")

    if args.list:
        for row in listing(todo):
            print(row)
        for note in inventory_notes(builds):   # stderr: the listing stays one line per job
            print(note, file=sys.stderr)
        return 0

    if args.index is not None:
        if not 1 <= args.index <= len(todo):
            raise SystemExit(f"ERROR: --index {args.index} outside 1-{len(todo)}.")
        todo = [todo[args.index - 1]]
    elif not (args.all or args.build or args.replica):
        raise SystemExit("ERROR: nothing selected. Pass --all, --build, --index or --list.")

    failures = 0
    for job in todo:
        out_dir = args.out / job.build.name
        target = out_dir / f"{job.replica}.json"
        if target.exists() and not args.force:
            print(f"skip {job.build.name}/{job.replica}: {target.name} exists (--force to redo)")
            continue
        print(f"--- {job.build.name} / {job.replica}: {job.dcd}", flush=True)
        try:
            summary = reduce_replica(job, out_dir, stride=args.stride, topology=args.topology)
        except Exception as exc:  # noqa: BLE001 -- see below
            # Deliberately broad: a replica killed at its wall-time leaves a TRUNCATED DCD, and
            # the readers raise whatever their decompressor or parser happens to raise (EOFError,
            # struct.error, ...). In a job array one unusable replica must not take the others
            # with it, so the type is reported rather than filtered, and the exit code still says
            # something failed.
            print(f"FAILED {job.build.name}/{job.replica}: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            failures += 1
            continue
        eq = summary["equilibration"]
        assert isinstance(eq, dict)
        print(f"    {summary['n_frames']} frames, {summary['length_ns']:.0f} ns; "
              f"equilibrated at {eq['t0_ns']:.0f} ns (n_eff {eq['n_eff']:.0f}); "
              f"{summary['wall_seconds']:.0f} s")
        for w in summary["warnings"]:  # type: ignore[union-attr]
            print(f"    WARNING: {w}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
