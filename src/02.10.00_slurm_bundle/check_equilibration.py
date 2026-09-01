#!/usr/bin/env python3
"""Decide whether an equilibrated membrane system is ready for production.

Runs after 02_equilibrate.py (step 3) and before the unrestrained pre-production leg (step 3.5).
The question it answers is not "did the job finish" but "has the system stopped changing" --
which for a PACKMOL-built bilayer is governed by the lipids, not the protein. The protein relaxes
in a few hundred ps; the membrane, built in an artificially ordered extended conformation, needs
tens of ns for area per lipid to converge. So a clean PASS here on a 2.25 ns ramp means "nothing
is broken", NOT "sampling can start" -- see the note printed at the end.

Reads what 02_equilibrate.py now writes: the state log (thermodynamics), the DCD (the periodic box
per frame, which is what area-per-lipid and thickness are computed from), and the stage manifest,
so convergence is judged on the FINAL, least-restrained stage rather than smeared over the ramp.

Usage: python check_equilibration.py --top system.prmtop --eq system_eq --receptor receptor.pdb

Exits non-zero only on FAIL (something is physically wrong). Non-convergence is a WARN: it is the
expected state after a short ramp and the reason the pre-production leg exists.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.rms import rmsd as rmsd_superposed
from MDAnalysis.lib.distances import capped_distance

# --- thresholds -------------------------------------------------------------------------------
# Ranges are for 9:1 POPC:cholesterol + OPC water + 0.15 M NaCl at 310 K (SPECIFICATION D-12).
DENSITY_RANGE = (0.98, 1.10)   # g/mL; < 0.98 means a void, i.e. the packing cell was too large
DENSITY_DRIFT = 0.003          # g/mL tolerated between halves of the final stage
TEMP_TOL_K = 2.0               # |mean T - 310| for a system of this size
BOX_DRIFT_A = 0.75             # per-axis A between halves of the final stage
APL_RANGE = (55.0, 70.0)       # A^2/lipid: pure POPC is ~64-68, cholesterol condenses it to ~60-64
THICK_RANGE = (34.0, 44.0)     # A, phosphate-to-phosphate; the build measured 38.6 A
CORE_WATER_MAX = 10            # waters in the tail region and NOT inside the protein
CA_RMSD_MAX = 3.0              # A, whole receptor vs the staged (OPM-oriented) receptor
TM_RMSD_MAX = 2.0              # A, membrane-embedded CA only -- the number that should be tight
REGISTRATION_DRIFT_A = 2.0     # A of vertical drift out of the OPM slab over the run
SG_RANGE = (1.90, 2.30)        # A, C142-C219 disulfide
SS_CUTOFF_A = 2.50             # SG-SG separation counted as bonded (cf. make_tleap.py)
LIG_RMSD_MAX = 3.0             # A, receptor-aligned ligand drift over the ramp
CORE_HALF_A = 8.0              # |z - midplane| defining the hydrophobic core
TM_HALF_A = 15.0               # |z - midplane| defining the membrane-embedded CA set

# Lipid21 is MODULAR: it does not have a POPC residue. Each phospholipid is split into a headgroup
# residue plus one residue per acyl chain -- POPC is PC + PA (palmitoyl) + OL (oleoyl) -- so a
# `resname POPC` selection matches nothing in the prmtop even though the packed PDB used that name
# (which is why the build-time checks, which read the PDB, are unaffected). Counting *residues*
# would also treat one POPC as three lipids and put the area per lipid out by 3x.
#
# So molecules are counted naming-agnostically: one phosphorus per phospholipid, plus sterol
# residues. Nothing else in this system contains P.
STEROL_RESN = {"CHL", "CHL1", "CLR", "CHOL", "CHO"}
WATER_RESN = {"WAT", "HOH", "OPC", "TIP3", "SOL"}
ION_RESN = {"NA", "CL", "K", "NA+", "CL-", "K+", "MG", "CA2", "ZN"}

# Atom NAMES are the least portable thing in a prmtop and have now bitten twice: Lipid21 calls the
# phosphatidylcholine phosphorus `P31`, not `P`, and MDAnalysis reports empty elements for this
# topology ("Unknown ATOMIC_NUMBER"), so neither `name P` nor `element P` finds it. Masses are
# always present and unambiguous, so select on those. HMR does not perturb this: it is applied when
# OpenMM builds the System, never written back to the prmtop, and phosphorus carries no hydrogens
# anyway. Windows are tight enough to exclude the neighbours -- Na 22.99, S 32.06, Cl 35.45.
MASS_P = (30.5, 31.5)
MASS_O = (15.5, 16.5)


class Report:
    """Collects named checks so the console table, the JSON and the exit code stay consistent."""

    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str, str]] = []
        self.values: dict[str, object] = {}
        self.failed = False
        self.warned = False

    def add(self, name, value, verdict, detail=""):
        self.rows.append((name, value, verdict, detail))
        if verdict == "FAIL":
            self.failed = True
        elif verdict == "WARN":
            self.warned = True

    def verdict_of(self, name):
        return next((r[2] for r in self.rows if r[0] == name), None)

    def judge_drift(self, name, value, tol, fmt="{:.2f}", unit="", soft=True, detail=""):
        """Judge |drift| against a tolerance but DISPLAY the sign.

        The sign is the diagnostically useful half: an area per lipid still falling means the
        membrane is condensing and needs longer, while the same magnitude oscillating about a
        plateau means it is done. Reporting abs() throws that away.
        """
        ok = abs(value) <= tol
        verdict = "PASS" if ok else ("WARN" if soft else "FAIL")
        self.values[name] = float(value)
        arrow = "" if ok else ("  (falling)" if value < 0 else "  (rising)")
        self.add(name, fmt.format(value) + unit + arrow, verdict,
                 detail or f"|drift| <= {tol}{unit} over the window")
        return ok

    def judge(self, name, value, lo, hi, fmt="{:.2f}", unit="", soft=False, detail=""):
        ok = lo <= value <= hi
        verdict = "PASS" if ok else ("WARN" if soft else "FAIL")
        self.values[name] = float(value)
        self.add(name, fmt.format(value) + unit, verdict,
                 detail or f"expected {fmt.format(lo)}-{fmt.format(hi)}{unit}")
        return ok


def read_state_log(path: Path) -> dict[str, np.ndarray]:
    """Parse an OpenMM StateDataReporter CSV by COLUMN NAME.

    The reporter emits columns in its own fixed order, not the order the kwargs were passed, so
    positional parsing silently mislabels energy as temperature the moment a field is added.
    """
    with open(path) as fh:
        rows = list(csv.reader(fh))
    if not rows:
        raise SystemExit(f"ERROR: {path} is empty.")
    header = [h.lstrip("#").strip().strip('"') for h in rows[0]]
    data = np.array([[float(x) for x in r] for r in rows[1:] if len(r) == len(header)])
    if data.size == 0:
        raise SystemExit(f"ERROR: {path} has a header but no data rows; did the run die early?")

    def col(*needles):
        for i, h in enumerate(header):
            if all(n.lower() in h.lower() for n in needles):
                return data[:, i]
        return None

    return {"step": col("Step"), "time_ps": col("Time"), "pe": col("Potential", "Energy"),
            "temp": col("Temperature"), "density": col("Density"), "volume": col("Box", "Volume")}


def halves_drift(x: np.ndarray) -> float:
    """Mean of the second half minus mean of the first: the plateau test."""
    if len(x) < 4:
        return float("nan")
    h = len(x) // 2
    return float(np.mean(x[h:]) - np.mean(x[:h]))


def residue_census(u, limit=14):
    """(resname, count) by descending count -- what the topology actually contains."""
    counts = {}
    for r in u.residues:
        counts[r.resname.strip()] = counts.get(r.resname.strip(), 0) + 1
    return sorted(counts.items(), key=lambda kv: -kv[1])[:limit]


def leaflet_counts(u, phos, midplane_z):
    """(upper, lower) lipid MOLECULE counts: one per phosphorus, plus sterol residues.

    Computed ONCE: lipids do not flip-flop on an equilibration timescale, and re-deriving it per
    frame means walking every residue in the system (mostly water) on every frame.
    """
    zp = phos.positions[:, 2]
    up, lo = int((zp >= midplane_z).sum()), int((zp < midplane_z).sum())
    for r in u.residues:
        if r.resname.strip().upper() in STEROL_RESN:
            if r.atoms.positions[:, 2].mean() >= midplane_z:
                up += 1
            else:
                lo += 1
    return up, lo


def protein_area(prot_pos, midplane_z) -> float | None:
    """xy convex-hull area of the membrane-embedded protein, to subtract from the box area.

    A hull overestimates a non-convex cross-section, so the area per lipid derived from it is a
    lower bound. Reported alongside the gross value rather than instead of it.
    """
    try:
        from scipy.spatial import ConvexHull
    except ImportError:
        return None
    sl = prot_pos[np.abs(prot_pos[:, 2] - midplane_z) < TM_HALF_A][:, :2]
    if len(sl) < 3:
        return None
    return float(ConvexHull(sl).volume)  # 'volume' is area for 2-D input


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", required=True, help="system.prmtop")
    ap.add_argument("--eq", default="system_eq",
                    help="equilibration output basename (.log/.dcd/_stages.json/.pdb)")
    ap.add_argument("--receptor", default="receptor.pdb", help="staged OPM-oriented receptor")
    ap.add_argument("--out", default="eq_qc")
    ap.add_argument("--lig", default=None,
                    help="ligand resname; defaults to system.json's. Skipped for an apo build")
    ap.add_argument("--skip-piercing", action="store_true")
    args = ap.parse_args()

    # system.json records what the build IS. Builds made before it existed have none, and there
    # the intent is genuinely unknown -- reporting "system.json expects a ligand" against a file
    # that does not exist sends the reader to the wrong place.
    meta_known = False
    if args.lig is None:
        try:
            args.lig = json.loads(Path("system.json").read_text()).get("ligand_resname") or "apo"
            meta_known = True
        except (OSError, ValueError):
            args.lig = "LIG"
    else:
        meta_known = True

    eq = Path(args.eq)
    log_p, dcd_p, stg_p = Path(f"{eq}.log"), Path(f"{eq}.dcd"), Path(f"{eq}_stages.json")
    # 02_equilibrate.py writes <eq>.pdb; 03_production.py (the pre-production leg) writes
    # <eq>_final.pdb. Either is a valid final frame for the piercing check.
    pdb_p = next((q for q in (Path(f"{eq}.pdb"), Path(f"{eq}_final.pdb")) if q.exists()),
                 Path(f"{eq}.pdb"))
    missing = [p for p in (log_p, dcd_p) if not p.exists()]
    if missing:
        print(f"ERROR: {', '.join(str(m) for m in missing)} not found.", file=sys.stderr)
        print("  Equilibrations run before the reporters were added to 02_equilibrate.py wrote only"
              "\n  the final .xml/.pdb, and cannot be assessed -- there is no time series. Re-run"
              "\n  step 3 with the current script.", file=sys.stderr)
        return 2

    rep = Report()
    log = read_state_log(log_p)

    # --- restrict the plateau tests to the final, least-restrained stage ------------------------
    final_lo, final_hi, n_stages = None, None, None
    if stg_p.exists():
        man = json.loads(stg_p.read_text())
        last = man["stages"][-1]
        final_lo, final_hi, n_stages = last["first_step"], last["last_step"], len(man["stages"])
        sel = (log["step"] >= final_lo) & (log["step"] <= final_hi)
        scope = f"stage {last['stage']}/{n_stages} (k_bb={last['k_backbone']:g})"
    else:
        sel = log["step"] >= log["step"][len(log["step"]) // 2]
        scope = "second half of the run (no stage manifest)"
    if sel.sum() < 4:
        sel = np.ones_like(log["step"], dtype=bool)
        scope += " -- too few samples, using the whole run"
    print(f"convergence scope: {scope}\n")

    # --- 1. thermodynamics ----------------------------------------------------------------------
    dens, temp, pe = log["density"][sel], log["temp"][sel], log["pe"][sel]
    rep.judge("density", float(dens.mean()), *DENSITY_RANGE, unit=" g/mL", fmt="{:.4f}")
    rep.judge_drift("density drift", halves_drift(dens), DENSITY_DRIFT,
                    unit=" g/mL", fmt="{:+.4f}")
    rep.judge("temperature", float(temp.mean()), 310 - TEMP_TOL_K, 310 + TEMP_TOL_K,
              unit=" K", fmt="{:.1f}")
    # Energy drift is scored against its own fluctuation: an absolute kJ/mol threshold is
    # meaningless across system sizes, but drift much larger than sigma means still-relaxing.
    pe_sigma = float(pe.std()) or 1.0
    rep.judge_drift("potential-energy drift", halves_drift(pe) / pe_sigma, 1.0,
                    unit=" sigma", fmt="{:+.2f}", detail="|drift| <= 1 sigma of its own fluctuation")

    # --- 2. box, membrane, protein (from the trajectory) ----------------------------------------
    u = mda.Universe(args.top, str(dcd_p))
    n_frames = len(u.trajectory)
    if n_frames < 2:
        print("ERROR: the DCD has fewer than 2 frames; nothing to assess.", file=sys.stderr)
        return 2

    ca = u.select_atoms("name CA")            # lipids/water carry no CA, so this is the protein
    # One phosphorus per phospholipid, found by mass rather than by name (see MASS_P above).
    phos = u.select_atoms(f"prop mass > {MASS_P[0]} and prop mass < {MASS_P[1]}")
    prot = u.select_atoms("protein")
    if not len(prot):
        # `protein` keys off a built-in resname list that may not know CYX/HID/ASH. Fall back to
        # "everything that is not solvent, ion or lipid" rather than to CA alone, which would
        # under-count the protein cross-section and the core-water proximity test.
        skip = WATER_RESN | ION_RESN | STEROL_RESN | {"PC", "PE", "PS", "PA", "OL", "ST", "MY",
                                                      "LAL", "DHA", "SA", "PGR", "POPC", "POPE"}
        prot = u.atoms[[a.index for a in u.atoms
                        if a.residue.resname.strip().upper() not in skip]]
    water = u.select_atoms(f"resname {' '.join(sorted(WATER_RESN))}")
    wat_o = water.select_atoms(f"prop mass > {MASS_O[0]} and prop mass < {MASS_O[1]}") \
        if len(water) else water
    sg = u.select_atoms("name SG")
    print(f"selections: {len(phos)} phosphorus, {len(ca)} CA, {len(prot)} protein atoms, "
          f"{len(wat_o)} waters, {len(sg)} SG")
    lig = u.select_atoms("") if args.lig.lower() == "apo" else u.select_atoms(f"resname {args.lig}")

    if not len(phos):
        print("ERROR: no phosphorus atoms found; is this a bilayer?", file=sys.stderr)
        print("  Selection is by mass "
              f"({MASS_P[0]}-{MASS_P[1]} amu), so this is not an atom-naming problem -- the "
              "topology genuinely has no P.", file=sys.stderr)
        print("  Residues present (name x count):", file=sys.stderr)
        for name, n in residue_census(u):
            print(f"    {name:<6} {n}", file=sys.stderr)
        return 2

    ref_ca = mda.Universe(args.receptor).select_atoms("name CA").positions.copy()
    ca_ok = len(ref_ca) == len(ca)
    if not ca_ok:
        rep.add("CA count vs receptor.pdb", f"{len(ca)} vs {len(ref_ca)}", "WARN",
                "cannot compute RMSD against the staged receptor")

    box, apl_g, apl_n, thick, rmsd_all, rmsd_tm, reg, sgd = ([] for _ in range(8))
    lig_rmsd, lig_ref = [], []
    tm_mask = None
    n_up, n_lo = leaflet_counts(u, phos, float(phos.positions[:, 2].mean()))
    per_leaflet = max((n_up + n_lo) / 2.0, 1.0)
    n_sterol = sum(1 for r in u.residues if r.resname.strip().upper() in STEROL_RESN)
    print(f"lipids: {n_up + n_lo} molecules ({len(phos)} phospholipid + {n_sterol} sterol; "
          f"{n_up} upper / {n_lo} lower), {per_leaflet:g} per leaflet")
    print(f"residues: " + ", ".join(f"{k} x{v}" for k, v in residue_census(u, 8)) + "\n")
    for ts in u.trajectory:
        lx, ly, lz = ts.dimensions[:3]
        box.append((lx, ly, lz))
        zp = phos.positions[:, 2]
        mid = float(zp.mean())
        up, lo = zp[zp >= mid], zp[zp < mid]
        thick.append(float(up.mean() - lo.mean()) if len(up) and len(lo) else np.nan)

        apl_g.append(lx * ly / per_leaflet)
        pa = protein_area(prot.positions, mid)
        apl_n.append((lx * ly - pa) / per_leaflet if pa is not None else np.nan)

        cap = ca.positions
        if tm_mask is None:
            tm_mask = np.abs(cap[:, 2] - mid) < TM_HALF_A
        reg.append(float(cap[tm_mask][:, 2].mean() - mid))
        if ca_ok:
            # SUPERPOSED. A raw coordinate difference is meaningless here: packmol-memgen
            # translates the solute when it packs (that shift is what check_placement.py
            # measures), so comparing absolute positions against the pre-packing receptor.pdb
            # reports the translation, not the conformational change -- it read 84.8 A on a
            # perfectly good system. Rigid-body motion is covered separately by the OPM
            # registration drift, which is measured against the bilayer midplane and so is
            # already translation-invariant.
            rmsd_all.append(float(rmsd_superposed(cap, ref_ca, superposition=True)))
            rmsd_tm.append(float(rmsd_superposed(cap[tm_mask], ref_ca[tm_mask],
                                                 superposition=True)))
        if len(sg) >= 2:
            # The construct has many cysteines (12 SG here), only one of which is the modelled
            # C142-C219 disulfide, so pick out bonded PAIRS by distance rather than assuming the
            # only two SG atoms in the system are the pair. Same approach make_tleap.py uses to
            # emit the  line in the first place.
            p_ = sg.positions
            dm = np.linalg.norm(p_[:, None, :] - p_[None, :, :], axis=2)
            iu = np.triu_indices(len(sg), k=1)
            close = [float(dm[i, j]) for i, j in zip(*iu) if dm[i, j] <= SS_CUTOFF_A]
            sgd.append(close)
        if len(lig):
            if not lig_ref:
                lig_ref.append(lig.positions.copy())
            # Receptor-aligned: the restraint ramp holds the protein, so drift here is the ligand
            # leaving the pocket rather than the whole system translating.
            lig_rmsd.append(float(np.sqrt(
                ((lig.positions - cap[tm_mask].mean(axis=0)
                  - (lig_ref[0] - ref_ca[tm_mask].mean(axis=0))) ** 2).sum(axis=1).mean()))
                if ca_ok else float(np.sqrt(((lig.positions - lig_ref[0]) ** 2).sum(axis=1).mean())))

    box = np.array(box)
    rep.judge_drift("box Lx drift", halves_drift(box[:, 0]), BOX_DRIFT_A, unit=" A", fmt="{:+.2f}")
    rep.judge_drift("box Lz drift", halves_drift(box[:, 2]), BOX_DRIFT_A, unit=" A", fmt="{:+.2f}")
    rep.values["box_final_A"] = [round(float(v), 2) for v in box[-1]]

    # Report BOTH estimates. The protein cross-section is subtracted using an xy convex hull,
    # which necessarily OVERestimates a non-convex 7TM bundle, so the net area per lipid is a
    # LOWER bound and the gross value an upper one. Judging on the net alone makes a marginal
    # reading look like a membrane problem when it may just be the hull.
    have_net = np.isfinite(apl_n).all()
    apl = np.array(apl_n if have_net else apl_g)
    apl_kind = "area/lipid" if have_net else "area/lipid (gross)"
    half = slice(-max(len(apl) // 2, 1), None)
    if have_net:
        rep.add("area/lipid (gross)", f"{float(np.array(apl_g)[half].mean()):.2f} A^2", "PASS",
                "no protein subtraction -- an upper bound; the true value lies between the two")
    rep.judge(apl_kind, float(apl[half].mean()), *APL_RANGE, unit=" A^2", soft=True,
              detail=f"expected {APL_RANGE[0]:.0f}-{APL_RANGE[1]:.0f} A^2; convex-hull protein "
                     "subtraction makes this a lower bound")
    # Convergence is judged on the GROSS drift, because that is box area alone. The net value
    # also moves when the protein's convex hull breathes, which is not membrane condensation: on
    # the first 100 ns leg the net drift read -1.29 A^2 of which only -0.65 was the box, the rest
    # being +64 A^2 of hull as the receptor relaxed. Judging on the net therefore overstates how
    # far the membrane still has to go.
    rep.judge_drift("area/lipid drift (gross)", halves_drift(np.array(apl_g)), 0.5,
                    unit=" A^2", fmt="{:+.2f}",
                    detail="|drift| <= 0.5 A^2 -- box area only; the membrane convergence signal")
    if have_net:
        rep.add("area/lipid drift (net)", f"{halves_drift(apl):+.2f} A^2", "PASS",
                "includes the protein hull breathing, so not a membrane convergence signal")
    rep.judge("bilayer thickness (P-P)", float(np.nanmean(thick[-len(thick) // 2:])),
              *THICK_RANGE, unit=" A", soft=True)

    # --- 3. water in the hydrophobic core -------------------------------------------------------
    # PACKMOL traps waters in the tail region. Waters INSIDE the protein are legitimate, so only
    # count those in the lipid phase: in the core slab and not within 6 A of any protein atom.
    zp = phos.positions[:, 2]
    mid = float(zp.mean())
    core = wat_o[np.abs(wat_o.positions[:, 2] - mid) < CORE_HALF_A] if len(wat_o) else wat_o
    n_core = 0
    if len(core):
        pairs = capped_distance(core.positions, prot.positions, max_cutoff=6.0,
                                box=u.trajectory.ts.dimensions, return_distances=False)
        near_protein = set(pairs[:, 0].tolist()) if len(pairs) else set()
        n_core = len(core) - len(near_protein)
    rep.judge("waters in the lipid core", n_core, 0, CORE_WATER_MAX, fmt="{:.0f}", soft=True,
              detail=f"<= {CORE_WATER_MAX}; excludes waters within 6 A of the protein")

    # --- 4. protein ------------------------------------------------------------------------------
    if ca_ok:
        rep.judge("CA RMSD, whole receptor", rmsd_all[-1], 0.0, CA_RMSD_MAX, unit=" A")
        rep.judge("CA RMSD, membrane-embedded", rmsd_tm[-1], 0.0, TM_RMSD_MAX, unit=" A")
    # Vertical drift out of the OPM slab. check_placement.py cannot be reused here: it requires the
    # receptor to be a RIGID translation of receptor.pdb (max residual 0.05 A) and refuses anything
    # that has relaxed. So measure the drift directly, against this run's own first frame.
    rep.judge("OPM registration drift", abs(reg[-1] - reg[0]), 0.0, REGISTRATION_DRIFT_A,
              unit=" A", detail=f"vertical shift vs frame 0 (now {reg[-1]:+.2f} A off the midplane)")
    if len(lig):
        rep.judge(f"ligand RMSD ({args.lig})", lig_rmsd[-1], 0.0, LIG_RMSD_MAX, unit=" A",
                  detail=f"<= {LIG_RMSD_MAX} A; the ligand should stay in the pocket")
    elif args.lig.lower() != "apo":
        if meta_known:
            rep.add(f"ligand ({args.lig})", "absent", "FAIL",
                    "a ligand was requested but the topology has none -- this built apo")
        else:
            rep.add("ligand", "none in topology", "WARN",
                    "no system.json (build predates it), so the intent is unknown; pass --lig apo "
                    "if this is the apo arm")
    bonded = sgd[-1] if sgd else []
    if len(bonded) == 1:
        rep.judge("disulfide SG-SG", bonded[0], *SG_RANGE, unit=" A",
                  detail=f"1 bonded pair among {len(sg)} SG atoms (C142-C219)")
    elif not bonded:
        rep.add("disulfide SG-SG", f"0 of {len(sg)} SG bonded", "FAIL",
                f"no SG-SG pair within {SS_CUTOFF_A} A -- the tleap bond did not take")
    else:
        rep.add("disulfide SG-SG", f"{len(bonded)} pairs", "WARN",
                f"expected 1 (C142-C219) among {len(sg)} SG atoms: "
                + ", ".join(f"{d:.2f}" for d in bonded) + " A")

    # --- 5. lipid ring piercing (topological: minimisation cannot undo it) ----------------------
    if not args.skip_piercing and pdb_p.exists() and Path("check_piercing.py").exists():
        r = subprocess.run([sys.executable, "check_piercing.py", str(pdb_p)],
                           capture_output=True, text=True)
        tail = [ln for ln in r.stdout.strip().splitlines() if ln.strip()]
        pierced = bool(r.returncode) or any("threaded" in ln and "no " not in ln for ln in tail)
        rep.add("ring piercing (final frame)", "none found" if not pierced else "PIERCING FOUND",
                "WARN" if pierced else "PASS", "see the check_piercing.py output above")
        if r.stdout.strip():
            print(r.stdout.strip() + "\n")

    # --- report ----------------------------------------------------------------------------------
    w = max(len(r[0]) for r in rep.rows) + 2
    print(f"{'check'.ljust(w)}{'value'.ljust(22)}{'verdict'.ljust(9)}note")
    print("-" * (w + 22 + 9 + 44))
    for name, value, verdict, detail in rep.rows:
        print(f"{name.ljust(w)}{value[:20].ljust(22)}{verdict.ljust(9)}{detail}")

    rep.values.update({"n_frames": n_frames, "scope": scope,
                       "verdicts": {r[0]: r[2] for r in rep.rows}})
    Path(f"{args.out}.json").write_text(json.dumps(rep.values, indent=2))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t = log["time_ps"] / 1000.0
        f = np.arange(n_frames)
        fig, ax = plt.subplots(2, 3, figsize=(15, 7))
        ax[0][0].plot(t, log["density"]); ax[0][0].set_title("density (g/mL)"); ax[0][0].set_xlabel("ns")
        ax[0][1].plot(t, log["temp"]);    ax[0][1].set_title("temperature (K)"); ax[0][1].set_xlabel("ns")
        ax[0][2].plot(t, log["pe"]);      ax[0][2].set_title("potential energy (kJ/mol)"); ax[0][2].set_xlabel("ns")
        ax[1][0].plot(f, box[:, 0], label="Lx"); ax[1][0].plot(f, box[:, 2], label="Lz")
        ax[1][0].legend(); ax[1][0].set_title("box (A)"); ax[1][0].set_xlabel("frame")
        ax[1][1].plot(f, apl); ax[1][1].set_title(apl_kind + " (A^2)"); ax[1][1].set_xlabel("frame")
        ax[1][2].plot(f, thick, label="P-P thickness")
        if ca_ok:
            ax[1][2].plot(f, rmsd_tm, label="TM CA RMSD")
        ax[1][2].legend(); ax[1][2].set_title("thickness / RMSD (A)"); ax[1][2].set_xlabel("frame")
        fig.tight_layout(); fig.savefig(f"{args.out}.png", dpi=150)
        print(f"\nwrote {args.out}.json and {args.out}.png")
    except ImportError:
        print(f"\nwrote {args.out}.json (matplotlib absent, no figure)")

    # --- verdict ---------------------------------------------------------------------------------
    print()
    if rep.failed:
        print("==> FAIL: something is physically wrong. Do NOT start production; read the table.")
        return 1
    if rep.warned:
        print("==> WARN: nothing is broken, but the system is still relaxing.")
    else:
        print("==> PASS: thermodynamics and geometry are sound.")

    # Tailor the next step to WHICH leg this is. The stage manifest only exists for the restrained
    # ramp, so its absence means this was an unrestrained run -- and telling someone to go and run
    # the pre-production leg they have just finished is worse than saying nothing.
    ramp = stg_p.exists()
    drifting = [n for n in ("area/lipid drift (gross)", "box Lz drift", "box Lx drift",
                            "density drift")
                if rep.verdict_of(n) == "WARN"]
    if ramp:
        print("    A short restrained ramp cannot equilibrate a PACKMOL-built bilayer: the lipids")
        print("    are packed in an artificially ordered, extended conformation and area per lipid")
        print("    takes tens of ns to converge. Run the unrestrained leg and re-check:")
        print("        ./submit.sh preprod")
    elif drifting:
        print(f"    Volume has settled, but {', '.join(drifting)} still moving: the membrane is")
        print("    condensing laterally, which outlasts volume equilibration. Read the SIGN above --")
        print("    consistently falling or rising means keep going, oscillating about a plateau means")
        print("    it is done. To extend, raise ZH_PREPROD_NS in cluster.env and resubmit; or start")
        print("    production and discard the first tens of ns as equilibration.")
    else:
        print("    Ready for production:  ./submit.sh prod")
    return 0


if __name__ == "__main__":
    sys.exit(main())
