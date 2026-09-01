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
CORE_HALF_A = 8.0              # |z - midplane| defining the hydrophobic core
TM_HALF_A = 15.0               # |z - midplane| defining the membrane-embedded CA set

LIPID_RESN = {"POPC", "CHL1", "CLR", "POPE", "POPS", "PSM"}
WATER_RESN = {"WAT", "HOH", "OPC", "TIP3", "SOL"}


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


def leaflet_counts(u, midplane_z):
    """(upper, lower) lipid counts, split on each lipid's centre of geometry.

    Computed ONCE: lipids do not flip-flop on an equilibration timescale, and re-deriving this
    per frame means walking every residue in the system (mostly water) on every frame.
    """
    lipids = [r for r in u.residues if r.resname in LIPID_RESN]
    if not lipids:
        return 0, 0
    zs = np.array([r.atoms.positions[:, 2].mean() for r in lipids])
    return int((zs >= midplane_z).sum()), int((zs < midplane_z).sum())


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
    ap.add_argument("--skip-piercing", action="store_true")
    args = ap.parse_args()

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
    rep.judge("density drift", abs(halves_drift(dens)), 0.0, DENSITY_DRIFT,
              unit=" g/mL", fmt="{:.4f}", soft=True, detail=f"<= {DENSITY_DRIFT} g/mL over the stage")
    rep.judge("temperature", float(temp.mean()), 310 - TEMP_TOL_K, 310 + TEMP_TOL_K,
              unit=" K", fmt="{:.1f}")
    # Energy drift is scored against its own fluctuation: an absolute kJ/mol threshold is
    # meaningless across system sizes, but drift much larger than sigma means still-relaxing.
    pe_sigma = float(pe.std()) or 1.0
    rep.judge("potential-energy drift", abs(halves_drift(pe)) / pe_sigma, 0.0, 1.0,
              unit=" sigma", fmt="{:.2f}", soft=True, detail="<= 1 sigma of its own fluctuation")

    # --- 2. box, membrane, protein (from the trajectory) ----------------------------------------
    u = mda.Universe(args.top, str(dcd_p))
    n_frames = len(u.trajectory)
    if n_frames < 2:
        print("ERROR: the DCD has fewer than 2 frames; nothing to assess.", file=sys.stderr)
        return 2

    ca = u.select_atoms("name CA")            # lipids/water carry no CA, so this is the protein
    phos = u.select_atoms("resname POPC and name P")
    prot = u.select_atoms("protein") if len(u.select_atoms("protein")) else ca
    wat_o = u.select_atoms(f"resname {' '.join(WATER_RESN)} and name O OW O1")
    sg = u.select_atoms("name SG")

    if not len(phos):
        print("ERROR: no POPC phosphorus atoms found; is this a bilayer?", file=sys.stderr)
        return 2

    ref_ca = mda.Universe(args.receptor).select_atoms("name CA").positions.copy()
    ca_ok = len(ref_ca) == len(ca)
    if not ca_ok:
        rep.add("CA count vs receptor.pdb", f"{len(ca)} vs {len(ref_ca)}", "WARN",
                "cannot compute RMSD against the staged receptor")

    box, apl_g, apl_n, thick, rmsd_all, rmsd_tm, reg, sgd = ([] for _ in range(8))
    tm_mask = None
    n_up, n_lo = leaflet_counts(u, float(phos.positions[:, 2].mean()))
    per_leaflet = max((n_up + n_lo) / 2.0, 1.0)
    print(f"lipids: {n_up + n_lo} ({n_up} upper / {n_lo} lower), {per_leaflet:g} per leaflet\n")
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
            rmsd_all.append(float(np.sqrt(((cap - ref_ca) ** 2).sum(axis=1).mean())))
            rmsd_tm.append(float(np.sqrt(
                ((cap[tm_mask] - ref_ca[tm_mask]) ** 2).sum(axis=1).mean())))
        if len(sg) >= 2:
            d = np.linalg.norm(sg.positions[0] - sg.positions[1])
            sgd.append(float(d) if len(sg) == 2 else np.nan)

    box = np.array(box)
    rep.judge("box Lx drift", abs(halves_drift(box[:, 0])), 0.0, BOX_DRIFT_A,
              unit=" A", soft=True, detail=f"<= {BOX_DRIFT_A} A over the stage")
    rep.judge("box Lz drift", abs(halves_drift(box[:, 2])), 0.0, BOX_DRIFT_A,
              unit=" A", soft=True, detail=f"<= {BOX_DRIFT_A} A over the stage")
    rep.values["box_final_A"] = [round(float(v), 2) for v in box[-1]]

    apl = np.array(apl_n if np.isfinite(apl_n).all() else apl_g)
    apl_kind = "area/lipid" if np.isfinite(apl_n).all() else "area/lipid (gross, no protein subtraction)"
    rep.judge(apl_kind, float(apl[-len(apl) // 2:].mean()), *APL_RANGE, unit=" A^2", soft=True)
    rep.judge("area/lipid drift", abs(halves_drift(apl)), 0.0, 0.5, unit=" A^2", soft=True,
              detail="<= 0.5 A^2 -- the slowest observable; expect this to fail on a short ramp")
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
    if len(sg) == 2:
        rep.judge("disulfide SG-SG", sgd[-1], *SG_RANGE, unit=" A")
    else:
        rep.add("disulfide SG-SG", f"{len(sg)} SG atoms", "WARN",
                "expected exactly 2 (C142/C219); cannot verify the tleap bond")

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
    print("    A short restrained ramp cannot equilibrate a PACKMOL-built bilayer: the lipids are")
    print("    packed in an artificially ordered, extended conformation and area per lipid takes")
    print("    tens of ns to converge. Run the unrestrained pre-production leg and re-check:")
    print("        ./submit.sh preprod    # then: python check_equilibration.py --eq preprod ...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
