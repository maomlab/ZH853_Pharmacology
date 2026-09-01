#!/usr/bin/env python3
"""6-stage restrained equilibration of the MOR-ZH853 membrane system (OpenMM).

Adapts the CHARMM-GUI membrane-equilibration schedule: strong positional restraints on protein
and lipid heavy atoms, released over six stages (NVT then NPT), before free production. Input is
an Amber system (prmtop/rst7) from PACKMOL-Memgen + tleap.

Normally submitted as step 3 via `sbatch submit_equilibrate.sbatch` from the build directory
(intermediate/02.10.00_build/<D250>_<timestamp>/), not invoked by hand. To run it directly, do so
from that same directory and keep the `system_eq` output basename -- submit_production.sbatch
reads back `${SYS}_eq.xml`, so `--out eq` would strand step 4:

    python 02_equilibrate.py --prmtop system.prmtop --inpcrd system.rst7 --out system_eq
"""

from __future__ import annotations

import argparse
import json

import openmm as mm
from openmm import app, unit

# (backbone_k, sidechain_k, lipid_k) in kJ/mol/nm^2, ensemble, steps  -- CHARMM-GUI-style ramp.
STAGES = [
    (4000.0, 2000.0, 1000.0, "NVT", 125000),
    (2000.0, 1000.0,  400.0, "NVT", 125000),
    (1000.0,  500.0,  400.0, "NPT", 125000),
    ( 500.0,  200.0,  200.0, "NPT", 250000),
    ( 100.0,   50.0,   40.0, "NPT", 250000),
    (  10.0,    0.0,    0.0, "NPT", 250000),
]
TEMP = 310 * unit.kelvin
DT = 2 * unit.femtoseconds  # 2 fs during equilibration (HMR/4 fs switched on in production)


def restraint_force(system, prmtop, positions, sel_backbone, sel_sidechain, sel_lipid):
    """Add a CustomExternalForce holding selected atoms; return (force, {group: [indices]})."""
    force = mm.CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addPerParticleParameter("k")
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    groups = {"backbone": [], "sidechain": [], "lipid": []}
    for atom in prmtop.topology.atoms():
        if atom.element is None or atom.element == app.element.hydrogen:
            continue
        i = atom.index
        res = atom.residue.name.strip().upper()
        if res in LIPID_RESN:
            grp = "lipid"
        elif atom.name in ("N", "CA", "C", "O"):
            grp = "backbone"
        elif res not in SOLVENT_RESN:
            grp = "sidechain"
        else:
            continue
        x, y, z = positions[i].value_in_unit(unit.nanometer)
        force.addParticle(i, [0.0, x, y, z])
        groups[grp].append(force.getNumParticles() - 1)
    system.addForce(force)
    return force, groups


# Lipid21 is MODULAR: one POPC is three residues -- PC (headgroup) + PA (palmitoyl) + OL (oleoyl)
# -- and cholesterol is CHL, not CHL1 (CHL1 is what packmol-memgen calls it on the way IN). A
# lipid residue missing from this set is not an error: it silently falls through to the
# `sidechain` group and follows the wrong restraint schedule. Hence the group census printed
# below, so a naming surprise is visible in the job log rather than inferred months later.
LIPID_RESN = {
    # whole-lipid names (as packed, and other force fields)
    "POPC", "POPE", "POPS", "PSM", "CHL1", "CLR", "CHOL",
    # Lipid21 headgroups
    "PC", "PE", "PS", "PGR", "PH-", "PSA",
    # Lipid21 acyl chains
    "PA", "OL", "ST", "MY", "LAL", "DHA", "SA",
    # Lipid21 cholesterol
    "CHL",
}
SOLVENT_RESN = {"WAT", "HOH", "OPC", "TIP3", "SOL", "NA", "CL", "K", "NA+", "CL-", "K+"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--inpcrd", required=True)
    ap.add_argument("--out", default="eq")
    ap.add_argument("--report-ps", type=float, default=5.0,
                    help="state-log cadence in ps (energy/T/density/volume)")
    ap.add_argument("--dcd-ps", type=float, default=20.0,
                    help="trajectory cadence in ps; frames carry the box, so this sets the "
                         "resolution of the area-per-lipid and thickness traces")
    args = ap.parse_args()

    prmtop = app.AmberPrmtopFile(args.prmtop)
    inpcrd = app.AmberInpcrdFile(args.inpcrd)
    system = prmtop.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds, rigidWater=True,
    )
    force, groups = restraint_force(system, prmtop, inpcrd.positions, "N CA C O", "sc", "lipid")
    # If a lipid residue name is not in LIPID_RESN it lands in `sidechain` and gets the wrong
    # restraint ramp, silently. Print the census (and the residue names behind it) so that shows up.
    resnames = sorted({a.residue.name.strip().upper() for a in prmtop.topology.atoms()})
    print("restraint groups: " + ", ".join(f"{g}={len(v)}" for g, v in groups.items()), flush=True)
    print("  lipid residues seen: "
          + (", ".join(r for r in resnames if r in LIPID_RESN) or "NONE"), flush=True)
    unknown = [r for r in resnames
               if r not in LIPID_RESN and r not in SOLVENT_RESN and len(r) != 3]
    if unknown:
        print(f"  NOTE: unrecognised residue names (check they are meant to be restrained as "
              f"protein): {', '.join(unknown)}", flush=True)

    integrator = mm.LangevinMiddleIntegrator(TEMP, 1 / unit.picosecond, DT)
    platform = mm.Platform.getPlatformByName("CUDA")
    sim = app.Simulation(prmtop.topology, system, integrator, platform,
                         {"Precision": "mixed"})
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    sim.minimizeEnergy(maxIterations=5000)

    # Equilibration is only assessable if it leaves a trace. The state log carries the
    # thermodynamic plateau (energy/T/density/volume); the DCD carries the periodic box per frame,
    # which is what area-per-lipid and bilayer thickness are computed from. check_equilibration.py
    # reads both, and the stage manifest below so it can restrict "has it converged" to the final,
    # least-restrained stage rather than averaging over the restraint ramp.
    log_every = max(1, int((args.report_ps * unit.picoseconds) / DT))
    dcd_every = max(1, int((args.dcd_ps * unit.picoseconds) / DT))
    sim.reporters.append(app.StateDataReporter(
        f"{args.out}.log", log_every, step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, temperature=True, density=True, volume=True, speed=True,
    ))
    sim.reporters.append(app.DCDReporter(f"{args.out}.dcd", dcd_every))

    manifest = {"dt_fs": DT.value_in_unit(unit.femtoseconds), "log_every_steps": log_every,
                "dcd_every_steps": dcd_every, "stages": []}
    barostat = None
    for n, (kb, ks, kl, ensemble, steps) in enumerate(STAGES, 1):
        for grp, k in (("backbone", kb), ("sidechain", ks), ("lipid", kl)):
            for idx in groups[grp]:
                params = force.getParticleParameters(idx)
                force.setParticleParameters(idx, params[0], [k, *params[1][1:]])
        force.updateParametersInContext(sim.context)
        if ensemble == "NPT" and barostat is None:
            barostat = mm.MonteCarloMembraneBarostat(
                1 * unit.bar, 0 * unit.bar * unit.nanometer, TEMP,
                mm.MonteCarloMembraneBarostat.XYIsotropic,
                mm.MonteCarloMembraneBarostat.ZFree, 100,
            )
            system.addForce(barostat)
            sim.context.reinitialize(preserveState=True)
        sim.context.setVelocitiesToTemperature(TEMP)
        print(f"[stage {n}/6] {ensemble} k=({kb},{ks},{kl}) steps={steps}", flush=True)
        first = sim.currentStep + 1
        sim.step(steps)
        manifest["stages"].append({
            "stage": n, "ensemble": ensemble, "k_backbone": kb, "k_sidechain": ks, "k_lipid": kl,
            "first_step": first, "last_step": sim.currentStep,
        })

    with open(f"{args.out}_stages.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    sim.saveState(f"{args.out}.xml")
    with open(f"{args.out}.pdb", "w") as fh:
        app.PDBFile.writeFile(prmtop.topology, sim.context.getState(getPositions=True).getPositions(), fh)
    print(f"Equilibration complete -> {args.out}.xml")
    print(f"QC next: python check_equilibration.py --top {args.prmtop} --eq {args.out} "
          f"--receptor receptor.pdb")


if __name__ == "__main__":
    main()
