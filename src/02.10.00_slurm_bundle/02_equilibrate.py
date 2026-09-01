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
        res = atom.residue.name
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


LIPID_RESN = {"POPC", "CHL1", "CLR", "POPE", "POPS", "PSM", "OL", "PC", "PA"}
SOLVENT_RESN = {"WAT", "HOH", "OPC", "NA", "CL", "K", "TIP3", "SOL"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--inpcrd", required=True)
    ap.add_argument("--out", default="eq")
    args = ap.parse_args()

    prmtop = app.AmberPrmtopFile(args.prmtop)
    inpcrd = app.AmberInpcrdFile(args.inpcrd)
    system = prmtop.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds, rigidWater=True,
    )
    force, groups = restraint_force(system, prmtop, inpcrd.positions, "N CA C O", "sc", "lipid")

    integrator = mm.LangevinMiddleIntegrator(TEMP, 1 / unit.picosecond, DT)
    platform = mm.Platform.getPlatformByName("CUDA")
    sim = app.Simulation(prmtop.topology, system, integrator, platform,
                         {"Precision": "mixed"})
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    sim.minimizeEnergy(maxIterations=5000)

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
        sim.step(steps)

    sim.saveState(f"{args.out}.xml")
    with open(f"{args.out}.pdb", "w") as fh:
        app.PDBFile.writeFile(prmtop.topology, sim.context.getState(getPositions=True).getPositions(), fh)
    print(f"Equilibration complete -> {args.out}.xml")


if __name__ == "__main__":
    main()
