#!/usr/bin/env python3
"""Production MD for the MOR-ZH853 membrane system (OpenMM).

Restraint-free NPT production with a semi-isotropic membrane barostat and hydrogen-mass
repartitioning (4 fs). Resumes from the equilibrated state; writes DCD + checkpoints + a
state-data log. Run one instance per replica (distinct --seed / --out).

Normally submitted as step 4 via `sbatch submit_production.sbatch` (a 1-3 job array, one replica
per task) from the build directory, after step 3 has written the equilibrated state. To run one
replica by hand, from that same directory:

    python 03_production.py --prmtop system.prmtop --state system_eq.xml --out prod_r1 \
                            --ns 500 --seed 1
"""

from __future__ import annotations

import argparse

import openmm as mm
from openmm import app, unit

TEMP = 310 * unit.kelvin
DT = 4 * unit.femtoseconds  # enabled by HMR (hydrogenMass below)
REPORT_PS = 100  # trajectory/log cadence in ps


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--state", required=True, help="equilibrated state .xml")
    ap.add_argument("--out", default="prod")
    ap.add_argument("--ns", type=float, default=500.0, help="production length (ns)")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    prmtop = app.AmberPrmtopFile(args.prmtop)
    system = prmtop.createSystem(
        nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds, rigidWater=True,
        hydrogenMass=4 * unit.amu,  # HMR -> stable 4 fs
    )
    system.addForce(mm.MonteCarloMembraneBarostat(
        1 * unit.bar, 0 * unit.bar * unit.nanometer, TEMP,
        mm.MonteCarloMembraneBarostat.XYIsotropic,
        mm.MonteCarloMembraneBarostat.ZFree, 100,
    ))
    integrator = mm.LangevinMiddleIntegrator(TEMP, 1 / unit.picosecond, DT)
    integrator.setRandomNumberSeed(args.seed)
    platform = mm.Platform.getPlatformByName("CUDA")
    sim = app.Simulation(prmtop.topology, system, integrator, platform, {"Precision": "mixed"})
    sim.loadState(args.state)
    sim.context.setVelocitiesToTemperature(TEMP, args.seed)

    steps = int((args.ns * unit.nanoseconds) / DT)
    report = int((REPORT_PS * unit.picoseconds) / DT)
    sim.reporters.append(app.DCDReporter(f"{args.out}.dcd", report))
    sim.reporters.append(app.CheckpointReporter(f"{args.out}.chk", report * 10))
    sim.reporters.append(app.StateDataReporter(
        f"{args.out}.log", report, step=True, time=True, potentialEnergy=True,
        temperature=True, density=True, volume=True, speed=True,
    ))
    print(f"Production: {args.ns} ns, {steps} steps @ {DT}, seed {args.seed}", flush=True)
    sim.step(steps)
    sim.saveState(f"{args.out}_final.xml")
    # Final frame as PDB: check_equilibration.py hands it to check_piercing.py, which needs
    # coordinates rather than a serialised State.
    with open(f"{args.out}_final.pdb", "w") as fh:
        app.PDBFile.writeFile(prmtop.topology,
                              sim.context.getState(getPositions=True).getPositions(), fh)
    print(f"Done -> {args.out}.dcd")


if __name__ == "__main__":
    main()
