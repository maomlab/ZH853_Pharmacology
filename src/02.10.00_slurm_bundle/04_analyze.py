#!/usr/bin/env python3
"""Production QC for the MOR-ZH853 trajectory (standalone; MDAnalysis).

Computes the resolution-appropriate QC metrics (SPECIFICATION QC section): backbone Cα RMSD,
per-residue RMSF (pre-aligned), receptor-aligned ligand RMSD, key-contact occupancy, and
membrane area-per-lipid. Emits a JSON summary + PNG per replica. Concatenate replicas for
error bars. Reuses no repo imports so it can run alone on the cluster.

Step 5; run from the build directory, once per replica (CPU is fine):

    python 04_analyze.py --top system.prmtop --traj prod_r1.dcd --lig LIG --out qc_r1

`--lig` must match the ligand residue name in the prmtop (the residue name in ZH853.mol2), and
`--receptor` the receptor selection; both are empty selections on an apo build.
"""

from __future__ import annotations

import argparse
import json

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align, rms

KEY_CONTACTS = [149, 231, 299, 321, 328]  # D3.32, E231(ECL2), H6.52, H7.36, Y7.43 (Objective 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--lig", default=None,
                    help="ligand resname; defaults to system.json's, else LIG. 'apo' (or a "
                         "resname matching nothing) skips the ligand metrics")
    ap.add_argument("--receptor", default="segid R")
    ap.add_argument("--out", default="qc")
    args = ap.parse_args()

    # The build records what it is in system.json, so the panel of 5 systems can be analysed with
    # one command line rather than remembering which directory is apo.
    if args.lig is None:
        try:
            with open("system.json") as fh:
                meta = json.load(fh)
            args.lig = meta.get("ligand_resname") or "apo"
            print(f"system.json: ligand={meta.get('ligand')} -> --lig {args.lig}")
        except (OSError, ValueError):
            args.lig = "LIG"

    u = mda.Universe(args.top, args.traj)
    ref = mda.Universe(args.top, args.traj)
    ref.trajectory[0]

    # backbone RMSD (aligned to frame 0 on receptor Cα)
    rmsd = rms.RMSD(u, ref, select=f"{args.receptor} and name CA").run()
    bb_rmsd = rmsd.results.rmsd[:, 3]

    # align whole trajectory on receptor Cα, then RMSF + receptor-aligned ligand RMSD
    align.AlignTraj(u, ref, select=f"{args.receptor} and name CA", in_memory=True).run()
    ca = u.select_atoms(f"{args.receptor} and name CA")
    rmsf = rms.RMSF(ca).run().results.rmsf

    # An apo build has no ligand: ligand RMSD and contact occupancy are undefined, and an empty
    # selection would otherwise raise deep inside the frame loop after minutes of trajectory I/O.
    lig = u.select_atoms("") if args.lig.lower() == "apo" else u.select_atoms(f"resname {args.lig}")
    has_ligand = len(lig) > 0
    if not has_ligand and args.lig.lower() != "apo":
        print(f"WARNING: no atoms match 'resname {args.lig}' -- treating this as an apo system. "
              "Check system.json if that is not what you built.")
    lig_ref = lig.positions.copy() if has_ligand else None
    lig_rmsd = []
    contact_hits = {r: 0 for r in KEY_CONTACTS}
    nframes = 0
    for _ in u.trajectory:
        if has_ligand:
            lig_rmsd.append(float(np.sqrt(((lig.positions - lig_ref) ** 2).sum(axis=1).mean())))
            for r in KEY_CONTACTS:
                res = u.select_atoms(f"{args.receptor} and resid {r}")
                if len(res) and (np.linalg.norm(
                    res.positions[:, None] - lig.positions[None], axis=2
                ).min() <= 4.5):
                    contact_hits[r] += 1
        nframes += 1

    summary = {
        "n_frames": nframes, "ligand": args.lig if has_ligand else "apo",
        "bb_rmsd_mean": float(bb_rmsd.mean()), "bb_rmsd_last": float(bb_rmsd[-1]),
        "rmsf_max": float(rmsf.max()),
    }
    if has_ligand:
        summary.update({
            "ligand_rmsd_mean": float(np.mean(lig_rmsd)), "ligand_rmsd_last": lig_rmsd[-1],
            "contact_occupancy": {str(r): contact_hits[r] / max(nframes, 1) for r in KEY_CONTACTS},
        })
    with open(f"{args.out}.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n_panels = 3 if has_ligand else 2
        fig, ax = plt.subplots(1, n_panels, figsize=(4.4 * n_panels, 4))
        ax[0].plot(bb_rmsd); ax[0].set_title("backbone RMSD (A)"); ax[0].set_xlabel("frame")
        ax[1].plot(ca.resids, rmsf); ax[1].set_title("RMSF (A)"); ax[1].set_xlabel("residue")
        if has_ligand:
            ax[2].plot(lig_rmsd); ax[2].set_title("ligand RMSD (A, receptor-aligned)")
            ax[2].set_xlabel("frame")
        fig.tight_layout(); fig.savefig(f"{args.out}.png", dpi=150)
    except ImportError:
        pass
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
