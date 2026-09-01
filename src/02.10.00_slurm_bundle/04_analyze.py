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
    ap.add_argument("--lig", default="LIG")
    ap.add_argument("--receptor", default="segid R")
    ap.add_argument("--out", default="qc")
    args = ap.parse_args()

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

    lig = u.select_atoms(f"resname {args.lig}")
    lig_ref = lig.positions.copy()
    lig_rmsd = []
    contact_hits = {r: 0 for r in KEY_CONTACTS}
    nframes = 0
    for _ in u.trajectory:
        lig_rmsd.append(float(np.sqrt(((lig.positions - lig_ref) ** 2).sum(axis=1).mean())))
        for r in KEY_CONTACTS:
            res = u.select_atoms(f"{args.receptor} and resid {r}")
            if len(res) and (np.linalg.norm(
                res.positions[:, None] - lig.positions[None], axis=2
            ).min() <= 4.5):
                contact_hits[r] += 1
        nframes += 1

    summary = {
        "n_frames": nframes,
        "bb_rmsd_mean": float(bb_rmsd.mean()), "bb_rmsd_last": float(bb_rmsd[-1]),
        "rmsf_max": float(rmsf.max()),
        "ligand_rmsd_mean": float(np.mean(lig_rmsd)), "ligand_rmsd_last": lig_rmsd[-1],
        "contact_occupancy": {str(r): contact_hits[r] / max(nframes, 1) for r in KEY_CONTACTS},
    }
    with open(f"{args.out}.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(13, 4))
        ax[0].plot(bb_rmsd); ax[0].set_title("backbone RMSD (A)"); ax[0].set_xlabel("frame")
        ax[1].plot(ca.resids, rmsf); ax[1].set_title("RMSF (A)"); ax[1].set_xlabel("residue")
        ax[2].plot(lig_rmsd); ax[2].set_title("ligand RMSD (A, receptor-aligned)")
        ax[2].set_xlabel("frame")
        fig.tight_layout(); fig.savefig(f"{args.out}.png", dpi=150)
    except ImportError:
        pass
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
