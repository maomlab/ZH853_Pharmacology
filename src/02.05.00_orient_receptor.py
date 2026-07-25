#!/usr/bin/env python3
"""Orient the rebuilt receptor to the membrane normal for PACKMOL-Memgen (Phase 2).

PACKMOL-Memgen's `--preoriented` assumes the protein's membrane normal is along +z and the
bilayer centre is at z=0. The deposited/rebuilt receptor is still in the cryo-EM frame, so we
rotate the TM-bundle principal axis onto z and centre the CA cloud at the origin.

This is a deterministic default; for production, PPM/OPM orientation (which uses the actual
membrane-insertion energetics) is preferred -- see the note in 01_build_system.sh.

Run (local analysis env): python src/02.05.00_orient_receptor.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import numpy as np  # noqa: E402

from zh853mor import paths, structure  # noqa: E402

IN_PDB = paths.INTERMEDIATE / "02.03.00_receptor" / "receptorR_fixed_heavy.pdb"


def rotation_onto_z(normal: np.ndarray) -> np.ndarray:
    """Rotation matrix mapping the unit vector ``normal`` onto +z (Rodrigues)."""
    n = normal / np.linalg.norm(normal)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(n, z)
    s = np.linalg.norm(v)
    c = float(np.dot(n, z))
    if s < 1e-8:  # already aligned (c=+1) or anti-aligned (c=-1)
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def main() -> int:
    if not IN_PDB.exists():
        print(f"ERROR: {IN_PDB} not found -- run `make prep-receptor` first.", file=sys.stderr)
        return 1
    u = structure.load(IN_PDB)
    ca = u.select_atoms("name CA")

    # membrane normal = TM-bundle principal axis
    x = ca.positions - ca.positions.mean(axis=0)
    normal = np.linalg.svd(x, full_matrices=False)[2][0]
    r = rotation_onto_z(normal)

    # centre CA cloud at origin, then rotate the whole receptor so normal -> +z
    all_atoms = u.atoms
    all_atoms.positions = (all_atoms.positions - ca.center_of_geometry()) @ r.T

    out_dir = paths.ensure_dir(paths.INTERMEDIATE / "02.05.00_oriented")
    out = out_dir / "receptorR_oriented.pdb"
    all_atoms.write(str(out))

    # report: TM span should now lie along z
    zspan = float(u.select_atoms("name CA").positions[:, 2].ptp())
    print(f"Oriented normal -> +z; wrote {out}")
    print(f"CA z-extent (membrane-spanning length) = {zspan:.1f} A")
    print("NOTE: --preoriented is now valid for this file. For production, prefer PPM/OPM.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
