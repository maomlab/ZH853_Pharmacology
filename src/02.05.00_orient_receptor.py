#!/usr/bin/env python3
"""Orient the rebuilt receptor to the membrane normal for PACKMOL-Memgen (Phase 2).

PACKMOL-Memgen's `--preoriented` assumes the membrane normal is along +z and the bilayer centre
is at z=0. We rotate the TM-bundle principal axis onto z and centre on the modeled cholesterols.

IMPORTANT: this is a quick structural proxy for a first-pass local build. For production the
membrane should be placed with the OPM/PPM transfer-energy method (the community standard; hydrophobic
thickness ~32 A for MOR, OPM 4DKL 32.0+/-1.0 A), NOT from the 3 cholesterols alone -- those bind at
site-specific motifs and fix the midplane only to ~2 A. src/02.06.00 validates this placement against
the Trp/Tyr aromatic girdle (~30 A) and experimental POPC thickness; see docs/METHODS_md_prep.md.

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

    # membrane normal = TM-bundle principal axis (rebuilt receptor is in the deposited frame)
    x = ca.positions - ca.positions.mean(axis=0)
    normal = np.linalg.svd(x, full_matrices=False)[2][0]
    r = rotation_onto_z(normal)

    # membrane MIDPLANE from the 3 modeled cholesterols (84 atoms) in the deposited structure (same frame),
    # NOT the Cα centroid -- the receptor is asymmetric along z, so the Cα centroid is offset from
    # the bilayer centre and mis-places packmol-memgen's head/tail planes.
    orig = structure.load(paths.CRYOEM_PDB)
    clr = orig.select_atoms("segid R and resname CLR")
    if not len(clr):
        print("ERROR: no cholesterol (CLR) in the deposited model to define the membrane centre.",
              file=sys.stderr)
        return 1

    # Rotate normal -> z, then centre: xy on the receptor, z=0 on the cholesterol (membrane) midplane.
    ca_xy = (ca.positions @ r.T)[:, :2].mean(axis=0)
    clr_z = float((clr.positions @ r.T)[:, 2].mean())
    shift = np.array([ca_xy[0], ca_xy[1], clr_z])

    def orient(atomgroup) -> None:
        atomgroup.positions = atomgroup.positions @ r.T - shift

    orient(u.atoms)

    out_dir = paths.ensure_dir(paths.INTERMEDIATE / "02.05.00_oriented")
    out = out_dir / "receptorR_oriented.pdb"
    u.atoms.write(str(out))

    # Also emit the oriented complex (receptor + cholesterol + ligand) for the membrane-placement
    # figure -- same transform. Requires prep-assess (02.01.00) to have written the components.
    reclig = paths.INTERMEDIATE / "02.01.00_components" / "receptor_ligand.pdb"
    if reclig.exists():
        cx = structure.load(reclig)
        orient(cx.atoms)
        cx.atoms.write(str(out_dir / "complex_oriented.pdb"))

    caz = u.select_atoms("name CA").positions[:, 2]
    clr_span = float(np.ptp((clr.positions @ r.T)[:, 2]))
    print(f"Oriented normal -> +z; membrane midplane (cholesterol) at z=0; wrote {out}")
    print(f"CA z-extent = {float(np.ptp(caz)):.1f} A; membrane extends {caz.min():.1f}..{caz.max():.1f} "
          f"around z=0 (cholesterol z-span {clr_span:.1f} A)")
    print("NOTE: --preoriented is now valid. For production, PPM/OPM is the rigorous alternative.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
