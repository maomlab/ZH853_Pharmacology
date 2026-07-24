#!/usr/bin/env python3
"""Compute a canonical-GPCR camera frame for the overview render.

Orients the MOR-Gi complex with the membrane normal vertical, the extracellular face UP and
the intracellular G-protein interface DOWN. Writes overview_camera.json (center, up, direction,
radius) consumed by render.js.

Membrane normal = the receptor (chain R) TM-bundle principal axis; its sign is set so it points
away from the G-protein (chains A/B/C) = toward extracellular. The horizontal viewing direction
is the least-spread in-plane axis (widest silhouette).
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PDB = REPO / "data" / "mor_gi_scfv16_refine-coot-22_real_space_refined_169_edit.pdb"


def main() -> int:
    u = mda.Universe(str(PDB))
    rec = u.select_atoms("segid R and protein and name CA")
    gp = u.select_atoms("segid A B C and name CA")
    lig = u.select_atoms("resname L01")
    allca = u.select_atoms("protein and name CA")

    # membrane normal = receptor principal axis
    r = rec.positions - rec.positions.mean(axis=0)
    normal = np.linalg.svd(r, full_matrices=False)[2][0]
    normal /= np.linalg.norm(normal)

    # sign so 'up' points away from the G protein (extracellular)
    to_gprotein = gp.positions.mean(axis=0) - rec.positions.mean(axis=0)
    up = -np.sign(np.dot(normal, to_gprotein)) * normal

    # Viewing direction: look along the horizontal offset between the receptor and G protein
    # so the two stack vertically in the image (receptor above, G protein below).
    rec_c = rec.positions.mean(axis=0)
    gp_c = gp.positions.mean(axis=0)
    offset = gp_c - rec_c
    view_dir = offset - np.dot(offset, up) * up  # remove the vertical part
    if np.linalg.norm(view_dir) < 1e-3:  # already well-stacked -> any horizontal axis
        c0 = allca.positions - allca.positions.mean(axis=0)
        horiz = c0 - np.outer(c0 @ up, up)
        view_dir = np.linalg.svd(horiz, full_matrices=False)[2][1]
        view_dir = view_dir - np.dot(view_dir, up) * up
    view_dir /= np.linalg.norm(view_dir)

    # center the receptor/G-protein stack; size to include the whole complex (scFv16 included)
    center = 0.5 * (rec_c + gp_c)
    c = allca.positions - center
    radius = float(np.linalg.norm(c, axis=1).max())

    # sanity: extracellular (ligand) should project higher than the G protein along up
    lig_h = float(np.dot(lig.positions.mean(axis=0) - center, up))
    gp_h = float(np.dot(gp.positions.mean(axis=0) - center, up))
    assert lig_h > gp_h, "orientation check failed: ligand not above G protein"

    cam = {
        "center": [float(x) for x in center],
        "up": [float(x) for x in up],
        "direction": [float(x) for x in view_dir],
        "radius": radius,
    }
    (HERE / "overview_camera.json").write_text(json.dumps(cam, indent=2))
    print(f"up={cam['up']}  view={cam['direction']}  radius={radius:.1f}")
    print(f"ligand height along up = {lig_h:.1f} > G-protein {gp_h:.1f}  (extracellular is up)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
