#!/usr/bin/env python3
"""Build a MolViewSpec (MVS) scene of the ZH853 orthosteric pocket + key interactions.

Emits pocket.mvsj (loaded by render.js via MolStar's loadMvsData): receptor cartoon (faded),
ZH853 + key pocket residues as sticks, and the salt-bridge/H-bond contacts drawn as labeled
distance lines. Coordinates are read from the deposited structure.
"""

from __future__ import annotations

import base64
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import molviewspec as mvs  # noqa: E402
import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PDB = REPO / "intermediate" / "02.01.00_components" / "receptor_ligand.pdb"

# Key pocket residues (auth numbering) -> label.
KEY_RESIDUES = {149: "D3.32", 231: "E231 (ECL2)", 299: "H6.52", 321: "H7.36",
                328: "Y7.43", 129: "N2.63", 150: "Y3.33", 153: "M3.36"}
# Interaction lines: (receptor resid, receptor atom, ligand atom, kind, color)
INTERACTIONS = [
    (149, ("OD1", "OD2"), "N09", "salt bridge", "#e8820c"),
    (328, ("OH",), "N09", "H-bond", "#1f77b4"),
    (231, ("OE1", "OE2"), "N51", "salt bridge", "#e8820c"),
    (129, ("OD1",), "N31", "H-bond", "#1f77b4"),
]


def coord(u, sel: str) -> list[float]:
    ag = u.select_atoms(sel)
    return [float(x) for x in ag.positions.mean(axis=0)]


def nearest_pair(u, resid: int, atoms: tuple[str, ...], lig_atom: str):
    rec = u.select_atoms(f"segid R and resid {resid} and name {' '.join(atoms)}")
    lig = u.select_atoms(f"resname L01 and name {lig_atom}")
    if not len(rec) or not len(lig):
        return None
    d = np.linalg.norm(rec.positions[:, None] - lig.positions[None], axis=2)
    i, j = np.unravel_index(d.argmin(), d.shape)
    return ([float(x) for x in rec.positions[i]], [float(x) for x in lig.positions[j]], float(d[i, j]))


def main() -> int:
    u = mda.Universe(str(PDB))
    data_url = "data:text/plain;base64," + base64.b64encode(PDB.read_bytes()).decode()

    b = mvs.create_builder()
    struct = b.download(url=data_url).parse(format="pdb").model_structure()

    # faded receptor cartoon
    (struct.component(selector="polymer")
        .representation(type="cartoon")
        .color(color="#c8c8c8")
        .opacity(opacity=0.32))

    # ligand ZH853 (chain E) as ball-and-stick, yellow carbons (unlabeled to keep it unobscured)
    lig = struct.component(selector=mvs.ComponentExpression(auth_asym_id="E"))
    lig.representation(type="ball_and_stick").color(color="#f2c811")

    # key pocket residues as sticks + labels
    for resid, lab in KEY_RESIDUES.items():
        comp = struct.component(selector=mvs.ComponentExpression(auth_asym_id="R", auth_seq_id=resid))
        comp.representation(type="ball_and_stick").color(color="#4a90c2")
        comp.label(text=lab)

    # interaction distance lines
    prims = struct.primitives()
    for resid, atoms, lig_atom, kind, color in INTERACTIONS:
        pair = nearest_pair(u, resid, atoms, lig_atom)
        if pair is None:
            continue
        start, end, dist = pair
        prims.distance(start=start, end=end, color=color, radius=0.08,
                       label_template=f"{kind} {dist:.1f} A", label_color=color)

    # focus camera on the ligand pocket
    lig.focus(radius_factor=2.2, direction=(0.2, -0.3, -1.0), up=(0, 1, 0))

    out = HERE / "pocket.mvsj"
    state = b.get_state()
    out.write_text(state.dumps() if hasattr(state, "dumps") else state.json())
    print(f"wrote {out} ({out.stat().st_size // 1024} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
