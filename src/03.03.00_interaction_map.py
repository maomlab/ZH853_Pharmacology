#!/usr/bin/env python3
"""PoseView-style 2D interaction map for ZH853 (Objective 1 / manuscript).

Draws the ZH853 2D skeleton (atom-colored) with each contacting receptor residue placed in the
direction of its real contact atom and connected by an interaction-typed line (salt bridge,
H-bond, pi-stack, cation-pi, hydrophobic). Vector PDF output for the manuscript.

Run: ``python src/03.03.00_interaction_map.py``  (or ``make interaction-map``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import rdCoordGen  # noqa: E402
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate  # noqa: E402

from zh853mor import chem, comparators, interactions, structure  # noqa: E402
from zh853mor.paths import PRODUCT, ensure_dir  # noqa: E402

RDLogger.DisableLog("rdApp.*")

ELEM_COLOR = {"N": "#2060d0", "O": "#d02020", "C": "#303030", "S": "#c0a000"}
RESTYPE_COLOR = {  # residue node fill by chemical class
    "acidic": "#f4a3a3", "basic": "#a3c4f4", "polar": "#b6e3b6",
    "aromatic": "#d9b3e6", "hydrophobic": "#e0ddd0",
}
INTERACTION_STYLE = {  # (color, linestyle, linewidth, label)
    "ionic": ("#e8820c", (0, (4, 2)), 2.4, "salt bridge"),
    "hbond": ("#1f77b4", (0, (2, 2)), 1.8, "H-bond"),
    "aromatic": ("#2ca02c", (0, (5, 2)), 1.8, "$\\pi$-stack"),
    "cation_pi": ("#7b3fa0", (0, (1, 1)), 1.8, "cation-$\\pi$"),
    "hydrophobic": ("#9a9a9a", (0, (1, 3)), 1.2, "hydrophobic"),
    "contact": ("#cccccc", (0, (1, 4)), 0.9, "contact"),
}
ACIDIC, BASIC = {"ASP", "GLU"}, {"LYS", "ARG", "HIS"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
AROMATIC = {"PHE", "TRP", "TYR", "HIS"}


def restype(resname: str) -> str:
    if resname in ACIDIC:
        return "acidic"
    if resname in BASIC:
        return "basic"
    if resname in AROMATIC:
        return "aromatic"
    if resname in POLAR:
        return "polar"
    return "hydrophobic"


def build_ligand_2d():
    """Return (coords2d Nx2, elements, bonds, pdb_atom_names, lig_positions3d)."""
    cx = comparators.load_complex("ZH853")
    lig = cx.ligand
    names = [a.name for a in lig.atoms]
    pos3d = lig.positions
    block = ["HETATM%5d %-4s L01 E   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
             % (i + 1, a.name, *a.position, a.name[0]) for i, a in enumerate(lig.atoms)]
    pdbmol = Chem.MolFromPDBBlock("\n".join(block) + "\nEND\n", sanitize=False, proximityBonding=True)
    mol = AssignBondOrdersFromTemplate(Chem.MolFromSmiles(chem.ANALOGS["ZH853"][0]), pdbmol)
    rdCoordGen.AddCoords(mol)
    conf = mol.GetConformer()
    coords2d = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                         for i in range(mol.GetNumAtoms())])
    elements = [a.GetSymbol() for a in mol.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()) for b in mol.GetBonds()]
    return coords2d, elements, bonds, names, pos3d


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    coords2d, elements, bonds, names, pos3d = build_ligand_2d()
    coords2d = coords2d - coords2d.mean(axis=0)
    scale = 1.0 / (np.abs(coords2d).max() + 1e-9)
    xy = coords2d * scale  # normalized to ~[-1,1]

    cx = comparators.load_complex("ZH853")
    fp = interactions.fingerprint(cx.receptor, cx.ligand, cx.offset_to_human)
    rec = cx.receptor

    fig, ax = plt.subplots(figsize=(9, 9))
    # ligand bonds
    for i, j, order in bonds:
        p, q = xy[i], xy[j]
        offs = [0.0]
        if order == 2:
            offs = [-0.012, 0.012]
        d = q - p
        perp = np.array([-d[1], d[0]])
        perp = perp / (np.linalg.norm(perp) + 1e-9)
        for o in offs:
            ax.plot([p[0] + perp[0] * o, q[0] + perp[0] * o],
                    [p[1] + perp[1] * o, q[1] + perp[1] * o], "-", color="#303030", lw=1.6, zorder=3)
    # heteroatom labels
    for i, el in enumerate(elements):
        if el != "C":
            ax.text(*xy[i], el, ha="center", va="center", fontsize=9, fontweight="bold",
                    color=ELEM_COLOR.get(el, "#303030"),
                    bbox={"boxstyle": "circle,pad=0.05", "fc": "white", "ec": "none"}, zorder=4)

    center = xy.mean(axis=0)
    placed: list[np.ndarray] = []
    for hid in sorted(fp):
        r = fp[hid]
        # nearest ligand atom (3D) -> its 2D position gives the contact direction
        rc = rec.select_atoms(f"resid {hid - cx.offset_to_human}").positions
        dists = np.linalg.norm(pos3d[:, None] - rc[None], axis=2).min(axis=1)
        atom_idx = int(dists.argmin())
        anchor = xy[atom_idx]
        direction = anchor - center
        direction = direction / (np.linalg.norm(direction) + 1e-9)
        node = anchor + direction * 0.55
        # de-overlap: nudge if too close to an existing node
        for _ in range(24):
            if all(np.linalg.norm(node - p) > 0.30 for p in placed):
                break
            node = node + direction * 0.06 + np.array([direction[1], -direction[0]]) * 0.05
        placed.append(node)

        itypes = sorted(r.interactions, key=lambda t: -interactions.__dict__.get("ARO_CUT", 0))
        primary = next((t for t in ("ionic", "cation_pi", "aromatic", "hbond", "hydrophobic")
                        if t in r.interactions), "contact")
        color, ls, lw, _ = INTERACTION_STYLE[primary]
        ax.plot([anchor[0], node[0]], [anchor[1], node[1]], linestyle=ls, color=color, lw=lw, zorder=2)

        fc = RESTYPE_COLOR[restype(r.resname)]
        label = f"{r.resname.title()}{hid}\n{structure.bw(hid)}"
        box = FancyBboxPatch((node[0] - 0.13, node[1] - 0.065), 0.26, 0.13,
                             boxstyle="round,pad=0.02,rounding_size=0.04",
                             fc=fc, ec="#555555", lw=0.8, zorder=5)
        ax.add_patch(box)
        ax.text(node[0], node[1], label, ha="center", va="center", fontsize=6.5, zorder=6)

    # legends
    inter_handles = [plt.Line2D([], [], color=c, linestyle=ls, lw=lw, label=lab)
                     for (c, ls, lw, lab) in INTERACTION_STYLE.values()]
    res_handles = [plt.Line2D([], [], marker="s", ls="", mfc=c, mec="#555", label=k)
                   for k, c in RESTYPE_COLOR.items()]
    leg1 = ax.legend(handles=inter_handles, loc="upper left", fontsize=8, title="interaction",
                     framealpha=0.9)
    ax.add_artist(leg1)
    ax.legend(handles=res_handles, loc="upper right", fontsize=8, title="residue type", framealpha=0.9)

    ax.set_title("ZH853 – MOR interaction map (from the 3.5 Å cryo-EM pose)", fontsize=12)
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.9, 1.9)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()

    ensure_dir(PRODUCT)
    pdf = PRODUCT / f"03.03.00_interaction_map_{today}.pdf"
    png = PRODUCT / f"03.03.00_interaction_map_{today}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=200)
    fig.savefig(ensure_dir(PRODUCT / "manuscript" / "figures") / "fig5_interaction_map.pdf")
    plt.close(fig)
    print(f"Wrote {pdf.name}, {png.name}, and manuscript fig5.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
