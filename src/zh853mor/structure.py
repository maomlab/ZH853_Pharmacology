"""Structure-analysis helpers for the MOR-Gi-scFv16-ZH853 cryo-EM model.

Thin, typed wrappers over MDAnalysis that capture the verified Phase-1 analyses:
chain inventory, receptor gaps, disulfides, and the ZH853 binding-pocket contact map.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)

import MDAnalysis as mda  # noqa: E402  (import after warning filter is intentional)

# Ligand ZH853 is deposited as HETATM residue L01 (chain E).
LIGAND_RESNAME = "L01"
# Receptor is chain/segid R in the deposited model.
RECEPTOR_SEGID = "R"

# Canonical human OPRM1 (P35372) orthosteric residues -> Ballesteros-Weinstein label.
# Used to confirm the construct's numbering scheme.
CANONICAL_POCKET: dict[int, tuple[str, str]] = {
    126: ("GLN", "2.60"),
    129: ("ASN", "2.63"),
    149: ("ASP", "3.32"),
    150: ("TYR", "3.33"),
    153: ("MET", "3.36"),
    235: ("LYS", "5.39"),
    295: ("TRP", "6.48"),
    298: ("ILE", "6.51"),
    299: ("HIS", "6.52"),
    302: ("VAL", "6.55"),
    320: ("TRP", "7.35"),
    324: ("ILE", "7.39"),
    327: ("GLY", "7.42"),
    328: ("TYR", "7.43"),
}

# Ballesteros-Weinstein label for pocket-lining residues (human OPRM1 numbering).
# Superset of CANONICAL_POCKET covering the full ZH853 contact shell. A trailing "?"
# marks a loop/less-certain generic-number assignment; region tags used where no clean
# BW number applies. Well-established TM positions carry no "?".
POCKET_BW: dict[int, str] = {
    77: "1.39?",
    124: "2.58",
    126: "2.60",
    129: "2.63",
    130: "2.64",
    135: "ECL1",
    145: "3.28",
    146: "3.29",
    149: "3.32",
    150: "3.33",
    153: "3.36",
    219: "45.50",  # ECL2 disulfide Cys (to C142/3.25)
    220: "45.51",
    221: "45.52",
    223: "45.54",
    231: "ECL2",
    234: "5.38",
    235: "5.39",
    238: "5.42",
    295: "6.48",
    298: "6.51",
    299: "6.52",
    302: "6.55",
    320: "7.35",
    321: "7.36",
    324: "7.39",
    327: "7.42",
    328: "7.43",
}


def bw(resid: int) -> str:
    """Return the BW/region label for a human-OPRM1 residue id, or '-' if unknown."""
    return POCKET_BW.get(resid, "-")


@dataclass(frozen=True)
class Contact:
    """A receptor residue in contact with the ligand."""

    resid: int
    resname: str
    n_atoms: int  # number of contacting heavy atoms within the cutoff


def load(pdb: Path) -> mda.Universe:
    """Load a PDB into an MDAnalysis Universe."""
    return mda.Universe(str(pdb))


def verify_numbering(u: mda.Universe) -> dict[int, bool]:
    """Check each canonical pocket residue against the model.

    Returns a mapping ``resid -> matches_expected_resname``. All-True implies the
    construct uses human OPRM1 (P35372) numbering.
    """
    receptor = u.select_atoms(f"segid {RECEPTOR_SEGID} and protein")
    by_id = {int(a.resid): str(a.resname) for a in receptor}
    return {
        resid: by_id.get(resid) == expected_resname
        for resid, (expected_resname, _bw) in CANONICAL_POCKET.items()
    }


def receptor_gaps(u: mda.Universe) -> list[tuple[int, int]]:
    """Return (before, after) residue-id pairs bracketing internal chain-R gaps."""
    ca = u.select_atoms(f"segid {RECEPTOR_SEGID} and protein and name CA")
    resids = sorted(int(r) for r in ca.resids)
    return [(a, b) for a, b in zip(resids, resids[1:], strict=False) if b - a > 1]


def disulfides(u: mda.Universe, cutoff: float = 2.5) -> list[tuple[int, int, float]]:
    """Return (resid_i, resid_j, distance) for chain-R Cys SG pairs within ``cutoff``."""
    sg = u.select_atoms(f"segid {RECEPTOR_SEGID} and resname CYS and name SG")
    out: list[tuple[int, int, float]] = []
    pos = sg.positions
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            d = float(np.linalg.norm(pos[i] - pos[j]))
            if d < cutoff:
                out.append((int(sg[i].resid), int(sg[j].resid), d))
    return out


def pocket_contacts(u: mda.Universe, cutoff: float = 4.5) -> list[Contact]:
    """Receptor residues with any heavy atom within ``cutoff`` A of ligand ZH853."""
    lig = u.select_atoms(f"resname {LIGAND_RESNAME}")
    near = u.select_atoms(f"(protein and around {cutoff} group lig)", lig=lig)
    counts: dict[tuple[int, str], int] = {}
    for atom in near:
        key = (int(atom.resid), str(atom.resname))
        counts[key] = counts.get(key, 0) + 1
    return sorted(
        (Contact(resid=r, resname=n, n_atoms=c) for (r, n), c in counts.items()),
        key=lambda c: c.resid,
    )
