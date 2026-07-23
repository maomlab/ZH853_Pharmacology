"""MD system-preparation helpers for the MOR-Gi-scFv16-ZH853 complex.

Structure assessment (missing atoms, termini, titratable/His residues, the Na+ pocket),
component splitting, and membrane-plane geometry from the modeled cholesterol. These inform
the documented prep decisions in docs/METHODS_md_prep.md; the heavy lifting (adding atoms,
membrane packing, parameterization) runs in the analysis/cluster envs.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)

import MDAnalysis as mda  # noqa: E402

# Standard heavy-atom counts per residue (backbone N,CA,C,O + sidechain; no OXT/H).
HEAVY_ATOMS = {
    "ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9, "GLY": 4,
    "HIS": 10, "ILE": 8, "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11, "PRO": 7, "SER": 6,
    "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7,
}

# Chain roles in the deposited model.
CHAIN_ROLES = {
    "A": "G-alpha-i", "B": "G-beta-1", "C": "G-gamma-2", "D": "scFv16",
    "E": "ZH853 (ligand L01)", "R": "MOR (OPRM1) + cholesterol",
}

# Key functional residues (human OPRM1 numbering; verified against the construct).
D250_SODIUM = 116  # D2.50 (ASP116), allosteric Na+ pocket -- ambiguous protonation
DRY_MOTIF = (166, 167, 168)  # D3.49-R3.50-Y3.51 (ASP166-ARG167-TYR168)


@dataclass
class ResidueIssue:
    chain: str
    resid: int
    resname: str
    kind: str  # "incomplete" | "his" | "titratable" | "terminus"
    detail: str


def incomplete_residues(u: mda.Universe) -> list[ResidueIssue]:
    """Protein residues modeled with fewer heavy atoms than the standard count."""
    issues = []
    for res in u.select_atoms("protein").residues:
        expected = HEAVY_ATOMS.get(res.resname)
        if expected is None:
            continue
        n = len(res.atoms.select_atoms("not name H*"))
        if n < expected:
            issues.append(ResidueIssue(res.segid, int(res.resid), res.resname,
                                       "incomplete", f"{n}/{expected} heavy atoms"))
    return issues


def chain_termini(u: mda.Universe) -> list[ResidueIssue]:
    """First/last modeled residue of each protein chain (need capping if truncated)."""
    out = []
    for seg in sorted({a.segid for a in u.select_atoms("protein")}):
        ca = u.select_atoms(f"segid {seg} and protein and name CA")
        if not len(ca):
            continue
        resids = sorted(int(r) for r in ca.resids)
        role = CHAIN_ROLES.get(seg, "?")
        out.append(ResidueIssue(seg, resids[0], _resname(u, seg, resids[0]), "terminus",
                                f"N-terminus of {role} (cap ACE if truncated)"))
        out.append(ResidueIssue(seg, resids[-1], _resname(u, seg, resids[-1]), "terminus",
                                f"C-terminus of {role} (cap NME if truncated)"))
    return out


def histidines(u: mda.Universe, segid: str = "R") -> list[ResidueIssue]:
    """His residues in a chain (each needs a HID/HIE/HIP tautomer assignment)."""
    his = u.select_atoms(f"segid {segid} and resname HIS and name CA")
    return [ResidueIssue(segid, int(a.resid), "HIS", "his", "assign HID/HIE/HIP tautomer")
            for a in his]


def _resname(u: mda.Universe, seg: str, resid: int) -> str:
    sel = u.select_atoms(f"segid {seg} and resid {resid} and name CA")
    return str(sel[0].resname) if len(sel) else "?"


def membrane_frame(u: mda.Universe, segid: str = "R") -> dict[str, object]:
    """Estimate the membrane normal and slab from the receptor TM bundle + cholesterol.

    The membrane normal is the principal axis of the receptor CA cloud (longest axis for a
    7TM bundle). Cholesterol atoms define the bilayer center and thickness along that normal.
    """
    ca = u.select_atoms(f"segid {segid} and protein and name CA")
    coords = ca.positions - ca.positions.mean(axis=0)
    # principal axes via SVD; first singular vector = membrane normal
    _, _, vt = np.linalg.svd(coords, full_matrices=False)
    normal = vt[0]
    clr = u.select_atoms(f"segid {segid} and resname CLR")
    result: dict[str, object] = {"normal": normal.tolist(), "n_cholesterol_atoms": len(clr)}
    if len(clr):
        center = clr.positions.mean(axis=0)
        proj = (clr.positions - center) @ normal
        result["membrane_center"] = center.tolist()
        result["cholesterol_span_along_normal"] = float(proj.max() - proj.min())
    # receptor extent along the normal (membrane-spanning length)
    rproj = coords @ normal
    result["receptor_span_along_normal"] = float(rproj.max() - rproj.min())
    return result
