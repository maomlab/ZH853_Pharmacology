"""Load MOR comparator complexes with a consistent receptor/ligand selection.

Each comparator PDB uses different chain IDs, ligand codes, and residue numbering
(human vs mouse OPRM1). This module isolates the orthosteric receptor chain and the
bound agonist, and detects the offset that maps the receptor onto human OPRM1 (P35372)
numbering so all downstream analyses share one coordinate frame.

Ligand-selection strategy (verified against RCSB, see docs/references.md):
  - Peptide agonists (DAMGO, endomorphin-1) are split across a short protein chain plus
    non-standard-residue HETATMs (DAL/MEA/ETA for DAMGO, NH2 cap for endomorphin-1).
  - Where a ligand appears in two receptor copies (asymmetric unit), we keep the copy
    nearest the receptor-R D3.32 anchor.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import paths, structure

# pdb id -> (receptor segid, ligand selection, human-readable ligand name, ligand class)
SPEC: dict[str, tuple[str, str, str, str]] = {
    "ZH853": ("R", f"resname {structure.LIGAND_RESNAME}", "ZH853", "cyclic-peptide"),
    "8F7R": ("R", "segid P Q or resname NH2", "endomorphin-1", "peptide"),
    "8EFQ": ("R", "segid P or resname DAL MEA ETA", "DAMGO", "peptide"),
    "6DDE": ("R", "segid D or resname DAL MEA ETA", "DAMGO", "peptide"),
    "8F7Q": ("R", "segid P", "beta-endorphin", "peptide"),
    "9WST": ("R", "segid P or resname DAL MEA ETA", "DAMGO (Gz)", "peptide"),
    "9WSV": ("R", "segid P or resname DAL MEA ETA", "DAMGO (arrestin)", "peptide"),
    "5C1M": ("A", "resname VF1", "BU72", "small-molecule"),
    "7T2G": ("R", "resname EIG", "mitragynine-PI", "small-molecule"),
    "8EF5": ("R", "resname 7V7", "fentanyl", "small-molecule"),
    "8EFB": ("R", "resname WH2", "oliceridine", "small-molecule"),
    "8EFL": ("R", "resname WH9", "SR-17018", "small-molecule"),
    "8EFO": ("R", "resname 8QY", "PZM21", "small-molecule"),
    "4DKL": ("A", "resname BF0", "beta-FNA", "antagonist"),
}


@dataclass
class Complex:
    """A loaded comparator: receptor and ligand atom groups + numbering offset."""

    pdb: str
    name: str
    ligand_class: str
    universe: object
    receptor: object  # MDAnalysis AtomGroup (protein, receptor chain)
    ligand: object  # MDAnalysis AtomGroup (agonist)
    offset_to_human: int  # add to construct resid to get human OPRM1 numbering


def _detect_offset(receptor) -> int:
    """Return the offset (0 human, +2 mouse) mapping the receptor to human OPRM1."""
    by_id = {int(a.resid): str(a.resname) for a in receptor}
    if by_id.get(149) == "ASP" and by_id.get(299) == "HIS" and by_id.get(328) == "TYR":
        return 0
    if by_id.get(147) == "ASP" and by_id.get(297) == "HIS" and by_id.get(326) == "TYR":
        return 2
    raise ValueError("Could not detect OPRM1 numbering scheme (no D3.32/H6.52/Y7.43 landmark)")


def _nearest_copy(ligand, receptor, offset: int):
    """If the ligand spans multiple chains, keep the copy nearest receptor D3.32."""
    segids = sorted(set(ligand.segids))
    if len(segids) <= 1:
        return ligand
    d332 = receptor.select_atoms(f"resid {149 - offset} and name OD1 OD2")
    if not len(d332):
        return ligand
    anchor = d332.positions.mean(axis=0)
    best, best_d = segids[0], np.inf
    for seg in segids:
        sub = ligand.select_atoms(f"segid {seg}")
        d = float(np.linalg.norm(sub.positions - anchor, axis=1).min())
        if d < best_d:
            best_d, best = d, seg
    return ligand.select_atoms(f"segid {best}")


def load_complex(pdb: str, pdb_path: Path | None = None) -> Complex:
    """Load a comparator (or 'ZH853') into a :class:`Complex`."""
    seg, ligsel, name, lclass = SPEC[pdb]
    path = paths.CRYOEM_PDB if pdb == "ZH853" else (pdb_path or paths.COMPARATORS / f"{pdb}.pdb")
    u = structure.load(path)
    receptor = u.select_atoms(f"segid {seg} and protein")
    offset = _detect_offset(receptor)
    ligand = _nearest_copy(u.select_atoms(ligsel), receptor, offset)
    if not len(ligand):
        raise ValueError(f"{pdb}: empty ligand selection {ligsel!r}")
    return Complex(pdb, name, lclass, u, receptor, ligand, offset)
