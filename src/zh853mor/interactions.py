"""Heavy-atom geometric protein-ligand interaction fingerprints.

Deliberately hydrogen-free: the cryo-EM model is 3.5 A with no modeled hydrogens, so
angle-based H-bond criteria would be over-precise. We classify interactions from
heavy-atom geometry, which is the appropriate resolution-matched approach and applies
uniformly to ZH853 and every comparator.

Interaction types per receptor residue:
  - ``ionic``       opposite-charge groups within IONIC_CUT (salt bridge)
  - ``hbond``       ligand N/O -- residue polar N/O within HBOND_CUT (donor/acceptor proxy)
  - ``hydrophobic`` ligand C -- residue hydrophobic sidechain C within HYDROPHOBIC_CUT
  - ``aromatic``    ligand aromatic ring centroid -- residue aromatic centroid within ARO_CUT
  - ``cation_pi``   residue cationic N -- ligand aromatic centroid within CATIONPI_CUT

Cutoffs follow common PLIP/ProLIF-style geometric conventions, relaxed slightly for 3.5 A.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

CONTACT_CUT = 4.5
IONIC_CUT = 4.0
HBOND_CUT = 3.5
HYDROPHOBIC_CUT = 4.0
ARO_CUT = 5.5
CATIONPI_CUT = 6.0

# Protein charged sidechain atoms.
ANIONIC = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
CATIONIC = {"LYS": {"NZ"}, "ARG": {"NH1", "NH2", "NE"}, "HIS": {"ND1", "NE2"}}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS", "TYR"}

# Aromatic ring atom names for standard residues (six-membered ring where applicable).
AROMATIC_RING_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # benzene ring of indole
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}


@dataclass
class ResidueFingerprint:
    """Interactions between one receptor residue and the ligand."""

    resid: int  # human OPRM1 numbering
    resname: str
    interactions: set[str] = field(default_factory=set)
    min_dist: float = np.inf


def _element(atom) -> str:
    """Best-effort element symbol from an MDAnalysis atom."""
    el = getattr(atom, "element", "") or ""
    if not el:
        name = atom.name.strip()
        el = name[0] if name and not name[0].isdigit() else (name[1:2] if len(name) > 1 else "")
    return el.upper()


def _ligand_aromatic_centroids(ligand) -> list[np.ndarray]:
    """Ring centroids of the ligand's aromatic rings.

    Standard aromatic amino-acid residues (in peptide ligands) use an atom-name template.
    Non-standard/small-molecule ligands fall back to RDKit bond perception.
    """
    centroids: list[np.ndarray] = []
    used_perception = False
    for res in ligand.residues:
        template = AROMATIC_RING_ATOMS.get(res.resname)
        if template:
            ring = res.atoms.select_atoms("name " + " ".join(template))
            if len(ring) >= 5:
                centroids.append(ring.positions.mean(axis=0))
        else:
            used_perception = True
    if used_perception:
        centroids.extend(_rdkit_aromatic_centroids(ligand))
    return centroids


def _rdkit_aromatic_centroids(ligand, planarity_tol: float = 0.15) -> list[np.ndarray]:
    """Perceive aromatic rings on a non-standard ligand by geometry.

    Template- and bond-order-free: build connectivity from atom proximity, then treat any
    5-/6-membered ring of C/N/O atoms that is planar (max out-of-plane deviation below
    ``planarity_tol`` A) as aromatic. This distinguishes flat aromatic rings from puckered
    aliphatic rings (e.g. the ZH853 macrocyclic lactam) without needing hydrogens or charges.
    """
    centroids: list[np.ndarray] = []
    try:
        from rdkit import Chem
    except ImportError:
        return centroids
    lines = []
    for i, a in enumerate(ligand.atoms):
        x, y, z = a.position
        el = _element(a).capitalize() or "C"
        lines.append(
            f"HETATM{i + 1:>5} {a.name:<4} LIG A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}"
        )
    block = "\n".join(lines) + "\nEND\n"
    mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True, proximityBonding=True)
    if mol is None:
        return centroids
    Chem.FastFindRings(mol)
    conf = mol.GetConformer()
    aromatic_elems = {"C", "N", "O"}
    for ring in mol.GetRingInfo().AtomRings():
        if not 5 <= len(ring) <= 6:
            continue
        if any(mol.GetAtomWithIdx(i).GetSymbol().upper() not in aromatic_elems for i in ring):
            continue
        pts = np.array([list(conf.GetAtomPosition(i)) for i in ring])
        centered = pts - pts.mean(axis=0)
        # smallest singular vector is the plane normal; its singular value ~ out-of-plane spread
        normal = np.linalg.svd(centered)[2][-1]
        if float(np.abs(centered @ normal).max()) < planarity_tol:
            centroids.append(pts.mean(axis=0))
    return centroids


def _protein_aromatic_centroid(residue) -> np.ndarray | None:
    template = AROMATIC_RING_ATOMS.get(residue.resname)
    if not template:
        return None
    ring = residue.atoms.select_atoms("name " + " ".join(template))
    return ring.positions.mean(axis=0) if len(ring) >= 5 else None


def fingerprint(receptor, ligand, offset_to_human: int = 0) -> dict[int, ResidueFingerprint]:
    """Compute the per-residue interaction fingerprint for one complex.

    ``receptor``/``ligand`` are MDAnalysis AtomGroups; ``offset_to_human`` is added to
    receptor resids so keys are always human OPRM1 numbers.
    """
    lig_pos = ligand.positions
    lig_el = np.array([_element(a) for a in ligand.atoms])
    lig_N = lig_pos[lig_el == "N"]
    lig_O = lig_pos[lig_el == "O"]
    lig_C = lig_pos[lig_el == "C"]
    lig_polar = lig_pos[(lig_el == "N") | (lig_el == "O")]
    lig_aro = _ligand_aromatic_centroids(ligand)

    # Restrict to receptor residues with any atom within CONTACT_CUT of the ligand.
    # Evaluate the geometric selection on the universe, then intersect with the receptor.
    near = receptor.universe.select_atoms(f"around {CONTACT_CUT} group lig", lig=ligand) & receptor
    out: dict[int, ResidueFingerprint] = {}

    for res in near.residues:
        hid = int(res.resid) + offset_to_human
        fp = ResidueFingerprint(resid=hid, resname=str(res.resname))
        atoms = res.atoms
        pos = atoms.positions
        names = atoms.names
        el = np.array([_element(a) for a in atoms])

        # min heavy-heavy distance
        fp.min_dist = float(np.linalg.norm(pos[:, None, :] - lig_pos[None, :, :], axis=2).min())

        # hydrophobic: residue sidechain C -- ligand C
        if res.resname in HYDROPHOBIC_RES and len(lig_C):
            c = pos[el == "C"]
            if len(c) and np.linalg.norm(c[:, None] - lig_C[None], axis=2).min() < HYDROPHOBIC_CUT:
                fp.interactions.add("hydrophobic")

        # hbond proxy: residue polar N/O -- ligand polar N/O
        res_polar = pos[(el == "N") | (el == "O")]
        if len(res_polar) and len(lig_polar) and (
            np.linalg.norm(res_polar[:, None] - lig_polar[None], axis=2).min() < HBOND_CUT
        ):
            fp.interactions.add("hbond")

        # ionic: residue anionic O -- ligand N ; residue cationic N -- ligand O
        ani = ANIONIC.get(res.resname)
        if ani:
            a_pos = np.array([p for p, n in zip(pos, names, strict=False) if n in ani])
            if len(a_pos) and len(lig_N) and np.linalg.norm(
                a_pos[:, None] - lig_N[None], axis=2
            ).min() < IONIC_CUT:
                fp.interactions.add("ionic")
        cat = CATIONIC.get(res.resname)
        cat_pos = None
        if cat:
            cat_pos = np.array([p for p, n in zip(pos, names, strict=False) if n in cat])
            if len(cat_pos) and len(lig_O) and np.linalg.norm(
                cat_pos[:, None] - lig_O[None], axis=2
            ).min() < IONIC_CUT:
                fp.interactions.add("ionic")

        # aromatic pi-stacking: residue ring centroid -- ligand ring centroid
        res_aro = _protein_aromatic_centroid(res)
        if res_aro is not None and lig_aro and (
            min(np.linalg.norm(res_aro - lc) for lc in lig_aro) < ARO_CUT
        ):
            fp.interactions.add("aromatic")

        # cation-pi: residue cationic N -- ligand aromatic centroid
        if cat_pos is not None and len(cat_pos) and lig_aro and (
            min(np.linalg.norm(cp - lc) for cp in cat_pos for lc in lig_aro) < CATIONPI_CUT
        ):
            fp.interactions.add("cation_pi")

        if fp.interactions or fp.min_dist <= CONTACT_CUT:
            out[hid] = fp

    return out


def summarize(fp: dict[int, ResidueFingerprint]) -> str:
    """One-line-per-residue text summary, sorted by residue id."""
    lines = []
    for hid in sorted(fp):
        r = fp[hid]
        types = ",".join(sorted(r.interactions)) or "contact"
        lines.append(f"{r.resname}{hid}: {types} ({r.min_dist:.2f} A)")
    return "\n".join(lines)
