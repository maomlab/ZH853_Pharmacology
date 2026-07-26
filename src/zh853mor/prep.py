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


# --- receptor finalisation: His tautomers and neutral terminal caps ---------------------------
# Both exist because a downstream default silently overrides an upstream decision. tleap maps a
# residue named HIS to HIE regardless of the tautomer that was actually determined, and a chain
# that starts/ends on a standard residue is built with charged termini. Neither raises an error,
# so both are fixed here, in the receptor, rather than trusted to the assembly step.

# cols: 1-6 record, 7-11 serial, 13-16 name, 17 altLoc, 18-20 resName, 22 chain, 23-26 resSeq,
# 31-54 xyz, 55-60 occupancy, 61-66 B, 77-78 element.
_PDB_FMT = ("ATOM  {serial:5d} {name:<4s} {resname:>3s} {chain:1s}{resid:4d}    "
            "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}  ")


def _pdb_atoms(path) -> list[dict]:
    out = []
    with open(path) as fh:
        lines = list(fh)
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            out.append({"name": line[12:16].strip(), "resname": line[17:20].strip(),
                        "chain": line[21], "resid": int(line[22:26]),
                        "xyz": np.array([float(line[30:38]), float(line[38:46]),
                                         float(line[46:54])]),
                        "elem": line[76:78].strip() or line[12:16].strip()[0]})
    return out


def his_tautomers(protonated_pdb) -> dict[int, str]:
    """Read back the HID/HIE/HIP choice OpenMM made, from which N carries a hydrogen.

    PDBFixer's ``addMissingHydrogens`` picks the tautomer from the local hydrogen-bond network
    (it estimates both HD1 and HE2 positions and counts nearby acceptors), which is the standard
    way to make this call. It does *not* rename the residue, so the decision is invisible unless
    it is recovered here and written into the residue name.
    """
    borne: dict[int, set[str]] = {}
    for a in _pdb_atoms(protonated_pdb):
        if a["resname"] in ("HIS", "HID", "HIE", "HIP"):
            borne.setdefault(a["resid"], set()).add(a["name"])
    out = {}
    for resid, names in borne.items():
        hd1, he2 = "HD1" in names, "HE2" in names
        out[resid] = "HIP" if hd1 and he2 else "HID" if hd1 else "HIE" if he2 else "HIS"
    return out


def _place(a: np.ndarray, b: np.ndarray, c: np.ndarray,
           bond: float, angle: float, dihedral: float) -> np.ndarray:
    """Natural-extension-reference-frame placement of a fourth atom from three known ones."""
    th, ph = np.radians(angle), np.radians(dihedral)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-bond * np.cos(th), bond * np.sin(th) * np.cos(ph), bond * np.sin(th) * np.sin(ph)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def cap_termini(pdb_in, pdb_out) -> dict:
    """Write a copy of a heavy-atom protein PDB with neutral ACE/NME caps on the termini.

    The construct is an internal fragment of full-length OPRM1 (69--349 of 400), so charged
    termini would put two formal charges where the real protein has peptide bonds --- one of them
    at the extracellular face. Cap geometry is ideal-internal-coordinate; the free backbone
    torsion is chosen by scanning it and taking the rotamer furthest from every existing atom, so
    the caps start clash-free rather than needing minimisation to escape a bad guess. The
    C-terminal OXT is dropped, since that carboxylate oxygen becomes the NME amide nitrogen.
    """
    atoms = _pdb_atoms(pdb_in)
    resids = sorted({a["resid"] for a in atoms})
    first, last = resids[0], resids[-1]
    chain = atoms[0]["chain"]
    xyz = {(a["resid"], a["name"]): a["xyz"] for a in atoms}

    def best(build, anchor: int) -> tuple[list[tuple[str, np.ndarray, str]], float, float]:
        """Scan the free backbone torsion; keep the rotamer furthest from every other atom.

        The anchor residue is excluded from the score: the cap is covalently bonded to it, so its
        1-2/1-3 distances are fixed by the ideal geometry and would otherwise floor the score at
        the C-N bond length (~1.34 A) and make the scan blind to real clashes elsewhere.
        """
        cloud = np.array([a["xyz"] for a in atoms if a["resid"] != anchor])
        top = (None, -1.0, 0.0)
        for tor in range(0, 360, 5):
            placed = build(float(tor))
            gap = min(float(np.linalg.norm(cloud - p, axis=1).min()) for _, p, _ in placed)
            if gap > top[1]:
                top = (placed, gap, float(tor))
        return top

    # ACE: C bonded to N(first); torsion scanned is C(ACE)-N-CA-C, i.e. phi of the first residue.
    def build_ace(tor: float):
        c = _place(xyz[(first, "C")], xyz[(first, "CA")], xyz[(first, "N")], 1.335, 121.7, tor)
        o = _place(xyz[(first, "CA")], xyz[(first, "N")], c, 1.229, 122.9, 0.0)
        ch3 = _place(xyz[(first, "CA")], xyz[(first, "N")], c, 1.522, 116.6, 180.0)
        return [("CH3", ch3, "C"), ("C", c, "C"), ("O", o, "O")]

    # NME: N bonded to C(last); torsion scanned is N(NME)-C-CA-N, i.e. psi of the last residue.
    def build_nme(tor: float):
        n = _place(xyz[(last, "N")], xyz[(last, "CA")], xyz[(last, "C")], 1.335, 116.6, tor)
        ch3 = _place(xyz[(last, "CA")], xyz[(last, "C")], n, 1.449, 121.7, 180.0)
        return [("N", n, "N"), ("CH3", ch3, "C")]

    ace, ace_gap, phi = best(build_ace, first)
    nme, nme_gap, psi = best(build_nme, last)

    out, serial = [], 0

    def emit(name, resname, resid, pos, elem):
        nonlocal serial
        serial += 1
        pdb_name = f" {name:<3s}" if len(name) < 4 else name
        out.append(_PDB_FMT.format(serial=serial, name=pdb_name, resname=resname, chain=chain,
                                   resid=resid, x=pos[0], y=pos[1], z=pos[2], elem=elem))

    for name, pos, elem in ace:
        emit(name, "ACE", first - 1, pos, elem)
    for a in atoms:
        if a["name"] == "OXT":     # replaced by the NME amide nitrogen
            continue
        emit(a["name"], a["resname"], a["resid"], a["xyz"], a["elem"])
    for name, pos, elem in nme:
        emit(name, "NME", last + 1, pos, elem)
    out.append("TER")
    out.append("END")
    with open(pdb_out, "w") as fh:
        fh.write("\n".join(out) + "\n")

    return {"ace_resid": first - 1, "nme_resid": last + 1, "phi": phi, "psi": psi,
            "min_contact": round(min(ace_gap, nme_gap), 2), "dropped_oxt": True}


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
