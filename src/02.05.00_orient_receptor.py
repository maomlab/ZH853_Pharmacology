#!/usr/bin/env python3
"""Orient the receptor in the membrane by superposing onto an OPM/PPM reference (Phase 2).

Best practice is the OPM/PPM transfer-energy method (community standard). Since our cryo-EM model
is not itself in OPM, we transfer the placement structurally: fetch the OPM-oriented MOR reference
(6DDF, MOR-Gi-DAMGO; membrane normal along z, midplane at z=0, boundaries marked by DUM pseudo-atoms
-> hydrophobic thickness 31.4 A), superpose our receptor onto it (Kabsch on matched Ca), and apply
the transform to the whole complex. PACKMOL-Memgen `--preoriented` is then valid and the membrane
thickness comes from OPM, not the bound cholesterols.

If the OPM reference cannot be fetched (offline), falls back to a principal-axis + cholesterol-midplane
proxy (clearly weaker; see docs/METHODS_md_prep.md and SPECIFICATION D-14). Either way, 03.04.00
validates the placement against the Trp/Tyr aromatic girdle and experimental POPC thickness.

Run (local analysis env): python src/02.05.00_orient_receptor.py
"""

from __future__ import annotations

import sys
import urllib.error
import urllib.request

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import numpy as np  # noqa: E402

from zh853mor import paths, structure  # noqa: E402

IN_PDB = paths.INTERMEDIATE / "02.03.00_receptor" / "receptorR_fixed_heavy.pdb"
COMPLEX_PDB = paths.INTERMEDIATE / "02.01.00_components" / "receptor_ligand.pdb"
OPM_REF = "6ddf"  # MOR-Gi-DAMGO, OPM-oriented; DUM z = +/-15.7 -> 31.4 A hydrophobic thickness
OPM_URL = f"https://opm-assets.storage.googleapis.com/pdb/{OPM_REF}.pdb"
OPM_CACHE = paths.DATA / "opm" / f"{OPM_REF}.pdb"


def fetch_opm() -> bool:
    if OPM_CACHE.exists():
        return True
    paths.ensure_dir(OPM_CACHE.parent)
    try:
        with urllib.request.urlopen(OPM_URL, timeout=30) as r:  # noqa: S310 (trusted OPM host)
            OPM_CACHE.write_bytes(r.read())
        return True
    except (urllib.error.URLError, TimeoutError) as exc:
        print(f"WARNING: could not fetch OPM {OPM_REF}: {exc}", file=sys.stderr)
        return False


def ref_offset(ref_rec) -> int:
    """Residue offset to add to OUR (human) resid to index the reference (0 human, -2 mouse)."""
    by_id = {int(a.resid): a.resname for a in ref_rec if a.name == "CA"}
    if by_id.get(149) == "ASP" and by_id.get(328) == "TYR":
        return 0
    if by_id.get(147) == "ASP" and by_id.get(326) == "TYR":
        return -2
    raise ValueError("reference numbering not recognized (no D3.32/Y7.43 landmark)")


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t mapping points p onto q (least-squares)."""
    pc, qc = p - p.mean(0), q - q.mean(0)
    u, _, vt = np.linalg.svd(pc.T @ qc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, q.mean(0) - r @ p.mean(0)


def opm_orient(u):
    """Superpose our receptor onto the OPM reference; return (R, t, thickness, rmsd) or None."""
    if not fetch_opm():
        return None
    ref = structure.load(OPM_CACHE)
    ref_rec = ref.select_atoms("segid R and protein")
    if not len(ref_rec):
        ref_rec = ref.select_atoms("chainID R and protein")
    off = ref_offset(ref_rec)
    dum = ref.select_atoms("resname DUM")
    thickness = float(np.ptp(dum.positions[:, 2])) if len(dum) else float("nan")

    our_ca = {int(a.resid): a.position for a in u.select_atoms("segid R and name CA")}
    ref_ca = {int(a.resid): a.position for a in ref_rec.select_atoms("name CA")}
    common = sorted(r for r in our_ca if (r + off) in ref_ca)
    if len(common) < 50:
        print(f"WARNING: only {len(common)} matched Ca; skipping OPM orient", file=sys.stderr)
        return None
    p = np.array([our_ca[r] for r in common])
    q = np.array([ref_ca[r + off] for r in common])
    r, t = kabsch(p, q)
    rmsd = float(np.sqrt((((p @ r.T + t) - q) ** 2).sum(axis=1).mean()))
    return r, t, thickness, rmsd


def principal_axis_fallback(u):
    """Weaker proxy: normal = TM principal axis; midplane = cholesterol centre."""
    ca = u.select_atoms("segid R and name CA")
    x = ca.positions - ca.positions.mean(0)
    normal = np.linalg.svd(x, full_matrices=False)[2][0]
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(normal, z)
    s, c = np.linalg.norm(v), float(np.dot(normal, z))
    if s < 1e-8:
        r = np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        r = np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))
    orig = structure.load(paths.CRYOEM_PDB)
    clr = orig.select_atoms("segid R and resname CLR")
    ca_xy = (ca.positions @ r.T)[:, :2].mean(0)
    clr_z = float((clr.positions @ r.T)[:, 2].mean())
    return r, -np.array([ca_xy[0], ca_xy[1], clr_z])   # t: xy on receptor, z=0 at cholesterol midplane


def main() -> int:
    if not IN_PDB.exists():
        print(f"ERROR: {IN_PDB} not found -- run `make prep-receptor` first.", file=sys.stderr)
        return 1
    u = structure.load(IN_PDB)

    result = opm_orient(u)
    if result is not None:
        r, t, thickness, rmsd = result
        method = f"OPM ({OPM_REF.upper()}); hydrophobic thickness {thickness:.1f} A; " \
                 f"superposition RMSD {rmsd:.2f} A over the receptor"
    else:
        r, t = principal_axis_fallback(u)
        thickness = float("nan")
        method = "FALLBACK principal-axis + cholesterol midplane (OPM unavailable; weaker)"

    def apply(ag) -> None:
        ag.positions = ag.positions @ r.T + t

    apply(u.atoms)
    out_dir = paths.ensure_dir(paths.INTERMEDIATE / "02.05.00_oriented")
    u.atoms.write(str(out_dir / "receptorR_oriented.pdb"))
    if COMPLEX_PDB.exists():
        cx = structure.load(COMPLEX_PDB)
        apply(cx.atoms)
        cx.atoms.write(str(out_dir / "complex_oriented.pdb"))

    caz = u.select_atoms("segid R and name CA").positions[:, 2]
    print(f"Oriented by {method}")
    print(f"Membrane normal -> z, midplane at z=0; receptor Ca z-extent {np.ptp(caz):.1f} A")
    print("Validate the thickness with `make membrane-plot` (03.04.00). For a definitive per-structure")
    print("value, submit this model to the PPM 3.0 server (opm.phar.umich.edu/ppm_server).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
