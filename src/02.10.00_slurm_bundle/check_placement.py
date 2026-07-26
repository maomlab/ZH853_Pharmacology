#!/usr/bin/env python3
"""Verify that PACKMOL-Memgen preserved the OPM membrane registration (SPECIFICATION D-14).

02.05.00 superposes the receptor onto the OPM reference (6DDF) so that the OPM midplane is at
z = 0. PACKMOL-Memgen, however, re-centres the solute on its own z bounding box before building
the bilayer -- and our receptor's bbox centre is ~5 A below the OPM midplane, because the
intracellular face (H8, ICL3, C-term) protrudes further than the extracellular face. If memgen
re-centres, the receptor ends up ~5 A too high in the membrane: the Trp/Tyr girdle falls out of
register with the interface and hydrophobic belt residues face the headgroups.

This script measures the offset directly in the built system instead of assuming either way:
  1. the receptor is a rigid body common to receptor.pdb (OPM frame) and bilayer_system.pdb, so
     the z shift memgen applied is just the mean CA displacement;
  2. the bilayer midplane is measured independently from the lipid phosphate planes;
  3. misregistration = receptor shift - bilayer midplane, i.e. how far the receptor sits above
     where OPM puts it, in the frame of the membrane that was actually built.

Exits non-zero if |misregistration| exceeds TOL_A so 01_build_system.sh fails loudly rather than
handing a silently mis-embedded system to tleap.

Run: python check_placement.py [built.pdb] [oriented_reference.pdb]
"""

from __future__ import annotations

import sys

import numpy as np

TOL_A = 1.5  # max tolerated receptor/bilayer misregistration

AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASH", "GLH",
    "LYN", "ACE", "NME",
}
TRP_RING = {"CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"}
TYR_RING = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"}
OPM_HALF = 15.7  # OPM 6DDF DUM boundary -> 31.4 A hydrophobic thickness (cf. 02.06.00)


def read_pdb(path):
    """Yield (resname, resseq, atomname, xyz) for every ATOM/HETATM record, in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                out.append((
                    line[17:20].strip(), line[22:27].strip(), line[12:16].strip(),
                    np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                ))
    return out


def ca_positions(atoms) -> np.ndarray:
    return np.array([xyz for rn, _, an, xyz in atoms if an == "CA" and rn in AA])


def ring_centroids(atoms) -> np.ndarray:
    """z of every Trp/Tyr ring centroid -- the interfacial 'aromatic girdle' marker."""
    byres, order = {}, []
    for rn, ri, an, xyz in atoms:
        if rn in ("TRP", "TYR"):
            key = (rn, ri)
            if key not in byres:
                byres[key] = []
                order.append(key)
            byres[key].append((an, xyz))
    zs = []
    for rn, ri in order:
        ring = TRP_RING if rn == "TRP" else TYR_RING
        pts = [xyz for an, xyz in byres[(rn, ri)] if an in ring]
        if len(pts) >= 5:
            zs.append(float(np.mean([p[2] for p in pts])))
    return np.array(zs)


def phosphate_planes(atoms):
    """(midplane z, P-P thickness, n) from the lipid phosphate atoms of both leaflets."""
    pz = np.array([xyz[2] for rn, _, an, xyz in atoms if an in ("P", "P31") and rn not in AA])
    if len(pz) < 4:
        return None
    lower, upper = pz[pz < pz.mean()], pz[pz >= pz.mean()]
    return float((lower.mean() + upper.mean()) / 2), float(upper.mean() - lower.mean()), len(pz)


def main() -> int:
    built = sys.argv[1] if len(sys.argv) > 1 else "bilayer_system.pdb"
    ref = sys.argv[2] if len(sys.argv) > 2 else "receptor.pdb"

    ca_built, ca_ref = ca_positions(read_pdb(built)), ca_positions(read_pdb(ref))
    if len(ca_built) != len(ca_ref):
        print(f"ERROR: {len(ca_built)} CA in {built} vs {len(ca_ref)} in {ref}; cannot compare.",
              file=sys.stderr)
        return 2

    delta = ca_built - ca_ref
    shift = delta.mean(axis=0)
    residual = float(np.abs(delta - shift).max())
    if residual > 0.05:
        print(f"ERROR: receptor is not a rigid translation of {ref} (max residual {residual:.2f} A)."
              " It was rotated or rebuilt; the OPM frame no longer applies.", file=sys.stderr)
        return 2

    atoms_built = read_pdb(built)
    planes = phosphate_planes(atoms_built)
    if planes is None:
        print(f"ERROR: fewer than 4 lipid phosphate atoms in {built}; is this a bilayer?",
              file=sys.stderr)
        return 2
    midplane, pp, n_p = planes

    # In the OPM frame the midplane is z = 0, so the receptor sits `shift[2]` above where OPM put
    # it, measured against the bilayer that was actually built.
    offset = float(shift[2]) - midplane

    zs = ring_centroids(atoms_built) - midplane
    upper = zs[(zs > 7) & (zs < 22)]
    lower = zs[(zs < -7) & (zs > -22)]

    print(f"receptor translation applied by packmol-memgen  = "
          f"({shift[0]:+.2f}, {shift[1]:+.2f}, {shift[2]:+.2f}) A")
    print(f"bilayer midplane from {n_p} lipid phosphates      = {midplane:+.2f} A "
          f"(P-P thickness {pp:.1f} A)")
    if len(upper) and len(lower):
        print(f"Trp/Tyr girdle about the built midplane        = "
              f"{lower.mean():+.1f} / {upper.mean():+.1f} A "
              f"(girdle thickness {upper.mean() - lower.mean():.1f} A; OPM slab "
              f"{-OPM_HALF:+.1f} / {OPM_HALF:+.1f})")
    print(f"OPM misregistration                            = {offset:+.2f} A "
          f"(tolerance {TOL_A:.1f} A)")

    if abs(offset) <= TOL_A:
        print("OK: the OPM registration from 02.05.00 survived the build.")
        return 0

    print()
    print(f"FAIL: the receptor sits {offset:+.2f} A off the OPM position in the built bilayer.")
    print("packmol-memgen re-centred the solute on its z bounding box, which is not the OPM")
    print("midplane. Remedies, in order of preference:")
    print(f"  1. Re-run packmol-memgen with an explicit z offset of {-offset:+.2f} A if your build")
    print("     supports it:  packmol-memgen --help | grep -iE 'translate|offset|center'")
    print("  2. Failing that, defeat the re-centring by handing memgen a receptor whose z bounding")
    print(f"     box is centred on the OPM midplane, then shift the packed bilayer back by "
          f"{offset:+.2f} A.")
    print("  3. Do NOT 'fix' this by translating the receptor inside the packed box -- the lipids")
    print("     were packed around its old position and would be driven through the protein.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
