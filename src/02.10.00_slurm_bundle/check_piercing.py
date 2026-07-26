#!/usr/bin/env python3
"""Find lipid tails threaded through rings in a packed membrane system.

PACKMOL enforces a minimum pairwise distance, which does not prevent a lipid tail from passing
*through* the middle of a ring -- every atom can be >2 A from every other atom while the tail is
threaded through a Phe/Tyr/Trp/His/Pro ring or one of cholesterol's four fused rings. This is the
failure mode packmol-memgen tries to catch and could not here ("Lipid piercing finder failed").

It matters more than an ordinary clash because it is topological: minimisation pushes atoms apart
along straight lines, and no such path unthreads a ring. The lipid stays trapped for the whole
trajectory, distorting the ring, the local packing and any observable that depends on them.

Method: rings are found as 5- and 6-cycles in each residue's own bond graph (bonds inferred by
distance, so this works for cholesterol and the ligand without hardcoded atom names). For every
bond in a *different* residue that crosses a ring's plane, the crossing point is computed and
compared with the ring radius -- the standard segment/disc intersection test.

Run: python check_piercing.py [built.pdb]   (exit 1 if any piercing is found)
"""

from __future__ import annotations

import sys
from collections import defaultdict

import numpy as np

BOND_A = 1.75      # max heavy-atom covalent bond length
NEAR_A = 9.0       # only test bonds whose midpoint is within this of a ring centre
MARGIN = 1.0       # allow a bond to cross this far outside the ring radius and still count

SOLVENT = {"WAT", "HOH", "TIP3", "SOL", "NA", "NA+", "CL", "CL-", "K", "K+"}


def read_heavy(path):
    """[(resid_key, resname, atomname, xyz)] for heavy atoms, skipping water and ions."""
    out = []
    with open(path) as fh:
        lines = list(fh)
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name, resname = line[12:16].strip(), line[17:20].strip()
        if resname in SOLVENT or name.startswith("H") or line[76:78].strip() == "H":
            continue
        key = (line[21], line[22:27])          # chain + resSeq/icode, as written
        out.append((key, resname, name,
                    np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return out


def bonds_within(idx: list[int], xyz: np.ndarray) -> list[tuple[int, int]]:
    """Distance-inferred covalent bonds among a residue's atoms."""
    pts = xyz[idx]
    d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
    i, j = np.where((d > 0.5) & (d < BOND_A))
    return [(idx[a], idx[b]) for a, b in zip(i, j, strict=True) if a < b]


def find_rings(idx: list[int], bond_list: list[tuple[int, int]]) -> list[list[int]]:
    """All simple 5- and 6-cycles in a small bond graph."""
    adj = defaultdict(set)
    for a, b in bond_list:
        adj[a].add(b)
        adj[b].add(a)
    rings, seen = [], set()

    def walk(start, node, path):
        if len(path) > 6:
            return
        for nxt in adj[node]:
            if nxt == start and len(path) in (5, 6):
                key = frozenset(path)
                if key not in seen:
                    seen.add(key)
                    rings.append(list(path))
            elif nxt not in path and nxt > start:
                walk(start, nxt, path + [nxt])

    for a in idx:
        walk(a, a, [a])
    return rings


def main() -> int:
    path = sys.argv[1] if len(sys.argv) > 1 else "bilayer_system.pdb"
    atoms = read_heavy(path)
    if not atoms:
        print(f"ERROR: no non-solvent heavy atoms read from {path}.", file=sys.stderr)
        return 2
    xyz = np.array([a[3] for a in atoms])

    by_res = defaultdict(list)
    for i, (key, _resname, _, _) in enumerate(atoms):
        by_res[key].append(i)

    rings, all_bonds = [], []
    for key, idx in by_res.items():
        bl = bonds_within(idx, xyz)
        all_bonds.extend((key, a, b) for a, b in bl)
        for ring in find_rings(idx, bl):
            pts = xyz[ring]
            centre = pts.mean(axis=0)
            # plane normal = smallest-variance direction of the ring atoms
            normal = np.linalg.svd(pts - centre)[2][-1]
            radius = float(np.linalg.norm(pts - centre, axis=1).max())
            rings.append((key, atoms[ring[0]][1], centre, normal, radius))

    print(f"{len(rings)} rings and {len(all_bonds)} bonds in {len(by_res)} non-solvent residues")

    hits = []
    bond_mid = np.array([(xyz[a] + xyz[b]) / 2 for _, a, b in all_bonds])
    for rkey, rname, centre, normal, radius in rings:
        near = np.where(np.linalg.norm(bond_mid - centre, axis=1) < NEAR_A)[0]
        for bi in near:
            bkey, a, b = all_bonds[bi]
            if bkey == rkey:
                continue
            da = float(np.dot(xyz[a] - centre, normal))
            db = float(np.dot(xyz[b] - centre, normal))
            if da * db >= 0:                      # both endpoints on the same side
                continue
            t = da / (da - db)
            cross = xyz[a] + t * (xyz[b] - xyz[a])
            offset = float(np.linalg.norm(cross - centre))
            if offset < radius + MARGIN:
                hits.append((offset, radius, rkey, rname, bkey, atoms[a][1],
                             atoms[a][2], atoms[b][2]))

    if not hits:
        print("OK: no lipid tail or sidechain is threaded through a ring.")
        return 0

    hits.sort()
    print(f"\nFAIL: {len(hits)} ring piercing(s). A threaded tail cannot unthread during "
          "minimisation or MD --")
    print("it stays trapped for the whole trajectory. Re-pack (a different packmol seed usually "
          "clears it)")
    print("rather than trying to minimise it out; verify visually before spending GPU time.\n")
    for offset, radius, rkey, rname, bkey, bname, an, bn in hits:
        print(f"  {rname} {rkey[0]}{rkey[1].strip():>5s} ring (r={radius:.2f} A) pierced by "
              f"{bname} {bkey[0]}{bkey[1].strip():>5s} bond {an}-{bn} "
              f"at {offset:.2f} A from the ring centre")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
