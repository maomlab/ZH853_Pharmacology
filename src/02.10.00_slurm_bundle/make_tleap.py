#!/usr/bin/env python3
"""Generate tleap_run.in from the tleap.in template and the actual packed system.

Two values in the tleap input can only be known after PACKMOL-Memgen has run:

  disulfides  loadpdb renumbers residues sequentially from 1 across the entire system, so the
              OPRM1 numbering used to write `bond sys.142.SG sys.219.SG` no longer exists. This
              finds SG-SG pairs geometrically and emits the bond in tleap's numbering. The bonded
              cysteines are also renamed CYS -> CYX, without which ff19SB builds an HG onto a
              sulfur that is about to get a second bond.

  box         the periodic cell must be the one PACKMOL packed into. Candidates are collected
              from the CRYST1 record, packmol-memgen.json, and the `inside box` constraints in
              packmol.inp, then each is checked against the packed coordinates -- the first that
              actually contains them wins. A cell smaller than its own contents makes atoms wrap
              onto their periodic images, which no amount of minimisation repairs.

Run (in the build directory, after packmol-memgen): python make_tleap.py
"""

from __future__ import annotations

import contextlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SS_CUTOFF_A = 2.5  # SG-SG distance below which a disulfide is called (bonded is ~2.05 A)
LIGAND_RESNAMES = ("ZH8", "L01", "LIG", "MOL", "UNL")


def residue_key(line: str):
    """Fields tleap starts a new residue on: chain, sequence number, insertion code, name."""
    return line[21], line[22:26], line[26], line[17:20]


def load(path: Path):
    """Return (all lines, [(seq_index, resname, [(atomname, xyz, line_index)])]) in file order."""
    lines = path.read_text().splitlines()
    residues, prev = [], None
    for i, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = residue_key(line)
        if key != prev:
            residues.append((len(residues) + 1, line[17:20].strip(), []))
            prev = key
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        residues[-1][2].append((line[12:16].strip(), xyz, i))
    return lines, residues


def find_disulfides(residues):
    """[(seq_i, seq_j, distance)] for every SG pair closer than SS_CUTOFF_A."""
    sgs = [(idx, xyz) for idx, rn, atoms in residues if rn in ("CYS", "CYX")
           for an, xyz, _ in atoms if an == "SG"]
    out = []
    for a in range(len(sgs)):
        for b in range(a + 1, len(sgs)):
            d = float(np.linalg.norm(sgs[a][1] - sgs[b][1]))
            if d < SS_CUTOFF_A:
                out.append((sgs[a][0], sgs[b][0], d))
    return out


def rename_cyx(lines, residues, bonded: set[int]) -> list[str]:
    out = list(lines)
    for idx, _, atoms in residues:
        if idx in bonded:
            for _, _, li in atoms:
                out[li] = out[li][:17] + "CYX" + out[li][20:]
    return out


def box_candidates(lines, packmol_inp: Path, json_path: Path) -> list[tuple[str, tuple]]:
    """Every periodic cell we can derive, best source first.

    These disagree, and the difference matters: a box smaller than the coordinates means atoms
    wrap on top of their own periodic images, which minimisation cannot fix.
    """
    out = []
    for line in lines:
        if line.startswith("CRYST1"):
            dims = (float(line[6:15]), float(line[15:24]), float(line[24:33]))
            if all(d > 1.0 for d in dims):
                out.append(("the CRYST1 record", dims))
            break

    if json_path.exists():
        try:
            blob = json.loads(json_path.read_text())
        except (ValueError, OSError):
            blob = None
        if isinstance(blob, dict):
            for key in ("box", "box_size", "boxsize", "dimensions", "cell", "pbc"):
                v = blob.get(key)
                if isinstance(v, (list, tuple)) and len(v) >= 3:
                    with contextlib.suppress(TypeError, ValueError):
                        out.append((f"packmol-memgen.json['{key}']",
                                    tuple(float(x) for x in v[:3])))
            if not any("json" in s for s, _ in out):
                print(f"note: {json_path.name} has no recognised box key; top-level keys are "
                      f"{sorted(blob)[:12]}", file=sys.stderr)

    if packmol_inp.exists():
        lo, hi = None, None
        for line in packmol_inp.read_text().splitlines():
            f = line.split()
            if len(f) == 8 and f[0] == "inside" and f[1] == "box":
                v = np.array([float(x) for x in f[2:]])
                lo = v[:3] if lo is None else np.minimum(lo, v[:3])
                hi = v[3:] if hi is None else np.maximum(hi, v[3:])
        if lo is not None:
            out.append((f"the `inside box` bounds in {packmol_inp}",
                        tuple(float(x) for x in (hi - lo))))
    return out


def verify_box(dims, residues) -> tuple[bool, str]:
    """Check the cell actually contains the packed coordinates.

    PACKMOL's `inside box` constrains the molecules it places, but nothing constrains the fixed
    solute, and the reported cell can be smaller than what was packed. If the span exceeds the
    cell, atoms near one face wrap onto atoms near the opposite face.
    """
    pts = np.array([xyz for _, _, atoms in residues for _, xyz, _ in atoms])
    span = pts.max(axis=0) - pts.min(axis=0)
    excess = span - np.array(dims)
    msg = (f"coordinate span {span[0]:.2f} x {span[1]:.2f} x {span[2]:.2f} A "
           f"vs cell {dims[0]:.2f} x {dims[1]:.2f} x {dims[2]:.2f} A")
    if (excess <= 0.5).all():
        return True, msg + "  (fits)"

    worst = []
    for ax, name in enumerate("xyz"):
        if excess[ax] > 0.5:
            half = dims[ax] / 2.0
            centre = (pts[:, ax].max() + pts[:, ax].min()) / 2.0
            outside = defaultdict(int)
            for _, resname, atoms in residues:
                for _, xyz, _ in atoms:
                    if abs(xyz[ax] - centre) > half:
                        outside[resname] += 1
            top = ", ".join(f"{k} x{v}" for k, v in
                            sorted(outside.items(), key=lambda kv: -kv[1])[:5])
            worst.append(f"    {name}: {excess[ax]:+.2f} A over; "
                         f"{sum(outside.values())} atoms outside, mostly {top}")
    return False, msg + "  (DOES NOT FIT)\n" + "\n".join(worst)


def main() -> int:
    build = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("bilayer_system.pdb")
    template = Path("tleap.in")
    if not build.exists() or not template.exists():
        print(f"ERROR: need both {build} and {template} in the current directory.", file=sys.stderr)
        return 1

    lines, residues = load(build)
    ss = find_disulfides(residues)
    if not ss:
        print("ERROR: no SG-SG pair within "
              f"{SS_CUTOFF_A} A. The conserved OPRM1 C142-C219 disulfide is missing from the packed"
              " system -- check that the receptor survived packmol-memgen's preprocessing.",
              file=sys.stderr)
        return 1

    bonded = {i for i, j, _ in ss} | {j for i, j, _ in ss}
    out_pdb = build.with_name(build.stem + "_ff.pdb")
    out_pdb.write_text("\n".join(rename_cyx(lines, residues, bonded)) + "\n")

    bonds = "\n".join(f"bond sys.{i}.SG sys.{j}.SG    # SG-SG {d:.2f} A" for i, j, d in ss)
    for i, j, d in ss:
        print(f"disulfide: residues {i}-{j} (SG-SG {d:.2f} A), renamed CYX")

    cands = box_candidates(lines, Path("packmol.inp"), Path("packmol-memgen.json"))
    if not cands:
        print("ERROR: could not determine the periodic box from CRYST1, packmol-memgen.json, "
              "or packmol.inp.", file=sys.stderr)
        return 1
    for src, d in cands:
        print(f"candidate cell {d[0]:.2f} x {d[1]:.2f} x {d[2]:.2f} A (from {src})")

    dims = None
    for src, d in cands:
        ok, msg = verify_box(d, residues)
        if ok:
            dims = d
            print(f"box {d[0]:.2f} x {d[1]:.2f} x {d[2]:.2f} A (from {src}); {msg}")
            break
    if dims is None:
        src, d = cands[0]
        _, msg = verify_box(d, residues)
        print(f"\nERROR: no candidate periodic cell contains the packed coordinates.\n{msg}\n",
              file=sys.stderr)
        print("Do not just enlarge the cell: that leaves a low-density gap at the faces that the\n"
              "membrane barostat will close by collapsing the bilayer. Re-pack instead, and check\n"
              "whether the overflowing residues are the fixed solute (not constrained by PACKMOL's\n"
              "`inside box`) or packed lipids/water (a PACKMOL convergence failure -- see the final\n"
              "'Maximum violation of the restraints' in packmol.log).", file=sys.stderr)
        return 1

    if not any(rn in LIGAND_RESNAMES for _, rn, _ in residues):
        print()
        print("WARNING: no ZH853 ligand residue found in the packed system. 01_build_system.sh")
        print("  stages the receptor-only PDB, so packmol-memgen packed an APO receptor and")
        print("  loadmol2 alone will not put the ligand in the box. This will build System B-apo,")
        print("  not System A. See README 'Systems to build'.")

    text = (template.read_text()
            .replace("@SYSTEM_PDB@", out_pdb.name)
            .replace("@DISULFIDES@", bonds)
            .replace("@SETBOX@", f"set sys box {{ {dims[0]:.3f} {dims[1]:.3f} {dims[2]:.3f} }}"))
    Path("tleap_run.in").write_text(text)
    print(f"Wrote tleap_run.in and {out_pdb.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
