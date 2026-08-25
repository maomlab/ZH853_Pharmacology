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
              onto their periodic images, which no amount of minimisation repairs. packmol-memgen
              2025.1.29 supplies none of the recorded sources reliably (no CRYST1 is ever written;
              the JSON appears only for --dry_run; its box summary is logged at DEBUG level), and
              PACKMOL leaves rim lipids/water ~1-2 A outside the `inside box` walls when the
              all-together runs end at maxiter -- so if every candidate fails, the packed extent
              + 2*CONTENT_MARGIN_A is used as a last resort rather than aborting the build.

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
# Per-face breathing room for the last-resort cell derived from the packed extent itself: keeps
# nearest periodic images >= 2.5 A apart instead of letting boundary atoms land on top of each
# other, without diluting the density PACKMOL actually realised (a much larger cell WOULD leave a
# vacuum gap for the membrane barostat to collapse).
CONTENT_MARGIN_A = 1.25
# Overflow past a recorded cell wall is only absorbable if it is rim squeeze (edge molecules'
# tails/waters poking a few A past the wall, most of each molecule still inside) rather than
# molecules parked well outside their packing region (holes inside the box, overlaps across the
# periodic seam -- a repack). Penetration DEPTH is the discriminator, not the fraction of a
# residue beyond the wall: memgen writes each POPC as three Lipid21 residues (PA/OL/PC), so a
# squeezed tail fragment is a "residue" that can be mostly outside while the lipid is fine.
# Depth needs no molecule grouping at all. PACKMOL's tolerance is 2.0 A and the systematic
# maxiter squeeze measured <= ~3.7 A on every build so far, so 6 A is a generous ceiling:
BENIGN_MAX_DEPTH_A = 6.0
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
        bounds = inside_box_bounds(packmol_inp)
        if bounds is not None:
            out.append((f"the `inside box` bounds in {packmol_inp}",
                        tuple(float(x) for x in (bounds[1] - bounds[0]))))
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


def content_cell(residues) -> tuple:
    """Last-resort periodic cell: the packed extent plus CONTENT_MARGIN_A on every face."""
    pts = np.array([xyz for _, _, atoms in residues for _, xyz, _ in atoms])
    return tuple(float(x) for x in pts.max(axis=0) - pts.min(axis=0) + 2 * CONTENT_MARGIN_A)


def inside_box_bounds(packmol_inp: Path):
    """(lo, hi) corner arrays of the union of the `inside box` regions, or None."""
    lo, hi = None, None
    if packmol_inp.exists():
        for line in packmol_inp.read_text().splitlines():
            f = line.split()
            if len(f) == 8 and f[0] == "inside" and f[1] == "box":
                v = np.array([float(x) for x in f[2:]])
                lo = v[:3] if lo is None else np.minimum(lo, v[:3])
                hi = v[3:] if hi is None else np.maximum(hi, v[3:])
    return (lo, hi) if lo is not None else None


def classify_overflow(residues, dims, bounds=None) -> tuple[bool, str]:
    """Decide whether overflow past a cell of size `dims` is safe to absorb.

    PACKMOL's `inside box` constrains every atom of the structure (memgen writes it at
    `structure` level), so atoms beyond a wall are real restraint violations from its maxiter
    stop. What matters is how DEEP they go past the wall. Edge molecules squeezed partway
    through -- tails/waters poking a few A out -- are 'rim squeeze': the interior stays dense
    and the overflow is absorbed by adopting the slightly larger cell the coordinates define.
    Atoms driven deeper than BENIGN_MAX_DEPTH_A mean misplaced molecules: holes inside the box,
    overlaps across the periodic seam. Repack.

    `bounds`, when available (the absolute `inside box` corners), anchors the walls so a lone
    far-out molecule cannot shift the reference frame and hide its own depth; otherwise the
    span midpoint is used, which is exact only for near-symmetric overflow.

    Returns (benign, report).
    """
    pts_all = np.array([xyz for _, _, atoms in residues for _, xyz, _ in atoms])
    span = pts_all.max(axis=0) - pts_all.min(axis=0)
    lines, all_shallow = [], True
    for ax, name in enumerate("xyz"):
        excess = span[ax] - dims[ax]
        if excess <= 0.5:
            continue
        n_res = n_at = 0
        worst_depth, worst_desc = 0.0, ""
        for seq, resname, atoms in residues:
            xs = np.array([xyz[ax] for _, xyz, _ in atoms])
            if bounds is not None:
                depth = np.maximum(bounds[0][ax] - xs, xs - bounds[1][ax])
            else:  # span-midpoint fallback
                centre = (pts_all[:, ax].max() + pts_all[:, ax].min()) / 2.0
                depth = np.abs(xs - centre) - dims[ax] / 2.0
            out = depth > 0
            if not out.any():
                continue
            n_res += 1
            n_at += int(out.sum())
            d_max = float(depth.max())
            if d_max > worst_depth:
                worst_depth, worst_desc = d_max, f"{resname} {seq}"
        bad = worst_depth > BENIGN_MAX_DEPTH_A
        all_shallow &= not bad
        anchor = "" if bounds is not None else " (span-midpoint reference)"
        lines.append(f"  {name}: +{excess:.2f} A over the wall{anchor}; {n_at} atoms in {n_res} "
                     f"residues poke out, deepest {worst_depth:.2f} A ({worst_desc})"
                     + ("   <-- deeper than rim squeeze; molecule(s) displaced" if bad else ""))
    return all_shallow, "\n".join(lines)


def diagnose_box(build: Path) -> int:
    """--diagnose: print every candidate cell and the overflow verdict; build nothing."""
    lines, residues = load(build)
    cands = box_candidates(lines, Path("packmol.inp"), Path("packmol-memgen.json"))
    pts_all = np.array([xyz for _, _, atoms in residues for _, xyz, _ in atoms])
    span = pts_all.max(axis=0) - pts_all.min(axis=0)
    print(f"{build}: coordinate span {span[0]:.2f} x {span[1]:.2f} x {span[2]:.2f} A, "
          f"{len(pts_all)} atoms in {len(residues)} residues")
    if not cands:
        print("no candidate cells found (CRYST1 / packmol-memgen.json / packmol.inp)")
        return 1
    for src, d in cands:
        ok, msg = verify_box(d, residues)
        print(f"\ncandidate {d[0]:.2f} x {d[1]:.2f} x {d[2]:.2f} A (from {src}): "
              f"{'fits' if ok else 'DOES NOT FIT'}")
        if not ok:
            bounds = inside_box_bounds(Path("packmol.inp")) if "inside box" in src else None
            benign, report = classify_overflow(residues, d, bounds)
            print(report)
            print("verdict: " + ("rim squeeze (absorbable)" if benign
                                 else "escaped molecules -- re-pack"))
    return 0


def main() -> int:
    args = [a for a in sys.argv[1:] if a != "--diagnose"]
    if "--diagnose" in sys.argv[1:]:
        return diagnose_box(Path(args[0]) if args else Path("bilayer_system.pdb"))
    build = Path(args[0]) if args else Path("bilayer_system.pdb")
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
        # Every recorded source failed. With packmol-memgen 2025.1.29 that is expected -- it
        # records its cell nowhere machine-readable (see box_candidates) -- but the coordinates
        # exceeding every candidate still means PACKMOL ended its maxiter runs with atoms past
        # the `inside box` walls. Absorb that with a larger cell only if it is rim squeeze;
        # escaped molecules mean holes inside the box and clashes across the seam: repack.
        ref_dims = max((d for _, d in cands), key=lambda d: d[0] * d[1] * d[2])
        ref_src = next(s for s, d in cands if d == ref_dims)
        bounds = inside_box_bounds(Path("packmol.inp")) if "inside box" in ref_src else None
        benign, report = classify_overflow(residues, ref_dims, bounds)
        print(f"\nOverflow diagnosis vs largest recorded candidate "
              f"{ref_dims[0]:.2f} x {ref_dims[1]:.2f} x {ref_dims[2]:.2f} A:")
        print(report)
        d = content_cell(residues)
        _, msg = verify_box(d, residues)
        if not benign:
            print("\nERROR: molecules sit wholly outside their packing region -- a PACKMOL",
                  file=sys.stderr)
            print("convergence failure, not rim squeeze. Re-pack in a fresh directory",
                  file=sys.stderr)
            print("(01_build_system.sh always builds pristine) and compare the final",
                  file=sys.stderr)
            print("'Maximum violation of the constraints' between the two packmol.log stages.",
                  file=sys.stderr)
            return 1
        print(f"\nWARNING: no recorded cell contains the packed coordinates; the overflow above")
        print(f"  is rim squeeze, absorbed with the packed extent + {2 * CONTENT_MARGIN_A:.2f} A.")
        print(f"box {d[0]:.2f} x {d[1]:.2f} x {d[2]:.2f} A (from packed coordinates); {msg}")
        dims = d

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
