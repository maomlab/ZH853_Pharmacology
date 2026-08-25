#!/usr/bin/env python3
"""Normalise the ACE/NME terminal caps in the packed system to tleap's templates.

packmol-memgen's preprocessing mangles the neutral caps on this stack in two ways:

  * reduce protonates the capping amides and some of its hydrogens survive memgen's
    H-stripping under non-template names -- an `HN2` rode along on the NME nitrogen;
  * the ff19SB libraries name the NME methyl carbon `C` (hydrogens H1/H2/H3), while the
    staged receptor from 02.03.00 follows the older pdbfixer convention `CH3`. loadpdb
    matches atoms by NAME, so the file's `CH3` arrives as an unknown atom ("Created a new
    atom named: CH3"), leap builds its own `C` at a default position, and check sys dies
    with "Atom .R<NME ...>.A<CH3> does not have a type". ACE (`CH3/C/O`) matches and loads
    cleanly -- but is audited anyway.

The fix is pure bookkeeping where possible: atoms whose identity is unambiguous keep their
packed coordinates and are renamed in place (`CH3` -> `C` moves nothing); stray hydrogens
are dropped because leap rebuilds every H from templates anyway; a heavy atom genuinely
absent from the file is filled from the reference receptor, rigidly translated by the
median offset over all paired CA atoms (the same offset check_placement measures).

The expected heavy-atom names come from the loaded library itself
($AMBERHOME/dat/leap/lib/aminoct12.lib), so a future force-field switch cannot silently
re-introduce the mismatch; a hardcoded fallback covers unreadable libraries.

Run (in the build directory, after packmol-memgen, before make_tleap.py):
    python fix_caps.py bilayer_system.pdb receptor.pdb
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np

CAP_RESNAMES = ("ACE", "NME")
# ff19SB-era heavy-atom conventions, used only when the library cannot be parsed:
FALLBACK_HEAVIES = {"ACE": ["CH3", "C", "O"], "NME": ["N", "C"]}
LIB_RELATIVE = Path("dat/leap/lib/aminoct12.lib")
ATOM_FMT = ("ATOM  {serial:5d} {name:>4s} {resname:<3s} {chain}{resid:4d}{icode:1s}   "
            "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}")


def template_heavies(resname: str) -> list[str]:
    """Heavy-atom names for a cap, read from the installed AmberTools library."""
    amberhome = os.environ.get("AMBERHOME")
    if amberhome:
        lib = Path(amberhome) / LIB_RELATIVE
        if lib.exists():
            in_table, names = False, []
            for line in lib.read_text().splitlines():
                if line.startswith(f"!entry.{resname}.unit.atoms table"):
                    in_table = True
                    continue
                if in_table:
                    if line.startswith("!"):  # next section: table complete
                        break
                    m = re.match(r'\s*"([^"]+)"\s+"([^"]+)"', line)
                    if m and not m.group(2).startswith("H"):  # types H*/HC/H1... are hydrogen
                        names.append(m.group(1))
            if names:
                return names
        print(f"note: {LIB_RELATIVE} unreadable under AMBERHOME={amberhome}; "
              f"using fallback template for {resname}", file=sys.stderr)
    return list(FALLBACK_HEAVIES[resname])


def element_of(line: str) -> str:
    """Element symbol: PDB columns 77-78 when populated, else inferred from the atom name."""
    elem = line[76:78].strip().upper()
    if elem:
        return elem
    name = line[12:16].strip()
    return name[0].upper() if not name[0].isdigit() else name[1].upper()


def cap_residues(lines: list[str]) -> list[dict]:
    """Group ATOM/HETATM lines into residues (as make_tleap.load does); return the caps."""
    residues, prev = [], None
    for i, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (line[21], line[22:26], line[26], line[17:20])
        if key != prev:
            residues.append({"first": line, "name": line[17:20].strip(), "atoms": []})
            prev = key
        residues[-1]["atoms"].append((line[12:16].strip(), i))
    return [r for r in residues if r["name"] in CAP_RESNAMES]


def ca_translation(packed: list[str], reference: list[str]) -> np.ndarray | None:
    """Median offset between the packed and staged receptor over paired CA atoms."""
    def cas(ls):
        pts = []
        for line in ls:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        return np.array(pts)
    p, r = cas(packed), cas(reference)
    if len(p) == 0 or len(p) != len(r):
        return None
    return np.median(p - r, axis=0)


def set_atom_name(line: str, new: str) -> str:
    return line[:12] + f"{new:>4s}" + line[16:]


def set_xyz(line: str, xyz: np.ndarray) -> str:
    return line[:30] + "".join(f"{v:8.3f}" for v in xyz) + line[54:]


def main() -> int:
    packed_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("bilayer_system.pdb")
    ref_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("receptor.pdb")
    if not packed_path.exists():
        print(f"ERROR: {packed_path} not found.", file=sys.stderr)
        return 1
    lines = packed_path.read_text().splitlines()

    caps = cap_residues(lines)
    if not caps:
        print("ERROR: no ACE/NME residues found in "
              f"{packed_path.name} -- the caps did not survive packing and must be restored, "
              "not repaired.", file=sys.stderr)
        return 1

    drops: set[int] = set()
    changes: list[str] = []
    appends: list[str] = []
    for cap in caps:
        want = template_heavies(cap["name"])
        have = cap["atoms"]
        heavies = [(n, li) for n, li in have if element_of(lines[li]) != "H"]
        stray_h = [(n, li) for n, li in have if element_of(lines[li]) == "H"]
        drops.update(li for _, li in stray_h)
        if stray_h:
            changes.append(f"{cap['name']}: dropped stray H ({', '.join(n for n, _ in stray_h)})"
                           " -- leap rebuilds them")

        used: set[int] = set()
        missing: list[str] = []
        for name in want:
            match = next(((n, li) for n, li in heavies if li not in used and n == name), None)
            if match is None:  # any unused heavy of the right element pairs with this slot
                match = next(((n, li) for n, li in heavies if li not in used
                              and element_of(lines[li]) == name[0]), None)
                if match is not None and match[0] != name:
                    changes.append(f"{cap['name']}: renamed {match[0]} -> {name}")
                    lines[match[1]] = set_atom_name(lines[match[1]], name)
            if match is None:
                missing.append(name)
            else:
                used.add(match[1])

        unknown = [(n, li) for n, li in heavies if li not in used]
        if unknown:
            print(f"ERROR: {cap['name']} contains heavy atoms matching no template slot: "
                  f"{[n for n, _ in unknown]} -- inspect manually.", file=sys.stderr)
            return 1

        if missing:
            if not ref_path.exists():
                print(f"ERROR: {cap['name']} is missing heavy atoms {missing} and "
                      f"{ref_path} is unavailable to restore them.", file=sys.stderr)
                return 1
            ref_lines = ref_path.read_text().splitlines()
            delta = ca_translation(lines, ref_lines)
            if delta is None:
                print("ERROR: cannot derive the packing translation (CA counts differ between "
                      "the packed and reference receptor); restore caps manually.",
                      file=sys.stderr)
                return 1
            # reference coordinates by element, from every reference copy of this cap
            ref_xyz: dict[str, np.ndarray] = {}
            for r in cap_residues(ref_lines):
                if r["name"] != cap["name"]:
                    continue
                for _, li in r["atoms"]:
                    line = ref_lines[li]
                    if element_of(line) == "H":
                        continue
                    ref_xyz.setdefault(element_of(line), np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
            first = cap["first"]
            for k, name in enumerate(missing):
                elem = name[0]
                if elem not in ref_xyz:
                    print(f"ERROR: no {elem}-element atom in the reference {cap['name']} to "
                          f"restore {name}.", file=sys.stderr)
                    return 1
                xyz = ref_xyz.pop(elem) + delta
                appends.append(ATOM_FMT.format(
                    serial=len(lines) + len(appends) + 1, name=name,
                    resname=cap["name"], chain=first[21], resid=int(first[22:26]),
                    icode=first[26], x=xyz[0], y=xyz[1], z=xyz[2], elem=elem))
                changes.append(f"{cap['name']}: restored {name} from reference "
                               f"(+{np.round(delta, 2).tolist()} translation)")

    if not changes:
        print(f"Caps already conform to the templates ({', '.join(r['name'] for r in caps)}); "
              "nothing to do.")
        return 0

    out = [l for i, l in enumerate(lines) if i not in drops]
    if appends:
        # restored atoms must land inside the coordinate block -- after TER/END they are
        # invisible to loadpdb
        term = next((k for k, l in enumerate(out) if l.startswith(("TER", "END"))), len(out))
        out = out[:term] + appends + out[term:]
    packed_path.write_text("\n".join(out) + "\n")
    for c in changes:
        print(c)
    print(f"Wrote normalised caps to {packed_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
