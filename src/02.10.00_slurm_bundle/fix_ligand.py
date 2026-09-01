#!/usr/bin/env python3
"""Reconcile the packed ligand residue with its mol2 template, before tleap.

Same class of problem fix_caps.py solves for ACE/NME, and it bites for the same reason: `loadpdb`
matches atoms to a unit template BY NAME. Two mismatches are guaranteed here.

  1. Residue name. The deposited ZH853 pose calls the ligand `L01`; antechamber is run with
     `-rn LIG`. A residue tleap has no unit for is loaded as untyped atoms and `saveamberparm`
     dies with "Atom .R<L01 ...>.A<C20> does not have a type".

  2. Atom names. The PDB carries the deposited names (C20, C30, C40 ...). SDF has no atom-name
     field, so antechamber GENERATES names (C1, C2, ...) when it builds the mol2 from
     <ligand>_prepared.sdf. Nothing makes those agree, and the failure is per-atom and total.

Both are fixed by mapping positionally rather than by name: the prepared SDF is written from the
deposited molecule with RDKit AddHs, which preserves heavy-atom order and appends hydrogens, so
heavy atom i of the packed residue is heavy atom i of the mol2. That assumption is CHECKED (count
and element sequence) and the script aborts rather than renaming atoms onto the wrong templates --
a silent mis-map would parameterize the ligand as a differently-connected molecule.

Hydrogens are not required in the PDB: tleap builds them from the template. The packed residue is
the 59 heavy atoms of the deposited pose.

Usage: python fix_ligand.py bilayer_system.pdb ZH853.mol2 --input-resname L01 [--resname LIG]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ELEM_COL = slice(76, 78)


def pdb_element(line: str) -> str:
    """Element for an ATOM/HETATM line: columns 77-78 if present, else inferred from the name."""
    e = line[ELEM_COL].strip()
    if e:
        return e.upper()
    name = line[12:16].strip()
    return (name[0] if not name[0].isdigit() else name[1:2]).upper()


def read_mol2_atoms(path: Path) -> list[tuple[str, str]]:
    """[(atom_name, element)] from a mol2 @<TRIPOS>ATOM block, in file order."""
    out, in_block = [], False
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            in_block = s.upper().startswith("@<TRIPOS>ATOM")
            continue
        if in_block and s:
            f = s.split()
            if len(f) >= 6:
                # SYBYL type is element[.hybridisation], e.g. C.3, N.am, O.co2
                out.append((f[1], f[5].split(".")[0].upper()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb", help="packed system PDB, edited in place")
    ap.add_argument("mol2", help="antechamber mol2 whose names/order are authoritative")
    ap.add_argument("--input-resname", required=True, help="what the ligand is called in the PDB")
    ap.add_argument("--resname", default="LIG", help="what to rename it to (must match the mol2 unit)")
    args = ap.parse_args()

    pdb = Path(args.pdb)
    lines = pdb.read_text().splitlines(keepends=True)
    mol2 = read_mol2_atoms(Path(args.mol2))
    if not mol2:
        print(f"ERROR: no @<TRIPOS>ATOM block in {args.mol2}.", file=sys.stderr)
        return 2

    want = args.input_resname.strip()
    idx = [i for i, l in enumerate(lines)
           if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == want]
    if not idx:
        present = sorted({l[17:20].strip() for l in lines if l.startswith(("ATOM", "HETATM"))})
        print(f"ERROR: no residue named '{want}' in {pdb}.", file=sys.stderr)
        print(f"  Residues present: {' '.join(present)}", file=sys.stderr)
        print("  The ligand did not survive packing -- packmol-memgen's reduce preprocessing can", file=sys.stderr)
        print("  drop HETATM records. Check that the staged complex still has it, and that the", file=sys.stderr)
        print("  build used the complex and not the receptor-only PDB.", file=sys.stderr)
        return 1

    # Heavy atoms only, on both sides: the deposited pose carries no hydrogens and tleap rebuilds
    # them from the template anyway.
    pdb_heavy = [i for i in idx if pdb_element(lines[i]) != "H"]
    mol2_heavy = [(n, e) for n, e in mol2 if e != "H"]

    if len(pdb_heavy) != len(mol2_heavy):
        print(f"ERROR: {len(pdb_heavy)} heavy atoms in the packed '{want}' residue but "
              f"{len(mol2_heavy)} in {args.mol2}.", file=sys.stderr)
        print("  These must be the same molecule in the same order. Either the mol2 was built from", file=sys.stderr)
        print("  a different SDF than the pose, or packing lost atoms. Not renaming anything.", file=sys.stderr)
        return 1

    mismatch = [(k, pdb_element(lines[i]), mol2_heavy[k][1])
                for k, i in enumerate(pdb_heavy) if pdb_element(lines[i]) != mol2_heavy[k][1]]
    if mismatch:
        print(f"ERROR: element sequence differs at {len(mismatch)} of {len(pdb_heavy)} positions; "
              "the positional mapping is not valid.", file=sys.stderr)
        for k, a, b in mismatch[:8]:
            print(f"  heavy atom {k + 1}: PDB {a} vs mol2 {b}", file=sys.stderr)
        print("  Refusing to rename: this would parameterize the ligand as a different molecule.",
              file=sys.stderr)
        return 1

    renamed = 0
    for k, i in enumerate(pdb_heavy):
        new = mol2_heavy[k][0]
        line = lines[i]
        if line[12:16].strip() != new:
            renamed += 1
        # PDB atom-name convention: 4-char field, 1-3 char names start at column 14
        field = f"{new:<4}" if len(new) >= 4 else f" {new:<3}"
        lines[i] = line[:12] + field[:4] + line[16:17] + f"{args.resname:>3}" + line[20:]

    # Any hydrogens that came through get dropped; tleap rebuilds them from the template and a
    # stray H with a non-template name is exactly what kills loadpdb.
    drop = {i for i in idx if pdb_element(lines[i]) == "H"}
    if drop:
        lines = [l for i, l in enumerate(lines) if i not in drop]

    pdb.write_text("".join(lines))
    print(f"fix_ligand: {len(pdb_heavy)} heavy atoms matched by element sequence; "
          f"{renamed} atom names rewritten from {args.mol2}, "
          f"residue {want} -> {args.resname}"
          + (f", {len(drop)} stray H dropped" if drop else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
