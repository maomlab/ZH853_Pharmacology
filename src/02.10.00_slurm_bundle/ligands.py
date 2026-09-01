#!/usr/bin/env python3
"""Registry of the systems this bundle can build: apo plus the four ZH cyclic peptides.

One source of truth, queried by 01_build_system.sh (`python ligands.py --field ...`) and imported
by the QC scripts, so the shell and the Python agree on what a ligand is called and where its
inputs live.

The built system always uses resname LIG regardless of what the input PDB called the ligand
(the deposited ZH853 pose calls it L01). fix_ligand.py does that rename, so `--lig LIG` is
correct for every non-apo build and per-ligand identity is recorded in system.json instead.

Usage:
    python ligands.py --list
    python ligands.py --field resname --ligand ZH853
    python ligands.py --check ZH853 --repo /path/to/repo
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

RESNAME = "LIG"  # resname in the built system, for every ligand

# net_charge is the formal charge antechamber is told to fit to (-nc); all four carry the
# protonated Tyr1 alpha-amine that salt-bridges D3.32 (Asp149), hence +1.
LIGANDS: dict[str, dict] = {
    "apo": {
        "description": "receptor only, no ligand",
        "net_charge": None, "input_resname": None, "sdf": None, "complex": None,
    },
    "ZH853": {
        "description": "Tyr-cyclo[D-Lys-Trp-Phe-Glu]-Gly-NH2 (deposited ligand)",
        "net_charge": 1, "input_resname": "L01",
        "sdf": "intermediate/02.04.00_ligand/ZH853_prepared.sdf",
        # the deposited pose predates the multi-ligand layout and has no name in its filename
        "complex": ["intermediate/02.05.00_oriented/complex_ZH853_oriented.pdb",
                    "intermediate/02.05.00_oriented/complex_oriented.pdb"],
    },
    "ZH850": {
        "description": "Tyr-cyclo[D-Lys-Trp-Phe-Glu]-NH2 (analog 1; ZH853 minus the Gly-NH2 tail)",
        "net_charge": 1, "input_resname": "L01",
        "sdf": "intermediate/02.04.00_ligand/ZH850_prepared.sdf",
        "complex": ["intermediate/02.05.00_oriented/complex_ZH850_oriented.pdb"],
    },
    "ZH831": {
        "description": "Tyr-cyclo[D-Glu-Phe-Phe-Lys]-NH2 (analog 2; Trp->Phe, reversed cycle)",
        "net_charge": 1, "input_resname": "L01",
        "sdf": "intermediate/02.04.00_ligand/ZH831_prepared.sdf",
        "complex": ["intermediate/02.05.00_oriented/complex_ZH831_oriented.pdb"],
    },
    "ZH809": {
        "description": "Tyr-cyclo[D-Lys-Trp-Phe-Asp]-NH2 (analog 3; Glu->Asp, one CH2 shorter)",
        "net_charge": 1, "input_resname": "L01",
        "sdf": "intermediate/02.04.00_ligand/ZH809_prepared.sdf",
        "complex": ["intermediate/02.05.00_oriented/complex_ZH809_oriented.pdb"],
    },
}

NAMES = list(LIGANDS)
LIGAND_NAMES = [n for n in NAMES if n != "apo"]


def get(name: str) -> dict:
    if name not in LIGANDS:
        raise SystemExit(f"ERROR: unknown ligand '{name}'. Known: {' '.join(NAMES)}")
    return LIGANDS[name]


def resolve_complex(name: str, repo: Path) -> Path | None:
    """First existing oriented receptor+ligand complex for `name`, or None."""
    for rel in get(name)["complex"] or []:
        p = repo / rel
        if p.exists():
            return p
    return None


def missing_inputs(name: str, repo: Path) -> list[str]:
    """Which prep artefacts a build of `name` still needs. Empty means ready."""
    spec = get(name)
    if name == "apo":
        return []
    missing = []
    if resolve_complex(name, repo) is None:
        missing.append("oriented receptor+ligand complex: " + " or ".join(spec["complex"]))
    if not (repo / spec["sdf"]).exists():
        missing.append(f"prepared ligand SDF: {spec['sdf']}")
    return missing


# The finalised receptor (ACE/NME caps + named His tautomers from 02.03.00) and the oriented
# complex are DIFFERENT vintages: complex_oriented.pdb was written from the raw deposited complex
# and still has bare HIS and no caps, so building from it directly trips the staleness guard --
# and would otherwise have produced charged termini and default tautomers. Both files carry the
# same OPM transform though, so the fix is to graft the ligand onto the finalised receptor rather
# than to re-run the orientation chain. Verified at build time, never assumed.
GRAFT_TOL_A = 0.10
RECEPTOR_REL = "intermediate/02.05.00_oriented/receptorR_oriented.pdb"


def _ca(path: Path) -> dict[int, tuple[float, float, float]]:
    out = {}
    for l in path.read_text().splitlines():
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            out[int(l[22:26])] = tuple(float(l[30 + 8 * i:38 + 8 * i]) for i in range(3))
    return out


def graft(name: str, repo: Path, out: Path) -> int:
    """Write receptor.pdb: the finalised receptor, plus the ligand for a non-apo build."""
    receptor = repo / RECEPTOR_REL
    if not receptor.exists():
        print(f"ERROR: {receptor} not found.", file=sys.stderr)
        return 1
    rec_lines = [l for l in receptor.read_text().splitlines() if not l.startswith("END")]
    if name == "apo":
        out.write_text("\n".join(rec_lines) + "\nEND\n")
        print(f"staged {out.name}: apo receptor, {sum(1 for l in rec_lines if l.startswith(('ATOM', 'HETATM')))} atoms")
        return 0

    cx = resolve_complex(name, repo)
    if cx is None:
        print(f"ERROR: no oriented complex for {name}.", file=sys.stderr)
        return 1

    # Same OPM frame? If the two files were regenerated out of step the ligand would be placed
    # against a receptor it was never docked to, and nothing downstream would notice.
    a, b = _ca(receptor), _ca(cx)
    shared = sorted(set(a) & set(b))
    if not shared:
        print(f"ERROR: {receptor.name} and {cx.name} share no CA residue numbering.", file=sys.stderr)
        return 1
    worst = max(sum((a[r][i] - b[r][i]) ** 2 for i in range(3)) ** 0.5 for r in shared)
    if worst > GRAFT_TOL_A:
        print(f"ERROR: {receptor.name} and {cx.name} are not in the same frame "
              f"(max CA displacement {worst:.3f} A over {len(shared)} residues, "
              f"tolerance {GRAFT_TOL_A} A).", file=sys.stderr)
        print("  The ligand pose cannot be transferred onto this receptor. Re-run the orientation", file=sys.stderr)
        print("  chain so both come from the same superposition:  make prep-receptor-rebuild prep-receptor-orient", file=sys.stderr)
        return 1

    resn = get(name)["input_resname"]
    lig = [l for l in cx.read_text().splitlines()
           if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == resn]
    if not lig:
        print(f"ERROR: no '{resn}' residue in {cx}.", file=sys.stderr)
        return 1
    out.write_text("\n".join(rec_lines + lig) + "\nEND\n")
    print(f"staged {out.name}: finalised receptor + {len(lig)} {resn} ligand atoms from "
          f"{cx.name} (frames agree to {worst:.3f} A over {len(shared)} CA)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="print the buildable system names")
    ap.add_argument("--describe", action="store_true", help="print name + description per line")
    ap.add_argument("--field", help="print one field for --ligand (resname/net_charge/sdf/...)")
    ap.add_argument("--ligand")
    ap.add_argument("--check", metavar="LIGAND", help="exit 1 listing what prep inputs are missing")
    ap.add_argument("--graft", metavar="LIGAND",
                    help="write the staged receptor(+ligand) PDB for LIGAND to --out")
    ap.add_argument("--out", default="receptor.pdb")
    ap.add_argument("--repo", default=".", help="repository root, for --check/--field complex")
    args = ap.parse_args()

    if args.list:
        print(" ".join(NAMES))
        return 0
    if args.describe:
        for n in NAMES:
            print(f"{n}\t{LIGANDS[n]['description']}")
        return 0
    if args.check:
        repo = Path(args.repo).resolve()
        missing = missing_inputs(args.check, repo)
        if missing:
            print(f"ERROR: cannot build '{args.check}' -- missing prep inputs:", file=sys.stderr)
            for m in missing:
                print(f"  - {m}", file=sys.stderr)
            print("  Only ZH853 has a deposited pose; the analogs need one generated before they",
                  file=sys.stderr)
            print("  can be simulated (see README 'Ligands and the apo system').", file=sys.stderr)
            return 1
        return 0
    if args.graft:
        return graft(args.graft, Path(args.repo).resolve(), Path(args.out))
    if args.field:
        spec = get(args.ligand or "apo")
        if args.field == "resname":
            print(RESNAME if (args.ligand or "apo") != "apo" else "")
        elif args.field == "complex":
            p = resolve_complex(args.ligand, Path(args.repo).resolve())
            print(p or "")
        else:
            v = spec.get(args.field)
            print("" if v is None else (json.dumps(v) if isinstance(v, list) else v))
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
