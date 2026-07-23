#!/usr/bin/env python3
"""Phase-1 QC of the deposited MOR-Gi-scFv16-ZH853 model.

Verifies residue numbering, reports receptor gaps, disulfides, and the ZH853
binding-pocket contact map. Writes a summary to product/. Reproduces the
interactive analyses that established the project's structural facts.

Run: ``python src/01.02.00_qc_structure.py``.
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from zh853mor import paths, structure  # noqa: E402


def main() -> int:
    u = structure.load(paths.CRYOEM_PDB)

    numbering = structure.verify_numbering(u)
    gaps = structure.receptor_gaps(u)
    ss = structure.disulfides(u)
    contacts = structure.pocket_contacts(u)

    lines: list[str] = []
    lines.append(f"# Structure QC -- {paths.CRYOEM_PDB.name}")
    lines.append(f"_Generated {date.today().isoformat()} by src/01.02.00_qc_structure.py_\n")

    n_ok = sum(numbering.values())
    lines.append(f"## Numbering: {n_ok}/{len(numbering)} canonical pocket residues match "
                 f"human OPRM1 (P35372)")
    for resid, ok in sorted(numbering.items()):
        resname, bw = structure.CANONICAL_POCKET[resid]
        lines.append(f"  - {resname}{resid} (BW {bw}): {'OK' if ok else 'MISMATCH'}")

    lines.append(f"\n## Receptor chain R gaps: {len(gaps) or 'none'}")
    for a, b in gaps:
        lines.append(f"  - {a} -> {b} (missing {b - a - 1})")

    lines.append("\n## Disulfides (chain R)")
    for i, j, d in ss:
        lines.append(f"  - CYS{i}-CYS{j}: {d:.2f} A")

    lines.append(f"\n## ZH853 pocket contacts (<=4.5 A): {len(contacts)} residues")
    for c in contacts:
        lines.append(f"  - {c.resname}{c.resid} ({c.n_atoms} atoms)")

    report = "\n".join(lines) + "\n"
    paths.ensure_dir(paths.PRODUCT)
    out = paths.PRODUCT / f"01.02.00_structure_qc_{date.today():%Y%m%d}.md"
    out.write_text(report)
    print(report)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
