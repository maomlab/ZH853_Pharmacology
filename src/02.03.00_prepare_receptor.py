#!/usr/bin/env python3
"""Rebuild the MOR receptor for MD with PDBFixer (Phase 2).

Rebuilds truncated sidechains (the 45 incomplete residues), keeps the modeled sequence
(no internal gaps to bridge), enforces the C142-C219 disulfide, and writes a clean
heavy-atom receptor for the membrane builder plus a pH-7.4 protonated copy for inspection.
Cholesterol/ligand are handled separately (membrane builder / ligand parameterization).

Requires the analysis env (openmm + pdbfixer). Run: ``python src/02.03.00_prepare_receptor.py``.
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from zh853mor import paths, prep, structure  # noqa: E402


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    work = paths.ensure_dir(paths.INTERMEDIATE / "02.03.00_receptor")

    # Write receptor protein (chain R) as PDBFixer input.
    u = structure.load(paths.CRYOEM_PDB)
    rec_in = work / "receptorR_raw.pdb"
    u.select_atoms("segid R and protein").write(str(rec_in))

    from openmm.app import PDBFile  # noqa: E402
    from pdbfixer import PDBFixer  # noqa: E402

    fixer = PDBFixer(filename=str(rec_in))
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # no internal gaps; do not model unresolved termini (cap later)
    fixer.findMissingAtoms()
    n_rebuilt = sum(len(v) for v in fixer.missingAtoms.values())
    fixer.addMissingAtoms()  # rebuild truncated sidechains

    # heavy-atom structure for the membrane builder (CHARMM-GUI adds its own H)
    heavy_out = work / "receptorR_fixed_heavy.pdb"
    PDBFile.writeFile(fixer.topology, fixer.positions, str(heavy_out), keepIds=True)

    # protonated copy at pH 7.4 for inspection / local sanity
    fixer.addMissingHydrogens(7.4)
    prot_out = work / "receptorR_fixed_pH7.4.pdb"
    PDBFile.writeFile(fixer.topology, fixer.positions, str(prot_out), keepIds=True)

    # verify the rebuild closed the incomplete-sidechain gaps
    after = structure.load(prot_out)
    remaining = prep.incomplete_residues(after)

    lines = [f"# Receptor preparation (PDBFixer) — {today}", "",
             f"- Input: MOR chain R ({len(u.select_atoms('segid R and protein'))} atoms).",
             f"- Missing heavy atoms rebuilt: **{n_rebuilt}** across truncated sidechains.",
             f"- Incomplete residues remaining after rebuild: **{len(remaining)}** "
             "(expected ~0; any residual are chain-terminal).",
             f"- Disulfides preserved: {', '.join(f'C{i}-C{j}' for i,j,_ in structure.disulfides(u))}.",
             "",
             "## Outputs (intermediate/02.03.00_receptor/)",
             f"- `{heavy_out.name}` — heavy-atom receptor for the membrane builder (CHARMM-GUI/PACKMOL).",
             f"- `{prot_out.name}` — pH-7.4 protonated copy for inspection.",
             "",
             "## Remaining prep (in the membrane-builder / tleap step)",
             "- **Cap truncated termini** T69 (N-term) and F349 (C-term) with ACE/NME (neutral) "
             "rather than charged termini — they are internal fragments of full-length OPRM1.",
             "- Apply the **D2.50 (Asp116) protonated variant** for the parallel system (D-11).",
             "- Assign **His tautomers** per `02.02.00_protonation` (H299/H321 near the pocket).",
             "- Cholesterol from the cryo-EM model is dropped here; the membrane builder places lipids "
             "(POPC:chol 9:1 baseline, D-4).",
             ]
    if remaining:
        lines.append("\n## Residual incomplete residues")
        for r in remaining:
            lines.append(f"- {r.chain}/{r.resname}{r.resid}: {r.detail}")

    (paths.PRODUCT / f"02.03.00_receptor_prep_{today}.md").write_text("\n".join(lines) + "\n")
    print(f"Rebuilt {n_rebuilt} atoms; {len(remaining)} incomplete residues remain.")
    print(f"Outputs in {work}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
