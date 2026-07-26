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

    # raw heavy-atom structure (pre-finalisation; kept for provenance)
    raw_heavy = work / "receptorR_fixed_heavy_raw.pdb"
    PDBFile.writeFile(fixer.topology, fixer.positions, str(raw_heavy), keepIds=True)

    # protonated copy at pH 7.4 -- also the source of the His tautomer assignment below
    fixer.addMissingHydrogens(7.4)
    prot_out = work / "receptorR_fixed_pH7.4.pdb"
    PDBFile.writeFile(fixer.topology, fixer.positions, str(prot_out), keepIds=True)

    # --- finalise the receptor so no downstream default can override these decisions ----------
    # 1. His tautomers: OpenMM already chose HID/HIE/HIP from the local H-bond network when it
    #    added hydrogens, but it does not rename the residue -- and tleap maps a residue still
    #    named HIS to HIE, silently. Recover the choice and write it into the residue name.
    taut = prep.his_tautomers(prot_out)
    named = work / "receptorR_fixed_heavy_his.pdb"
    with open(raw_heavy) as fh, open(named, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "HIS":
                line = line[:17] + f"{taut.get(int(line[22:26]), 'HIS'):>3s}" + line[20:]
            out.write(line)

    # 2. Neutral ACE/NME caps on the truncated termini (D-15).
    heavy_out = work / "receptorR_fixed_heavy.pdb"
    caps = prep.cap_termini(named, heavy_out)

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
             "## His tautomers (assigned here, not left to a downstream default)",
             "OpenMM/PDBFixer picks HID/HIE/HIP from the local hydrogen-bond network while adding",
             "hydrogens, but leaves the residue named `HIS` — and tleap maps `HIS` to **HIE**",
             "regardless. The choice is therefore written into the residue name here.",
             "",
             "| Residue | BW | assigned |",
             "|---|---|---|",
             *[f"| HIS{r} | {({299: 'H6.52', 321: 'H7.36'}).get(r, '—')} | **{t}** |"
               for r, t in sorted(taut.items())],
             "",
             "H299 (H6.52) and H321 (H7.36) line the orthosteric pocket, so verify these two against",
             "the ZH853 contact map before production (`02.02.00_protonation`).",
             "",
             "## Terminal caps",
             f"- ACE {caps['ace_resid']} / NME {caps['nme_resid']} added; C-terminal OXT dropped.",
             f"- Backbone torsions chosen by clash scan: phi {caps['phi']:.0f}°, psi {caps['psi']:.0f}°; "
             f"closest cap-to-protein contact **{caps['min_contact']} Å**.",
             "- Rationale: 69–349 is an internal fragment of full-length OPRM1 (400 aa), so charged "
             "termini would add two formal charges the real receptor does not have (D-15).",
             "",
             "## Outputs (intermediate/02.03.00_receptor/)",
             f"- `{heavy_out.name}` — **finalised** heavy-atom receptor (tautomers named, termini "
             "capped) for orientation (02.05.00) and the membrane builder.",
             f"- `{raw_heavy.name}` — pre-finalisation copy, for provenance.",
             f"- `{prot_out.name}` — pH-7.4 protonated copy; source of the tautomer assignment.",
             "",
             "## Remaining prep (in the build step)",
             "- The **D2.50 (Asp116) protonated variant** is a build-time residue rename "
             "(`D250=ASH ./01_build_system.sh`), not a separate prep run — the geometry is identical "
             "(D-11).",
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
