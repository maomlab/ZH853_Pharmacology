#!/usr/bin/env python3
"""Assess prep needs and split the complex into components (Phase 2).

Reports missing atoms, chain termini, His tautomer sites, the Na+ pocket, and disulfides;
writes per-component PDBs (receptor+ligand, Gi, scFv16) to intermediate/ for downstream prep.

Run: ``python src/02.01.00_assess_and_split.py``  (or ``make prep-complex-split``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from zh853mor import paths, prep, structure  # noqa: E402


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    u = structure.load(paths.CRYOEM_PDB)

    incomplete = prep.incomplete_residues(u)
    termini = prep.chain_termini(u)
    his = prep.histidines(u, "R")
    ss = structure.disulfides(u)
    mem = prep.membrane_frame(u, "R")

    lines = [f"# MD prep assessment ({today})", "",
             f"Source: `{paths.CRYOEM_PDB.name}` (3.5 A cryo-EM).", ""]

    lines.append("## Chain inventory")
    for seg, role in prep.CHAIN_ROLES.items():
        ca = u.select_atoms(f"segid {seg} and protein and name CA")
        n = len(ca)
        span = f"{int(ca.resids.min())}-{int(ca.resids.max())}" if n else "ligand/het"
        lines.append(f"- **{seg}** — {role}: {n} residues ({span})")

    lines.append(f"\n## Incomplete residues (truncated sidechains): {len(incomplete)}")
    if incomplete:
        lines.append("Need sidechain rebuild (pdbfixer/Modeller). Common at 3.5 A on the surface.")
        for i in incomplete[:60]:
            lines.append(f"- {i.chain}/{i.resname}{i.resid}: {i.detail}")
        if len(incomplete) > 60:
            lines.append(f"- … and {len(incomplete) - 60} more")
    else:
        lines.append("None — all modeled residues have complete heavy-atom sets.")

    lines.append("\n## Chain termini (cap ACE/NME if truncated vs the full sequence)")
    for t in termini:
        lines.append(f"- {t.chain}/{t.resname}{t.resid}: {t.detail}")

    lines.append(f"\n## Histidine tautomers to assign (chain R): {len(his)}")
    lines.append("Assign HID/HIE/HIP from local H-bonding; **H299 (6.52)** contacts the ligand and "
                 "**H321 (7.36)** is a ZH853-distinctive contact (Objective 1) — assign these carefully.")
    for h in his:
        lines.append(f"- HIS{h.resid}")

    lines.append("\n## Disulfides (enforce explicitly)")
    for i, j, d in ss:
        lines.append(f"- CYS{i}-CYS{j} ({d:.2f} A) — conserved ECL2-TM3 bridge")

    lines.append("\n## Na+ pocket / protonation-sensitive sites")
    d250 = u.select_atoms(f"segid R and resid {prep.D250_SODIUM} and name CA")
    d250_id = f"{d250[0].resname}{prep.D250_SODIUM}" if len(d250) else "absent"
    lines.append(f"- **D2.50 = {d250_id}**: allosteric sodium site — protonation ambiguous. "
                 "Build **parallel protonated/deprotonated systems** (SPECIFICATION D-11).")
    dry = u.select_atoms(f"segid R and resid {prep.DRY_MOTIF[0]}-{prep.DRY_MOTIF[2]} and name CA")
    dry_s = "-".join(f"{a.resname}{a.resid}" for a in dry)
    lines.append(f"- DRY motif (3.49-3.51): {dry_s} — active-state ionic-lock region; check E/D3.49.")

    lines.append("\n## Membrane frame (from TM bundle + modeled cholesterol)")
    lines.append(f"- Membrane normal (receptor principal axis): "
                 f"[{', '.join(f'{x:.2f}' for x in mem['normal'])}]")
    lines.append(f"- Modeled cholesterol atoms: {mem['n_cholesterol_atoms']}")
    if "cholesterol_span_along_normal" in mem:
        lines.append(f"- Cholesterol span along normal: "
                     f"{mem['cholesterol_span_along_normal']:.1f} A (bilayer hydrophobic thickness cue)")
    lines.append(f"- Receptor span along normal: {mem['receptor_span_along_normal']:.1f} A")
    lines.append("- Use these to pre-orient the receptor (PPM/OPM) before membrane building.")

    paths.ensure_dir(paths.PRODUCT)
    (paths.PRODUCT / f"02.01.00_prep_assessment_{today}.md").write_text("\n".join(lines) + "\n")

    # ---- split into components ----
    inter = paths.ensure_dir(paths.INTERMEDIATE / "02.01.00_components")
    comps = {
        "receptor_ligand": "segid R E",       # MOR + cholesterol + ZH853
        "gi_heterotrimer": "segid A B C",
        "scfv16": "segid D",
    }
    for name, sel in comps.items():
        ag = u.select_atoms(sel)
        out = inter / f"{name}.pdb"
        ag.write(str(out))
        lines_written = len(ag)
        print(f"  wrote {out.name}: {lines_written} atoms")

    print(f"\nAssessment: {len(incomplete)} incomplete, {len(his)} His, {len(ss)} disulfides.")
    print(f"Report: product/02.01.00_prep_assessment_{today}.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
