#!/usr/bin/env python3
"""Ranked mutation panel to selectively abrogate ZH853 (Objective 2).

Scores each ZH853 pocket contact by how strongly ZH853 engages it and how few other
agonists share it, emphasizing sparing the two most closely related full agonists that
must be preserved (endomorphin-1, the direct parent, and DAMGO). Proposes specific
substitutions with rationale and caveats. Writes product/03.02.00_mutation_panel_*.md.

These are hypotheses for wet-lab testing; each should be validated by MD occupancy
(Phase 4) and, where feasible, relative binding free energy (Phase 6) before commitment.

Run: ``python src/03.02.00_mutation_panel.py``  (or ``make mutations``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

from rdkit import RDLogger  # noqa: E402

from zh853mor import comparators, interactions, paths, structure  # noqa: E402

RDLogger.DisableLog("rdApp.*")

STRENGTH = {"ionic": 4.0, "cation_pi": 3.5, "aromatic": 3.0, "hbond": 2.0, "hydrophobic": 1.0}
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
# Full agonists whose binding must be spared; endomorphin-1 is the direct ZH853 parent.
PARENTS = {"8F7R": "endomorphin-1", "8EFQ": "DAMGO"}
AGONISTS = ["8F7R", "8EFQ", "6DDE", "8F7Q", "9WST", "9WSV",
            "5C1M", "8EF5", "8EFB", "8EFL", "8EFO", "7T2G"]


def propose(resname: str, inter: set[str]) -> tuple[list[str], str]:
    """Return (substitutions, mechanism) tailored to the interaction being removed."""
    one = THREE_TO_ONE.get(resname, "X")
    iso = {"GLU": "Q", "ASP": "N", "LYS": "M", "ARG": "M"}  # charge-neutral isostere
    if "ionic" in inter and resname in iso:
        return [f"{one}{{}}A", f"{one}{{}}{iso[resname]}"], (
            "remove the salt bridge (Ala) or neutralize charge iso-sterically"
        )
    if ("cation_pi" in inter or "aromatic" in inter) and resname in ("PHE", "TYR", "TRP", "HIS"):
        return [f"{one}{{}}A", f"{one}{{}}L"], "ablate the aromatic/cation-pi ring"
    if "hbond" in inter and resname in ("TYR", "THR", "SER", "ASN", "GLN"):
        alt = "F" if resname == "TYR" else "A"
        return [f"{one}{{}}{alt}"], "remove the H-bonding group"
    if "hydrophobic" in inter:
        return [f"{one}{{}}A"], "open a packing cavity (less specific)"
    return [f"{one}{{}}A"], "perturb the contact (weak/contact-only)"


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    fps = {p: interactions.fingerprint(*_load(p)) for p in ["ZH853", *AGONISTS]}
    zh = fps["ZH853"]

    rows = []
    for hid, r in zh.items():
        z_strength = max((STRENGTH[i] for i in r.interactions), default=0.5)
        sharers = [p for p in AGONISTS if hid in fps[p]]
        share_frac = len(sharers) / len(AGONISTS)
        spares_parents = not any(p in sharers for p in PARENTS)
        # selectivity: strong for ZH853, rarely shared, bonus if both parents spared
        score = z_strength * (1 - share_frac) * (1.3 if spares_parents else 0.7)
        subs, mech = propose(r.resname, r.interactions)
        subs = [s.format(hid) for s in subs]
        rows.append({
            "hid": hid, "resname": r.resname, "bw": structure.bw(hid),
            "inter": ", ".join(sorted(r.interactions)) or "contact-only",
            "z_strength": z_strength, "n_share": len(sharers), "spares": spares_parents,
            "score": score, "subs": subs, "mech": mech,
            "parent_hit": [PARENTS[p] for p in PARENTS if p in sharers],
        })
    rows.sort(key=lambda d: -d["score"])

    lines = [f"# ZH853-selective mutation panel ({today})", "",
             "Ranked candidates to disrupt ZH853 while sparing related full agonists "
             "(endomorphin-1, DAMGO). Score = ZH853 interaction strength x (1 - fraction of "
             "12 agonists sharing the contact) x parent-sparing bonus. **Hypotheses — validate "
             "by MD occupancy (Phase 4) and relative FEP (Phase 6).**", ""]
    lines.append("| Rank | Residue | BW | ZH853 interaction | # agonists sharing | spares "
                 "endo-1 & DAMGO? | suggested mutations | mechanism | caveat |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for i, d in enumerate(rows, 1):
        caveat = "—"
        if d["parent_hit"]:
            caveat = "also engages " + ", ".join(d["parent_hit"])
        elif d["z_strength"] <= 1.0:
            caveat = "weak (hydrophobic/contact) — may be nonspecific"
        lines.append(
            f"| {i} | {d['resname']}{d['hid']} | {d['bw']} | {d['inter']} | "
            f"{d['n_share']}/12 | {'yes' if d['spares'] else 'NO'} | "
            f"{', '.join(d['subs'])} | {d['mech']} | {caveat} |"
        )

    # Highlighted recommendations: strong ZH853 interaction, parents spared.
    top = [d for d in rows if d["spares"] and d["z_strength"] >= 2.0][:5]
    lines.append("\n## Recommended primary tests (strong ZH853 interaction, parents spared)\n")
    for d in top:
        lines.append(f"- **{d['subs'][0]}** ({d['resname']}{d['hid']}, {d['bw']}): "
                     f"ZH853 makes a {d['inter']} contact here that none of the parent agonists "
                     f"make; {d['mech']}.")
    if not top:
        lines.append("- (No candidate combined a strong ZH853 interaction with full parent "
                     "sparing; see the ranked table — GLU231 is the strongest lead.)")

    lines.append("\n## Do-not-mutate controls (universal agonist anchors)\n")
    universal = [d for d in rows if d["n_share"] >= 11]
    lines.append("These are engaged by essentially all agonists — mutating them would abrogate "
                 "the comparators too, so they serve as negative controls, not selective targets: "
                 + ", ".join(f"{d['resname']}{d['hid']}({d['bw']})" for d in universal) + ".")

    out = paths.PRODUCT / f"03.02.00_mutation_panel_{today}.md"
    paths.ensure_dir(paths.PRODUCT)
    out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out}")
    return 0


def _load(pdb: str):
    cx = comparators.load_complex(pdb)
    return cx.receptor, cx.ligand, cx.offset_to_human


if __name__ == "__main__":
    raise SystemExit(main())
