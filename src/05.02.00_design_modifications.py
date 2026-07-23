#!/usr/bin/env python3
"""Structure-guided drug-likeness modification design for ZH853 (Objective 3).

Combines (a) a structure-based burial analysis of the bound ZH853 ligand — separating the
buried pharmacophore from solvent-exposed, derivatizable positions — with (b) computational
enumeration of modification strategies (N-methylation, lipidation, halogenation), reporting
predicted property changes on the two PK axes (passive permeability vs plasma half-life).

Writes product/05.02.00_design_report_*.md and a property-shift figure.

Run: ``python src/05.02.00_design_modifications.py``  (or ``make design``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from zh853mor import chem, comparators, paths  # noqa: E402


def burial_analysis() -> tuple[int, int, list[str]]:
    """Per-atom burial of the bound ZH853: return (n_buried, n_exposed, exposed_atom_names)."""
    cx = comparators.load_complex("ZH853")
    lig = cx.ligand
    rec = cx.receptor
    rpos = rec.positions
    buried, exposed = 0, []
    for a in lig.atoms:
        d = float(np.linalg.norm(rpos - a.position, axis=1).min())
        if d <= 4.5:
            buried += 1
        else:
            exposed.append(f"{a.name}({d:.1f})")
    return buried, len(exposed), exposed


def variant_row(name: str, m, axis: str, risk: str, note: str) -> dict:
    p = chem.descriptors_2d(name, m)
    return {"name": name, "p": p, "axis": axis, "risk": risk, "note": note}


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    parent = chem.mol("ZH853")
    p0 = chem.descriptors_2d("ZH853", parent)

    n_buried, n_exposed, exposed = burial_analysis()

    rows = [variant_row("ZH853 (parent)", parent, "-", "-", "deposited ligand")]

    # ---- N-methylation (permeability) ----
    sites = chem.n_methyl_sites(parent)
    per_site = []
    for i, s in enumerate(sites):
        v = chem.add_n_methyls(parent, [s])
        pv = chem.descriptors_2d(f"NMe#{i + 1}", v)
        per_site.append((pv.tpsa, pv.hbd))
    rows.append(variant_row("tri-N-methyl (3 backbone NMe)", chem.add_n_methyls(parent, sites[:3]),
                            "permeability", "moderate (backbone H-bonds)",
                            "masks 3 donors; rigidifies; boosts passive permeability & protease "
                            "stability (cyclosporine-like)"))
    rows.append(variant_row("hexa-N-methyl (all backbone NMe)", chem.add_n_methyls(parent, sites),
                            "permeability", "high (may perturb bound conformation)",
                            "maximal HBD/TPSA reduction; likely potency trade-off — test subsets"))

    # ---- Halogenation (metabolic stability) ----
    f_smiles = chem.ANALOGS["ZH853"][0].replace("Cc4ccccc4", "Cc4ccc(F)cc4")
    rows.append(variant_row("4-F-Phe (para-fluoro-Phe)", chem.mol(f_smiles),
                            "metabolic", "low (conservative)",
                            "blocks Phe para-hydroxylation; minimal property change; potency-neutral"))

    # ---- Lipidation (half-life) ----
    for lname, (lsmi, desc) in chem.LINKERS.items():
        try:
            conj = chem.conjugate_at_primary_amide(parent, lsmi)
            rows.append(variant_row(f"C-term {lname}", conj, "half-life",
                                    "site-dependent (use exposed cap)", desc))
        except ValueError as exc:
            print(f"  {lname}: skipped ({exc})")

    # ---- report ----
    md = [f"# ZH853 structure-guided modification design ({today})", "",
          "## Structure-based derivatization map",
          f"Of {n_buried + n_exposed} ZH853 heavy atoms, **{n_buried} are buried** (<=4.5 A from the "
          f"receptor — the pharmacophore: the Tyr1 amine/D149 salt bridge and the aromatic pocket "
          f"contacts) and **{n_exposed} are solvent-exposed** and therefore candidate derivatization "
          "handles for half-life-extending conjugation without disrupting binding.",
          f"\nExposed atoms (name, min dist to receptor A): {', '.join(exposed)}.",
          "\n> Design rule: keep the buried Tyr1 amine + aromatic message intact; place lipid/PEG "
          "conjugation on the exposed C-terminal cap region. Permeability edits (N-methylation) should "
          "avoid backbone amides that H-bond the receptor (D149/Y328/E231 network from Objective 1).",
          "", "## Two PK axes (distinct goals)",
          "1. **Passive permeability / oral-CNS exposure** — lower TPSA & HBD: N-methylation, "
          "H-bond-donor reduction, halogenation. Moves the molecule toward the chameleonic bRo5 window.",
          "2. **Plasma half-life / duration** — albumin-binding **lipidation** (GLP-1/semaglutide "
          "strategy): a fatty di-acid + gamma-Glu + PEG linker. This *raises* polarity/size (keeps it "
          "injectable) but extends half-life from minutes to ~days. Orthogonal to axis 1.",
          "", "## Enumerated modifications (predicted properties)", "",
          "| Variant | PK axis | MW | TPSA | HBD | cLogP | dTPSA | dHBD | pharmacophore risk | rationale |",
          "|" + "---|" * 10]
    for r in rows:
        p = r["p"]
        md.append(f"| {r['name']} | {r['axis']} | {p.mw:.0f} | {p.tpsa:.0f} | {p.hbd} | "
                  f"{p.clogp:.2f} | {p.tpsa - p0.tpsa:+.0f} | {p.hbd - p0.hbd:+d} | "
                  f"{r['risk']} | {r['note']} |")

    per_site_dtpsa = round(sum(p0.tpsa - t for t, _ in per_site) / len(per_site))
    md += ["", "## Findings & recommendations",
           f"- **N-methylation is the highest-value permeability lever.** A single backbone N-methyl "
           f"removes one donor and ~{per_site_dtpsa} A^2 TPSA (avg over the {len(per_site)} sites); "
           "three together bring TPSA/HBD toward the chameleonic oral-bRo5 window. Because "
           "several backbone amides H-bond the receptor (Objective 1), **N-methylate solvent-facing "
           "amides only** — prioritize a 2-3 site subset and validate potency by relative FEP (Phase 6).",
           "- **Lipidation (semaglutide-style) is the half-life strategy**, placed on the exposed "
           "C-terminal cap identified above; expect injectable, long-acting PK rather than improved "
           "permeability. A C16 palmitoyl is the minimal variant; the gamma-Glu-2xAEEA-C18-diacid is the "
           "validated albumin-binder.",
           "- **4-F-Phe** is a low-risk metabolic-stability tweak (blocks para-hydroxylation), "
           "essentially property-neutral.",
           "- **Start from ZH831** (least liable analog) for permeability-focused series; keep ZH853 for "
           "potency-focused/half-life series.",
           "", "> All predictions are 2D/3D descriptor estimates. Permeability of macrocyclic peptides is "
           "conformation-dependent ('molecular chameleons'); confirm with 3D-PSA over conformer "
           "ensembles and, ideally, PAMPA/Caco-2. Any pharmacophore-adjacent edit must be checked for "
           "potency by relative FEP (Objective 4)."]
    out = paths.PRODUCT / f"05.02.00_design_report_{today}.md"
    out.write_text("\n".join(md) + "\n")

    # ---- figure: property shifts on the TPSA-vs-cLogP plane ----
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    colors = {"-": "k", "permeability": "tab:blue", "metabolic": "tab:green",
              "half-life": "tab:red"}
    for r in rows:
        p = r["p"]
        ax.scatter(p.clogp, p.tpsa, s=90, c=colors[r["axis"]], edgecolors="k", zorder=5)
        ax.annotate(r["name"].split(" (")[0], (p.clogp, p.tpsa),
                    textcoords="offset points", xytext=(6, 4), fontsize=7)
    ax.axhline(140, color="green", ls="--", lw=1)
    ax.axhline(250, color="orange", ls=":", lw=1)
    ax.set_xlabel("cLogP")
    ax.set_ylabel("TPSA (A^2)")
    ax.set_title("ZH853 modification design: property shifts by PK strategy")
    handles = [plt.Line2D([], [], marker="o", ls="", mfc=c, mec="k", label=k)
               for k, c in colors.items() if k != "-"]
    ax.legend(handles=handles, fontsize=8, title="PK axis")
    fig.tight_layout()
    png = paths.PRODUCT / f"05.02.00_design_property_shifts_{today}.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"Buried {n_buried} / exposed {n_exposed} ligand atoms.")
    print(f"Wrote {out} and {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
