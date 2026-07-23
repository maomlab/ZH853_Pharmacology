#!/usr/bin/env python3
"""Physicochemical property panel for ZH853 and analogs (Objective 3).

Quantifies the drug-likeness liabilities of the four cyclic-peptide analogs (2D + best-effort
3D descriptors, beyond-rule-of-5 classification) and plots them in Ro5 / bRo5 property space.
Writes product/05.01.00_analog_properties_*.{csv,md,png}.

Run: ``python src/05.01.00_analog_properties.py``  (or ``make analogs``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from zh853mor import chem, paths  # noqa: E402


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    props = []
    for name in chem.ANALOGS:
        m = chem.mol(name)
        p = chem.descriptors_2d(name, m)
        p = chem.add_descriptors_3d(p, m)
        props.append(p)
        print(f"  {name}: MW={p.mw:.0f} TPSA={p.tpsa:.0f} HBD={p.hbd} bRo5={p.bro5_flags}")

    # ---- CSV + MD table ----
    cols = ["name", "mw", "tpsa", "clogp", "hbd", "hba", "rotb", "aromatic_rings",
            "fsp3", "formal_charge", "rgyr", "psa3d_frac", "intramol_hbonds", "bro5_flags"]
    paths.ensure_dir(paths.PRODUCT)
    csv = paths.PRODUCT / f"05.01.00_analog_properties_{today}.csv"
    lines = [",".join(cols)]
    for p in props:
        lines.append(",".join([
            p.name, f"{p.mw:.1f}", f"{p.tpsa:.0f}", f"{p.clogp:.2f}", str(p.hbd), str(p.hba),
            str(p.rotb), str(p.aromatic_rings), f"{p.fsp3:.2f}", str(p.formal_charge),
            f"{p.rgyr:.2f}", f"{p.psa3d_frac:.2f}", f"{p.intramol_hbonds:.0f}",
            ";".join(p.bro5_flags),
        ]))
    csv.write_text("\n".join(lines) + "\n")

    md = [f"# ZH853 analog property panel ({today})", "",
          "Reference thresholds: Lipinski Ro5 (MW<=500, HBD<=5, HBA<=10, cLogP<=5, TPSA<=140); "
          "oral bRo5 space (Doak/Kihlberg: MW up to ~1000, TPSA up to ~250, HBD up to ~10 for "
          "*chameleonic* macrocycles). intramol_hbonds = shielded backbone donors (higher = better "
          "passive permeability potential).", "",
          "| Analog | MW | TPSA | cLogP | HBD | HBA | RotB | ArRings | Fsp3 | Rgyr | polarSASA% | "
          "intraHB | bRo5 flags |", "|" + "---|" * 13]
    for p in props:
        md.append(f"| {p.name} | {p.mw:.0f} | {p.tpsa:.0f} | {p.clogp:.2f} | {p.hbd} | {p.hba} | "
                  f"{p.rotb} | {p.aromatic_rings} | {p.fsp3:.2f} | {p.rgyr:.2f} | "
                  f"{100 * p.psa3d_frac:.0f} | {p.intramol_hbonds:.0f} | "
                  f"{', '.join(p.bro5_flags)} |")
    md += ["", "## Diagnosis",
           "All four analogs sit in **beyond-rule-of-5** space: MW 714-810, TPSA 235-280, HBD 8-10 "
           "all exceed oral Ro5 limits. The dominant permeability liabilities are **high TPSA and "
           "HBD count** (exposed backbone/side-chain donors). cLogP is low-to-moderate (-0.15 to "
           "0.73). Cyclization and D-residues already confer protease resistance; the remaining gap "
           "is **passive permeability** (and, separately, **plasma half-life**). See "
           "`05.02.00_design_report` for modification strategies.",
           "", "ZH831 (MW 714, TPSA 235, HBD 8) is the least-liable starting point; ZH853 (the most "
           "potent, deposited ligand) is the most liable (extra Gly cap adds HBD/TPSA)."]
    (paths.PRODUCT / f"05.01.00_analog_properties_{today}.md").write_text("\n".join(md) + "\n")

    # ---- Ro5 / bRo5 property-space plot ----
    fig, ax = plt.subplots(figsize=(7, 5.5))
    xs = [p.mw for p in props]
    ys = [p.tpsa for p in props]
    ax.axvspan(0, 500, ymin=0, ymax=140 / 500, color="green", alpha=0.06)
    ax.axhline(140, color="green", ls="--", lw=1, label="Ro5 limit (TPSA 140)")
    ax.axvline(500, color="green", ls="--", lw=1, label="Ro5 limit (MW 500)")
    ax.axhline(250, color="orange", ls=":", lw=1, label="oral bRo5 (TPSA ~250)")
    ax.axvline(1000, color="orange", ls=":", lw=1, label="oral bRo5 (MW ~1000)")
    ax.scatter(xs, ys, s=90, c="crimson", zorder=5, edgecolors="k")
    for p in props:
        ax.annotate(p.name, (p.mw, p.tpsa), textcoords="offset points", xytext=(6, 5), fontsize=9)
    ax.set_xlabel("Molecular weight (Da)")
    ax.set_ylabel("TPSA (A^2)")
    ax.set_xlim(400, 1100)
    ax.set_ylim(100, 320)
    ax.set_title("ZH853 analogs in Ro5 / beyond-Ro5 property space")
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    png = paths.PRODUCT / f"05.01.00_analog_property_space_{today}.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"\nWrote {csv}, .md, and {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
