#!/usr/bin/env python3
"""Membrane-placement determination figure (Phase 2 / manuscript).

Side view (lateral axis vs membrane normal z) of the oriented complex showing how the modeled
cholesterols fix the bilayer: their mean z sets the midplane (z=0) and their span sets the
hydrophobic thickness. Overlays receptor Cα and ZH853. Vector PDF for the manuscript.

Run: python src/02.06.00_membrane_placement.py   (after `make prep-orient`).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from zh853mor import paths, structure  # noqa: E402

ORIENTED = paths.INTERMEDIATE / "02.05.00_oriented" / "complex_oriented.pdb"


def main() -> int:
    if not ORIENTED.exists():
        print(f"ERROR: {ORIENTED} not found -- run `make prep-orient` first.", file=sys.stderr)
        return 1
    u = structure.load(ORIENTED)
    rec = u.select_atoms("protein and name CA")
    clr = u.select_atoms("resname CLR")
    lig = u.select_atoms("resname L01")

    half = float(np.ptp(clr.positions[:, 2])) / 2.0  # hydrophobic half-thickness (cholesterol span/2)
    today = f"{date.today():%Y%m%d}"

    fig, (ax, axh) = plt.subplots(1, 2, figsize=(8.5, 6), gridspec_kw={"width_ratios": [3, 1]},
                                  sharey=True)

    # side view: lateral x vs membrane normal z
    ax.scatter(rec.positions[:, 0], rec.positions[:, 2], s=10, c="#9aa0a6", label="receptor Cα",
               zorder=2)
    ax.scatter(clr.positions[:, 0], clr.positions[:, 2], s=14, c="#e8820c", label="cholesterol",
               zorder=3)
    ax.scatter(lig.positions[:, 0], lig.positions[:, 2], s=16, c="#2ca02c", label="ZH853", zorder=4)
    ax.axhspan(-half, half, color="#f2c811", alpha=0.15, zorder=0)
    ax.axhline(0, color="#333", lw=1.2, ls="-", zorder=1)
    ax.axhline(half, color="#e8820c", lw=1, ls="--", zorder=1)
    ax.axhline(-half, color="#e8820c", lw=1, ls="--", zorder=1)
    ax.text(ax.get_xlim()[1], 1.5, " midplane (z=0)", va="bottom", ha="right", fontsize=8)
    ax.text(ax.get_xlim()[1], half + 0.5, f" ±{half:.0f} Å (hydrophobic)", va="bottom", ha="right",
            fontsize=8, color="#e8820c")
    ax.annotate("extracellular", xy=(ax.get_xlim()[0], 26), fontsize=8, style="italic")
    ax.annotate("intracellular", xy=(ax.get_xlim()[0], -34), fontsize=8, style="italic")
    ax.set_xlabel("lateral position (Å)")
    ax.set_ylabel("position along membrane normal, z (Å)")
    ax.set_title("Oriented MOR–ZH853 with modeled cholesterol")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)

    # marginal z-histogram: cholesterol distribution defines the midplane
    bins = np.linspace(-40, 32, 37)
    axh.hist(rec.positions[:, 2], bins=bins, orientation="horizontal", color="#9aa0a6", alpha=0.6,
             label="Cα")
    axh.hist(clr.positions[:, 2], bins=bins, orientation="horizontal", color="#e8820c", alpha=0.85,
             label="CHL")
    axh.axhline(0, color="#333", lw=1.2)
    axh.axhline(half, color="#e8820c", lw=1, ls="--")
    axh.axhline(-half, color="#e8820c", lw=1, ls="--")
    axh.set_xlabel("count")
    axh.set_title("z-distribution")

    fig.tight_layout()
    paths.ensure_dir(paths.PRODUCT)
    pdf = paths.PRODUCT / f"02.06.00_membrane_placement_{today}.pdf"
    fig.savefig(pdf)
    fig.savefig(paths.ensure_dir(paths.PRODUCT / "manuscript" / "figures")
                / "fig_membrane_determination.pdf")
    plt.close(fig)
    print(f"midplane z=0 (cholesterol mean); hydrophobic half-thickness ±{half:.1f} Å")
    print(f"Wrote {pdf} and manuscript fig_membrane_determination.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
