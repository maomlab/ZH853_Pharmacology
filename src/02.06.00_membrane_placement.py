#!/usr/bin/env python3
"""Membrane-placement determination + validation figure (Phase 2 / manuscript).

Side view (lateral axis vs membrane normal z) of the oriented complex. The membrane position is
cross-checked three ways and compared to experiment/OPM:
  - modeled cholesterols (3 molecules, 84 atoms) -- a weak, site-specific marker;
  - the Trp/Tyr aromatic "girdle" (interfacial snorkeling residues) -- a robust structural marker;
  - literature: POPC hydrocarbon core 2Dc = 28.8 A (Kucerka 2011); OPM hydrophobic thickness for
    MOR = 32.0 +/- 1.0 A (4DKL); class-A GPCRs 31-35 A.
Best practice is OPM/PPM energy minimization (Lomize 2012, 2022); the cholesterol/aromatic checks
here confirm the placement is reasonable and set the build slab to ~31-32 A.

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
TRP_RING = ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"]
TYR_RING = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]
OPM_MOR = 32.0            # OPM hydrophobic thickness for MOR (4DKL), +/-1.0 A
POPC_CORE = 28.8         # POPC hydrocarbon core 2Dc (Kucerka 2011)


def aromatic_belt(u) -> tuple[np.ndarray, float, float]:
    """Return (ring-centroid z of Trp+Tyr, lower interface z, upper interface z)."""
    zc = []
    for resn, ring in (("TRP", TRP_RING), ("TYR", TYR_RING)):
        for res in u.select_atoms(f"segid R and resname {resn}").residues:
            rr = res.atoms.select_atoms("name " + " ".join(ring))
            if len(rr) >= 5:
                zc.append(float(rr.positions[:, 2].mean()))
    z = np.array(zc)
    upper = float(z[(z > 7) & (z < 22)].mean())    # extracellular interfacial band
    lower = float(z[(z < -7) & (z > -22)].mean())  # intracellular interfacial band
    return z, lower, upper


def main() -> int:
    if not ORIENTED.exists():
        print(f"ERROR: {ORIENTED} not found -- run `make prep-orient` first.", file=sys.stderr)
        return 1
    u = structure.load(ORIENTED)
    rec = u.select_atoms("protein and name CA")
    clr = u.select_atoms("resname CLR")
    lig = u.select_atoms("resname L01")

    chol_half = float(np.ptp(clr.positions[:, 2])) / 2.0
    arom_z, lo, up = aromatic_belt(u)
    arom_thick = up - lo
    today = f"{date.today():%Y%m%d}"

    fig, (ax, axh) = plt.subplots(1, 2, figsize=(9, 6), gridspec_kw={"width_ratios": [3.2, 1]},
                                  sharey=True)

    # side view: lateral x vs membrane normal z
    ax.scatter(rec.positions[:, 0], rec.positions[:, 2], s=9, c="#9aa0a6", label="receptor Cα", zorder=2)
    ax.scatter(clr.positions[:, 0], clr.positions[:, 2], s=13, c="#e8820c", label="cholesterol", zorder=3)
    ax.scatter(lig.positions[:, 0], lig.positions[:, 2], s=15, c="#2ca02c", label="ZH853", zorder=4)
    # aromatic girdle markers along the far side
    xr = ax.get_xlim() if False else (rec.positions[:, 0].min(), rec.positions[:, 0].max())
    ax.scatter(np.full(arom_z.shape, xr[1] + 2), arom_z, marker="D", s=22, c="#7b3fa0",
               label="Trp/Tyr ring", zorder=5, clip_on=False)

    # cholesterol slab (thin) and aromatic/OPM slab (recommended)
    ax.axhspan(-chol_half, chol_half, color="#f2c811", alpha=0.13, zorder=0)
    ax.axhline(0, color="#333", lw=1.0, ls="-", zorder=1)
    for zc in (chol_half, -chol_half):
        ax.axhline(zc, color="#e8820c", lw=1, ls="--", zorder=1)
    for zc in (up, lo):
        ax.axhline(zc, color="#7b3fa0", lw=1.2, ls=":", zorder=1)
    ax.annotate("extracellular", xy=(xr[0], 27), fontsize=8, style="italic")
    ax.annotate("intracellular", xy=(xr[0], -35), fontsize=8, style="italic")
    ax.set_xlabel("lateral position (Å)")
    ax.set_ylabel("position along membrane normal, z (Å)")
    ax.set_title("Membrane placement: cholesterol vs aromatic girdle")
    ax.legend(loc="lower left", fontsize=7.5, framealpha=0.9)

    # thickness comparison box
    txt = ("hydrophobic thickness\n"
           f"cholesterol span:  {2 * chol_half:.0f} Å\n"
           f"Trp/Tyr girdle:    {arom_thick:.0f} Å\n"
           f"OPM MOR (4DKL):    {OPM_MOR:.0f}±1 Å\n"
           f"POPC 2Dc (expt):   {POPC_CORE:.0f} Å")
    ax.text(0.98, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            family="monospace", bbox={"boxstyle": "round", "fc": "white", "ec": "#999", "alpha": 0.9})

    # marginal z-histogram
    bins = np.linspace(-40, 32, 37)
    axh.hist(rec.positions[:, 2], bins=bins, orientation="horizontal", color="#9aa0a6", alpha=0.6)
    axh.hist(clr.positions[:, 2], bins=bins, orientation="horizontal", color="#e8820c", alpha=0.85)
    axh.axhline(0, color="#333", lw=1.0)
    for zc in (up, lo):
        axh.axhline(zc, color="#7b3fa0", lw=1.2, ls=":")
    axh.set_xlabel("count")
    axh.set_title("z-distribution")

    fig.tight_layout()
    paths.ensure_dir(paths.PRODUCT)
    pdf = paths.PRODUCT / f"02.06.00_membrane_placement_{today}.pdf"
    fig.savefig(pdf)
    fig.savefig(paths.ensure_dir(paths.PRODUCT / "manuscript" / "figures")
                / "fig_membrane_determination.pdf")
    plt.close(fig)
    print(f"cholesterol span {2 * chol_half:.1f} A (midplane {clr.positions[:, 2].mean():.1f}); "
          f"aromatic girdle {arom_thick:.1f} A (interfaces {lo:.1f}, {up:.1f})")
    print(f"OPM MOR reference {OPM_MOR:.0f}+/-1 A; POPC 2Dc {POPC_CORE:.0f} A -> build slab ~31-32 A")
    print(f"Wrote {pdf} and manuscript fig_membrane_determination.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
