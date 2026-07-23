#!/usr/bin/env python3
"""2D atom-colored vector depictions of ZH853 and analogs (Objective 3 / manuscript).

Draws ZH853 and the ZH850/831/809 analogs as a scaffold-aligned 2D grid with standard atom
coloring, using RDKit's SVG (vector) backend, and converts to PDF for LaTeX inclusion.

Run: ``python src/05.03.00_ligand_depictions.py``  (or ``make depictions``).
"""

from __future__ import annotations

import sys
from datetime import date

sys.path.insert(0, __file__.rsplit("/src/", 1)[0] + "/src")

import cairosvg  # noqa: E402
from rdkit import RDLogger  # noqa: E402
from rdkit.Chem import rdCoordGen, rdFMCS  # noqa: E402
from rdkit.Chem.AllChem import GenerateDepictionMatching2DStructure  # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402

from zh853mor import chem, paths  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def aligned_mols() -> tuple[list, list[str]]:
    """Return analog mols with 2D coords aligned to the ZH853 shared scaffold, + legends."""
    names = list(chem.ANALOGS)
    mols = [chem.mol(n) for n in names]
    for m in mols:
        rdCoordGen.AddCoords(m)
    ref = mols[0]  # ZH853
    mcs = rdFMCS.FindMCS(mols, completeRingsOnly=True, ringMatchesRingOnly=True, timeout=30)
    core = mcs.queryMol
    if ref.HasSubstructMatch(core):
        for m in mols[1:]:
            if m.HasSubstructMatch(core):
                try:
                    GenerateDepictionMatching2DStructure(m, ref, refPatt=core)
                except (ValueError, RuntimeError):
                    pass
    legends = [f"{n}  ({chem.descriptors_2d(n, m).mw:.0f} Da)" for n, m in zip(names, mols, strict=False)]
    return mols, legends


def main() -> int:
    today = f"{date.today():%Y%m%d}"
    mols, legends = aligned_mols()

    d = rdMolDraw2D.MolDraw2DSVG(820, 760, 410, 380)  # total W,H; panel W,H -> 2x2 grid
    opts = d.drawOptions()
    opts.addStereoAnnotation = True
    opts.legendFontSize = 20
    opts.bondLineWidth = 2
    d.DrawMolecules(mols, legends=legends)
    d.FinishDrawing()
    svg = d.GetDrawingText()

    paths.ensure_dir(paths.PRODUCT)
    svg_path = paths.PRODUCT / f"05.03.00_ligand_depictions_{today}.svg"
    pdf_path = paths.PRODUCT / f"05.03.00_ligand_depictions_{today}.pdf"
    svg_path.write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(pdf_path))

    # copy to manuscript figures (stable name)
    fig = paths.ensure_dir(paths.PRODUCT / "manuscript" / "figures") / "fig4_ligand_depictions.pdf"
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(fig))

    print(f"Drew {len(mols)} analogs -> {svg_path.name}, {pdf_path.name}, {fig.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
