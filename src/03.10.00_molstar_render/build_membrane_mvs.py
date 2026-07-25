#!/usr/bin/env python3
"""MolViewSpec scene of the determined membrane placement (manuscript figure panel A).

Side view of the oriented MOR-ZH853 complex with the modeled cholesterols (which define the
bilayer) and a semi-transparent slab drawn at the hydrophobic core (z = +/- cholesterol span/2,
centered on the cholesterol midplane at z=0). Emits membrane.mvsj for render.js.
"""

from __future__ import annotations

import base64
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import molviewspec as mvs  # noqa: E402
import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PDB = REPO / "intermediate" / "02.05.00_oriented" / "complex_oriented.pdb"


def main() -> int:
    u = mda.Universe(str(PDB))
    clr = u.select_atoms("resname CLR")
    half = float(np.ptp(clr.positions[:, 2])) / 2.0            # hydrophobic half-thickness (report only)

    data_url = "data:text/plain;base64," + base64.b64encode(PDB.read_bytes()).decode()
    b = mvs.create_builder()
    struct = b.download(url=data_url).parse(format="pdb").model_structure()

    # Receptor cartoon; the 3 modeled cholesterols (84 atoms; orange) mark the bilayer directly -- explicit
    # slab/midplane is quantified in the companion determination plot (src/02.06.00). A MolStar box
    # primitive was tried for the slab but does not render legibly at usable transparency.
    (struct.component(selector="polymer").representation(type="cartoon")
        .color(color="#9aa0a6").opacity(opacity=0.7))
    (struct.component(selector=mvs.ComponentExpression(auth_comp_id="CLR"))
        .representation(type="ball_and_stick").color(color="#e8820c"))
    (struct.component(selector=mvs.ComponentExpression(auth_asym_id="E"))
        .representation(type="ball_and_stick").color(color="#2ca02c"))

    # side view: look along +y with z (membrane normal) up
    struct.component(selector="all").focus(direction=(0, 1, 0), up=(0, 0, 1), radius_factor=1.1)

    out = HERE / "membrane.mvsj"
    state = b.get_state()
    out.write_text(state.dumps() if hasattr(state, "dumps") else state.json())
    print(f"wrote {out} (half-thickness {half:.1f} A; slab {xext:.0f}x{yext:.0f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
