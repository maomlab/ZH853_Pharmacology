#!/usr/bin/env python3
"""MolViewSpec scene of the OPM-placed membrane (manuscript figure panel A).

Side view of the oriented MOR-ZH853 complex (receptor + cholesterols + ZH853). The bilayer band
is drawn by trim_figures.py behind the protein at the OPM hydrophobic thickness; here we also emit
`membrane_calib.json` (OPM half-thickness from the reference DUM atoms + the cholesterol z-range) so
the band can be calibrated from pixels to Angstroms in the render.
"""

from __future__ import annotations

import base64
import json
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import molviewspec as mvs  # noqa: E402
import MDAnalysis as mda  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PDB = REPO / "intermediate" / "02.05.00_oriented" / "complex_oriented.pdb"
OPM_REF = REPO / "data" / "opm" / "6ddf.pdb"


def opm_half_thickness() -> float:
    """Half of the OPM hydrophobic thickness from the reference DUM boundary atoms (fallback 15.7)."""
    if OPM_REF.exists():
        dum = mda.Universe(str(OPM_REF)).select_atoms("resname DUM")
        if len(dum):
            return float(np.ptp(dum.positions[:, 2])) / 2.0
    return 15.7


def main() -> int:
    u = mda.Universe(str(PDB))
    clr = u.select_atoms("resname CLR")
    opm_half = opm_half_thickness()

    # calibration for trim_figures: map cholesterol pixel extent -> z, then draw band at +/- opm_half
    (HERE / "membrane_calib.json").write_text(json.dumps({
        "opm_half": opm_half,
        "chol_z_min": float(clr.positions[:, 2].min()),
        "chol_z_max": float(clr.positions[:, 2].max()),
    }))

    data_url = "data:text/plain;base64," + base64.b64encode(PDB.read_bytes()).decode()
    b = mvs.create_builder()
    struct = b.download(url=data_url).parse(format="pdb").model_structure()

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
    print(f"wrote {out.name}; OPM hydrophobic thickness {2 * opm_half:.1f} A (band drawn at trim)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
